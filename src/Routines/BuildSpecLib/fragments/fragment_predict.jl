# Copyright (C) 2024 Nathan Wamsley
#
# This file is part of Pioneer.jl
#
# Pioneer.jl is free software: you can redistribute it and/or modify
# it under the terms of the GNU Affero General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU Affero General Public License for more details.
#
# You should have received a copy of the GNU Affero General Public License
# along with this program. If not, see <https://www.gnu.org/licenses/>.

# src/fragments/fragment_predict.jl

using JSON

"""
    load_fragment_mod_dictionaries(config_path::String)

Load structural and isotope modification dictionaries from a configuration JSON file and
return the fragment filtering parameters required for decoy cloning.
"""
function load_fragment_mod_dictionaries(config_path::String)
    params = JSON.parsefile(config_path)

    structural_mod_to_mass = Dict{String, Float32}()
    for (mass, name) in zip(
        get(params["fixed_mods"], "mass", Float64[]),
        get(params["fixed_mods"], "name", String[])
    )
        structural_mod_to_mass[name] = Float32(mass)
    end

    if haskey(params, "variable_mods")
        for (mass, name) in zip(
            get(params["variable_mods"], "mass", Float64[]),
            get(params["variable_mods"], "name", String[])
        )
            structural_mod_to_mass[name] = Float32(mass)
        end
    end

    iso_mods_dict = Dict{String, Dict{String, Float32}}()
    for mod_group in get(params, "isotope_mod_groups", Any[])
        name = mod_group["name"]
        iso_mods_dict[name] = Dict{String, Float32}()
        for channel in get(mod_group, "channels", Any[])
            iso_mods_dict[name][channel["channel"]] = Float32(channel["mass"])
        end
    end

    if get(params, "channel_decoys", false)
        for mod_group in get(params, "decoy_isotope_mod_groups", Any[])
            name = mod_group["name"]
            dict = get!(iso_mods_dict, name, Dict{String, Float32}())
            for channel in get(mod_group, "channels", Any[])
                dict[channel["channel"]] = Float32(channel["mass"])
            end
        end
    end

    mods_to_sulfur_diff = Dict{String, Int8}()
    for mod_group in get(params, "sulfur_mod_groups", Any[])
        mods_to_sulfur_diff[mod_group["name"]] = Int8(mod_group["sulfur_count"])
    end

    library_params = get(params, "library_params", Dict{String, Any}())
    include_immonium = get(library_params, "include_immonium", true)
    max_frag_rank = Int(get(library_params, "max_frag_rank", 255))
    length_to_frag_count_multiple = Float32(get(library_params, "length_to_frag_count_multiple", 255))

    return structural_mod_to_mass,
           iso_mods_dict,
           mods_to_sulfur_diff,
           include_immonium,
           max_frag_rank,
           length_to_frag_count_multiple
end

function get_fragment_annotation_info(
    annotation,
    model::KoinaModelType,
    ion_dictionary::Union{Nothing, Dict{Int32, String}},
    cache::Dict{Any, PioneerFragAnnotation},
    immonium_to_sulfur_count::Dict{String, Int8} = Dict{String, Int8}()
)
    if haskey(cache, annotation)
        return cache[annotation]
    end

    frag_annotation = if model isa SplineCoefficientModel
        if ion_dictionary === nothing
            error("Ion dictionary required for spline coefficient models")
        end
        ion_name = ion_dictionary[Int32(annotation)]
        UniSpecFragAnnotation(ion_name)
    elseif model isa InstrumentSpecificModel
        UniSpecFragAnnotation(String(annotation))
    else
        GenericFragAnnotation(String(annotation))
    end

    info = parse_fragment_annotation(
        frag_annotation;
        immonium_to_sulfur_count=immonium_to_sulfur_count,
    )
    cache[annotation] = info
    return info
end

function clone_decoy_fragments(
    peptides_df::DataFrame,
    target_fragments::DataFrame,
    model_type::KoinaModelType,
    config_path::String
)
    isempty(target_fragments) && return DataFrame()

    pair_to_target = Dict{UInt32, UInt32}()
    for (idx, row) in enumerate(eachrow(peptides_df))
        pair_to_target[row.pair_id] = UInt32(idx)
    end

    target_groups = Dict{UInt32, DataFrame}()
    if !isempty(target_fragments)
        for subdf in groupby(target_fragments, :precursor_idx)
            pid = first(subdf.precursor_idx)
            target_groups[pid] = DataFrame(subdf)
        end
    end

    structural_mod_to_mass,
    iso_mods_dict,
    _,
    include_immonium,
    max_frag_rank,
    length_to_frag_count_multiple = load_fragment_mod_dictionaries(config_path)
    immonium_to_sulfur_count = get_immonium_sulfur_dict(asset_path("immonium.txt"))

    ion_dictionary = nothing
    if model_type isa SplineCoefficientModel
        ion_dictionary = get_altimeter_ion_dict(joinpath(@__DIR__, "..", "..", "..", "..", "assets", "ion_dictionary.txt"))
    end

    aa_masses = zeros(Float32, 255)
    structural_mod_masses = zeros(Float32, 255)
    iso_mod_masses = zeros(Float32, 255)
    annotation_cache = Dict{Any, PioneerFragAnnotation}()

    mods_column = hasproperty(peptides_df, :mods) ? :mods : :structural_mods
    iso_column = hasproperty(peptides_df, :isotope_mods) ? :isotope_mods : :isotopic_mods

    decoy_tables = DataFrame[]
    for (idx, row) in enumerate(eachrow(peptides_df))
        target_idx = pair_to_target[row.pair_id]

        frag_df = copy(target_groups[target_idx])
        frag_df[!, :precursor_idx] .= UInt32(idx)

        filter_fragments!(frag_df, model_type)

        sequence = row.sequence
        struct_mods = if row[mods_column] !== missing
            String(row[mods_column])
        else
            ""
        end
        iso_mods = if row[iso_column] !== missing
            String(row[iso_column])
        else
            ""
        end

        get_aa_masses!(aa_masses, sequence)
        get_structural_mod_masses!(structural_mod_masses, struct_mods, structural_mod_to_mass)
        getIsoModMasses!(iso_mod_masses, struct_mods, iso_mods, iso_mods_dict)
        seq_length = UInt8(length(sequence))

        keep_rows = Int[]
        frag_infos = PioneerFragAnnotation[]
        for (row_idx, frag_row) in enumerate(eachrow(frag_df))
            info = get_fragment_annotation_info(
                frag_row.annotation,
                model_type,
                ion_dictionary,
                annotation_cache,
                immonium_to_sulfur_count,
            )
            if info.immonium && !include_immonium
                continue
            end
            push!(keep_rows, row_idx)
            push!(frag_infos, info)
        end

        isempty(keep_rows) && continue

        frag_df = frag_df[keep_rows, :]

        max_allowed = min(
            max_frag_rank,
            round(Int, seq_length * length_to_frag_count_multiple) + 1,
        )

        if max_allowed <= 0
            continue
        end

        if nrow(frag_df) > max_allowed
            order = if hasproperty(frag_df, :intensities)
                sortperm(frag_df.intensities; rev=true)
            elseif hasproperty(frag_df, :ranking)
                sortperm(frag_df.ranking)
            elseif hasproperty(frag_df, :rank)
                sortperm(frag_df.rank)
            else
                collect(1:nrow(frag_df))
            end
            order = order[1:max_allowed]
            frag_df = frag_df[order, :]
            frag_infos = frag_infos[order]
        end

        for (frag_row, info) in zip(eachrow(frag_df), frag_infos)
            start_idx, stop_idx = get_fragment_indices(info.base_type, info.frag_index, seq_length)
            frag_row.mz = get_fragment_mz(
                start_idx,
                stop_idx,
                info.base_type,
                info.charge,
                aa_masses,
                structural_mod_masses,
                iso_mod_masses
            )
        end

        push!(decoy_tables, frag_df)
    end

    return isempty(decoy_tables) ? DataFrame() : vcat(decoy_tables...)
end

"""
    predict_fragments(
        peptide_table_path::String,
        frags_out_path::String,
        model_type::KoinaModelType,
        instrument_type::String,
        max_koina_batches::Int,
        batch_size::Int,
        model_name::String;
        intensity_threshold::Float32 = 0.001f0
    )

General fragment prediction dispatcher that routes to appropriate method based on model type.
"""
function predict_fragments(
    peptide_table_path::String,
    frags_out_path::String,
    model_type::KoinaModelType,
    instrument_type::String,
    max_koina_batches::Int,
    batch_size::Int,
    model_name::String;
    intensity_threshold::Float32 = 0.001f0
)
    if !haskey(KOINA_URLS, model_name)
        error("Invalid model name: $model_name. Valid options: $(join(keys(KOINA_URLS), ", "))")
    end

    peptides_df = DataFrame(Arrow.Table(peptide_table_path))
    decoy_mask = hasproperty(peptides_df, :decoy) ? peptides_df.decoy : falses(nrow(peptides_df))
    target_rows = findall(!, decoy_mask)

    if isempty(target_rows)
        @warn "No target precursors available for fragment prediction."
        Arrow.write(frags_out_path, DataFrame())
        return
    end

    target_df = peptides_df[target_rows, :]
    target_precursor_indices = UInt32.(target_rows)

    koina_pool_size = max_koina_batches * 5
    batch_size = min(batch_size, 1000)
    n_targets = length(target_rows)
    batch_start_idxs = collect(1:batch_size*koina_pool_size:n_targets)

    rm(frags_out_path, force=true)

    target_batches = DataFrame[]
    for start_idx in ProgressBar(batch_start_idxs)
        stop_idx = min(start_idx + batch_size*koina_pool_size - 1, n_targets)
        batch_df = target_df[start_idx:stop_idx, :]
        batch_indices = target_precursor_indices[start_idx:stop_idx]

        frags_out = predict_fragments_batch(
            batch_df,
            model_type,
            instrument_type,
            batch_size,
            max_koina_batches,
            batch_indices,
        )

        push!(target_batches, frags_out)
    end

    target_fragments = isempty(target_batches) ? DataFrame() : vcat(target_batches...)

    config_path = joinpath(dirname(frags_out_path), "config.json")
    decoy_fragments = clone_decoy_fragments(peptides_df, target_fragments, model_type, config_path)

    fragments_df = isempty(decoy_fragments) ? target_fragments : vcat(target_fragments, decoy_fragments)
    sort_fragments!(fragments_df)

    Arrow.write(frags_out_path, fragments_df; file=false)
end

"""
Fragment prediction for instrument-specific models (e.g., UniSpec, AlphaPeptDeep).
"""
function predict_fragments_batch(
    peptides_df::DataFrame,
    model::InstrumentSpecificModel,
    instrument_type::String,
    batch_size::Int,
    concurrent_koina_requests::Int,
    precursor_indices::AbstractVector{UInt32}
)::DataFrame
    # Verify instrument compatibility
    if instrument_type ∉ MODEL_CONFIGS[model.name].instruments
        error("Invalid instrument: $instrument_type for model $(model.name). Valid names are ", MODEL_CONFIGS[model.name].instruments)
    end

    # Prepare batches
    json_batches = prepare_koina_batch(
        model,
        peptides_df,
        instrument_type,
        batch_size=batch_size
    )
    # Request predictions
    responses = make_koina_batch_requests(json_batches, KOINA_URLS[model.name]; concurrency=concurrent_koina_requests)
    # Process responses
    batch_dfs = []
    for (i, response) in enumerate(responses)
        batch_result = parse_koina_batch(model, response)
        start_idx = (i-1) * batch_size + 1
        stop_idx = min(i * batch_size, length(precursor_indices))
        batch_precursor_idxs = precursor_indices[start_idx:stop_idx]

        batch_df = batch_result.fragments
        batch_df[!, :precursor_idx] = repeat(batch_precursor_idxs, inner=batch_result.frags_per_precursor)
        # Filter and sort fragments
        filter_fragments!(batch_df, model)
        push!(batch_dfs, batch_df)
    end

    # Combine and filter results
    fragments_df = vcat(batch_dfs...)
    
    sort_fragments!(fragments_df)

    return fragments_df
end

"""
Fragment prediction for instrument-agnostic models (e.g., Prosit).
"""
function predict_fragments_batch(
    peptides_df::DataFrame,
    model::InstrumentAgnosticModel,
    _::String,  # instrument type not used
    batch_size::Int,
    concurrent_koina_requests::Int,
    precursor_indices::AbstractVector{UInt32}
)::DataFrame
    # Prepare batches (no instrument type needed)
    json_batches = prepare_koina_batch(
        model,
        peptides_df,
        batch_size=batch_size
    )

    # Request predictions
    responses = make_koina_batch_requests(json_batches, KOINA_URLS[model.name]; concurrency=concurrent_koina_requests)

    # Process responses
    batch_dfs = []
    for (i, response) in enumerate(responses)
        batch_result = parse_koina_batch(model, response)
        start_idx = (i-1) * batch_size + 1
        stop_idx = min(i * batch_size, length(precursor_indices))
        batch_precursor_idxs = precursor_indices[start_idx:stop_idx]

        batch_df = batch_result.fragments
        batch_df[!, :precursor_idx] = repeat(batch_precursor_idxs, inner=batch_result.frags_per_precursor)
        # Filter and sort fragments
        filter_fragments!(batch_df, model)
        push!(batch_dfs, batch_df)
    end

    fragments_df = vcat(batch_dfs...)
    sort_fragments!(fragments_df)

    return fragments_df
end

"""
Fragment prediction for spline coefficient models (e.g., Altimeter).
"""
function predict_fragments_batch(
    peptides_df::DataFrame,
    model::SplineCoefficientModel,
    instrument_type::String,
    batch_size::Int,
    concurrent_koina_requests::Int,
    precursor_indices::AbstractVector{UInt32}
)::DataFrame
    # Similar to InstrumentSpecificModel but handles spline coefficients
    json_batches = prepare_koina_batch(
        model,
        peptides_df,
        instrument_type,
        batch_size=batch_size
    )

    responses = make_koina_batch_requests(json_batches, KOINA_URLS[model.name]; concurrency=concurrent_koina_requests)

    batch_dfs = []
    knot_vectors = []
    
    for (i, response) in enumerate(responses)
        batch_result = parse_koina_batch(model, response)
        start_idx = (i-1) * batch_size + 1
        stop_idx = min(i * batch_size, length(precursor_indices))
        batch_precursor_idxs = precursor_indices[start_idx:stop_idx]

        batch_df = batch_result.fragments
        batch_df[!, :precursor_idx] = repeat(batch_precursor_idxs, inner=batch_result.frags_per_precursor)
        filter_fragments!(batch_df, model)
        push!(batch_dfs, batch_df)
        push!(knot_vectors, batch_result.extra_data)  # Store knot vectors
    end

    fragments_df = vcat(batch_dfs...)
    
    # Verify knot vectors are consistent
    if !all(k == first(knot_vectors) for k in knot_vectors)
        error("Inconsistent knot vectors across batches")
    end
    
    # Store knot vector with the data
    fragments_df[!, :knot_vector] .= Ref(first(knot_vectors))
    #For altimeter fragments are already sorted 
    #sort_fragments!(f)

    return fragments_df
end

"""
Filter fragments based on intensity and other criteria.
"""
function filter_fragments!(df::DataFrame, model::KoinaModelType)
    # Basic filtering common to all models
    filter!(:intensities => x -> x > 0.001f0, df)  # Remove very low intensity
    filter!(:mz => x -> x > 0, df)  # Remove invalid m/z
    
    # Model-specific filtering
    if model isa InstrumentSpecificModel
        filter!(row -> !occursin('i', row.annotation), df)  # Remove isotope peaks
    end
end

"""
Filter fragments based on intensity and other criteria.
"""
function filter_fragments!(df::DataFrame, model::SplineCoefficientModel)
    # Basic filtering common to all models
    #filter!(:coefficients => x -> x > zero(Float32), df)  # Remove very low intensity
    filter!(:mz => x -> x > 0, df)  # Remove invalid m/z
    
    # Model-specific filtering
    if model isa InstrumentSpecificModel
        filter!(row -> !occursin('i', row.annotation), df)  # Remove isotope peaks
    end
end


"""
Sort fragments by intensity within each precursor group.
"""
function sort_fragments!(df::DataFrame)
    nrow(df) <= 1 && return df

    if hasproperty(df, :intensities)
        sort!(df, [:precursor_idx, order(:intensities, rev=true)])
        return df
    elseif hasproperty(df, :ranking)
        sort!(df, [:precursor_idx, :ranking])
        return df
    elseif hasproperty(df, :rank)
        sort!(df, [:precursor_idx, :rank])
        return df
    end

    tmp_col = gensym(:frag_order)
    df[!, tmp_col] = collect(1:nrow(df))
    sort!(df, [:precursor_idx, tmp_col])
    select!(df, filter(col -> col != tmp_col, names(df)))
    return df
end

