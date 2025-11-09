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
    # Verify model configuration
    if !haskey(KOINA_URLS, model_name)
        error("Invalid model name: $model_name. Valid options: $(join(keys(KOINA_URLS), ", "))")
    end

    # Load data and split targets/decoys
    peptides_df = DataFrame(Arrow.Table(peptide_table_path))

    if :precursor_idx ∉ names(peptides_df)
        error("Peptide table $(peptide_table_path) is missing required :precursor_idx column")
    end

    if any(ismissing, peptides_df[!, :precursor_idx])
        error("Peptide table $(peptide_table_path) contains missing precursor_idx values")
    end

    peptides_df[!, :precursor_idx] = UInt32.(peptides_df[!, :precursor_idx])

    target_idxs = findall(.!peptides_df.decoy)
    decoy_idxs = findall(peptides_df.decoy)

    # Process in batches for targets only
    koina_pool_size = max_koina_batches * 5
    n_targets = length(target_idxs)
    batch_size = min(batch_size, 1000)
    target_batches = DataFrame[]
    if n_targets > 0
        batch_start_positions = 1:batch_size*koina_pool_size:n_targets
        rm(frags_out_path, force=true)

        for start_pos in ProgressBar(batch_start_positions)
            stop_pos = min(start_pos + batch_size*koina_pool_size - 1, n_targets)
            batch_indices = target_idxs[start_pos:stop_pos]
            batch_df = peptides_df[batch_indices, :]

            frags_out = predict_fragments_batch(
                batch_df,
                model_type,
                instrument_type,
                batch_size,
                max_koina_batches,
                UInt32.(peptides_df[batch_indices, :precursor_idx]),
            )

            push!(target_batches, frags_out)
        end
    end

    target_fragments_df = isempty(target_batches) ? DataFrame() : vcat(target_batches...)

    target_lookup = Dict{UInt32, SubDataFrame}()
    if :precursor_idx ∈ names(target_fragments_df)
        for group in groupby(target_fragments_df, :precursor_idx)
            target_lookup[UInt32(first(group.precursor_idx))] = group
        end
    end

    # Duplicate target predictions for decoys using partner mapping
    all_fragments_df = duplicate_decoy_fragments(
        target_fragments_df,
        target_lookup,
        peptides_df,
        target_idxs,
        decoy_idxs,
    )

    Arrow.write(frags_out_path, all_fragments_df)
end

function duplicate_decoy_fragments(
    target_fragments_df::DataFrame,
    target_lookup::Dict{UInt32, SubDataFrame},
    peptides_df::DataFrame,
    target_indices::AbstractVector{<:Integer},
    decoy_indices::AbstractVector{<:Integer},
)
    # Fast path: no decoys or no targets
    if isempty(decoy_indices) || isempty(target_fragments_df)
        return target_fragments_df
    end

    if :precursor_idx ∉ names(peptides_df)
        error("Peptide metadata missing required :precursor_idx column")
    end

    # Build lookup from (pair_id, charge) -> precursor index for targets
    partner_lookup = Dict{Tuple{UInt32, UInt8}, UInt32}()
    for idx in target_indices
        row_idx = Int(idx)
        pair_val = peptides_df.pair_id[row_idx]
        charge_val = peptides_df.precursor_charge[row_idx]
        precursor_val = peptides_df.precursor_idx[row_idx]
        if ismissing(pair_val) || ismissing(charge_val) || ismissing(precursor_val)
            continue
        end
        partner_lookup[(UInt32(pair_val), UInt8(charge_val))] = UInt32(precursor_val)
    end

    decoy_fragments = DataFrame[]
    for decoy_idx in ProgressBar(decoy_indices)
        row_idx = Int(decoy_idx)
        pair_val = peptides_df.pair_id[row_idx]
        charge_val = peptides_df.precursor_charge[row_idx]
        precursor_val = peptides_df.precursor_idx[row_idx]
        if ismissing(pair_val) || ismissing(charge_val)
            @warn "Decoy precursor $decoy_idx missing pair_id or precursor_charge metadata"
            continue
        end
        if ismissing(precursor_val)
            @warn "Decoy precursor $decoy_idx missing precursor_idx metadata"
            continue
        end
        partner_key = (UInt32(pair_val), UInt8(charge_val))

        if !haskey(partner_lookup, partner_key)
            @warn "No target partner found for decoy precursor $decoy_idx (pair_id=$(pair_val), charge=$(charge_val))"
            continue
        end

        target_precursor = partner_lookup[partner_key]
        subdf = get(target_lookup, target_precursor, nothing)
        if subdf === nothing
            @warn "No fragment predictions found for target precursor $target_precursor when duplicating decoy $decoy_idx"
            continue
        end

        decoy_df = copy(subdf)
        decoy_df[!, :precursor_idx] .= UInt32(precursor_val)
        push!(decoy_fragments, decoy_df)
    end

    if isempty(decoy_fragments)
        return target_fragments_df
    end

    combined = vcat(target_fragments_df, decoy_fragments...)
    sort!(combined, :precursor_idx)
    return combined
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
    precursor_indices::Vector{UInt32}
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
    batch_dfs = Vector{DataFrame}()
    for (i, response) in enumerate(responses)
        batch_result = parse_koina_batch(model, response)
        batch_start = (i-1) * batch_size + 1
        batch_stop = min(i * batch_size, length(precursor_indices))
        current_indices = precursor_indices[batch_start:batch_stop]

        batch_df = batch_result.fragments
        n_precursors_in_batch = length(current_indices)
        batch_df[!, :precursor_idx] = repeat(current_indices, inner=batch_result.frags_per_precursor)
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
    precursor_indices::Vector{UInt32}
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
    batch_dfs = Vector{DataFrame}()
    for (i, response) in enumerate(responses)
        batch_result = parse_koina_batch(model, response)
        batch_start = (i-1) * batch_size + 1
        batch_stop = min(i * batch_size, length(precursor_indices))
        current_indices = precursor_indices[batch_start:batch_stop]

        batch_df = batch_result.fragments
        n_precursors_in_batch = length(current_indices)
        batch_df[!, :precursor_idx] = repeat(current_indices, inner=batch_result.frags_per_precursor)
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
    precursor_indices::Vector{UInt32}
)::DataFrame
    # Similar to InstrumentSpecificModel but handles spline coefficients
    json_batches = prepare_koina_batch(
        model,
        peptides_df,
        instrument_type,
        batch_size=batch_size
    )

    responses = make_koina_batch_requests(json_batches, KOINA_URLS[model.name]; concurrency=concurrent_koina_requests)

    batch_dfs = Vector{DataFrame}()
    knot_vectors = []

    for (i, response) in enumerate(responses)
        batch_result = parse_koina_batch(model, response)
        batch_start = (i-1) * batch_size + 1
        batch_stop = min(i * batch_size, length(precursor_indices))
        current_indices = precursor_indices[batch_start:batch_stop]

        batch_df = batch_result.fragments
        n_precursors_in_batch = length(current_indices)
        batch_df[!, :precursor_idx] = repeat(current_indices, inner=batch_result.frags_per_precursor)
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
    sort!(df, [:precursor_idx, order(:intensities, rev=true)])
end

