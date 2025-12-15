"""
    buildPionLib(spec_lib_path::String,
                y_start_index::UInt8,
                y_start::UInt8,
                b_start_index::UInt8,
                b_start::UInt8,
                include_p_index::Bool,
                include_p::Bool,
                include_isotope::Bool,
                include_immonium::Bool,
                include_internal::Bool,
                include_neutral_diff::Bool,
                max_frag_charge::UInt8,
                max_frag_rank::UInt8,
                min_frag_intensity::Float32,
                rank_to_score::Vector{UInt8},
                rank_to_score_mode::String,
                frag_bounds::FragBoundModel,
                frag_bin_tol_ppm::Float32,
                rt_bin_tol_ppm::Float32,
                model_type::KoinaModelType = InstrumentSpecificModel("default"))

Build a Pioneer spectral library from preprocessed fragment and precursor data.

# Parameters
- `spec_lib_path`: Path to the directory containing input files and where output files will be written
- `y_start_index`: Minimum index for y-ions to include in the fragment index
- `y_start`: Minimum index for y-ions to include in detailed fragments
- `b_start_index`: Minimum index for b-ions to include in the fragment index
- `b_start`: Minimum index for b-ions to include in detailed fragments
- `include_p_index`: Whether to include precursor ions in the fragment index
- `include_p`: Whether to include precursor ions in detailed fragments
- `include_isotope`: Whether to include isotope peaks 
- `include_immonium`: Whether to include immonium ions
- `include_internal`: Whether to include internal fragment ions
- `include_neutral_diff`: Whether to include fragments with neutral losses
- `max_frag_charge`: Maximum fragment charge state to include
- `max_frag_rank`: Maximum number of fragments per precursor (ranked by intensity)
- `min_frag_intensity`: Minimum relative intensity threshold for fragments
- `rank_to_score`: Vector mapping intensity rank to scoring value
- `frag_bounds`: Model defining minimum and maximum fragment m/z bounds
- `frag_bin_tol_ppm`: Fragment binning tolerance in parts per million
- `rt_bin_tol_ppm`: Retention time binning tolerance in parts per million
- `model_type`: Type of prediction model used (default: InstrumentSpecificModel("default"))

# Returns
- `nothing`: Results are written to files in `spec_lib_path`

# Input Files
Requires the following files in `spec_lib_path`:
- fragments_table.arrow: Table of fragment data
- prec_to_frag.arrow: Table mapping precursors to fragments
- precursors_table.arrow: Table of precursor data

# Output Files
Creates the following files in `spec_lib_path`:
- f_index_fragments.arrow: Fragment index
- f_index_rt_bins.arrow: Retention time bins
- f_index_fragment_bins.arrow: Fragment m/z bins
- presearch_f_index_fragments.arrow: Presearch fragment index
- presearch_f_index_rt_bins.arrow: Presearch retention time bins
- presearch_f_index_fragment_bins.arrow: Presearch fragment m/z bins
- detailed_fragments.jld2: Detailed fragment information
- precursor_to_fragment_indices.jld2: Mapping of precursors to fragment indices
"""
function buildPionLib(spec_lib_path::String,
                      y_start_index::UInt8,
                      y_start::UInt8,
                      b_start_index::UInt8,
                      b_start::UInt8,
                      include_p_index::Bool,
                      include_p::Bool,
                      include_isotope::Bool,
                      include_immonium::Bool,
                      include_internal::Bool,
                      include_neutral_diff::Bool,
                      max_frag_charge::UInt8,
                      max_frag_rank::UInt8,
                      length_to_frag_count_multiple::AbstractFloat,
                      min_frag_intensity::Float32,
                      rank_to_score::Vector{UInt8},
                      rank_to_score_mode::String,
                      frag_bounds::FragBoundModel,
                      frag_bin_tol_ppm::Float32,
                      rt_bin_tol_ppm::Float32,
                      model_type::KoinaModelType,
                      )
    fragments_table, prec_to_frag, precursors_table = nothing, nothing, nothing
    try
        fragments_table = Arrow.Table(joinpath(spec_lib_path,"fragments_table.arrow"));
        prec_to_frag = Arrow.Table(joinpath(spec_lib_path,"prec_to_frag.arrow"));
        precursors_table = Arrow.Table(joinpath(spec_lib_path,"precursors_table.arrow"));
    catch e 
        throw(e)
        @error "could not find library..."
        return nothing
    end

    proportional_scores = nothing
    proportional_ranks = nothing
    effective_rank_mode = rank_to_score_mode
    if rank_to_score_mode == "intensity_proportional"
        @warn "Intensity proportional scoring requires spline coefficients; falling back to fixed ranks"
        effective_rank_mode = "fixed"
    end

    #Simple fragments that go into the fragment index
    #println("Get index fragments...")
    simple_frags = getSimpleFrags(
        fragments_table[:mz],
        fragments_table[:is_y],
        fragments_table[:is_b],
        fragments_table[:is_p],
        fragments_table[:fragment_index],
        fragments_table[:charge],
        fragments_table[:isotope],
        fragments_table[:is_internal],
        fragments_table[:is_immonium],
        fragments_table[:has_neutral_diff],
        precursors_table[:mz],
        precursors_table[:irt],
        precursors_table[:prec_charge],#precursor_charge],
        prec_to_frag[:start_idx],
        y_start_index,
        b_start_index,
        include_p_index,
        include_isotope,
        include_immonium,
        include_internal,
        include_neutral_diff,
        max_frag_charge,
        frag_bounds,
        rank_to_score,
        effective_rank_mode,
        proportional_scores,
        proportional_ranks
    );

    #println("Build fragment index...")
    ##########
    #Builds fragment indexes and saves them for the spec_lib_path
    sort!(simple_frags, by = x->x.prec_irt)
    buildFragmentIndex!(
        spec_lib_path,
        simple_frags,
        frag_bin_tol_ppm, #frag_tol_ppm
        rt_bin_tol_ppm,  #irt_bin_width
        index_name = ""
    );
    #println("Build presearch fragment index...")
    sort!(simple_frags, by = x->x.prec_irt)
    buildFragmentIndex!(
        spec_lib_path,
        simple_frags,
        frag_bin_tol_ppm, #frag_tol_ppm
        typemax(Float32),  #irt_tol
        index_name = "presearch_"
    );
    simple_frags = nothing
    GC.gc()

    #println("Get full fragments list...")
    detailed_frags, pid_to_fid = getDetailedFrags(
    fragments_table[:mz],
    fragments_table[:intensity],
    fragments_table[:is_y],
    fragments_table[:is_b],
    fragments_table[:is_p],
    fragments_table[:fragment_index],
    fragments_table[:charge],
    fragments_table[:sulfur_count],
    fragments_table[:ion_type],
    fragments_table[:isotope],
    fragments_table[:is_internal],
    fragments_table[:is_immonium],
    fragments_table[:has_neutral_diff],
    precursors_table[:mz],
    precursors_table[:prec_charge],#:precursor_charge],
    precursors_table[:length],
    prec_to_frag[:start_idx],
    y_start,
    b_start,
    include_p,
    include_isotope,
    include_immonium,
    include_internal,
    include_neutral_diff,
    max_frag_charge,
    frag_bounds,
    max_frag_rank,
    length_to_frag_count_multiple,
    min_frag_intensity,
    model_type
    ,
    effective_rank_mode,
    proportional_ranks
    );
    
    save_detailed_frags(
        joinpath(spec_lib_path, "detailed_fragments.jld2"),
        detailed_frags
    )

    jldsave(
        joinpath(spec_lib_path, "precursor_to_fragment_indices.jld2");
        pid_to_fid
    )

    return nothing
end


"""
    buildPionLib(spec_lib_path::String,
                y_start_index::UInt8,
                y_start::UInt8,
                b_start_index::UInt8,
                b_start::UInt8,
                include_p_index::Bool,
                include_p::Bool,
                include_isotope::Bool,
                include_immonium::Bool,
                include_internal::Bool,
                include_neutral_diff::Bool,
                max_frag_charge::UInt8,
                max_frag_rank::UInt8,
                min_frag_intensity::Float32,
                rank_to_score::Vector{UInt8},
                frag_bounds::FragBoundModel,
                frag_bin_tol_ppm::Float32,
                rt_bin_tol_ppm::Float32,
                model_type::KoinaModelType = InstrumentSpecificModel("default"))

Build a Pioneer spectral library from preprocessed fragment and precursor data.

# Parameters
- `spec_lib_path`: Path to the directory containing input files and where output files will be written
- `y_start_index`: Minimum index for y-ions to include in the fragment index
- `y_start`: Minimum index for y-ions to include in detailed fragments
- `b_start_index`: Minimum index for b-ions to include in the fragment index
- `b_start`: Minimum index for b-ions to include in detailed fragments
- `include_p_index`: Whether to include precursor ions in the fragment index
- `include_p`: Whether to include precursor ions in detailed fragments
- `include_isotope`: Whether to include isotope peaks 
- `include_immonium`: Whether to include immonium ions
- `include_internal`: Whether to include internal fragment ions
- `include_neutral_diff`: Whether to include fragments with neutral losses
- `max_frag_charge`: Maximum fragment charge state to include
- `max_frag_rank`: Maximum number of fragments per precursor (ranked by intensity)
- `min_frag_intensity`: Minimum relative intensity threshold for fragments
- `rank_to_score`: Vector mapping intensity rank to scoring value
- `frag_bounds`: Model defining minimum and maximum fragment m/z bounds
- `frag_bin_tol_ppm`: Fragment binning tolerance in parts per million
- `rt_bin_tol_ppm`: Retention time binning tolerance in parts per million
- `model_type`: Type of prediction model used (default: InstrumentSpecificModel("default"))

# Returns
- `nothing`: Results are written to files in `spec_lib_path`

# Input Files
Requires the following files in `spec_lib_path`:
- fragments_table.arrow: Table of fragment data
- prec_to_frag.arrow: Table mapping precursors to fragments
- precursors_table.arrow: Table of precursor data

# Output Files
Creates the following files in `spec_lib_path`:
- f_index_fragments.arrow: Fragment index
- f_index_rt_bins.arrow: Retention time bins
- f_index_fragment_bins.arrow: Fragment m/z bins
- presearch_f_index_fragments.arrow: Presearch fragment index
- presearch_f_index_rt_bins.arrow: Presearch retention time bins
- presearch_f_index_fragment_bins.arrow: Presearch fragment m/z bins
- detailed_fragments.jld2: Detailed fragment information
- precursor_to_fragment_indices.jld2: Mapping of precursors to fragment indices
"""
function buildPionLib(spec_lib_path::String,
                      y_start_index::UInt8,
                      y_start::UInt8,
                      b_start_index::UInt8,
                      b_start::UInt8,
                      include_p_index::Bool,
                      include_p::Bool,
                      include_isotope::Bool,
                      include_immonium::Bool,
                      include_internal::Bool,
                      include_neutral_diff::Bool,
                      max_frag_charge::UInt8,
                      max_frag_rank::UInt8,
                      length_to_frag_count_multiple::AbstractFloat,
                      min_frag_intensity::Float32,
                      rank_to_score::Vector{UInt8},
                      rank_to_score_mode::String,
                      frag_bounds::FragBoundModel,
                      frag_bin_tol_ppm::Float32,
                      rt_bin_tol_ppm::Float32,
                      model_type::SplineCoefficientModel,
                      )
    fragments_table, prec_to_frag, precursors_table = nothing, nothing, nothing
    try
        fragments_table = Arrow.Table(joinpath(spec_lib_path,"fragments_table.arrow"));
        prec_to_frag = Arrow.Table(joinpath(spec_lib_path,"prec_to_frag.arrow"));
        precursors_table = Arrow.Table(joinpath(spec_lib_path,"precursors_table.arrow"));
    catch e 
        @error "could not find library..."
        return nothing
    end

    proportional_scores = nothing
    proportional_ranks = nothing
    effective_rank_mode = rank_to_score_mode
    if rank_to_score_mode == "intensity_proportional"
        knots = if hasproperty(fragments_table, :knot_vector)
            Tuple(fragments_table[:knot_vector][1])
        elseif isfile(joinpath(spec_lib_path, "spline_knots.jld2"))
            try
                Tuple(load(joinpath(spec_lib_path, "spline_knots.jld2"))["spl_knots"])
            catch
                nothing
            end
        else
            nothing
        end

        proportional_scores, proportional_ranks = compute_proportional_fragment_scores(
            fragments_table[:coefficients],
            fragments_table[:intensity],
            fragments_table[:mz],
            fragments_table[:is_y],
            fragments_table[:is_b],
            fragments_table[:is_p],
            fragments_table[:fragment_index],
            fragments_table[:charge],
            fragments_table[:isotope],
            fragments_table[:is_internal],
            fragments_table[:is_immonium],
            fragments_table[:has_neutral_diff],
            precursors_table[:mz],
            precursors_table[:prec_charge],
            precursors_table[:length],
            prec_to_frag[:start_idx],
            y_start,
            b_start,
            include_p,
            include_isotope,
            include_immonium,
            include_internal,
            include_neutral_diff,
            max_frag_charge,
            frag_bounds,
            max_frag_rank,
            length_to_frag_count_multiple,
            min_frag_intensity,
            rank_to_score,
            knots,
        )

        if proportional_ranks === nothing
            effective_rank_mode = "fixed"
        end
    end

    #Simple fragments that go into the fragment index
    #println("Get index fragments...")
    simple_frags = getSimpleFrags(
        fragments_table[:mz],
        fragments_table[:is_y],
        fragments_table[:is_b],
        fragments_table[:is_p],
        fragments_table[:fragment_index],
        fragments_table[:charge],
        fragments_table[:isotope],
        fragments_table[:is_internal],
        fragments_table[:is_immonium],
        fragments_table[:has_neutral_diff],
        precursors_table[:mz],
        precursors_table[:irt],
        precursors_table[:prec_charge],#precursor_charge],
        prec_to_frag[:start_idx],
        y_start_index,
        b_start_index,
        include_p_index,
        include_isotope,
        include_immonium,
        include_internal,
        include_neutral_diff,
        max_frag_charge,
        frag_bounds,
        rank_to_score,
        effective_rank_mode,
        proportional_scores,
        proportional_ranks
    );

    #println("Build fragment index...")
    ##########
    #Builds fragment indexes and saves them for the spec_lib_path
    sort!(simple_frags, by = x->x.prec_irt)
    buildFragmentIndex!(
        spec_lib_path,
        simple_frags,
        frag_bin_tol_ppm, #frag_tol_ppm
        rt_bin_tol_ppm,  #irt_bin_width
        index_name = ""
    );
    #println("Build presearch fragment index...")
    sort!(simple_frags, by = x->x.prec_irt)
    buildFragmentIndex!(
        spec_lib_path,
        simple_frags,
        frag_bin_tol_ppm, #frag_tol_ppm
        typemax(Float32),  #irt_tol
        index_name = "presearch_"
    );
    simple_frags = nothing
    GC.gc()

    #println("Get full fragments list...")
    detailed_frags, pid_to_fid = getDetailedFrags(
    fragments_table[:mz],
    fragments_table[:coefficients],
    fragments_table[:intensity],
    fragments_table[:is_y],
    fragments_table[:is_b],
    fragments_table[:is_p],
    fragments_table[:fragment_index],
    fragments_table[:charge],
    fragments_table[:sulfur_count],
    fragments_table[:ion_type],
    fragments_table[:isotope],
    fragments_table[:is_internal],
    fragments_table[:is_immonium],
    fragments_table[:has_neutral_diff],
    precursors_table[:mz],
    precursors_table[:prec_charge],#:precursor_charge],
    precursors_table[:length],
    prec_to_frag[:start_idx],
    y_start,
    b_start,
    include_p,
    include_isotope,
    include_immonium,
    include_internal,
    include_neutral_diff,
    max_frag_charge,
    frag_bounds,
    max_frag_rank,
    length_to_frag_count_multiple,
    min_frag_intensity,
    model_type,
    effective_rank_mode,
    proportional_ranks
    );
    
    save_detailed_frags(
        joinpath(spec_lib_path, "detailed_fragments.jld2"),
        detailed_frags
    )

    jldsave(
        joinpath(spec_lib_path, "precursor_to_fragment_indices.jld2");
        pid_to_fid
    )

    return nothing
end

"""
    cleanUpLibrary(spec_lib_path::String)

Remove intermediate files after library building is complete.

# Parameters
- `spec_lib_path::String`: Path to the library directory

# Effects
Removes the following files if they exist:
- fragments_table.arrow
- prec_to_frag.arrow
- precursors.arrow

# Returns
- `nothing`
"""
function cleanUpLibrary(spec_lib_path::String)
    GC.gc()
    for fname in ["fragments_table.arrow", "prec_to_frag.arrow", "precursors.arrow"]
        fpath = joinpath(spec_lib_path, fname)
        if isfile(fpath)
            rm(fpath, force=true)
        end
    end
end

"""
    fragFilter(
        frag_is_y::Bool,
        frag_is_b::Bool,
        frag_is_p::Bool,
        frag_index::UInt8,
        frag_charge::UInt8,
        frag_isotope::UInt8,
        frag_internal::Bool,
        frag_immonium::Bool,
        frag_neutral_diff::Bool,
        frag_mz::Float32,
        frag_bounds::FragBoundModel,
        prec_mz::Float32,
        y_start::UInt8,
        b_start::UInt8,
        include_p::Bool,
        include_isotope::Bool,
        include_immonium::Bool,
        include_internal::Bool,
        include_neutral_diff::Bool,
        max_frag_charge::UInt8)::Bool

Filter fragments based on type, index, charge, and additional criteria.

# Parameters
- `frag_is_y`: Whether the fragment is a y-ion
- `frag_is_b`: Whether the fragment is a b-ion
- `frag_is_p`: Whether the fragment is a precursor ion
- `frag_index`: Index of the fragment in the peptide sequence
- `frag_charge`: Charge state of the fragment
- `frag_isotope`: Isotope state of the fragment
- `frag_internal`: Whether the fragment is an internal fragment
- `frag_immonium`: Whether the fragment is an immonium ion
- `frag_neutral_diff`: Whether the fragment has neutral losses
- `frag_mz`: m/z value of the fragment
- `frag_bounds`: Model defining valid m/z range based on precursor m/z
- `prec_mz`: m/z value of the precursor
- `y_start`: Minimum index for y-ions to include
- `b_start`: Minimum index for b-ions to include
- `include_p`: Whether to include precursor ions
- `include_isotope`: Whether to include isotope peaks
- `include_immonium`: Whether to include immonium ions
- `include_internal`: Whether to include internal fragment ions
- `include_neutral_diff`: Whether to include fragments with neutral losses
- `max_frag_charge`: Maximum fragment charge state to include

# Returns
- `Bool`: true if fragment passes all filters, false otherwise
"""
function fragFilter(
    frag_is_y::Bool,
    frag_is_b::Bool,
    frag_is_p::Bool,
    frag_index::UInt8,
    frag_charge::UInt8,
    frag_isotope::UInt8,
    frag_internal::Bool,
    frag_immonium::Bool,
    frag_neutral_diff::Bool,
    frag_mz::Float32,
    frag_bounds::FragBoundModel,
    prec_mz::Float32,
    y_start::UInt8,
    b_start::UInt8,
    include_p::Bool,
    include_isotope::Bool,
    include_immonium::Bool,
    include_internal::Bool,
    include_neutral_diff::Bool,
    max_frag_charge::UInt8)
    
    #println("frag_neutral_diff $frag_neutral_diff")
    min_frag_mz, max_frag_mz = frag_bounds(prec_mz)
    if (frag_mz < min_frag_mz) | (frag_mz > max_frag_mz)
        return false
    end
    if frag_is_y
        if frag_index < y_start
            return false
        end
    end
    if frag_is_b
        if frag_index < b_start
            return false
        end
    end
    if frag_is_p
        if !include_p
            return false
        end
    end
    if frag_immonium
        if !include_immonium
            return false
        end
    end
    if frag_internal
        if !include_internal
            return false
        end
    end
    if frag_neutral_diff
        if !include_neutral_diff
            return false
        end
    end
    if !include_isotope
        if !iszero(frag_isotope)
            return false
        end
    end
    if frag_charge > max_frag_charge
        return false
    end
    return true
end

"""
    getSimpleFrags(
        frag_mz::AbstractVector{Float32},
        frag_is_y::AbstractVector{Bool},
        frag_is_b::AbstractVector{Bool},
        frag_is_p::AbstractVector{Bool},
        frag_index::AbstractVector{UInt8},
        frag_charge::AbstractVector{UInt8},
        frag_isotope::AbstractVector{UInt8},
        frag_internal::AbstractVector{Bool},
        frag_immonium::AbstractVector{Bool},
        frag_neutral_diff::AbstractVector{Bool},
        precursor_mz::AbstractVector{Float32},
        precursor_irt::AbstractVector{Float32},
        precursor_charge::AbstractVector{UInt8},
        prec_to_frag_idx::AbstractVector{UInt64},
        y_start::UInt8,
        b_start::UInt8,
        include_p::Bool,
        include_isotope::Bool,
        include_immonium::Bool,
        include_internal::Bool,
        include_neutral_diff::Bool,
        max_frag_charge::UInt8,
        frag_bounds::FragBoundModel,
        rank_to_score::Vector{UInt8}
    )::Vector{SimpleFrag{Float32}}

Extract fragments for the fragment index from raw fragment data.

# Parameters
- `frag_mz`: Fragment m/z values
- `frag_is_y`: Whether each fragment is a y-ion
- `frag_is_b`: Whether each fragment is a b-ion  
- `frag_is_p`: Whether each fragment is a precursor ion
- `frag_index`: Index of each fragment in its peptide sequence
- `frag_charge`: Charge state of each fragment
- `frag_isotope`: Isotope state of each fragment
- `frag_internal`: Whether each fragment is an internal fragment
- `frag_immonium`: Whether each fragment is an immonium ion
- `frag_neutral_diff`: Whether each fragment has neutral losses
- `precursor_mz`: m/z values of precursors
- `precursor_irt`: Retention time values of precursors
- `precursor_charge`: Charge states of precursors
- `prec_to_frag_idx`: Index mapping precursors to their fragments
- `y_start`: Minimum index for y-ions to include
- `b_start`: Minimum index for b-ions to include
- `include_p`: Whether to include precursor ions
- `include_isotope`: Whether to include isotope peaks
- `include_immonium`: Whether to include immonium ions
- `include_internal`: Whether to include internal fragment ions
- `include_neutral_diff`: Whether to include fragments with neutral losses
- `max_frag_charge`: Maximum fragment charge state to include
- `frag_bounds`: Model defining valid m/z range based on precursor m/z
- `rank_to_score`: Vector mapping intensity rank to scoring value

# Returns
- Vector of SimpleFrag objects, containing filtered fragments for the index
"""
function compute_proportional_fragment_scores(
    frag_coef::AbstractVector{NTuple{N, Float32}},
    frag_intensity::AbstractVector{Float16},
    frag_mz::AbstractVector{Float32},
    frag_is_y::AbstractVector{Bool},
    frag_is_b::AbstractVector{Bool},
    frag_is_p::AbstractVector{Bool},
    frag_index::AbstractVector{UInt8},
    frag_charge::AbstractVector{UInt8},
    frag_isotope::AbstractVector{UInt8},
    frag_internal::AbstractVector{Bool},
    frag_immonium::AbstractVector{Bool},
    frag_neutral_diff::AbstractVector{Bool},
    precursor_mz::AbstractVector{Float32},
    precursor_charge::AbstractVector{UInt8},
    prec_len::AbstractVector{UInt8},
    prec_to_frag_idx::AbstractVector{UInt64},
    y_start::UInt8,
    b_start::UInt8,
    include_p::Bool,
    include_isotope::Bool,
    include_immonium::Bool,
    include_internal::Bool,
    include_neutral_diff::Bool,
    max_frag_charge::UInt8,
    frag_bounds::FragBoundModel,
    max_frag_rank::UInt8,
    length_to_frag_count_multiple::AbstractFloat,
    min_frag_intensity::AbstractFloat,
    rank_to_score::Vector{UInt8},
    knots::Union{Nothing, NTuple{K, Float32}},
    ) where {N,K}
    if knots === nothing
        return nothing, nothing
    end
    proportional_scores = zeros(UInt8, length(frag_coef))
    proportional_ranks = zeros(UInt8, length(frag_coef))

    degree = length(knots) - length(first(frag_coef)) - 1
    gqx, gqw = getSplineQuadrature(Float32, first(knots), last(knots))
    n_precursors = length(precursor_mz)
    pbar = ProgressBar(total=n_precursors)
    progress_lock = Threads.SpinLock()

    Threads.@threads for pid_int in range(one(Int), n_precursors)
        try
            prec_mz = precursor_mz[pid_int]
            allowed_max = min(max_frag_rank, UInt8(round((prec_len[pid_int]) * length_to_frag_count_multiple) + 1))
            frag_start_idx, frag_stop_idx = prec_to_frag_idx[pid_int], prec_to_frag_idx[pid_int+1] - 1
            candidate_indices = Int[]
            for frag_idx in range(frag_start_idx, frag_stop_idx)
                if !fragFilter(
                        frag_is_y[frag_idx],
                        frag_is_b[frag_idx],
                        frag_is_p[frag_idx],
                        frag_index[frag_idx],
                        frag_charge[frag_idx],
                        frag_isotope[frag_idx],
                        frag_internal[frag_idx],
                        frag_immonium[frag_idx],
                        frag_neutral_diff[frag_idx],
                        frag_mz[frag_idx],
                        frag_bounds,
                        prec_mz,
                        y_start,
                        b_start,
                        include_p,
                        include_isotope,
                        include_immonium,
                        include_internal,
                        include_neutral_diff,
                        max_frag_charge)
                    continue
                end
                if min_frag_intensity > frag_intensity[frag_idx]
                    continue
                end
                push!(candidate_indices, frag_idx)
            end
            if isempty(candidate_indices)
                continue
            end

            aucs = map(candidate_indices) do frag_idx
                splint(knots, frag_coef[frag_idx], degree, gqx, gqw)
            end
            ranked_perm = sortperm(aucs, rev=true)
            top_perm = ranked_perm[1:min(length(ranked_perm), Int(allowed_max))]
            top_indices = candidate_indices[top_perm]
            top_aucs = copy(aucs[top_perm])
            auc_sum = sum(top_aucs)
            top_count = length(top_indices)
            if auc_sum == 0
                top_aucs .= 1 / top_count
                auc_sum = sum(top_aucs)
            end
            weights = top_aucs ./ auc_sum
            target_count = min(top_count, length(rank_to_score))
            target_total = sum(rank_to_score[1:target_count])
            scaled = round.(Int, weights .* target_total)
            diff = target_total - sum(scaled)
            if diff != 0
                adjust_order = sortperm(weights, rev=true)
                for i in 1:abs(diff)
                    adj_idx = adjust_order[((i - 1) % length(adjust_order)) + 1]
                    scaled[adj_idx] += sign(diff)
                end
            end

            for (local_rank, frag_idx) in enumerate(top_indices)
                proportional_scores[frag_idx] = UInt8(clamp(scaled[local_rank], 0, typemax(UInt8)))
                proportional_ranks[frag_idx] = UInt8(local_rank)
            end
        finally
            Threads.lock(progress_lock) do
                update(pbar)
            end
        end
    end
    return proportional_scores, proportional_ranks
end

function getSimpleFrags(
    frag_mz::AbstractVector{Float32},
    frag_is_y::AbstractVector{Bool},
    frag_is_b::AbstractVector{Bool},
    frag_is_p::AbstractVector{Bool},
    frag_index::AbstractVector{UInt8},
    frag_charge::AbstractVector{UInt8},
    frag_isotope::AbstractVector{UInt8},
    frag_internal::AbstractVector{Bool},
    frag_immonium::AbstractVector{Bool},
    frag_neutral_diff::AbstractVector{Bool},
    precursor_mz::AbstractVector{Float32},
    precursor_irt::AbstractVector{Float32},
    precursor_charge::AbstractVector{UInt8},
    prec_to_frag_idx::AbstractVector{UInt64},
    y_start::UInt8,
    b_start::UInt8,
    include_p::Bool,
    include_isotope::Bool,
    include_immonium::Bool,
    include_internal::Bool,
    include_neutral_diff::Bool,
    max_frag_charge::UInt8,
    frag_bounds::FragBoundModel,
    rank_to_score::Vector{UInt8},
    rank_to_score_mode::String = "fixed",
    proportional_scores::Union{Nothing, AbstractVector{UInt8}} = nothing,
    proportional_ranks::Union{Nothing, AbstractVector{UInt8}} = nothing,
    )
    if (length(prec_to_frag_idx) - 1) != (length(precursor_mz))
        #println("mistake")
    end
    #Maximum ranked fragment that can be included in the fragment index
    max_rank_index = length(rank_to_score)
    #Number of precursors
    n_precursors = UInt32(length(precursor_mz))
    simple_frags = Vector{SimpleFrag{Float32}}(undef, n_precursors*max_rank_index)
    simple_frag_idx = 0
    for pid in range(one(UInt32), n_precursors)
        prec_mz = precursor_mz[pid]
        frag_start_idx, frag_stop_idx = prec_to_frag_idx[pid], prec_to_frag_idx[pid+1] - 1
        candidate_indices = Int[]
        for frag_idx in range(frag_start_idx, frag_stop_idx)
            if fragFilter(
                    frag_is_y[frag_idx],
                    frag_is_b[frag_idx],
                    frag_is_p[frag_idx],
                    frag_index[frag_idx],
                    frag_charge[frag_idx],
                    frag_isotope[frag_idx],
                    frag_internal[frag_idx],
                    frag_immonium[frag_idx],
                    frag_neutral_diff[frag_idx],
                    frag_mz[frag_idx],
                    frag_bounds,
                    prec_mz,
                    y_start,
                    b_start,
                    include_p,
                    include_isotope,
                    include_immonium,
                    include_internal,
                    include_neutral_diff,
                    max_frag_charge)==false
                continue
            end
            push!(candidate_indices, frag_idx)
        end

        if rank_to_score_mode == "intensity_proportional" && proportional_ranks !== nothing && proportional_scores !== nothing
            filter!(idx -> proportional_ranks[idx] > 0, candidate_indices)
            sort!(candidate_indices, by = idx -> proportional_ranks[idx])
        end

        rank = 1
        for frag_idx in candidate_indices
            if rank > max_rank_index
                break
            end
            score = rank_to_score_mode == "intensity_proportional" && proportional_scores !== nothing && proportional_scores[frag_idx] > 0 ? proportional_scores[frag_idx] : rank_to_score[rank]
            simple_frag_idx += 1
            simple_frags[simple_frag_idx] = SimpleFrag(
                frag_mz[frag_idx],
                pid,
                precursor_mz[pid],
                precursor_irt[pid],
                precursor_charge[pid],
                score
            )
            rank += 1
        end

    end
    return simple_frags[1:simple_frag_idx]
end

"""
    buildFragmentIndex!(
        folder_out::String,
        frag_ions::Vector{SimpleFrag{T}}, 
        frag_bin_tol_ppm::AbstractFloat, 
        rt_bin_tol::AbstractFloat;
        index_name::String = ""
    ) where {T<:AbstractFloat}

Build a hierarchical fragment index from a list of fragments.

# Parameters
- `folder_out`: Path to the directory where index files will be written
- `frag_ions`: Vector of SimpleFrag objects to index
- `frag_bin_tol_ppm`: Fragment binning tolerance in parts per million
- `rt_bin_tol`: Retention time bin width
- `index_name`: Prefix for output file names (default: "")

# Effects
Creates and writes the following files to `folder_out`:
- `index_name`f_index_fragments.arrow: Fragment entries sorted by retention time and m/z
- `index_name`f_index_rt_bins.arrow: Two-level binning structure for retention time
- `index_name`f_index_fragment_bins.arrow: Fragment m/z bins within each RT bin

# Returns
- `nothing`

# Implementation Notes
The indexing is hierarchical:
1. First bins by retention time
2. Within each RT bin, fragments are sorted by m/z
3. Then bins by m/z within each RT bin
This enables efficient fragment lookup during search.
"""
function buildFragmentIndex!(
                            folder_out::String,
                            frag_ions::Vector{SimpleFrag{T}}, 
                            frag_bin_tol_ppm::AbstractFloat, 
                            rt_bin_tol::AbstractFloat;
                            index_name::String = ""
                            ) where {T<:AbstractFloat}

    function buildFragIndex!(
        index_fragments::Vector{IndexFragment{T}}, 
        rt_bins::Vector{FragIndexBin{T}},
        frag_bins::Vector{FragIndexBin{T}},
        frag_ions::Vector{SimpleFrag{T}}, 
        frag_bin_tol_ppm::AbstractFloat,
        rt_bin_tol::AbstractFloat) where {T<:AbstractFloat}

        start_irt, stop_irt = getIRT(first(frag_ions)), getIRT(first(frag_ions))
        start_idx, stop_idx = 1, 1
        rt_bin_idx = 1
        frag_bin_idx = 1
        #Within each iRT bin (defined by `rt_size`) sort the fragments by precursor_mz

        for (i, frag_ion) in enumerate(frag_ions)
            stop_irt = getIRT(frag_ion)
            stop_idx = i
            #diff_mz = stop_fragmz - start_fragmz
            #mean_mz = (stop_fragmz + start_fragmz)/2
            #Is the difference between the first and last fragment greater than the bin_ppm?
            if ((stop_irt - start_irt) > rt_bin_tol) & (stop_idx != start_idx)#(diff_mz/(mean_mz/1e6) > frag_bin_tol_ppm) & (stop_idx != start_idx)
                stop_idx = i - 1 #i - 1 fragment is the last that should be incluced in the bin
                stop_irt = getIRT(frag_ions[stop_idx])
                #Within the fragment bin Sort by IRT
                sort!(@view(frag_ions[start_idx:stop_idx]), by = x->getMZ(x))#Try stable sorting algorithm for now, alg=QuickSort)
                #Add new fragbin
                #Build fragment bins for the retention time bin 
                first_frag_bin_idx = frag_bin_idx
                frag_bin_idx = buildFragBins!(index_fragments,
                            frag_bins,
                            frag_bin_idx,
                            frag_ions,
                            start_idx,stop_idx,
                            frag_bin_tol_ppm)
                rt_bins[rt_bin_idx] = FragIndexBin(start_irt, 
                                                    stop_irt, #important that stop_idx is i - 1 and not i
                                                        UInt32(first_frag_bin_idx),
                                                        UInt32(frag_bin_idx-1)
                                                        ) #-1 is critical
                rt_bin_idx += 1

                start_idx, stop_idx = i, i
                start_irt = getIRT(frag_ion)
            end
        end


        #Last bin is special case 
        if start_idx != length(frag_ions)
            stop_irt =  getIRT(frag_ions[stop_idx])
            sort!(@view(frag_ions[start_idx:stop_idx]), by = x->getMZ(x))# Try stable sorting for now, alg=QuickSort)
            #Add new fragbin
            #Build RT bins for the frag bin
            first_frag_bin_idx = frag_bin_idx
            frag_bin_idx = buildFragBins!(index_fragments,
                        frag_bins,
                        frag_bin_idx,
                        frag_ions,
                        start_idx,stop_idx,frag_bin_tol_ppm)

            rt_bins[rt_bin_idx] = FragIndexBin(start_irt, 
                                                    stop_irt, #important that stop_idx is i - 1 and not i
                                                    UInt32(first_frag_bin_idx),
                                                    UInt32(frag_bin_idx-1)) #-1 is critical
            rt_bin_idx += 1
        else
            first_frag_bin_idx = frag_bin_idx

            frag_bin_idx = buildFragBins!(index_fragments,
                        frag_bins,
                        frag_bin_idx,
                        frag_ions,
                        start_idx,stop_idx,frag_bin_tol_ppm)

            #Add new fragbin
            rt_bins[rt_bin_idx] = FragIndexBin(start_irt, 
                                                    getIRT(frag_ions[stop_idx]), #important that stop_idx is i - 1 and not i
                                                    UInt32(first_frag_bin_idx),
                                                    UInt32(frag_bin_idx-1)) #-1 is critical
            rt_bin_idx += 1
        end

        return frag_bin_idx, rt_bin_idx
    end

    function buildFragBins!(index_fragments::Vector{IndexFragment{T}}, 
                            frag_bins::Vector{FragIndexBin{T}},
                            frag_bin_idx::Int64,
                            frag_ions::Vector{SimpleFrag{T}}, 
                            start::Int64, 
                            stop::Int64, 
                            frag_bin_tol_ppm::AbstractFloat) where {T<:AbstractFloat}
        start_idx, stop_idx = start, start
        start_fragmz, stop_fragmz = getMZ(frag_ions[start]), getMZ(frag_ions[start])
        for i in range(start, stop)
            frag_ion = frag_ions[i]
            stop_fragmz = getMZ(frag_ion)
            stop_idx = i

            diff_mz = stop_fragmz - start_fragmz
            mean_mz = (stop_fragmz + start_fragmz)/2
            if (diff_mz/(mean_mz/1e6) > frag_bin_tol_ppm) & (stop_idx != start_idx)
                stop_idx = i - 1
                stop_fragmz =  getMZ(frag_ions[stop_idx]) #Need to set before sorting 
                sort!(@view(frag_ions[start_idx:stop_idx]), by = x->getPrecMZ(x))# Try stable sorting for now, alg=QuickSort)
                #Add new rt bin
                frag_bins[frag_bin_idx] = FragIndexBin(start_fragmz, 
                                                        stop_fragmz, #important that stop_idx is i - 1 and not i
                                                    UInt32(start_idx),
                                                    UInt32(stop_idx)
                                                )
                frag_bin_idx += 1
                for idx in range(start_idx, stop_idx)
                    index_fragments[idx] = IndexFragment(
                                                    getPrecID(frag_ions[idx]),
                                                    getPrecMZ(frag_ions[idx]),
                                                    getScore(frag_ions[idx]),
                                                    getPrecCharge(frag_ions[idx])
                                                    )
                end
                start_idx, stop_idx = i, i
                start_fragmz = getMZ(frag_ions[stop_idx])
            end
        end

        #Last bin is special case 
        if start_idx != stop
            stop_fragmz = getMZ(frag_ions[stop_idx])
            sort!(@view(frag_ions[start_idx:stop_idx]), by = x->getPrecMZ(x))# Try stable sorting for now, alg=QuickSort), alg=QuickSort)
            #Add new fragbin
            frag_bins[frag_bin_idx] = FragIndexBin(start_fragmz, 
                        stop_fragmz, #important that stop_idx is i - 1 and not i
                        UInt32(start_idx),
                        UInt32(stop_idx)
                    )
            frag_bin_idx += 1
            for idx in range(start_idx, stop_idx)
                index_fragments[idx] = IndexFragment(
                                                getPrecID(frag_ions[idx]),
                                                getPrecMZ(frag_ions[idx]),
                                                getScore(frag_ions[idx]),
                                                getPrecCharge(frag_ions[idx])
                                                )
            end
        else
            frag_bins[frag_bin_idx] = FragIndexBin(start_fragmz, 
                        getMZ(frag_ions[stop]), #important that stop_idx is i - 1 and not i
                        UInt32(start_idx),
                        UInt32(stop_idx)
                    )
            frag_bin_idx += 1
            index_fragments[stop_idx] = IndexFragment(
                getPrecID(frag_ions[stop_idx]),
                getPrecMZ(frag_ions[stop_idx]),
                getScore(frag_ions[stop_idx]),
                getPrecCharge(frag_ions[stop_idx])
                )
        end
        return frag_bin_idx
    end

    sort!(frag_ions, by = x->getIRT(x))# Try stable sorting for now, alg=QuickSort), alg=QuickSort)
    #diff = getPPM(getMZ(frag_ions[start]), bin_ppm) #ppm tolerance of the current fragment bin
    #Get smallest iRT in the library
    index_fragments = Vector{IndexFragment{T}}(undef, length(frag_ions))
    rt_bins = Vector{FragIndexBin{T}}(undef, length(frag_ions))
    frag_bins = Vector{FragIndexBin{T}}(undef, length(frag_ions))
    #println("building fragment index...")
    frag_bin_idx, rt_bin_idx = buildFragIndex!(index_fragments,
                    rt_bins,
                    frag_bins,
                    frag_ions,
                    frag_bin_tol_ppm,
                    rt_bin_tol)
                    
    fragments = (IndexFragment = index_fragments,)
    rt_bins  = (FragIndexBin = rt_bins[1:rt_bin_idx-1],)
    frag_bins = (FragIndexBin = frag_bins[1:frag_bin_idx-1],)
    #println("writing tables...")
    if !isdir(folder_out)
        mkdir(folder_out)
    end
    Arrow.write(joinpath(folder_out, index_name*"f_index_fragments.arrow"), fragments)
    Arrow.write(joinpath(folder_out, index_name*"f_index_rt_bins.arrow"), rt_bins)
    Arrow.write(joinpath(folder_out, index_name*"f_index_fragment_bins.arrow"), frag_bins)

    return 
end

"""
    Helper function to decide if the fragment should be filtered out. Filters based on maximum rank or
        a constant multiple of the peptide length, whichever is least 
"""
function filterFrag(rank::Int64, prec_len::UInt8, max_frag_rank::UInt8, length_to_frag_count_multiple::AbstractFloat)
    return rank > min(max_frag_rank, round((prec_len)*length_to_frag_count_multiple)+1)
end

"""
    getDetailedFrags(
        frag_mz::AbstractVector{Float32},
        frag_intensity::AbstractVector{Float16},
        frag_is_y::AbstractVector{Bool},
        frag_is_b::AbstractVector{Bool},
        frag_is_p::AbstractVector{Bool},
        frag_index::AbstractVector{UInt8},
        frag_charge::AbstractVector{UInt8},
        frag_sulfur_count::AbstractVector{UInt8},
        frag_ion_type::AbstractVector{UInt16},
        frag_isotope::AbstractVector{UInt8},
        frag_internal::AbstractVector{Bool},
        frag_immonium::AbstractVector{Bool},
        frag_neutral_diff::AbstractVector{Bool},
        precursor_mz::AbstractVector{Float32},
        precursor_charge::AbstractVector{UInt8},
        prec_to_frag_idx::AbstractVector{UInt64},
        y_start::UInt8,
        b_start::UInt8,
        include_p::Bool,
        include_isotope::Bool,
        include_immonium::Bool,
        include_internal::Bool,
        include_neutral_diff::Bool,
        max_frag_charge::UInt8,
        frag_bounds::FragBoundModel,
        max_frag_rank::UInt8,
        min_frag_intensity::AbstractFloat,
        koina_model::KoinaModelType
    )::Tuple{Vector{DetailedFrag{Float32}}, Vector{UInt64}}

Extract detailed fragment information for peptide scoring.

# Parameters
- `frag_mz`: Fragment m/z values
- `frag_intensity`: Fragment intensities
- `frag_is_y`: Whether each fragment is a y-ion
- `frag_is_b`: Whether each fragment is a b-ion
- `frag_is_p`: Whether each fragment is a precursor ion
- `frag_index`: Index of each fragment in its peptide sequence
- `frag_charge`: Charge state of each fragment
- `frag_sulfur_count`: Sulfur count for each fragment
- `frag_ion_type`: Ion type identifier for each fragment
- `frag_isotope`: Isotope state of each fragment
- `frag_internal`: Whether each fragment is an internal fragment
- `frag_immonium`: Whether each fragment is an immonium ion
- `frag_neutral_diff`: Whether each fragment has neutral losses
- `precursor_mz`: m/z values of precursors
- `precursor_charge`: Charge states of precursors
- `prec_to_frag_idx`: Index mapping precursors to their fragments
- `y_start`: Minimum index for y-ions to include
- `b_start`: Minimum index for b-ions to include
- `include_p`: Whether to include precursor ions
- `include_isotope`: Whether to include isotope peaks
- `include_immonium`: Whether to include immonium ions
- `include_internal`: Whether to include internal fragment ions
- `include_neutral_diff`: Whether to include fragments with neutral losses
- `max_frag_charge`: Maximum fragment charge state to include
- `frag_bounds`: Model defining valid m/z range based on precursor m/z
- `max_frag_rank`: Maximum number of fragments per precursor (ranked by intensity)
- `min_frag_intensity`: Minimum relative intensity threshold for fragments
- `koina_model`: Type of prediction model used

# Returns
- Tuple containing:
  - Vector of DetailedFrag objects with comprehensive fragment information
  - Vector mapping precursor indices to fragment indices
"""
function getDetailedFrags(
    frag_mz::AbstractVector{Float32},
    frag_intensity::AbstractVector{Float16},
    frag_is_y::AbstractVector{Bool},
    frag_is_b::AbstractVector{Bool},
    frag_is_p::AbstractVector{Bool},
    frag_index::AbstractVector{UInt8},
    frag_charge::AbstractVector{UInt8},
    frag_sulfur_count::AbstractVector{UInt8},
    frag_ion_type::AbstractVector{UInt16},
    frag_isotope::AbstractVector{UInt8},
    frag_internal::AbstractVector{Bool},
    frag_immonium::AbstractVector{Bool},
    frag_neutral_diff::AbstractVector{Bool},
    precursor_mz::AbstractVector{Float32},
    precursor_charge::AbstractVector{UInt8},
    prec_len::AbstractVector{UInt8},
    prec_to_frag_idx::AbstractVector{UInt64},
    y_start::UInt8,
    b_start::UInt8,
    include_p::Bool,
    include_isotope::Bool,
    include_immonium::Bool,
    include_internal::Bool,
    include_neutral_diff::Bool,
    max_frag_charge::UInt8,
    frag_bounds::FragBoundModel,
    max_frag_rank::UInt8,
    length_to_frag_count_multiple::AbstractFloat,
    min_frag_intensity::AbstractFloat,
    koina_model::KoinaModelType,
    rank_to_score_mode::String = "fixed",
    proportional_ranks::Union{Nothing, AbstractVector{UInt8}} = nothing)

    if (length(prec_to_frag_idx) - 1) != (length(precursor_mz))
        #println("mistake")
    end

    #Maximum ranked fragment that can be included in the fragment index
    max_frag_rank = UInt8(min(255, max_frag_rank))
    #Number of precursors 
    n_precursors = UInt32(length(precursor_mz))
    #Keep track of number of fragments to allocate 
    n_frags = zero(UInt64)
    precursor_fragments = Vector{Vector{Int}}(undef, Int(n_precursors))
    #println("counting fragments...")
    for pid in ProgressBar(range(one(UInt32), n_precursors))
        prec_mz = precursor_mz[pid] #Filter on precursor mass
        frag_start_idx, frag_stop_idx = prec_to_frag_idx[pid], prec_to_frag_idx[pid+1] - 1
        candidate_indices = Int[]
        for frag_idx in range(frag_start_idx, frag_stop_idx)
            if !fragFilter(
                    frag_is_y[frag_idx],
                    frag_is_b[frag_idx],
                    frag_is_p[frag_idx],
                    frag_index[frag_idx],
                    frag_charge[frag_idx],
                    frag_isotope[frag_idx],
                    frag_internal[frag_idx],
                    frag_immonium[frag_idx],
                    frag_neutral_diff[frag_idx],
                    frag_mz[frag_idx],
                    frag_bounds,
                    prec_mz,
                    y_start,
                    b_start,
                    include_p,
                    include_isotope,
                    include_immonium,
                    include_internal,
                    include_neutral_diff,
                    max_frag_charge)
                continue
            end
            if min_frag_intensity > frag_intensity[frag_idx]
                continue
            end
            push!(candidate_indices, frag_idx)
        end
        if rank_to_score_mode == "intensity_proportional" && proportional_ranks !== nothing
            filter!(idx -> proportional_ranks[idx] > 0, candidate_indices)
            sort!(candidate_indices, by = idx -> proportional_ranks[idx])
        end
        if isempty(candidate_indices)
            precursor_fragments[Int(pid)] = Int[]
            continue
        end
        if isempty(candidate_indices)
            precursor_fragments[Int(pid)] = Int[]
            continue
        end
        allowed_max = min(Int(max_frag_rank), Int(round((prec_len[pid]) * length_to_frag_count_multiple) + 1))
        precursor_fragments[Int(pid)] = candidate_indices[1:min(allowed_max, length(candidate_indices))]
        n_frags += UInt64(length(precursor_fragments[Int(pid)]))
    end
    detailed_frags = Vector{DetailedFrag{Float32}}(
                                    undef, 
                                    n_frags)   
    prec_to_frag_idx_new = Vector{UInt64}(undef, n_precursors + 1)
    detailed_frag_idx = 1
    #println("writing fragments...")
    for pid in ProgressBar(range(one(UInt32), n_precursors))
        prec_mz = precursor_mz[pid]
        #Index of the first fragment for the precursor 
        prec_to_frag_idx_new[pid] = UInt64(detailed_frag_idx)
        rank = 1
        for frag_idx in precursor_fragments[Int(pid)]
            is_y, is_internal, is_immonium = frag_is_y[frag_idx], frag_internal[frag_idx], frag_immonium[frag_idx]
            is_b, is_p = frag_is_b[frag_idx], frag_is_p[frag_idx]
            detailed_frags[detailed_frag_idx] = DetailedFrag(
                UInt32(pid),

                frag_mz[frag_idx],
                frag_intensity[frag_idx],

                frag_ion_type[frag_idx],
                is_y,
                is_b,
                is_p,
                frag_isotope[frag_idx]>0,

                frag_charge[frag_idx],
                frag_index[frag_idx],
                precursor_charge[pid],
                UInt8(rank),
                frag_sulfur_count[frag_idx]
            )
            detailed_frag_idx += 1
            rank += 1
        end

    end
    prec_to_frag_idx_new[end] = UInt64(detailed_frag_idx)
    return detailed_frags[1:detailed_frag_idx-1], prec_to_frag_idx_new
end

"""
    getDetailedFrags(
        frag_mz::AbstractVector{Float32},
        frag_coef::AbstractVector{NTuple{N, Float32}},
        frag_intensity::AbstractVector{Float16},
        frag_is_y::AbstractVector{Bool},
        frag_is_b::AbstractVector{Bool},
        frag_is_p::AbstractVector{Bool},
        frag_index::AbstractVector{UInt8},
        frag_charge::AbstractVector{UInt8},
        frag_sulfur_count::AbstractVector{UInt8},
        frag_ion_type::AbstractVector{UInt16},
        frag_isotope::AbstractVector{UInt8},
        frag_internal::AbstractVector{Bool},
        frag_immonium::AbstractVector{Bool},
        frag_neutral_diff::AbstractVector{Bool},
        precursor_mz::AbstractVector{Float32},
        precursor_charge::AbstractVector{UInt8},
        prec_to_frag_idx::AbstractVector{UInt64},
        y_start::UInt8,
        b_start::UInt8,
        include_p::Bool,
        include_isotope::Bool,
        include_immonium::Bool,
        include_internal::Bool,
        include_neutral_diff::Bool,
        max_frag_charge::UInt8,
        frag_bounds::FragBoundModel,
        max_frag_rank::UInt8,
        min_frag_intensity::AbstractFloat,
        koina_model::SplineCoefficientModel
    )::Tuple{Vector{SplineDetailedFrag{N, Float32}}, Vector{UInt64}} where {N}

Extract detailed fragment information for peptide scoring, including spline coefficients.

# Parameters
- `frag_mz`: Fragment m/z values
- `frag_coef`: Spline coefficients for each fragment
- `frag_intensity`: Fragment intensities
- `frag_is_y`: Whether each fragment is a y-ion
- `frag_is_b`: Whether each fragment is a b-ion
- `frag_is_p`: Whether each fragment is a precursor ion
- `frag_index`: Index of each fragment in its peptide sequence
- `frag_charge`: Charge state of each fragment
- `frag_sulfur_count`: Sulfur count for each fragment
- `frag_ion_type`: Ion type identifier for each fragment
- `frag_isotope`: Isotope state of each fragment
- `frag_internal`: Whether each fragment is an internal fragment
- `frag_immonium`: Whether each fragment is an immonium ion
- `frag_neutral_diff`: Whether each fragment has neutral losses
- `precursor_mz`: m/z values of precursors
- `precursor_charge`: Charge states of precursors
- `prec_to_frag_idx`: Index mapping precursors to their fragments
- `y_start`: Minimum index for y-ions to include
- `b_start`: Minimum index for b-ions to include
- `include_p`: Whether to include precursor ions
- `include_isotope`: Whether to include isotope peaks
- `include_immonium`: Whether to include immonium ions
- `include_internal`: Whether to include internal fragment ions
- `include_neutral_diff`: Whether to include fragments with neutral losses
- `max_frag_charge`: Maximum fragment charge state to include
- `frag_bounds`: Model defining valid m/z range based on precursor m/z
- `max_frag_rank`: Maximum number of fragments per precursor (ranked by intensity)
- `min_frag_intensity`: Minimum relative intensity threshold for fragments
- `koina_model`: Type of spline coefficient model used

# Returns
- Tuple containing:
  - Vector of SplineDetailedFrag objects with fragment info including spline coefficients
  - Vector mapping precursor indices to fragment indices

# Type Parameters
- `N`: Number of spline coefficients per fragment
"""
function getDetailedFrags(
    frag_mz::AbstractVector{Float32},
    frag_coef::AbstractVector{NTuple{N, Float32}},
    frag_intensity::AbstractVector{Float16},
    frag_is_y::AbstractVector{Bool},
    frag_is_b::AbstractVector{Bool},
    frag_is_p::AbstractVector{Bool},
    frag_index::AbstractVector{UInt8},
    frag_charge::AbstractVector{UInt8},
    frag_sulfur_count::AbstractVector{UInt8},
    frag_ion_type::AbstractVector{UInt16},
    frag_isotope::AbstractVector{UInt8},
    frag_internal::AbstractVector{Bool},
    frag_immonium::AbstractVector{Bool},
    frag_neutral_diff::AbstractVector{Bool},
    precursor_mz::AbstractVector{Float32},
    precursor_charge::AbstractVector{UInt8},
    prec_len::AbstractVector{UInt8},
    prec_to_frag_idx::AbstractVector{UInt64},
    y_start::UInt8,
    b_start::UInt8,
    include_p::Bool,
    include_isotope::Bool,
    include_immonium::Bool,
    include_internal::Bool,
    include_neutral_diff::Bool,
    max_frag_charge::UInt8,
    frag_bounds::FragBoundModel,
    max_frag_rank::UInt8,
    length_to_frag_count_multiple::AbstractFloat,
    min_frag_intensity::AbstractFloat,
    koina_model::SplineCoefficientModel,
    rank_to_score_mode::String = "fixed",
    proportional_ranks::Union{Nothing, AbstractVector{UInt8}} = nothing) where {N}

    if (length(prec_to_frag_idx) - 1) != (length(precursor_mz))
        #println("mistake")
    end

    #Maximum ranked fragment that can be included in the fragment index
    max_frag_rank = UInt8(min(255, max_frag_rank))
    #Number of precursors 
    n_precursors = UInt32(length(precursor_mz))
    #Keep track of number of fragments to allocate 
    n_frags = zero(UInt64)
    precursor_fragments = Vector{Vector{Int}}(undef, Int(n_precursors))
    #println("counting fragments...")
    for pid in ProgressBar(range(one(UInt32), n_precursors))
        prec_mz = precursor_mz[pid] #Filter on precursor mass
        frag_start_idx, frag_stop_idx = prec_to_frag_idx[pid], prec_to_frag_idx[pid+1] - 1
        candidate_indices = Int[]
        for frag_idx in range(frag_start_idx, frag_stop_idx)
            if !fragFilter(
                    frag_is_y[frag_idx],
                    frag_is_b[frag_idx],
                    frag_is_p[frag_idx],
                    frag_index[frag_idx],
                    frag_charge[frag_idx],
                    frag_isotope[frag_idx],
                    frag_internal[frag_idx],
                    frag_immonium[frag_idx],
                    frag_neutral_diff[frag_idx],
                    frag_mz[frag_idx],
                    frag_bounds,
                    prec_mz,
                    y_start,
                    b_start,
                    include_p,
                    include_isotope,
                    include_immonium,
                    include_internal,
                    include_neutral_diff,
                    max_frag_charge)
                continue
            end
            push!(candidate_indices, frag_idx)
        end
        if rank_to_score_mode == "intensity_proportional" && proportional_ranks !== nothing
            filter!(idx -> proportional_ranks[idx] > 0, candidate_indices)
            sort!(candidate_indices, by = idx -> proportional_ranks[idx])
        end
        allowed_max = min(Int(max_frag_rank), Int(round((prec_len[pid]) * length_to_frag_count_multiple) + 1))
        precursor_fragments[Int(pid)] = candidate_indices[1:min(allowed_max, length(candidate_indices))]
        n_frags += UInt64(length(precursor_fragments[Int(pid)]))
    end
    n_tuple_p = typeof(first(frag_coef))
    n_tuple_size = length(n_tuple_p.parameters)
    n_tuple_type = eltype(n_tuple_p)
    detailed_frags = Vector{SplineDetailedFrag{n_tuple_size, n_tuple_type}}(
                                    undef, 
                                    n_frags)   
    prec_to_frag_idx_new = Vector{UInt64}(undef, n_precursors + 1)
    detailed_frag_idx = 1
    #println("writing fragments...")
    for pid in ProgressBar(range(one(UInt32), n_precursors))
        prec_mz = precursor_mz[pid]
        #Index of the first fragment for the precursor 
        prec_to_frag_idx_new[pid] = UInt64(detailed_frag_idx)
        rank = 1
        for frag_idx in precursor_fragments[Int(pid)]
            is_y, is_internal, is_immonium = frag_is_y[frag_idx], frag_internal[frag_idx], frag_immonium[frag_idx]
            is_b, is_p = frag_is_b[frag_idx], frag_is_p[frag_idx]
            detailed_frags[detailed_frag_idx] = SplineDetailedFrag(
                UInt32(pid),

                frag_mz[frag_idx],
                frag_coef[frag_idx],

                frag_ion_type[frag_idx],
                is_y,
                is_b,
                is_p,
                frag_isotope[frag_idx]>0,

                frag_charge[frag_idx],
                frag_index[frag_idx],
                precursor_charge[pid],
                UInt8(rank),
                frag_sulfur_count[frag_idx]
            )
            detailed_frag_idx += 1
            rank += 1
        end

    end
    prec_to_frag_idx_new[end] = UInt64(detailed_frag_idx)
    return detailed_frags[1:detailed_frag_idx-1], prec_to_frag_idx_new
end
