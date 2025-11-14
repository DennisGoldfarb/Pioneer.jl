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

"""
    FirstPassSearch

Initial search method to identify PSMs and establish retention time calibration.

This search:
1. Performs initial PSM identification with learned scoring
2. Calculates retention time indices and FWHM statistics
3. Maps between library and empirical retention times
4. Generates RT calibration curves for subsequent searches

# Example Implementation
```julia
# Define search parameters
params = Dict(
    :isotope_err_bounds => (1, 0),
    :first_search_params => Dict(
        "n_train_rounds_probit" => 10,
        "max_iter_probit" => 100,
        "max_q_value_probit_rescore" => 0.01,
        "min_index_search_score" => 3,
        "min_frag_count" => 3,
        "min_spectral_contrast" => 0.1,
        "min_log2_matched_ratio" => -3.0,
        "min_topn_of_m" => (3, 5),
        "max_best_rank" => 3,
        "max_precursors_passing" => 5000,
        "abreviate_precursor_calc" => false
    ),
    :summarize_first_search_params => Dict(
        "min_inference_points" => 1000,
        "max_q_val_for_irt" => 0.01,
        "max_precursors" => 10000,
        "max_irt_bin_size" => 0.1,
        "max_prob_to_impute" => 0.99
    ),
    :irt_mapping_params => Dict(
        "min_prob" => 0.9
    )
)

# Execute search
results = execute_search(FirstPassSearch(), search_context, params)
```
"""
struct FirstPassSearch <: SearchMethod end

#==========================================================
Type Definitions
==========================================================#


"""
Results container for first pass search.
Holds FWHM statistics, PSM file paths, and model updates.
"""
struct FirstPassSearchResults <: SearchResults
    fwhms::Dictionary{Int64, @NamedTuple{median_fwhm::Float32,mad_fwhm::Float32}}
    psms::Base.Ref{DataFrame}
    ms1_mass_err_model::Base.Ref{<:MassErrorModel}
    ms1_ppm_errs::Vector{Float32}
    ms1_mass_plots::Vector{Plots.Plot}
    qc_plots_folder_path::String
    probit_precursor_candidates::Set{UInt32}
end

"""
Parameters for first pass search.
Configures PSM identification, scoring, and RT calibration.
"""
struct FirstPassSearchParameters{P<:PrecEstimation} <: FragmentIndexSearchParameters
    # Core parameters
    isotope_err_bounds::Tuple{UInt8, UInt8}
    min_fraction_transmitted::Float32
    frag_tol_ppm::Float32
    ms1_tol_ppm::Float32
    frag_err_quantile::Float32
    min_index_search_score::UInt8
    min_frag_count::Int64
    min_spectral_contrast::Float32
    min_log2_matched_ratio::Float32
    min_topn_of_m::Tuple{Int64, Int64}
    max_best_rank::UInt8
    n_frag_isotopes::Int64
    max_frag_rank::UInt8
    spec_order::Set{Int64}
    match_between_runs::Bool
    relative_improvement_threshold::Float32
    
    # Scoring parameters
    n_train_rounds_probit::Int64
    max_iter_probit::Int64
    max_q_value_probit_rescore::Float32
    max_PEP::Float32
    
    # RT parameters
    min_inference_points::Int64
    max_q_val_for_irt::Float32
    min_prob_for_irt_mapping::Float32
    max_irt_bin_size::Float32
    max_prob_to_impute::Float32
    fwhm_nstd::Float32
    irt_nstd::Float32
    plot_rt_alignment::Bool
    use_robust_fitting::Bool
    prec_estimation::P

    function FirstPassSearchParameters(params::PioneerParameters)
        # Extract relevant parameter groups
        global_params = params.global_settings
        first_params = params.first_search
        frag_params = first_params.fragment_settings
        score_params = first_params.scoring_settings
        rt_params = params.rt_alignment
        irt_mapping_params = first_params.irt_mapping
        # Convert isotope error bounds
        isotope_bounds = global_params.isotope_settings.err_bounds_first_pass
        # Determine precursor estimation strategy
        prec_estimation = global_params.isotope_settings.partial_capture ? PartialPrecCapture() : FullPrecCapture()
        
        new{typeof(prec_estimation)}(
            (UInt8(first(isotope_bounds)), UInt8(last(isotope_bounds))),
            0.0f0,  # No transmission threshold for first pass
            0.0f0,  # No fragment tolerance for first pass
            Float32(params.parameter_tuning.iteration_settings.ms1_tol_ppm),  # MS1 tolerance from config
            Float32(params.parameter_tuning.search_settings.frag_err_quantile),
            # Handle min_score as either single value or array (use first value if array)
            begin
                min_score_raw = frag_params.min_score
                if min_score_raw isa Vector
                    UInt8(first(min_score_raw))
                else
                    UInt8(min_score_raw)
                end
            end,
            Int64(frag_params.min_count),
            Float32(frag_params.min_spectral_contrast),
            Float32(frag_params.min_log2_ratio),
            (Int64(first(frag_params.min_top_n)), Int64(last(frag_params.min_top_n))),
            UInt8(1), # max_best_rank
            Int64(frag_params.n_isotopes),
            UInt8(frag_params.max_rank),
            Set{Int64}([2]),
            global_params.match_between_runs,
            Float32(frag_params.relative_improvement_threshold),
            
            Int64(score_params.n_train_rounds),
            Int64(score_params.max_iterations),
            Float32(score_params.max_q_value_probit_rescore),
            Float32(score_params.max_PEP),
            
            Int64(1000), # Default min_inference_points
            Float32(rt_params.min_probability),
            Float32(rt_params.min_probability),
            Float32(0.1), # Default max_irt_bin_size
            Float32(irt_mapping_params.max_prob_to_impute_irt),  # Default max_prob_to_impute
            Float32(irt_mapping_params.fwhm_nstd),   # Default fwhm_nstd
            Float32(irt_mapping_params.irt_nstd),   # Default irt_nstd
            Bool(hasproperty(irt_mapping_params, :plot_rt_alignment) ? irt_mapping_params.plot_rt_alignment : false),
            Bool(hasproperty(irt_mapping_params, :use_robust_fitting) ? irt_mapping_params.use_robust_fitting : true),
            prec_estimation
        )
    end
end


#==========================================================
Interface Implementation
==========================================================#

get_parameters(::FirstPassSearch, params::Any) = FirstPassSearchParameters(params)
getMs1MassErrorModel(ptsr::FirstPassSearchResults) = ptsr.ms1_mass_err_model[]
getMs1TolPpm(params::FirstPassSearchParameters) = params.ms1_tol_ppm

function init_search_results(
    ::FirstPassSearchParameters,
    search_context::SearchContext
)
    temp_folder = joinpath(getDataOutDir(search_context), "temp_data", "first_pass_psms")
    !isdir(temp_folder) && mkdir(temp_folder)
    out_dir = getDataOutDir(search_context)
    qc_dir = joinpath(out_dir, "qc_plots")
    ms1_mass_error_plots = joinpath(qc_dir, "ms1_mass_error_plots")
    !isdir(ms1_mass_error_plots ) && mkdir(ms1_mass_error_plots )
    return FirstPassSearchResults(
        Dictionary{Int64, NamedTuple{(:median_fwhm, :mad_fwhm), Tuple{Float32, Float32}}}(),
        Base.Ref{DataFrame}(),
        Base.Ref{MassErrorModel}(),
        Vector{Float32}(),
        Plots.Plot[],
        qc_dir,
        Set{UInt32}()
    )
end

#==========================================================
Core Processing Methods
==========================================================#

"""
Process a single MS file in the first pass search.
"""
function process_file!(
    results::FirstPassSearchResults,
    params::P, 
    search_context::SearchContext,
    ms_file_idx::Int64,
    spectra::MassSpecData
) where {P<:FirstPassSearchParameters}

    """
    Perform library search with current parameters.
    """
    function perform_library_search(
        spectra::MassSpecData,
        search_context::SearchContext,
        params::FirstPassSearchParameters,
        ms_file_idx::Int64)
        return library_search(spectra, search_context, params, ms_file_idx)
    end

    """
    Process PSMs from library search.
    """
    function process_psms!(
        psms::DataFrame,
        spectra::MassSpecData,
        search_context::SearchContext,
        params::FirstPassSearchParameters,
        ms_file_idx::Int64)

        file_label = try
            getParsedFileName(search_context, ms_file_idx)
        catch
            "file_$(ms_file_idx)"
        end

        function to_isolation_key(value)
            if value === missing
                return typemin(Int64)
            end
            v = Float64(value)
            if isnan(v)
                return typemin(Int64)
            end
            return round(Int64, v * 10_000)
        end

        function compute_psm_stats(
            psms::DataFrame,
            mask::AbstractVector{Bool}
        )
            idxs = findall(mask)
            n_entries = length(idxs)
            if n_entries == 0
                return (spectra = 0, precursor_windows = 0, precursors = 0)
            end
            combos = Set{Tuple{UInt32, Int64, Int64}}()
            precursors = Set{UInt32}()
            scan_idxs = psms[!, :scan_idx]
            precursor_idxs = psms[!, :precursor_idx]
            @inbounds for row_idx in idxs
                prec = precursor_idxs[row_idx]
                push!(precursors, prec)
                scan_idx_val = Int(scan_idxs[row_idx])
                center_key = to_isolation_key(getCenterMz(spectra, scan_idx_val))
                width_key = to_isolation_key(getIsolationWidthMz(spectra, scan_idx_val))
                push!(combos, (prec, center_key, width_key))
            end
            return (
                spectra = n_entries,
                precursor_windows = length(combos),
                precursors = length(precursors)
            )
        end

        """
        Select best PSMs based on criteria.
        """
        function select_best_psms!(
            psms::DataFrame,
            precursor_mzs::AbstractVector,
            params::FirstPassSearchParameters,
            search_context::SearchContext
        )
            fdr_scale_factor = getLibraryFdrScaleFactor(search_context)
            get_best_psms!(
                psms,
                precursor_mzs,
                max_PEP=params.max_PEP,
                fdr_scale_factor=fdr_scale_factor
            )
        end

        rt_model = getRtIrtModel(search_context, ms_file_idx)
        # Add columns
        add_psm_columns!(psms, spectra, search_context, rt_model, ms_file_idx)

        if !isempty(psms)
            target_mask_enter = psms[!, :target]
            decoy_mask_enter = .!target_mask_enter
            target_stats_enter = compute_psm_stats(psms, target_mask_enter)
            decoy_stats_enter = compute_psm_stats(psms, decoy_mask_enter)
            @user_info "FirstPassSearch file $(file_label) entering probit regression: " *
                "targets=$(target_stats_enter.spectra) spectra (" *
                "$(target_stats_enter.precursor_windows) precursor+isolation combos, " *
                "$(target_stats_enter.precursors) unique precursors); " *
                "decoys=$(decoy_stats_enter.spectra) spectra (" *
                "$(decoy_stats_enter.precursor_windows) precursor+isolation combos, " *
                "$(decoy_stats_enter.precursors) unique precursors)"
        end
        
        fallback_to_all_precursors = false
        chrom_candidates = compute_chromatogram_candidates(psms, spectra, search_context, params, ms_file_idx)

        if chrom_candidates === nothing
            log_chrom_precursor_counts(file_label, 0, 0)
            fallback_to_all_precursors = true
        elseif nrow(chrom_candidates) == 0
            log_chrom_precursor_counts(file_label, 0, 0)
            fallback_to_all_precursors = true
        elseif train_chromatogram_probit!(chrom_candidates, params, search_context, file_label)
            update_precursor_candidates_from_chromatograms!(
                results,
                chrom_candidates,
                params,
                file_label,
            )
        else
            log_chrom_precursor_counts(file_label, 0, 0)
            fallback_to_all_precursors = true
        end

        if fallback_to_all_precursors
            for prec in psms[!, :precursor_idx]
                push!(results.probit_precursor_candidates, UInt32(prec))
            end
        end

        # Score PSMs
        score_psms!(psms, params, search_context, spectra)

        # Get best PSMs
        select_best_psms!(
            psms,
            getMz(getPrecursors(getSpecLib(search_context))),#[:mz],
            params,
            search_context
        )


        if !isempty(psms)
            target_mask_pep = psms[!, :target]
            decoy_mask_pep = .!target_mask_pep
            target_stats_pep = compute_psm_stats(psms, target_mask_pep)
            decoy_stats_pep = compute_psm_stats(psms, decoy_mask_pep)
            @user_info "FirstPassSearch file $(file_label) passed PEP threshold (PEP <= $(Float64(params.max_PEP))): " *
                "targets=$(target_stats_pep.spectra) spectra (" *
                "$(target_stats_pep.precursor_windows) precursor+isolation combos, " *
                "$(target_stats_pep.precursors) unique precursors); " *
                "decoys=$(decoy_stats_pep.spectra) spectra (" *
                "$(decoy_stats_pep.precursor_windows) precursor+isolation combos, " *
                "$(decoy_stats_pep.precursors) unique precursors)\n"
        end

        return psms
    end

    """
    Add necessary columns to PSM DataFrame.
    """
    function add_psm_columns!(
        psms::DataFrame,
        spectra::MassSpecData,
        search_context::SearchContext,
        rt_model::RtConversionModel,
        ms_file_idx::Int64)
        add_main_search_columns!(
            psms,
            getModel(rt_model),
            getStructuralMods(getPrecursors(getSpecLib(search_context))),
            getMissedCleavages(getPrecursors(getSpecLib(search_context))),
            getIsDecoy(getPrecursors(getSpecLib(search_context))),
            getIrt(getPrecursors(getSpecLib(search_context))),
            getCharge(getPrecursors(getSpecLib(search_context))),
            getRetentionTimes(spectra),
            getTICs(spectra),
            getMzArrays(spectra)
        )
        # Calculate RT values
        psms[!, :irt_observed] = rt_model.(psms[!, :rt])
        psms[!, :irt_error] = Float16.(abs.(psms[!, :irt_observed] .- psms[!, :irt_predicted]))
        psms[!, :charge2] = UInt8.(psms[!, :charge] .== 2)
        psms[!, :ms_file_idx] .= UInt32(ms_file_idx)
    end

    function compute_chromatogram_candidates(
        psms::DataFrame,
        spectra::MassSpecData,
        search_context::SearchContext,
        _params::FirstPassSearchParameters,
        ms_file_idx::Int64,
    )
        required_cols = [
            :precursor_idx,
            :scan_idx,
            :rt,
            :weight,
            :gof,
            :matched_ratio,
            :fitted_manhattan_distance,
            :fitted_spectral_contrast,
            :scribe,
            :y_count,
            :target,
        ]

        if any(!hasproperty(psms, col) for col in required_cols)
            return nothing
        end

        chrom_source = select(psms, required_cols; copycols = true)
        if isempty(chrom_source)
            return DataFrame(
                precursor_idx = UInt32[],
                isotopes_captured = Tuple{Int8, Int8}[],
                target = Bool[],
                max_weight = Float32[],
                max_gof = Float32[],
                max_matched_ratio = Float32[],
                max_fitted_manhattan_distance = Float32[],
                max_fitted_spectral_contrast = Float32[],
                max_scribe = Float32[],
                y_ions_sum = Float32[],
                max_y_ions = Float32[],
                num_scans = Float32[],
                smoothness = Float32[],
                precursor_fraction_transmitted = Float32[],
                intercept = Float32[],
            )
        end

        chrom_source[!, :scan_idx] = UInt32.(chrom_source[!, :scan_idx])
        sort!(chrom_source, [:precursor_idx, :rt])

        precursors = getPrecursors(getSpecLib(search_context))
        get_isotopes_captured!(
            chrom_source,
            SeperateTraces(),
            getQuadTransmissionModel(search_context, ms_file_idx),
            getSearchData(search_context),
            chrom_source[!, :scan_idx],
            getCharge(precursors),
            getMz(precursors),
            getSulfurCount(precursors),
            getCenterMzs(spectra),
            getIsolationWidthMzs(spectra),
        )

        function safe_max32(iter)
            best = -Inf32
            for v in iter
                fv = Float32(v)
                if isfinite(fv) && fv > best
                    best = fv
                end
            end
            return best == -Inf32 ? 0f0 : best
        end

        function safe_sum32(iter)
            total = 0f0
            for v in iter
                fv = Float32(v)
                if isfinite(fv)
                    total += fv
                end
            end
            return total
        end

        function compute_smoothness(weights::Vector{Float32}, rts::Vector{Float32}, apex_idx::Int)
            apex_weight = weights[apex_idx]
            if apex_weight == 0f0
                return 0f0
            end
            n = length(weights)
            if n == 1
                return ((-2f0 * weights[1]) / apex_weight)^2
            end

            smoothness = 0f0
            for i in 1:n
                if i == 1
                    Δrt = rts[i + 1] - rts[i]
                    if Δrt != 0f0
                        term = ((weights[i + 1] - weights[i]) / Δrt + (-weights[i]) / Δrt) / apex_weight
                        smoothness += term^2
                    end
                elseif i == n
                    Δrt = rts[i] - rts[i - 1]
                    if Δrt != 0f0
                        term = ((weights[i - 1] - weights[i]) / Δrt + (-weights[i]) / Δrt) / apex_weight
                        smoothness += term^2
                    end
                else
                    back_rt = rts[i] - rts[i - 1]
                    fwd_rt = rts[i + 1] - rts[i]
                    term = 0f0
                    if back_rt != 0f0
                        term += (weights[i - 1] - weights[i]) / back_rt
                    end
                    if fwd_rt != 0f0
                        term += (weights[i + 1] - weights[i]) / fwd_rt
                    end
                    if term != 0f0
                        smoothness += (term / apex_weight)^2
                    end
                end
            end
            return smoothness
        end

        rows = NamedTuple[]
        grouped = groupby(chrom_source, [:precursor_idx, :isotopes_captured])
        for sub in grouped
            n = nrow(sub)
            if n == 0
                continue
            end

            weights = Float32.(sub[!, :weight])
            if all(iszero, weights)
                continue
            end

            rts = Float32.(sub[!, :rt])
            _, apex_idx = findmax(weights)
            window_start = max(1, apex_idx - 2)
            window_stop = min(n, apex_idx + 2)
            window = window_start:window_stop

            max_weight = safe_max32(weights)
            max_gof = safe_max32(sub[i, :gof] for i in window)
            max_matched_ratio = safe_max32(sub[i, :matched_ratio] for i in window)
            max_manhattan = safe_max32(sub[i, :fitted_manhattan_distance] for i in window)
            max_spec = safe_max32(sub[i, :fitted_spectral_contrast] for i in window)
            max_scribe = safe_max32(sub[i, :scribe] for i in window)
            y_sum = safe_sum32(sub[i, :y_count] for i in window)
            max_y = safe_max32(sub[i, :y_count] for i in window)
            smoothness = compute_smoothness(weights, rts, apex_idx)
            pct_trans = Float32(sub[apex_idx, :precursor_fraction_transmitted])
            if !isfinite(pct_trans)
                pct_trans = 0f0
            end

            target_val = Bool(sub[apex_idx, :target])
            prec_idx_val = UInt32(sub[apex_idx, :precursor_idx])
            iso_val = sub[apex_idx, :isotopes_captured]

            push!(
                rows,
                (
                    precursor_idx = prec_idx_val,
                    isotopes_captured = iso_val,
                    target = target_val,
                    max_weight = max_weight,
                    max_gof = max_gof,
                    max_matched_ratio = max_matched_ratio,
                    max_fitted_manhattan_distance = max_manhattan,
                    max_fitted_spectral_contrast = max_spec,
                    max_scribe = max_scribe,
                    y_ions_sum = y_sum,
                    max_y_ions = max_y,
                    num_scans = Float32(n),
                    smoothness = smoothness,
                    precursor_fraction_transmitted = pct_trans,
                    intercept = 1.0f0,
                ),
            )
        end

        return DataFrame(rows)
    end

    function train_chromatogram_probit!(
        chrom_df::DataFrame,
        params::FirstPassSearchParameters,
        search_context::SearchContext,
        file_label::String,
    )
        n = nrow(chrom_df)
        if n == 0
            return false
        end

        n_targets = count(t -> t, chrom_df[!, :target])
        n_decoys = n - n_targets
        if n_targets == 0 || n_decoys == 0
            @user_warn "FirstPassSearch file $(file_label) has insufficient chromatogram targets/decoys for probit training"
            return false
        end

        feature_cols = [
            :max_weight,
            :max_gof,
            :max_matched_ratio,
            :max_fitted_manhattan_distance,
            :max_fitted_spectral_contrast,
            :max_scribe,
            :y_ions_sum,
            :max_y_ions,
            :num_scans,
            :smoothness,
            :precursor_fraction_transmitted,
            :intercept,
        ]

        chrom_df[!, :score] = zeros(Float32, n)
        chrom_df[!, :prob] = zeros(Float32, n)
        chrom_df[!, :PEP] = zeros(Float32, n)

        tasks_per_thread = 10
        chunk_size = max(1, n ÷ (tasks_per_thread * Threads.nthreads()))
        data_chunks = partition(1:n, chunk_size)
        feature_df = select(chrom_df, feature_cols; copycols = false)

        try
            β = zeros(Float64, length(feature_cols))
            β = ProbitRegression(
                β,
                feature_df,
                chrom_df[!, :target],
                data_chunks,
                max_iter = params.max_iter_probit,
            )

            ModelPredict!(chrom_df[!, :score], feature_df, β, data_chunks)
            chrom_df[!, :score] = Float32.(chrom_df[!, :score])
            get_probs!(chrom_df, chrom_df[!, :score])

            fdr_scale_factor = getLibraryFdrScaleFactor(search_context)
            get_PEP!(
                chrom_df[!, :score],
                chrom_df[!, :target],
                chrom_df[!, :PEP];
                doSort = false,
                fdr_scale_factor = fdr_scale_factor,
            )
            return true
        catch e
            @user_warn "FirstPassSearch chromatogram probit failed for file $(file_label): $(e)"
            return false
        end
    end

    log_chrom_precursor_counts(file_label::String, target_count::Int, decoy_count::Int) =
        @user_info "FirstPassSearch file $(file_label) chromatogram probit passed $(target_count) target precursors and $(decoy_count) decoy precursors"

    function update_precursor_candidates_from_chromatograms!(
        results::FirstPassSearchResults,
        chrom_df::DataFrame,
        params::FirstPassSearchParameters,
        file_label::String,
    )
        pass_mask = chrom_df[!, :PEP] .<= Float32(params.max_PEP)
        target_precursors = Set{UInt32}()
        decoy_precursors = Set{UInt32}()

        for (idx, pass) in enumerate(pass_mask)
            if !pass
                continue
            end
            prec = UInt32(chrom_df[idx, :precursor_idx])
            if chrom_df[idx, :target]
                push!(target_precursors, prec)
            else
                push!(decoy_precursors, prec)
            end
        end

        log_chrom_precursor_counts(
            file_label,
            length(target_precursors),
            length(decoy_precursors),
        )

        for prec in target_precursors
            push!(results.probit_precursor_candidates, prec)
        end
        for prec in decoy_precursors
            push!(results.probit_precursor_candidates, prec)
        end

        return !isempty(target_precursors) || !isempty(decoy_precursors)
    end

    """
    Score PSMs using probit model.
    """
    function score_psms!(
        psms::DataFrame,
        params::FirstPassSearchParameters,
        search_context::SearchContext,
        spectra::MassSpecData)
        column_names = [
            :spectral_contrast, :city_block, :entropy_score, :scribe, :percent_theoretical_ignored,
            :charge2, :poisson, :irt_error, 
            :missed_cleavage, 
            :Mox,
            #:charge, Only works with charge 2 if at least 3 charge states presence. otherwise singular error
            #:b_count, might be good for non-tryptic enzymes
            :TIC, :y_count, :err_norm, :spectrum_peak_count, :intercept
        ]

        # Avoid singular error if no peaks were ignored
        if maximum(psms.percent_theoretical_ignored) == 0
            deleteat!(column_names, findfirst(==(:percent_theoretical_ignored), column_names))
        end


        # Select scoring columns
        select!(psms, vcat(column_names, [:ms_file_idx, :score, :precursor_idx, :scan_idx,
            :q_value, :log2_summed_intensity, :irt, :rt, :irt_predicted, :target]))
        sort!(psms, [:rt, :precursor_idx])
        # Score PSMs
        fdr_scale_factor = getLibraryFdrScaleFactor(search_context)
        try
            score_main_search_psms!(
                psms,
                column_names,
                n_train_rounds=params.n_train_rounds_probit,
                max_iter_per_round=params.max_iter_probit,
                max_q_value=Float64(params.max_q_value_probit_rescore),
                fdr_scale_factor=fdr_scale_factor
            )
        catch
            column_names = [
            :spectral_contrast, :city_block, :entropy_score, :scribe,
            :charge2, :poisson, :irt_error, :TIC, :y_count, :err_norm, :spectrum_peak_count, :intercept
            ]
            score_main_search_psms!(
                psms,
                column_names,
                n_train_rounds=params.n_train_rounds_probit,
                max_iter_per_round=params.max_iter_probit,
                max_q_value=Float64(params.max_q_value_probit_rescore),
                fdr_scale_factor=fdr_scale_factor
            )
        end
        # Process scores
       
        select!(psms, [:ms_file_idx, :score, :precursor_idx, :scan_idx,
            :q_value, :log2_summed_intensity, :irt, :rt, :irt_predicted, :target])
        get_probs!(psms, psms[!,:score])
    end

    try
        # Get models and update fragment lookup table
        psms = perform_library_search(spectra, search_context, params, ms_file_idx)
        results.psms[] = process_psms!(psms, spectra, search_context, params, ms_file_idx)

        temp_psms = results.psms[] 
        temp_psms = temp_psms[temp_psms[!,:q_value].<=0.001,:]

        # Check if we have any PSMs left after filtering
        if nrow(temp_psms) > 0
            most_intense = sortperm(temp_psms[!,:log2_summed_intensity], rev = true)
            ms1_errs = vcat(
                mass_error_search(
                    spectra,
                    temp_psms[most_intense[1:(min(3000, length(most_intense)))],:scan_idx],
                    temp_psms[most_intense[1:(min(3000, length(most_intense)))],:precursor_idx],
                    UInt32(ms_file_idx),
                    getSpecLib(search_context),
                    getSearchData(search_context),
                    MassErrorModel(
                    0.0f0,
                    (getMs1TolPpm(params), getMs1TolPpm(params))  # Use MS1 tolerance from JSON config
                    ),
                    params,
                    getNceModel(search_context, ms_file_idx),
                    MS1CHROM()
                )...
            )
        else
            ms1_errs = Float32[]
        end

        if length(ms1_errs) > 1
            mad_dev = mad(ms1_errs; normalize=true)
            med_errs = median(ms1_errs)
            low_bound, high_bound = med_errs - mad_dev*7, med_errs + mad_dev*7
            filter!(x->(low_bound<x)&(high_bound>x), ms1_errs)

            # Check again after filtering if we still have data
            if length(ms1_errs) > 0
                ms1_mass_err_model, ms1_ppm_errs = mass_err_ms1(ms1_errs, params)
                results.ms1_mass_err_model[] = ms1_mass_err_model
                append!(results.ms1_ppm_errs, ms1_ppm_errs)
            else
                # Filtering removed all data points, use defaults
                results.ms1_mass_err_model[] = getMassErrorModel(search_context, ms_file_idx)
                append!(results.ms1_ppm_errs, Float32[])
            end
            select!(results.psms[] , Not(:log2_summed_intensity))
        else
            #ms1_mass_err_model, ms1_ppm_errs = mass_err_ms1(ms1_errs, params)
            #Default to MS2 pattern
            results.ms1_mass_err_model[] = getMassErrorModel(search_context, ms_file_idx)
            append!(results.ms1_ppm_errs, Float32[])
            select!(results.psms[] , Not(:log2_summed_intensity))
        end
    catch e
        # Get file name for debugging
        file_name = try
            getFileIdToName(getMSData(search_context), ms_file_idx)
        catch
            "file_$ms_file_idx"
        end

        reason = "FirstPassSearch failed: $e"
        markFileFailed!(search_context, ms_file_idx, reason)
        # Also explicitly set the failed indicator for downstream tracking
        setFailedIndicator!(getMSData(search_context), ms_file_idx, true)
        @user_warn "First pass search failed for MS data file: $file_name. Error: $e. Creating empty results to continue pipeline."

        # Log full stack trace for debugging
        try
            bt = catch_backtrace()
            @user_error "Full error details:\n" * sprint(showerror, e, bt)
        catch
        end
        
        # Create an empty but properly structured DataFrame to avoid downstream errors
        empty_psms = DataFrame(
            ms_file_idx = UInt32[],
            scan_idx = UInt32[], 
            precursor_idx = UInt32[],
            rt = Float32[],
            irt_predicted = Float32[],
            q_value = Float32[],
            score = Float32[], 
            prob = Float32[],
            scan_count = UInt32[],
            fwhm = Float32[]  # Add this to prevent missing column error
        )
        results.psms[] = empty_psms
        
        # Set default mass error model
        results.ms1_mass_err_model[] = getMassErrorModel(search_context, ms_file_idx)
        #rethrow(e)
    end

    return results
end

"""
Initial file processing complete, no additional processing needed.
"""
function process_search_results!(
    results::FirstPassSearchResults,
    params::P,
    search_context::SearchContext,
    ms_file_idx::Int64,
    ::MassSpecData
) where {P<:FirstPassSearchParameters}
    psms = results.psms[]
    fwhms = skipmissing(psms[!, :fwhm])
    fwhm_points = count(!ismissing, fwhms)
    if fwhm_points >= 1#params.min_inference_points
        insert!(results.fwhms, ms_file_idx, (
            median_fwhm = median(fwhms),
            mad_fwhm = mad(fwhms, normalize=true)))
    else
        @user_warn "Insuficient fwhm_points to estimate for $ms_file_idx"
        insert!(results.fwhms, ms_file_idx, (
            median_fwhm = 0.2f0,
            mad_fwhm = 0.2f0))
    end
    parsed_fname = getParsedFileName(search_context, ms_file_idx)
    temp_path = joinpath(getDataOutDir(search_context), "temp_data", "first_pass_psms", parsed_fname * ".arrow")
    psms[!, :ms_file_idx] .= UInt32(ms_file_idx)
    Arrow.write(
        temp_path,
        select!(psms, [:ms_file_idx, :scan_idx, :precursor_idx, :rt,
            :irt_predicted, :q_value, :score, :prob, :scan_count])
    )
    setFirstPassPsms!(getMSData(search_context), ms_file_idx, temp_path)

    #####
    #MS1 mass tolerance 
    ms1_mass_error_folder = getMs1MassErrPlotFolder(search_context)
    parsed_fname = getParsedFileName(search_context, ms_file_idx)
    # Generate mass error plot
    push!(results.ms1_mass_plots, generate_ms1_mass_error_plot(results, parsed_fname))
    # Update models in search context
    setMs1MassErrorModel!(search_context, ms_file_idx, getMs1MassErrorModel(results))
end

"""
No cleanup needed between files.
"""
function reset_results!(results::FirstPassSearchResults)
    empty!(results.psms[])
    resize!(results.ms1_ppm_errs, 0)
    return nothing
end

"""
Summarize results across all files.
"""
function summarize_results!(
    results::FirstPassSearchResults,
    params::P,
    search_context::SearchContext
) where {P<:FirstPassSearchParameters}
    
    """
    Process precursors and calculate iRT errors.
    """
    function get_best_precursors_accross_runs!(
        search_context::SearchContext,
        results::FirstPassSearchResults,
        params::FirstPassSearchParameters
    )

        # Filter out failed files
        valid_indices = get_valid_file_indices(search_context)
        all_psms_paths = getFirstPassPsms(getMSData(search_context))
        valid_psms_paths = [all_psms_paths[i] for i in valid_indices]
        
        # Create RT-IRT map for valid files only
        all_rt_irt = getRtIrtModel(search_context)
        valid_rt_irt = Dict{Int64, RtConversionModel}(i => all_rt_irt[i] for i in valid_indices if haskey(all_rt_irt, i))
        
        if isempty(valid_psms_paths)
            @user_warn "No valid files for cross-run precursor analysis"
            return Dictionary{UInt32, @NamedTuple{best_prob::Float32, best_ms_file_idx::UInt32, best_scan_idx::UInt32, best_irt::Float32, mean_irt::Union{Missing, Float32}, var_irt::Union{Missing, Float32}, n::Union{Missing, UInt16}, mz::Float32}}()
        end
        
        # Get best precursors from valid files only
        return get_best_precursors_accross_runs(
            valid_psms_paths,
            getMz(getPrecursors(getSpecLib(search_context))),#[:mz],
            valid_rt_irt,
            max_q_val=params.max_q_val_for_irt
        )
    end
    # Map retention times and update iRT observations
    map_retention_times!(search_context, results, params)
    correct_first_pass_irt_errors!(search_context)
    # Process precursors
    precursor_dict = get_best_precursors_accross_runs!(search_context, results, params)

    # Ensure all charge states of identified peptides are available for downstream searches
    precursors = getPrecursors(getSpecLib(search_context))
    sequences = getSequence(precursors)
    structural_mods = getStructuralMods(precursors)
    isotopic_mods = getIsotopicMods(precursors)
    is_decoy = getIsDecoy(precursors)
    mz_vals = getMz(precursors)
    peptide_to_precursors = Dict{Tuple{Bool, Any, Any, Any}, Vector{UInt32}}()

    for prec_idx in 1:length(mz_vals)
        key = (
            is_decoy[prec_idx],
            sequences[prec_idx],
            structural_mods[prec_idx],
            isotopic_mods[prec_idx],
        )
        push!(get!(peptide_to_precursors, key, UInt32[]), UInt32(prec_idx))
    end

    passing_peptides = Set{Tuple{Bool, Any, Any, Any}}()
    for prec_idx in keys(precursor_dict)
        idx = Int(prec_idx)
        push!(passing_peptides, (
            is_decoy[idx],
            sequences[idx],
            structural_mods[idx],
            isotopic_mods[idx],
        ))
    end

    probit_candidates = results.probit_precursor_candidates

    for key in passing_peptides
        precursors_for_peptide = get(peptide_to_precursors, key, UInt32[])

        best_existing_entry = nothing
        best_existing_prob = typemin(Float32)
        for prec_idx in precursors_for_peptide
            if haskey(precursor_dict, prec_idx)
                entry = precursor_dict[prec_idx]
                if best_existing_entry === nothing || entry.best_prob > best_existing_prob
                    best_existing_entry = entry
                    best_existing_prob = entry.best_prob
                end
            end
        end

        if best_existing_entry !== nothing
            best_irt_seed = best_existing_entry.best_irt
            mean_irt_seed = best_existing_entry.mean_irt
            var_irt_seed = best_existing_entry.var_irt

            for prec_idx in precursors_for_peptide
                if !haskey(precursor_dict, prec_idx) && in(prec_idx, probit_candidates)
                    idx = Int(prec_idx)
                    mz_value = mz_vals[idx]

                    insert!(
                        precursor_dict,
                        prec_idx,
                        (
                            best_prob = 0f0,
                            best_ms_file_idx = zero(UInt32),
                            best_scan_idx = zero(UInt32),
                            best_irt = best_irt_seed,
                            mean_irt = mean_irt_seed,
                            var_irt = var_irt_seed,
                            n = zero(UInt16),
                            mz = Float32(coalesce(mz_value, 0f0)),
                        ),
                    )
                end
            end

            for prec_idx in precursors_for_peptide
                if haskey(precursor_dict, prec_idx)
                    entry = precursor_dict[prec_idx]
                    precursor_dict[prec_idx] = (
                        best_prob = entry.best_prob,
                        best_ms_file_idx = entry.best_ms_file_idx,
                        best_scan_idx = entry.best_scan_idx,
                        best_irt = best_irt_seed,
                        mean_irt = mean_irt_seed,
                        var_irt = var_irt_seed,
                        n = entry.n,
                        mz = entry.mz,
                    )
                end
            end
        end
    end

    setPrecursorDict!(search_context, precursor_dict)

    precursors = getPrecursors(getSpecLib(search_context))
    is_decoy = getIsDecoy(precursors)
    target_precursor_count = 0
    decoy_precursor_count = 0
    for pid in keys(precursor_dict)
        if is_decoy[pid]
            decoy_precursor_count += 1
        else
            target_precursor_count += 1
        end
    end
    @user_info "FirstPassSearch shared library for second search contains $(target_precursor_count) target precursors and $(decoy_precursor_count) decoy precursors"

    # Calculate RT indices
    create_rt_indices!(search_context, results, precursor_dict, params)
    
    # Merge mass error plots
    ms1_mass_error_folder = getMs1MassErrPlotFolder(search_context)
    output_path = joinpath(ms1_mass_error_folder, "ms1_mass_error_plots.pdf")
    try
        if isfile(output_path)
            rm(output_path)
        end
    catch e
        @user_warn "Could not clear existing file: $e"
    end

    if !isempty(results.ms1_mass_plots)
        save_multipage_pdf(results.ms1_mass_plots, output_path)
        empty!(results.ms1_mass_plots)
    end
end

