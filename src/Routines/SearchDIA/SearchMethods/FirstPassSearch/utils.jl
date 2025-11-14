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
    add_main_search_columns!(psms::DataFrame, rt_irt::UniformSpline,
                           structural_mods::AbstractVector{Union{Missing, String}},
                           prec_missed_cleavages::Arrow.Primitive{UInt8, Vector{UInt8}},
                           prec_is_decoy::Arrow.BoolVector{Bool},
                           prec_irt::Arrow.Primitive{T, Vector{T}},
                           prec_charge::Arrow.Primitive{UInt8, Vector{UInt8}},
                           scan_retention_time::AbstractVector{Float32},
                           tic::AbstractVector{Float32},
                           masses::AbstractArray) where {T<:AbstractFloat}

Adds essential columns to PSM DataFrame for scoring and analysis.

# Arguments
- `psms`: DataFrame to modify
- `rt_irt`: Spline for RT to iRT conversion
- `structural_mods`: Vector of structural modifications
- `prec_missed_cleavages`: Vector of missed cleavage counts
- `prec_is_decoy`: Vector indicating decoy status
- `prec_irt`: Vector of iRT values
- `prec_charge`: Vector of precursor charges
- `scan_retention_time`: Vector of scan retention times
- `tic`: Vector of total ion currents
- `masses`: Array of mass spectra

# Added Columns
- Basic metrics: RT, iRT, charge, TIC
- Modification info: missed_cleavage, Mox
- Scoring columns: score, q_value
- Analysis columns: target, spectrum_peak_count, err_norm
"""
function add_main_search_columns!(psms::DataFrame,
                                rt_irt::T,
                                structural_mods::AbstractVector{Union{Missing, String}},
                                prec_missed_cleavages::Arrow.Primitive{UInt8, Vector{UInt8}},
                                prec_is_decoy::Arrow.BoolVector{Bool},
                                prec_irt::Arrow.Primitive{U, Vector{U}},
                                prec_charge::Arrow.Primitive{UInt8, Vector{UInt8}},
    scan_retention_time::AbstractVector{Float32},
    tic::AbstractVector{Float32},
    masses::AbstractArray;
) where {T,U<:AbstractFloat}
    
    ###########################
    #Allocate new columns
    N = size(psms, 1)
    missed_cleavage = zeros(UInt8, N);
    Mox = zeros(UInt8, N);
    irt_pred = zeros(Float32, N);
    rt = zeros(Float32, N);
    irt = zeros(Float32, N);
    TIC = zeros(Float16, N);
    charge = zeros(UInt8, N);
    err_norm = zeros(Float16, N);
    targets = zeros(Bool, N);
    spectrum_peak_count = zeros(UInt32, N);
    scan_idx::Vector{UInt32} = psms[!,:scan_idx]
    precursor_idx::Vector{UInt32} = psms[!,:precursor_idx] 
    error::Vector{Float32} = psms[!,:error]
    total_ions::Vector{UInt8} = psms[!,:y_count] .+ psms[!,:b_count]
    function countMOX(seq::String)
        return UInt8(count("Unimod:35", seq))
    end

    #Split data into chunk ranges
    tasks_per_thread = 10
    chunk_size = max(1, size(psms, 1) ÷ (tasks_per_thread * Threads.nthreads()))
    data_chunks = partition(1:size(psms, 1), chunk_size) # partition your data into chunks that

    tasks = map(data_chunks) do chunk
        Threads.@spawn begin
            for i in chunk
                prec_idx = precursor_idx[i]

                targets[i] = prec_is_decoy[prec_idx] == false;
                missed_cleavage[i] = prec_missed_cleavages[prec_idx]
                Mox[i] = countMOX(coalesce(structural_mods[prec_idx], ""))::UInt8 #UInt8(length(collect(eachmatch(r"ox",  precursors[precursor_idx[i]].sequence))))
                irt_pred[i] = Float32(prec_irt[prec_idx]);
                rt[i] = Float32(scan_retention_time[scan_idx[i]]);
                irt[i] = rt_irt(rt[i])
                TIC[i] = Float16(log2(tic[scan_idx[i]]));
                charge[i] = UInt8(prec_charge[prec_idx]);
                err_norm[i] = min(Float16.(abs(error[i]/total_ions[i])), Float16(6e4))#Float16(min(abs((error[i])/(total_ions[i])), 6e4))
                spectrum_peak_count[i] = UInt32(length(masses[scan_idx[i]]))
            end
        end
    end
    fetch.(tasks)
    psms[!,:irt_predicted] = irt_pred
    psms[!,:rt] = rt
    psms[!,:irt] = irt
    psms[!,:TIC] = TIC
    psms[!,:total_ions] = total_ions
    psms[!,:err_norm] = err_norm
    psms[!,:target] = targets
    psms[!,:missed_cleavage] = missed_cleavage
    psms[!,:Mox] = Mox
    psms[!,:charge] = charge
    psms[!,:spectrum_peak_count] = spectrum_peak_count
    psms[!,:score] = zeros(Float32, N);
    psms[!,:q_value] = zeros(Float16, N);
    psms[!,:intercept] = ones(Float16, N)
end

const FIRST_PASS_PROBIT_FEATURE_COLUMNS = (
    :spectral_contrast,
    :city_block,
    :entropy_score,
    :scribe,
    :percent_theoretical_ignored,
    :charge2,
    :poisson,
    :irt_error,
    :missed_cleavage,
    :Mox,
    :TIC,
    :y_count,
    :err_norm,
    :spectrum_peak_count,
)

function prepare_chromatogram_feature_source(psms::DataFrame)
    if nrow(psms) == 0
        return DataFrame()
    end

    chrom_df = DataFrame()
    mandatory_cols = (:precursor_idx, :scan_idx, :rt, :target)
    optional_cols = (
        FIRST_PASS_PROBIT_FEATURE_COLUMNS...,
    )

    missing_mandatory = Symbol[]
    for col in mandatory_cols
        if hasproperty(psms, col)
            chrom_df[!, col] = copy(psms[!, col])
        else
            push!(missing_mandatory, col)
            default_value = if col === :target
                false
            elseif col === :rt
                0.0f0
            else
                UInt32(0)
            end
            chrom_df[!, col] = fill(default_value, nrow(psms))
        end
    end

    missing_optional = Symbol[]
    for col in optional_cols
        if hasproperty(psms, col)
            chrom_df[!, col] = copy(psms[!, col])
        else
            push!(missing_optional, col)
        end
    end

    if !isempty(missing_mandatory)
        @user_warn "Chromatogram feature source missing mandatory columns: $(join(string.(missing_mandatory), ", "))"
    end
    if !isempty(missing_optional)
        @user_info "Chromatogram feature source missing optional columns: $(join(string.(missing_optional), ", "))"
    end

    return chrom_df
end

function assign_pair_indices!(chrom_df::DataFrame, search_context::SearchContext)
    if nrow(chrom_df) == 0 || hasproperty(chrom_df, :pair_idx)
        return chrom_df
    end
    pair_ids = getPairIdx(getPrecursors(getSpecLib(search_context)))
    chrom_df[!, :pair_idx] = Vector{UInt32}(undef, nrow(chrom_df))
    @inbounds for (idx, prec_idx) in enumerate(chrom_df[!, :precursor_idx])
        chrom_df[idx, :pair_idx] = extract_pair_idx(pair_ids, prec_idx)
    end
    return chrom_df
end

function annotate_isotopes!(chrom_df::DataFrame,
                            search_context::SearchContext,
                            spectra::MassSpecData,
                            ms_file_idx::Int)
    if nrow(chrom_df) == 0
        return chrom_df
    end
    quad_model = getQuadTransmissionModel(search_context, ms_file_idx)
    precursors = getPrecursors(getSpecLib(search_context))
    get_isotopes_captured!(
        chrom_df,
        SeperateTraces(),
        quad_model,
        getSearchData(search_context),
        chrom_df[!, :scan_idx],
        getCharge(precursors),
        getMz(precursors),
        getSulfurCount(precursors),
        getCenterMzs(spectra),
        getIsolationWidthMzs(spectra)
    )
    return chrom_df
end

function summarize_chromatograms(chrom_df::DataFrame)
    if nrow(chrom_df) == 0 || !hasproperty(chrom_df, :isotopes_captured)
        return DataFrame()
    end

    grouped = groupby(chrom_df, [:precursor_idx, :pair_idx, :isotopes_captured])
    n_groups = length(grouped)
    if n_groups == 0
        return DataFrame()
    end

    precursor_type = eltype(chrom_df[!, :precursor_idx])
    pair_type = hasproperty(chrom_df, :pair_idx) ? eltype(chrom_df[!, :pair_idx]) : UInt32
    isotope_type = eltype(chrom_df[!, :isotopes_captured])

    summary_cols = Dict{Symbol, AbstractVector}()
    summary_cols[:precursor_idx] = Vector{precursor_type}(undef, n_groups)
    summary_cols[:pair_idx] = Vector{pair_type}(undef, n_groups)
    summary_cols[:isotopes_captured] = Vector{isotope_type}(undef, n_groups)
    summary_cols[:target] = Vector{Bool}(undef, n_groups)
    summary_cols[:num_scans] = Vector{Float32}(undef, n_groups)

    has_spectral_contrast = hasproperty(chrom_df, :spectral_contrast)
    has_city_block = hasproperty(chrom_df, :city_block)
    has_entropy_score = hasproperty(chrom_df, :entropy_score)
    has_scribe = hasproperty(chrom_df, :scribe)
    has_percent_theoretical_ignored = hasproperty(chrom_df, :percent_theoretical_ignored)
    has_charge2 = hasproperty(chrom_df, :charge2)
    has_poisson = hasproperty(chrom_df, :poisson)
    has_irt_error = hasproperty(chrom_df, :irt_error)
    has_missed_cleavage = hasproperty(chrom_df, :missed_cleavage)
    has_Mox = hasproperty(chrom_df, :Mox)
    has_TIC = hasproperty(chrom_df, :TIC)
    has_y_count = hasproperty(chrom_df, :y_count)
    has_err_norm = hasproperty(chrom_df, :err_norm)
    has_spectrum_peaks = hasproperty(chrom_df, :spectrum_peak_count)

    has_spectral_contrast && begin
        summary_cols[:spectral_contrast] = Vector{Float32}(undef, n_groups)
        summary_cols[:max_spectral_contrast] = Vector{Float32}(undef, n_groups)
    end
    has_city_block && begin
        summary_cols[:city_block] = Vector{Float32}(undef, n_groups)
        summary_cols[:max_city_block] = Vector{Float32}(undef, n_groups)
    end
    has_entropy_score && begin
        summary_cols[:entropy_score] = Vector{Float32}(undef, n_groups)
        summary_cols[:max_entropy_score] = Vector{Float32}(undef, n_groups)
    end
    has_scribe && begin
        summary_cols[:scribe] = Vector{Float32}(undef, n_groups)
        summary_cols[:max_scribe] = Vector{Float32}(undef, n_groups)
    end
    has_percent_theoretical_ignored && begin
        summary_cols[:percent_theoretical_ignored] = Vector{Float32}(undef, n_groups)
        summary_cols[:max_percent_theoretical_ignored] = Vector{Float32}(undef, n_groups)
    end
    has_charge2 && begin
        summary_cols[:charge2] = Vector{Float32}(undef, n_groups)
        summary_cols[:max_charge2] = Vector{Float32}(undef, n_groups)
    end
    has_poisson && begin
        summary_cols[:poisson] = Vector{Float32}(undef, n_groups)
        summary_cols[:max_poisson] = Vector{Float32}(undef, n_groups)
    end
    has_irt_error && begin
        summary_cols[:irt_error] = Vector{Float32}(undef, n_groups)
        summary_cols[:min_irt_error] = Vector{Float32}(undef, n_groups)
    end
    has_missed_cleavage && begin
        summary_cols[:missed_cleavage] = Vector{Float32}(undef, n_groups)
        summary_cols[:max_missed_cleavage] = Vector{Float32}(undef, n_groups)
    end
    has_Mox && begin
        summary_cols[:Mox] = Vector{Float32}(undef, n_groups)
        summary_cols[:max_Mox] = Vector{Float32}(undef, n_groups)
    end
    has_TIC && begin
        summary_cols[:TIC] = Vector{Float32}(undef, n_groups)
        summary_cols[:max_TIC] = Vector{Float32}(undef, n_groups)
    end
    has_y_count && begin
        summary_cols[:y_count] = Vector{Float32}(undef, n_groups)
        summary_cols[:max_y_count] = Vector{Float32}(undef, n_groups)
        summary_cols[:sum_y_count] = Vector{Float32}(undef, n_groups)
    end
    has_err_norm && begin
        summary_cols[:err_norm] = Vector{Float32}(undef, n_groups)
        summary_cols[:max_err_norm] = Vector{Float32}(undef, n_groups)
    end
    has_spectrum_peaks && begin
        summary_cols[:spectrum_peak_count] = Vector{Float32}(undef, n_groups)
        summary_cols[:max_spectrum_peak_count] = Vector{Float32}(undef, n_groups)
    end

    max_or_zero(values) = begin
        data = collect(skipmissing(values))
        return isempty(data) ? 0.0f0 : Float32(maximum(data))
    end
    sum_or_zero(values) = begin
        data = collect(skipmissing(values))
        return isempty(data) ? 0.0f0 : Float32(sum(Float64.(data)))
    end
    min_or_zero(values) = begin
        data = collect(skipmissing(values))
        return isempty(data) ? 0.0f0 : Float32(minimum(data))
    end

    for (i, group) in enumerate(grouped)
        summary_cols[:precursor_idx][i] = first(group[!, :precursor_idx])
        summary_cols[:pair_idx][i] = first(group[!, :pair_idx])
        summary_cols[:isotopes_captured][i] = first(group[!, :isotopes_captured])
        summary_cols[:target][i] = Bool(first(coalesce.(group[!, :target], false)))
        summary_cols[:num_scans][i] = Float32(nrow(group))
        if has_spectral_contrast
            val = max_or_zero(group[!, :spectral_contrast])
            summary_cols[:spectral_contrast][i] = val
            summary_cols[:max_spectral_contrast][i] = val
        end
        if has_city_block
            val = max_or_zero(group[!, :city_block])
            summary_cols[:city_block][i] = val
            summary_cols[:max_city_block][i] = val
        end
        if has_entropy_score
            val = max_or_zero(group[!, :entropy_score])
            summary_cols[:entropy_score][i] = val
            summary_cols[:max_entropy_score][i] = val
        end
        if has_scribe
            val = max_or_zero(group[!, :scribe])
            summary_cols[:scribe][i] = val
            summary_cols[:max_scribe][i] = val
        end
        if has_percent_theoretical_ignored
            val = max_or_zero(group[!, :percent_theoretical_ignored])
            summary_cols[:percent_theoretical_ignored][i] = val
            summary_cols[:max_percent_theoretical_ignored][i] = val
        end
        if has_charge2
            val = max_or_zero(group[!, :charge2])
            summary_cols[:charge2][i] = val
            summary_cols[:max_charge2][i] = val
        end
        if has_poisson
            val = max_or_zero(group[!, :poisson])
            summary_cols[:poisson][i] = val
            summary_cols[:max_poisson][i] = val
        end
        if has_irt_error
            val = min_or_zero(group[!, :irt_error])
            summary_cols[:irt_error][i] = val
            summary_cols[:min_irt_error][i] = val
        end
        if has_missed_cleavage
            val = max_or_zero(group[!, :missed_cleavage])
            summary_cols[:missed_cleavage][i] = val
            summary_cols[:max_missed_cleavage][i] = val
        end
        if has_Mox
            val = max_or_zero(group[!, :Mox])
            summary_cols[:Mox][i] = val
            summary_cols[:max_Mox][i] = val
        end
        if has_TIC
            val = max_or_zero(group[!, :TIC])
            summary_cols[:TIC][i] = val
            summary_cols[:max_TIC][i] = val
        end
        if has_y_count
            val = max_or_zero(group[!, :y_count])
            summary_cols[:y_count][i] = val
            summary_cols[:max_y_count][i] = val
            summary_cols[:sum_y_count][i] = sum_or_zero(group[!, :y_count])
        end
        if has_err_norm
            val = max_or_zero(group[!, :err_norm])
            summary_cols[:err_norm][i] = val
            summary_cols[:max_err_norm][i] = val
        end
        if has_spectrum_peaks
            val = max_or_zero(group[!, :spectrum_peak_count])
            summary_cols[:spectrum_peak_count][i] = val
            summary_cols[:max_spectrum_peak_count][i] = val
        end
    end

    return DataFrame(summary_cols)
end

function score_chromatogram_features!(chrom_summary::DataFrame,
                                      params::FirstPassSearchParameters,
                                      search_context::SearchContext)
    n = nrow(chrom_summary)
    if n == 0
        return chrom_summary
    end

    chrom_summary[!, :intercept] = ones(Float32, n)
    chrom_summary[!, :score] = zeros(Float32, n)
    chrom_summary[!, :q_value] = ones(Float16, n)
    chrom_summary[!, :PEP] = ones(Float16, n)

    feature_candidates = Symbol[]
    has_num_scans = hasproperty(chrom_summary, :num_scans)
    has_num_scans && push!(feature_candidates, :num_scans)

    aggregate_feature_order = (
        :max_spectral_contrast,
        :max_city_block,
        :max_entropy_score,
        :max_scribe,
        :max_percent_theoretical_ignored,
        :max_charge2,
        :max_poisson,
        :min_irt_error,
        :max_missed_cleavage,
        :max_Mox,
        :max_TIC,
        :max_y_count,
        :sum_y_count,
        :max_err_norm,
        :max_spectrum_peak_count,
    )

    missing_for_scoring = Symbol[]
    for col in aggregate_feature_order
        if hasproperty(chrom_summary, col)
            push!(feature_candidates, col)
        else
            push!(missing_for_scoring, col)
        end
    end

    if !isempty(missing_for_scoring)
        @user_info "Skipping unavailable chromatogram-level features: $(join(string.(missing_for_scoring), ", "))"
    end

    push!(feature_candidates, :intercept)
    feature_columns = feature_candidates

    zero_var_cols = Symbol[]
    for col in feature_columns
        col === :intercept && continue
        min_val, max_val = extrema(chrom_summary[!, col])
        if min_val == max_val
            push!(zero_var_cols, col)
        end
    end
    if !isempty(zero_var_cols)
        dropped = join(string.(zero_var_cols), ", ")
        @user_warn "Dropping $(length(zero_var_cols)) zero-variance chromatogram features: $dropped"
        filter!(col -> !(col in zero_var_cols), feature_columns)
    end
    if length(feature_columns) <= 1
        @user_warn "Skipping chromatogram-level probit scoring because no varying features remain after filtering"
        return chrom_summary
    end

    fdr_scale_factor = getLibraryFdrScaleFactor(search_context)
    score_main_search_psms!(
        chrom_summary,
        feature_columns,
        n_train_rounds = params.n_train_rounds_probit,
        max_iter_per_round = params.max_iter_probit,
        max_q_value = Float64(params.max_q_value_probit_rescore),
        fdr_scale_factor = fdr_scale_factor
    )
    get_probs!(chrom_summary, chrom_summary[!, :score])
    get_PEP!(chrom_summary[!, :score], chrom_summary[!, :target], chrom_summary[!, :PEP]; doSort=false, fdr_scale_factor=fdr_scale_factor)
    return chrom_summary
end

function chromatogram_precursor_candidates!(chrom_source::DataFrame,
                                            scored_psms::DataFrame,
                                            params::FirstPassSearchParameters,
                                            search_context::SearchContext,
                                            spectra::MassSpecData,
                                            ms_file_idx::Int)
    if nrow(chrom_source) == 0 || nrow(scored_psms) == 0
        return (Set{UInt32}(), 0, 0)
    end
    scores_df = select(scored_psms, [:precursor_idx, :scan_idx, :score, :prob])
    chrom_data = leftjoin(chrom_source, scores_df, on=[:precursor_idx, :scan_idx]; makeunique=true)
    assign_pair_indices!(chrom_data, search_context)
    annotate_isotopes!(chrom_data, search_context, spectra, ms_file_idx)
    chrom_summary = summarize_chromatograms(chrom_data)
    chrom_summary = score_chromatogram_features!(chrom_summary, params, search_context)
    if nrow(chrom_summary) == 0
        return (Set{UInt32}(), 0, 0)
    end
    passing_mask = chrom_summary[!, :PEP] .<= Float16(params.max_PEP)
    passing = chrom_summary[passing_mask, :]
    targets = Set{UInt32}(passing[passing[!, :target], :precursor_idx])
    decoys = Set{UInt32}(passing[.!passing[!, :target], :precursor_idx])
    return (Set{UInt32}(passing[!, :precursor_idx]), length(targets), length(decoys))
end

"""
    score_main_search_psms!(psms::DataFrame, column_names::Vector{Symbol};
                           n_train_rounds::Int64=2,
                           max_iter_per_round::Int64=50,
                           max_q_value::Float64=0.01,
                           fdr_scale_factor::Float32=1.0f0)

Performs iterative probit regression scoring of PSMs.

# Arguments
- `psms`: DataFrame containing PSMs to score
- `column_names`: Column names to use for scoring
- `n_train_rounds`: Number of training iterations
- `max_iter_per_round`: Maximum iterations per training round
- `max_q_value`: Maximum q-value threshold for filtering
- `fdr_scale_factor`: Scale factor to correct for library target/decoy ratio

# Process
1. First round: Trains on all data
2. Subsequent rounds: Trains on decoys and high-scoring targets
3. Updates scores and q-values after each round
"""
function score_main_search_psms!(psms::DataFrame, column_names::Vector{Symbol};
                             n_train_rounds::Int64 = 2,
                             max_iter_per_round::Int64 = 50,
                             max_q_value::Float64 = 0.01,
                             fdr_scale_factor::Float32 = 1.0f0
)

    β = zeros(Float64, length(column_names));
    best_psms = nothing#ones(Bool, size(psms, 1))

    tasks_per_thread = 10
    M = size(psms, 1)
    chunk_size = max(1, M ÷ (tasks_per_thread * Threads.nthreads()))
    data_chunks = partition(1:M, chunk_size) # partition your data into chunks

    for i in range(1, n_train_rounds)
        if i < 2 #Train on top scribe data during first round
            get_qvalues!(psms[!,:scribe], psms[!,:target], psms[!,:q_value]; fdr_scale_factor = fdr_scale_factor)
        end

        if i < n_train_rounds #Get Data to train on
            best_psms = ((psms[!,:q_value].<=max_q_value).&(psms[!,:target])) .| (psms[!,:target].==false);
        end

        psms_targets = psms[best_psms,:target]
        M = size(psms_targets, 1)
        sub_chunk_size = max(1, M ÷ (tasks_per_thread * Threads.nthreads()))
        sub_data_chunks = partition(1:M, sub_chunk_size) # partition your data into chunks that
        β = ProbitRegression(β, psms[best_psms,column_names], psms_targets, sub_data_chunks, max_iter = max_iter_per_round);
        ModelPredict!(psms[!,:score], psms[!,column_names], β, data_chunks); #Get Z-scores 
    end

    return 
end

"""
    get_probs!(psms::DataFrame, psms_scores::Vector{Float32})

Calculates probability scores for PSMs based on their Z-scores using error function.

# Arguments
- `psms`: DataFrame to modify
- `psms_scores`: Vector of PSM Z-scores

Adds 'prob' column to psms DataFrame containing probability scores.
Uses parallel processing for efficiency.
"""
function get_probs!(psms::DataFrame,
                    psms_scores::Vector{Float32}
)

    psms_probs = zeros(Float32, size(psms, 1))
    tasks_per_thread = 10
    M = size(psms, 1)
    chunk_size = max(1, M ÷ (tasks_per_thread * Threads.nthreads()))
    data_chunks = partition(1:M, chunk_size) # partition your data into chunks that
    tasks = map(data_chunks) do chunk
        Threads.@spawn begin
            for i in chunk
                psms_probs[i] = (1 + SpecialFunctions.erf(psms_scores[i]/sqrt(2)))/2
            end
        end
    end
    fetch.(tasks)
    psms[!,:prob] = psms_probs
    return nothing
end


"""
    get_best_psms!(psms::DataFrame, prec_mz::Arrow.Primitive{T, Vector{T}}; 
                   max_PEP::Float32=0.9f0, fdr_scale_factor::Float32=1.0f0) where {T<:AbstractFloat}

Processes PSMs to identify the best matches and calculate peak characteristics.

# Arguments
- `psms`: DataFrame containing peptide-spectrum matches
- `prec_mz`: Vector of precursor m/z values
- `max_PEP`: Maximum local FDR threshold for filtering PSMs
- `fdr_scale_factor`: Scale factor to correct for library target/decoy ratio

# Modifies PSMs DataFrame to add:
- `best_psm`: Boolean indicating highest scoring PSM for each precursor
- `fwhm`: Full width at half maximum of chromatographic peak
- `scan_count`: Number of scans below q-value threshold
- `prec_mz`: Precursor m/z values

Note: Assumes psms is sorted by retention time in ascending order.
"""
function get_best_psms!(psms::DataFrame,
                        prec_mz::Arrow.Primitive{T, Vector{T}};
                        max_PEP::Float32 = 0.9f0,
                        fdr_scale_factor::Float32 = 1.0f0
) where {T<:AbstractFloat}

    #highest scoring psm for a given precursor
    psms[!,:PEP] = zeros(Float16, size(psms, 1))
    #highest scoring psm for a given precursor
    psms[!,:best_psm] = zeros(Bool, size(psms, 1))
    #fwhm estimate of the precursor
    psms[!,:fwhm] = zeros(Union{Missing, Float32}, size(psms, 1))
    #number of scans below the q value threshold for hte precursor
    psms[!,:scan_count] = zeros(UInt16,size(psms, 1))

    
    #Get best psm for each precursor 
    #ASSUMES psms IS SOrtED BY rt IN ASCENDING ORDER
    gpsms = groupby(psms,:precursor_idx)
    for (precursor_idx, prec_psms) in pairs(gpsms)

        #Get the best scoring psm, and the its row index. 
        #Get the maximum intensity psm (under the q_value threshold) and its row index. 
        max_irt, min_irt = missing, missing
        max_log2_intensity, max_idx = missing, one(Int64)
        best_psm_score, best_psm_idx = zero(Float32), one(Int64)
        scan_count = one(UInt8)
        for i in range(1, size(prec_psms, 1))
            if coalesce(max_log2_intensity, zero(Float32)) < prec_psms[i,:log2_summed_intensity]
                max_log2_intensity = prec_psms[i,:log2_summed_intensity]
                max_idx = i
            end
            if prec_psms[i,:score]>best_psm_score
                best_psm_idx = i
                best_psm_score = prec_psms[i,:score]
            end
        end
        #Mark the best psm 
        prec_psms[best_psm_idx,:best_psm] = true

        #Try to estimate the fwhm. 
        i = max_idx - 1
        while i > 0
            #Is the i'th psm above half the maximum 
            if (prec_psms[i,:log2_summed_intensity] > (max_log2_intensity - 1.0))
                scan_count += 1
                min_irt = prec_psms[i,:irt]
            else
                break
            end
            i -= 1
        end
        
        i = max_idx + 1
        while i <= size(prec_psms, 1)
            #Is the i'th psm above half the maximum 
            if (prec_psms[i,:log2_summed_intensity] > (max_log2_intensity - 1.0))
                scan_count += 1
                max_irt = prec_psms[i,:irt]
            else
                break
            end
            i += 1
        end

        prec_psms[best_psm_idx,:fwhm] = max_irt - min_irt
        prec_psms[best_psm_idx,:scan_count] = scan_count
    end

    filter!(x->x.best_psm, psms);
    sort!(psms,:score, rev = true)
    # Will use PEP for final filter
    get_PEP!(psms[!,:score], psms[!,:target], psms[!,:PEP]; doSort=false, fdr_scale_factor=fdr_scale_factor);

    n = size(psms, 1)
    select!(psms, [:precursor_idx,:log2_summed_intensity,:rt,:irt_predicted,:q_value,:score,:prob,:fwhm,:scan_count,:scan_idx,:PEP,:target])

    first_fail = searchsortedfirst(psms[!,:PEP], Float16(max_PEP))
    if first_fail <= n
        deleteat!(psms, first_fail:n)
    end
    #println("unique IDs prefilter: ", n, " ", first_fail, "\n\n")

    mz = zeros(T, size(psms, 1));
    precursor_idx = psms[!,:precursor_idx]::Vector{UInt32}
    Threads.@threads for i in range(1, size(psms, 1))
        mz[i] = prec_mz[precursor_idx[i]];
    end
    psms[!,:prec_mz] = mz

    return
end

"""
    map_retention_times!(search_context::SearchContext, results::FirstPassSearchResults, params::FirstPassSearchParameters)

Maps retention times between library and empirical scales for each MS file. For files with sufficient high-confidence PSMs 
(probability > min_prob_for_irt_mapping), creates spline-based conversion models between retention time (RT) and indexed 
retention time (iRT) scales.

# Arguments
- `search_context`: Contains MS data and mapping information
- `results`: FirstPassSearch results container
- `params`: Search parameters including minimum probability threshold

Creates and stores:
- RT to iRT mapping splines
- iRT to RT mapping splines

Throws an error if insufficient high-confidence PSMs are found for mapping.
"""
function map_retention_times!(
    search_context::SearchContext,
    results::FirstPassSearchResults,
    params::FirstPassSearchParameters
)

    # Only process valid (non-failed) files
    valid_files = get_valid_file_indices(search_context)
    all_psms_paths = getFirstPassPsms(getMSData(search_context))

    for ms_file_idx in valid_files
        # Skip files that have been marked as failed
        if is_file_failed(search_context, ms_file_idx)
            continue
        end
        psms_path = all_psms_paths[ms_file_idx]
        psms = Arrow.Table(psms_path)
        best_hits = psms[:prob].>params.min_prob_for_irt_mapping#Map rts using only the best psms
        try#if sum(best_hits) > 100
            if params.use_robust_fitting
                # Use robust fitting with RANSAC/penalty
                best_psms_df = DataFrame(
                    rt = psms[:rt][best_hits],
                    irt_predicted = psms[:irt_predicted][best_hits]
                )

                # Fit RT to iRT model using shared robust fitting
                rt_model, valid_rt, valid_irt, irt_mad = Pioneer.fit_irt_model(
                    best_psms_df;
                    lambda_penalty = Float32(0.1),  # Default penalty
                    ransac_threshold = 1000,        # Default RANSAC threshold
                    min_psms = 10,                  # Default minimum PSMs
                    spline_degree = 3,              # Default spline degree
                    max_knots = 7,                  # Default max knots
                    outlier_threshold = Float32(5.0)  # Default outlier threshold
                )

                # Store RT to iRT model
                setRtIrtMap!(search_context, rt_model, ms_file_idx)

                # Fit inverse model (iRT to RT) - swap input/output
                irt_to_rt_df = DataFrame(
                    rt = valid_irt,           # Use iRT as input
                    irt_predicted = valid_rt  # Use RT as output
                )
                irt_model, _, _, _ = Pioneer.fit_irt_model(
                    irt_to_rt_df;
                    lambda_penalty = Float32(0.1),
                    ransac_threshold = 1000,
                    min_psms = 10,
                    spline_degree = 3,
                    max_knots = 7,
                    outlier_threshold = Float32(5.0)
                )
                setIrtRtMap!(search_context, irt_model, ms_file_idx)

                # Optionally generate plots
                if params.plot_rt_alignment
                    plot_rt_alignment_firstpass(
                        valid_rt,
                        valid_irt,
                        rt_model,
                        ms_file_idx,
                        getDataOutDir(search_context)
                    )
                end
            else
                # Use simple UniformSpline (legacy behavior)
                best_rts = psms[:rt][best_hits]
                best_irts = psms[:irt_predicted][best_hits]

                irt_to_rt_spline = UniformSpline(
                    best_rts,
                    best_irts,
                    3,
                    5
                )
                rt_to_irt_spline = UniformSpline(
                    best_irts,
                    best_rts,
                    3,
                    5
                )

                #Build rt=>irt and irt=> rt mappings for the file and add to the dictionaries
                setRtIrtMap!(search_context, SplineRtConversionModel(rt_to_irt_spline), ms_file_idx)
                setIrtRtMap!(search_context, SplineRtConversionModel(irt_to_rt_spline), ms_file_idx)
            end
        catch e
            # Get file name for debugging
            file_name = try
                getFileIdToName(getMSData(search_context), ms_file_idx)
            catch
                "file_$ms_file_idx"
            end
            
            # Safely compute PSM count to avoid excessive output
            n_good_psms = try
                sum(best_hits)
            catch
                0  # Default to 0 if calculation fails
            end
            
            @user_warn "RT mapping failed for MS data file: $file_name ($n_good_psms good PSMs found, need >100 for spline). Using identity RT model."
            
            # Use identity mapping as fallback - no RT to iRT conversion
            identity_model = IdentityModel()
            setRtIrtMap!(search_context, identity_model, ms_file_idx)
            setIrtRtMap!(search_context, identity_model, ms_file_idx)
        end
    end
    
    # Set identity models for failed files
    ms_data = getMassSpecData(search_context)
    for failed_idx in 1:length(ms_data.file_paths)
        if getFailedIndicator(ms_data, failed_idx)
            file_name = try
                getFileIdToName(getMSData(search_context), failed_idx)
            catch
                "file_$failed_idx"
            end
            @user_warn "Setting identity RT models for failed file: $file_name"
            setRtIrtMap!(search_context, IdentityModel(), failed_idx)
            setIrtRtMap!(search_context, IdentityModel(), failed_idx)
        end
    end
    
    return nothing
end

"""
    plot_rt_alignment_firstpass(rt, irt, rt_model, ms_file_idx, output_dir)

Generate RT alignment diagnostic plot for FirstPassSearch.

# Arguments
- `rt`: Observed retention times
- `irt`: Predicted indexed retention times
- `rt_model`: Fitted RT conversion model
- `ms_file_idx`: MS file index
- `output_dir`: Output directory path

Creates a scatter plot of RT vs iRT with the fitted model curve overlay.
Saves to FirstPass RT alignment folder in QC plots directory.
"""
function plot_rt_alignment_firstpass(
    rt::Vector{Float32},
    irt::Vector{Float32},
    rt_model::RtConversionModel,
    ms_file_idx::Int,
    output_dir::String
)
    n = length(rt)

    # Create scatter plot
    p = plot(
        rt,
        irt,
        seriestype = :scatter,
        title = "FirstPass RT Alignment (File $ms_file_idx)\nn = $n",
        xlabel = "Observed RT (min)",
        ylabel = "Predicted iRT",
        label = nothing,
        alpha = 0.3,
        markersize = 2,
        size = (800, 600)
    )

    # Add fitted model line
    rt_range = LinRange(minimum(rt), maximum(rt), 100)
    fitted_irt = [rt_model(r) for r in rt_range]
    plot!(p, rt_range, fitted_irt, color = :red, linewidth = 2, label = nothing)

    # Create output directory if needed
    rt_plot_folder = joinpath(output_dir, "qc_plots", "rt_alignment_plots", "firstpass")
    !isdir(rt_plot_folder) && mkpath(rt_plot_folder)

    # Save plot
    plot_path = joinpath(rt_plot_folder, "file_$(ms_file_idx)_rt_alignment.pdf")
    savefig(p, plot_path)

    return nothing
end


PrecToIrtType = Dictionary{UInt32, 
    NamedTuple{
        (:best_prob, :best_ms_file_idx, :best_scan_idx, :best_irt, :mean_irt, :var_irt, :n, :mz), 
        Tuple{Float32, UInt32, UInt32, Float32, Union{Missing, Float32}, Union{Missing, Float32}, Union{Missing, UInt16}, Float32}
    }
}

"""
    create_rt_indices!(search_context::SearchContext, results::FirstPassSearchResults,
                      precursor_dict::PrecToIrtType, params::FirstPassSearchParameters)

Creates retention time indices for improved search efficiency.

# Arguments
- `search_context`: Search context containing MS data
- `results`: FirstPassSearch results
- `precursor_dict`: Dictionary mapping precursors to iRT data
- `params`: Search parameters

# Process
1. Calculates iRT errors using peak width and variation
2. Creates precursor to iRT mapping
3. Generates RT indices for each MS file
4. Stores indices in search context
"""
function create_rt_indices!(
    search_context::SearchContext,
    results::FirstPassSearchResults,
    precursor_dict::PrecToIrtType,
    params::FirstPassSearchParameters
)


    # Calculate iRT errors
    irt_errs = get_irt_errs(results.fwhms, precursor_dict, params)


    setIrtErrors!(search_context, irt_errs)

    # Create precursor to iRT mapping
    prec_to_irt = map(x -> (irt=x[:best_irt], mz=x[:mz]), 
                      precursor_dict)

    # Set up indices folder
    rt_indices_folder = joinpath(getDataOutDir(search_context), "temp_data", "rt_indices")
    !isdir(rt_indices_folder) && mkdir(rt_indices_folder)

    # Get PSM paths and RT models
    all_psm_paths = getFirstPassPsms(getMSData(search_context))
    all_rt_models = getRtIrtMap(search_context)

    # Filter to get valid file indices and their paths
    valid_indexed_data = []
    for (idx, path) in enumerate(all_psm_paths)
        if !isempty(path) && haskey(all_rt_models, idx)
            push!(valid_indexed_data, (idx, path))
        end
    end

    # If no valid paths, skip RT index creation
    if isempty(valid_indexed_data)
        @user_warn "No valid PSM files for RT index creation - all files may have failed"
        return nothing
    end

    # Extract just the paths and indices
    valid_indices = [idx for (idx, _) in valid_indexed_data]
    valid_psm_paths = [path for (_, path) in valid_indexed_data]
    valid_rt_models = Dict(i => all_rt_models[idx] for (i, idx) in enumerate(valid_indices))

    # Make RT indices only for valid files
    rt_index_paths = makeRTIndices(
        rt_indices_folder,
        valid_psm_paths,
        prec_to_irt,
        valid_rt_models,
        min_prob=params.max_prob_to_impute
    )

    # Map the results back to original indices
    for (new_idx, orig_idx) in enumerate(valid_indices)
        if new_idx <= length(rt_index_paths)
            setRtIndex!(getMSData(search_context), orig_idx, rt_index_paths[new_idx])
        end
    end
end


"""
    get_irt_errs(fwhms::Dictionary, prec_to_irt::Dictionary, params::FirstPassSearchParameters)

Calculates iRT error tolerances based on peak widths and cross-run variation.

# Arguments
- `fwhms`: Dictionary of FWHM statistics per file
- `prec_to_irt`: Dictionary of precursor iRT data
- `params`: Parameters including FWHM and iRT standard deviation multipliers

# Returns
Dictionary mapping file indices to iRT tolerances, combining:
- Peak width variation (FWHM + n*MAD)
- Cross-run iRT variation
"""
function get_irt_errs(
    fwhms::Dictionary{Int64, 
                        @NamedTuple{
                            median_fwhm::Float32,
                            mad_fwhm::Float32
                        }},
    prec_to_irt::Dictionary{UInt32, 
    @NamedTuple{best_prob::Float32, 
                best_ms_file_idx::UInt32, 
                best_scan_idx::UInt32, 
                best_irt::Float32, 
                mean_irt::Union{Missing, Float32}, 
                var_irt::Union{Missing, Float32}, 
                n::Union{Missing, UInt16}, 
                mz::Float32}}
    ,
    params::FirstPassSearchParameters
)
    #Get upper bound on peak fwhm. Use median + n*standard_deviation
    #estimate standard deviation by the median absolute deviation. 
    #n is a user-defined paramter. 
    fwhms = map(x->x[:median_fwhm] + params.fwhm_nstd*x[:mad_fwhm],
    fwhms)
    #Get variance in irt of apex accross runs. Only consider precursor identified below q-value threshold
    #in more than two runs .
    irt_std = nothing
    variance_  = collect(skipmissing(map(x-> (x[:n] > 2) ? sqrt(x[:var_irt]/(x[:n] - 1)) : missing, prec_to_irt)))
    if !iszero(length(variance_))
        irt_std = median(variance_)
    else
        #This could happen if only two files are being searched 
        variance_  = collect(skipmissing(map(x-> (x[:n] == 2) ? sqrt(x[:var_irt]) : missing, prec_to_irt)))
        if iszero(length(variance_)) #only searching one file so 
            irt_std = 0.0f0
        else
            irt_std = median(variance_)
        end
    end
    #Number of standard deviations to cover
    irt_std *= params.irt_nstd
    #dictionary maping file name to irt tolerance.
    return map(x->Float32((x+irt_std))::Float32, fwhms)::Dictionary{Int64, Float32}
end

"""
    mass_err_ms1(ppm_errs::Vector{Float32}, params::P) where {P<:FragmentIndexSearchParameters}

Estimate MS1 mass error model using exponential distribution fitting (consistent with MS2 approach).

Fits separate exponential distributions to positive and negative mass error tails,
then calculates tolerance bounds at the specified quantile level.

# Arguments
- `ppm_errs`: Vector of PPM mass errors from MS1 precursor matches
- `params`: Search parameters containing fragment error quantile

# Returns
- `(MassErrorModel, Vector{Float32})`: Fitted model and bias-corrected errors

# Details
Uses the same exponential distribution fitting approach as MS2 fragment tolerance estimation:
1. Calculates median bias and removes it from errors
2. Separates positive and negative errors for asymmetric modeling
3. Fits exponential distributions to each tail: Exponential(mean_of_tail)
4. Calculates tolerance bounds at 1-frag_err_quantile level
5. Returns zero tolerances if exponential fitting fails (consistent with MS2)
"""
function mass_err_ms1(ppm_errs::Vector{Float32},
    params::P) where {P<:FragmentIndexSearchParameters}
    mass_err = median(ppm_errs)
    pmp_errs_corrected = ppm_errs .- mass_err

    # Calculate error bounds using exponential distribution fitting (matching MS2 approach)
    frag_err_quantile = getFragErrQuantile(params)

    # Separate positive and negative errors for asymmetric modeling
    r_errs = abs.(pmp_errs_corrected[pmp_errs_corrected .> 0.0f0])
    l_errs = abs.(pmp_errs_corrected[pmp_errs_corrected .< 0.0f0])

    r_bound, l_bound = nothing, nothing
    try
        # Fit exponential distribution to right tail (positive errors)
        # Exponential(θ) where θ is the scale parameter (mean)
        r_bound = quantile(Exponential(
            sum(r_errs)/length(r_errs)  # Mean of positive errors
        ), 1 - frag_err_quantile)

        # Fit exponential distribution to left tail (negative errors)
        l_bound = quantile(Exponential(
            sum(l_errs)/length(l_errs)  # Mean of negative errors
        ), 1 - frag_err_quantile)
    catch e
        @debug_l1 "MS1 exponential fitting failed: length(r_errs)=$(length(r_errs)) length(l_errs)=$(length(l_errs)) sum(r_errs)=$(sum(r_errs)) sum(l_errs)=$(sum(l_errs))"
        # Return zero tolerances on fitting failure (consistent with MS2 approach)
        return MassErrorModel(
            Float32(mass_err),
            (zero(Float32), zero(Float32))
        ), pmp_errs_corrected
    end

    return MassErrorModel(
        Float32(mass_err),
        (Float32(abs(l_bound)), Float32(abs(r_bound)))
    ), pmp_errs_corrected
end



