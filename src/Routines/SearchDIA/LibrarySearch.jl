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

function searchFragmentIndex(
                    scan_to_prec_idx::Vector{Union{Missing, UnitRange{Int64}}},
                    frag_index::FragmentIndex{Float32},
                    spectra::MassSpecData,
                    thread_task::Vector{Int64},
                    search_data::S,
                    params::P,
                    qtm::Q,
                    mem::M,
                    rt_to_irt_spline::Any,
                    irt_tol::AbstractFloat
                    ) where {M<:MassErrorModel, Q<:QuadTransmissionModel, S<:SearchDataStructures, P<:FragmentIndexSearchParameters}

    prec_id = 0
    precursors_passed_scoring = Vector{UInt32}(undef, 250000)
    rt_bin_idx = 1

    base_counter = getPrecursorScores(search_data)
    counter_size = length(base_counter.counts)
    raw_counters = Counter{UInt32, UInt8}[base_counter, Counter(UInt32, UInt8, counter_size), Counter(UInt32, UInt8, counter_size)]
    free_counters = Counter{UInt32, UInt8}[raw_counters...]
    smoothed_counter = Counter(UInt32, UInt8, counter_size)

    @inline uses_smoothed_scores(::FragmentIndexSearchParameters) = true
    @inline uses_smoothed_scores(::ParameterTuningSearchParameters) = false
    @inline uses_smoothed_scores(::NceTuningSearchParameters) = false
    @inline uses_smoothed_scores(::QuadTuningSearchParameters) = false

    local_score_threshold(min_score::UInt8) = UInt8(max(Int(min_score) - 1, 0))

    acquire_counter!() = pop!(free_counters)

    function release_counter!(counter::Counter{UInt32, UInt8})
        reset!(counter)
        push!(free_counters, counter)
        return nothing
    end

    @inline function normalize_float(val)
        return ismissing(val) ? val : Float32(val)
    end

    @inline function same_isolation_window(center_a, width_a, center_b, width_b)
        return !(ismissing(center_a) || ismissing(width_a) || ismissing(center_b) || ismissing(width_b)) &&
               center_a == center_b && width_a == width_b
    end

    function add_scores!(dest::Counter{UInt32, UInt8}, source::Counter{UInt32, UInt8})
        @inbounds @fastmath for i in 1:(getSize(source) - 1)
            id = getID(source, i)
            or_inc!(dest, id, getCount(source, id))
        end
        return nothing
    end

    function smooth_scores!(dest::Counter{UInt32, UInt8}, target, left, right)
        reset!(dest)
        if left !== nothing && same_isolation_window(target.center_mz, target.isolation_width, left.center_mz, left.isolation_width)
            add_scores!(dest, left.counter)
        end
        add_scores!(dest, target.counter)
        if right !== nothing && same_isolation_window(target.center_mz, target.isolation_width, right.center_mz, right.isolation_width)
            add_scores!(dest, right.counter)
        end
        return nothing
    end

    function record_matches!(counter::Counter{UInt32, UInt8}, scan_idx::Int64)
        if getID(counter, 1) > 0
            start_idx = prec_id + 1
            n = 1
            while n <= counter.matches
                prec_id += 1
                if prec_id > length(precursors_passed_scoring)
                    append!(precursors_passed_scoring,
                            Vector{eltype(precursors_passed_scoring)}(undef, length(precursors_passed_scoring))
                            )
                end
                precursors_passed_scoring[prec_id] = getID(counter, n)
                n += 1
            end
            scan_to_prec_idx[scan_idx] = start_idx:prec_id#stop_idx
        else
            scan_to_prec_idx[scan_idx] = missing
        end
        return nothing
    end

    function finalize_scan!(target, left, right)
        min_score = getMinIndexSearchScore(params)
        if uses_smoothed_scores(params)
            smooth_scores!(smoothed_counter, target, left, right)
            filterPrecursorMatches!(
                smoothed_counter,
                target.counter,
                min_score,
                local_score_threshold(min_score)
            )
            record_matches!(smoothed_counter, target.scan_idx)
        else
            filterPrecursorMatches!(target.counter, min_score)
            record_matches!(target.counter, target.scan_idx)
        end
        return nothing
    end

    prev_scan = nothing
    curr_scan = nothing

    for scan_idx in thread_task
        (scan_idx <= 0 || scan_idx > length(spectra)) && continue
        getMsOrder(spectra, scan_idx) ∉ getSpecOrder(params) && continue

        raw_counter = acquire_counter!()

        center_mz = normalize_float(getCenterMz(spectra, scan_idx))
        isolation_width = normalize_float(getIsolationWidthMz(spectra, scan_idx))

        # Update RT bin index based on iRT window
        irt_lo, irt_hi = getRTWindow(rt_to_irt_spline(getRetentionTime(spectra, scan_idx)), irt_tol)
        while rt_bin_idx < length(getRTBins(frag_index)) && getHigh(getRTBin(frag_index, rt_bin_idx)) < irt_lo
            rt_bin_idx += 1
        end
        while rt_bin_idx > 1 && getLow(getRTBin(frag_index, rt_bin_idx)) > irt_lo
            rt_bin_idx -= 1
        end

        # Fragment index search for matching precursors
        searchScan!(
            raw_counter,
            getRTBins(frag_index),
            getFragBins(frag_index),
            getFragments(frag_index),
            getMzArray(spectra, scan_idx),
            getIntensityArray(spectra, scan_idx),
            rt_bin_idx,
            irt_hi,
            mem,
            getQuadTransmissionFunction(qtm, center_mz, isolation_width),
            getIsotopeErrBounds(params)
        )

        new_scan = (;
            counter = raw_counter,
            center_mz = center_mz,
            isolation_width = isolation_width,
            scan_idx = scan_idx,
        )

        if curr_scan !== nothing
            finalize_scan!(curr_scan, prev_scan, new_scan)
            if prev_scan !== nothing
                release_counter!(prev_scan.counter)
            end
        end

        prev_scan = curr_scan
        curr_scan = new_scan
    end

    if curr_scan !== nothing
        finalize_scan!(curr_scan, prev_scan, nothing)
        if prev_scan !== nothing
            release_counter!(prev_scan.counter)
        end
        release_counter!(curr_scan.counter)
    end

    return precursors_passed_scoring[1:prec_id]
end

function getPSMS(
    ms_file_idx::UInt32,
    spectra::MassSpecData,
    thread_task::Vector{Int64},
    precursors::LibraryPrecursors,
    ion_list::LibraryFragmentLookup,
    nce_model::NceModel{Float32},
    scan_to_prec_idx::Vector{Union{Missing, UnitRange{Int64}}},
    precursors_passed_scoring::Vector{UInt32},
    search_data::S,
    params::P,
    qtm::Q,
    mem::M,
    rt_to_irt_spline::Any,
    irt_tol::AbstractFloat) where {M<:MassErrorModel, Q<:QuadTransmissionModel, S<:SearchDataStructures, P<:SearchParameters}

    msms_counts = Dict{Int64, Int64}()
    last_val = 0
    Hs = SparseArray(UInt32(5000))
    isotopes = zeros(Float32, 5)
    precursor_transmission = zeros(Float32, 5)
    for scan_idx in thread_task
        (scan_idx == 0 || scan_idx > length(spectra)) && continue
        ismissing(scan_to_prec_idx[scan_idx]) && continue

        # Scan Filtering
        msn = getMsOrder(spectra, scan_idx)
        msn ∉ getSpecOrder(params) && continue
        msn ∈ keys(msms_counts) ? msms_counts[msn] += 1 : msms_counts[msn] = 1

        # Ion Template Selection
        ion_idx, _ = selectTransitions!(
            getIonTemplates(search_data),
            StandardTransitionSelection(),
            getPrecEstimation(params),
            ion_list,
            nce_model,
            scan_to_prec_idx[scan_idx], precursors_passed_scoring,
            getMz(precursors),
            getCharge(precursors),
            getSulfurCount(precursors),
            getIrt(precursors),
            getIsoSplines(search_data),
            getQuadTransmissionFunction(qtm, getCenterMz(spectra, scan_idx), getIsolationWidthMz(spectra, scan_idx)),
            precursor_transmission, isotopes, getNFragIsotopes(params),
            getMaxFragRank(params),
            Float32(rt_to_irt_spline(getRetentionTime(spectra, scan_idx))),
            Float32(irt_tol),
            (getLowMz(spectra, scan_idx), getHighMz(spectra, scan_idx));
            isotope_err_bounds = getIsotopeErrBounds(params)
        )

        ion_idx < 2 && continue


        # Match peaks
        nmatches, nmisses = matchPeaks!(
            getIonMatches(search_data), 
            getIonMisses(search_data), 
            getIonTemplates(search_data), 
            ion_idx, 
            getMzArray(spectra, scan_idx), 
            getIntensityArray(spectra, scan_idx), 
            mem,
            getHighMz(spectra, scan_idx),
            UInt32(scan_idx), 
            ms_file_idx
        )
        
        sort!(@view(getIonMatches(search_data)[1:nmatches]), by = x->(x.peak_ind, x.prec_id), alg=QuickSort)
        # Process matches
        if nmatches > 2
            buildDesignMatrix!(Hs, getIonMatches(search_data), getIonMisses(search_data), 
                             nmatches, nmisses, getIdToCol(search_data))

            if getIdToCol(search_data).size > length(getSpectralScores(search_data))
                new_entries = getIdToCol(search_data).size - length(getSpectralScores(search_data)) + 1000 
                append!(getSpectralScores(search_data), Vector{eltype(getSpectralScores(search_data))}(undef, new_entries))
                append!(getUnscoredPsms(search_data), [eltype(getUnscoredPsms(search_data))() for _ in 1:new_entries])
            end 

            getDistanceMetrics(Hs, getSpectralScores(search_data); params.relative_improvement_threshold)

            ScoreFragmentMatches!(
                getUnscoredPsms(search_data),
                getIdToCol(search_data),
                getIonMatches(search_data), 
                nmatches, 
                mem,
                last(getMinTopNofM(params))
            )

            last_val = Score!(
                getScoredPsms(search_data), 
                getUnscoredPsms(search_data),
                getSpectralScores(search_data),
                getIdToCol(search_data),
                nmatches/(nmatches + nmisses),
                last_val,
                Hs.n,
                Float32(sum(getIntensityArray(spectra, scan_idx))), 
                scan_idx,
                min_spectral_contrast = getMinSpectralContrast(params),
                min_log2_matched_ratio = getMinLog2MatchedRatio(params),
                min_frag_count = getMinFragCount(params),
                max_best_rank = getMaxBestRank(params),
                min_topn = first(getMinTopNofM(params)),
                block_size = 500000
            )
        end
        # Reset arrays
        for scan_idx in range(1, Hs.n)
            getUnscoredPsms(search_data)[scan_idx] = eltype(getUnscoredPsms(search_data))()
        end
        reset!(getIdToCol(search_data))
        reset!(Hs)
    end
    return DataFrame(@view(getScoredPsms(search_data)[1:last_val]))
end

function library_search(spectra::MassSpecData, search_context::SearchContext, search_parameters::P, ms_file_idx::Int64) where {P<:ParameterTuningSearchParameters}
    return vcat(LibrarySearch(
                    spectra,
                    UInt32(ms_file_idx),
                    getPresearchFragmentIndex(getSpecLib(search_context)),
                    getSpecLib(search_context),
                    getSearchData(search_context),
                    getQuadTransmissionModel(search_context, ms_file_idx),
                    getMassErrorModel(search_context, ms_file_idx),
                    getRtIrtModel(search_context, ms_file_idx),
                    search_parameters,
                    getNceModel(search_context, ms_file_idx),
                    getIRTTol(search_parameters),
                )...)
end

function library_search(spectra::MassSpecData, search_context::SearchContext, search_parameters::P, ms_file_idx::Int64) where {P<:SearchParameters}

    return vcat(LibrarySearch(
                    spectra,
                    UInt32(ms_file_idx),
                    getFragmentIndex(getSpecLib(search_context)),
                    getSpecLib(search_context),
                    getSearchData(search_context),
                    getQuadTransmissionModel(search_context, ms_file_idx),
                    getMassErrorModel(search_context, ms_file_idx),
                    getRtIrtModel(search_context, ms_file_idx),
                    search_parameters,
                    getNceModel(search_context, ms_file_idx),
                    getIrtErrors(search_context)[ms_file_idx]
                )...)
end

function library_search(
    spectra::MassSpecData,
    search_context::SearchContext,
    search_parameters::P,
    ms_file_idx::Int64) where {P<:NceTuningSearchParameters}

    return LibrarySearchNceTuning(
        spectra,
        UInt32(ms_file_idx),
        getFragmentIndex(getSpecLib(search_context)),
        getSpecLib(search_context),
        getSearchData(search_context),
        getQuadTransmissionModel(search_context, ms_file_idx),
        getMassErrorModel(search_context, ms_file_idx),
        getRtIrtModel(search_context, ms_file_idx),
        search_parameters,
        search_parameters.nce_grid,
        getIrtErrors(search_context)[ms_file_idx]
    )
end

function LibrarySearch(
    spectra::MassSpecData,
    ms_file_idx::UInt32,
    fragment_index::FragmentIndex{Float32},
    spec_lib::SpectralLibrary,
    search_data::AbstractVector{S},
    qtm::Q,
    mem::M,
    rt_to_irt_spline::Any,
    params::P,
    nce_model::NceModel{Float32},
    irt_tol::AbstractFloat) where {
        M<:MassErrorModel,
        Q<:QuadTransmissionModel,
        S<:SearchDataStructures,
        P<:FragmentIndexSearchParameters}
    thread_tasks = partition_scans(spectra, Threads.nthreads())

    scan_to_prec_idx = Vector{Union{Missing, UnitRange{Int64}}}(undef, length(spectra))
    tasks = map(thread_tasks) do thread_task
        Threads.@spawn begin 
            thread_id = first(thread_task)
                return searchFragmentIndex(
                        scan_to_prec_idx,
                        fragment_index,
                        spectra,
                        last(thread_task),
                        search_data[thread_id],
                        params,
                        qtm,
                        mem,
                        rt_to_irt_spline,
                        irt_tol
                    )
        end
    end
    
    precursors_passed_scoring = fetch.(tasks)

    tasks = map(thread_tasks) do thread_task
        Threads.@spawn begin 
            thread_id = first(thread_task)
            return getPSMS(
                                ms_file_idx,
                                spectra,
                                last(thread_task), #getRange(thread_task),
                                getPrecursors(spec_lib),
                                getFragmentLookupTable(spec_lib),
                                nce_model,
                                scan_to_prec_idx,
                                precursors_passed_scoring[thread_id],
                                search_data[thread_id],
                                params,
                                qtm,
                                mem,
                                rt_to_irt_spline,
                                irt_tol)
        end
    end

    return fetch.(tasks)
end

function LibrarySearchNceTuning(
    spectra::MassSpecData,
    ms_file_idx::UInt32,
    fragment_index::FragmentIndex{Float32},
    spec_lib::SpectralLibrary,
    search_data::AbstractVector{S},
    qtm::Q,
    mem::M,
    rt_to_irt_spline::Any,
    params::P,
    nce_grid::AbstractVector{<:AbstractFloat},
    irt_tol::AbstractFloat) where {
        M<:MassErrorModel,
        Q<:QuadTransmissionModel,
        S<:SearchDataStructures,
        P<:FragmentIndexSearchParameters}

    # Check for valid inputs
    if isempty(nce_grid)
        @user_warn "LibrarySearchNceTuning: Empty NCE grid provided"
        return DataFrame()
    end

    @debug_l2 "LibrarySearchNceTuning: Starting with $(length(spectra)) scans, NCE grid: $nce_grid"

    thread_tasks = partition_scans(spectra, Threads.nthreads())

    @debug_l2 "LibrarySearchNceTuning: Created $(length(thread_tasks)) thread tasks"
    for (i, task) in enumerate(thread_tasks)
        task_range = last(task)
        n_scans = length(task_range)
        @debug_l2 "  Thread $i: $n_scans scans"
    end
    scan_to_prec_idx = Vector{Union{Missing, UnitRange{Int64}}}(undef, length(spectra))
    # Do fragment index search once
    tasks = map(thread_tasks) do thread_task
        Threads.@spawn begin
            thread_id = first(thread_task)
            try
                return searchFragmentIndex(
                    scan_to_prec_idx,
                    fragment_index,
                    spectra,
                    last(thread_task),
                    search_data[thread_id],
                    params,
                    qtm,
                    mem,
                    rt_to_irt_spline,
                    irt_tol
                )
            catch e
                @user_warn "Fragment index search failed on thread $thread_id: $e"
                return UInt32[]  # Return empty result on error
            end
        end
    end

    precursors_passed_scoring = fetch.(tasks)

    # For each NCE value, run getPSMS using the same fragment index results
    all_results = map(nce_grid) do nce
        # Create NCE model for this grid point
        nce_model = PiecewiseNceModel(nce)

        # Run getPSMS with this NCE model
        tasks = map(thread_tasks) do thread_task
            Threads.@spawn begin
                thread_id = first(thread_task)
                try
                    psms = getPSMS(
                        ms_file_idx,
                        spectra,
                        last(thread_task),
                        getPrecursors(spec_lib),
                        getFragmentLookupTable(spec_lib),
                        nce_model,
                        scan_to_prec_idx,
                        precursors_passed_scoring[thread_id],
                        search_data[thread_id],
                        params,
                        qtm,
                        mem,
                        rt_to_irt_spline,
                        irt_tol
                    )
                    if !isempty(psms)
                        psms[!, :nce] .= nce
                        return psms
                    else
                        return missing
                    end
                catch e
                    @user_warn "PSM search failed on thread $thread_id for NCE $nce: $e"
                    return missing
                end
            end
        end
        # Ensure we always return a DataFrame, not Vector{Any}
        fetched = collect(skipmissing(fetch.(tasks)))
        isempty(fetched) ? DataFrame() : vcat(fetched...)
    end

    # filter out empty DFs (which are actually Vectors instead of DataFrames)
    nonempty_dfs = filter(df -> df isa DataFrame, all_results)
    # Ensure we always return a DataFrame, not Vector{Any}
    return isempty(nonempty_dfs) ? DataFrame() : vcat(nonempty_dfs...)
end

function getRTWindow(irt::U, irt_tol::T) where {T,U<:AbstractFloat}
    return Float32(irt - irt_tol), Float32(irt + irt_tol)
end
