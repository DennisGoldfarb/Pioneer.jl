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
    get_best_precursors_accross_runs(psms_paths::Vector{String},
                                    prec_mzs::AbstractVector{Float32},
                                    rt_to_library_irt::Dict{Int64, RtConversionModel};
                                    max_q_val::Float32=0.01f0,
                                    mad_trim_multiplier::Float32=6.0f0)
    -> Dictionary{UInt32, NamedTuple}

Identify and collect best precursor matches across multiple runs for retention time calibration.

# Arguments
- `psms_paths`: Paths to PSM files from first pass search
- `prec_mzs`: Vector of precursor m/z values
- `rt_to_library_irt`: Dictionary mapping file indices to RT→library_iRT conversion models
- `max_q_val`: Maximum q-value threshold for considering PSMs
- `mad_trim_multiplier`: MAD multiplier used to trim outlier iRTs before consensus calculation

# Returns
Dictionary mapping precursor indices to NamedTuple containing:
- `best_prob`: Highest probability score
- `best_ms_file_idx`: File index with best match
- `best_scan_idx`: Scan index of best match
- `best_library_irt`: Library iRT value of best match
- `consensus_library_irt`: Trimmed mean library iRT across qualifying matches
- `median_library_irt`: Median library iRT across qualifying matches
- `mad_library_irt`: Median absolute deviation of library iRT across qualifying matches
- `n`: Number of retained matches after trimming
- `mz`: Precursor m/z value

# Process
1. First pass: Collects best matches and records library iRT observations for each precursor across the runs
2. Filters to top N precursors by probability
3. Calculates median/MAD statistics and trims outlier runs prior to consensus iRT calculation
"""

function get_best_precursors_accross_runs(
                         psms_paths::Vector{String},
                         prec_mzs::AbstractVector{Float32},
                         rt_to_library_irt::Dict{Int64, RtConversionModel};
                         max_q_val::Float32 = 0.01f0,
                         mad_trim_multiplier::Float32 = 6.0f0
                         )

    precursor_irt_values = Dictionary{UInt32, Vector{Float32}}()

    function readPSMs!(
        prec_to_best_prob::Dictionary{UInt32, @NamedTuple{ best_prob::Float32,
                                                    best_ms_file_idx::UInt32,
                                                    best_scan_idx::UInt32,
                                                    best_library_irt::Float32,
                                                    consensus_library_irt::Float32,
                                                    median_library_irt::Float32,
                                                    mad_library_irt::Float32,
                                                    n::UInt16,
                                                    mz::Float32}},
        precursor_idxs::AbstractVector{UInt32},
        q_values::AbstractVector{Float16},
        probs::AbstractVector{Float32},
        rts::AbstractVector{Float32},
        scan_idxs::AbstractVector{UInt32},
        ms_file_idxs::AbstractVector{UInt32},
        rt_to_library_irt::RtConversionModel,
        max_q_val::Float32)

        for row in eachindex(precursor_idxs)
            # Extract current PSM information
            q_value = q_values[row]
            precursor_idx = precursor_idxs[row]
            prob = probs[row]
            library_irt = rt_to_library_irt(rts[row])
            scan_idx = UInt32(scan_idxs[row])
            ms_file_idx = UInt32(ms_file_idxs[row])

            # Initialize running statistics
            passed_q_val = (q_value <= max_q_val)
            n = passed_q_val ? one(UInt16) : zero(UInt16)
            mz = prec_mzs[precursor_idx]

            #Has the precursor been encountered in a previous raw file?
            if haskey(prec_to_best_prob, precursor_idx)
                # Update existing precursor entry
                best_prob, best_ms_file_idx, best_scan_idx, best_library_irt, consensus_library_irt, median_library_irt, mad_library_irt, old_n, mz = prec_to_best_prob[precursor_idx]

                # Update best match if current is better
                if (best_prob < prob)
                    best_prob = prob
                    best_library_irt = library_irt
                    best_scan_idx = scan_idx
                    best_ms_file_idx = ms_file_idx
                end

                # Update running statistics
                n += old_n
                prec_to_best_prob[precursor_idx] = (
                                                best_prob = best_prob,
                                                best_ms_file_idx = best_ms_file_idx,
                                                best_scan_idx = best_scan_idx,
                                                best_library_irt = best_library_irt,
                                                consensus_library_irt = consensus_library_irt,
                                                median_library_irt = median_library_irt,
                                                mad_library_irt = mad_library_irt,
                                                n = n,
                                                mz = mz)
            else
                # Create new precursor entry
                val = (best_prob = prob,
                        best_ms_file_idx = ms_file_idx,
                        best_scan_idx = scan_idx,
                        best_library_irt = library_irt,
                        consensus_library_irt = library_irt,
                        median_library_irt = library_irt,
                        mad_library_irt = 0.0f0,
                        n = n,
                        mz = mz)
                insert!(prec_to_best_prob, precursor_idx, val)
            end

            if passed_q_val
                push!(get!(precursor_irt_values, precursor_idx, Float32[]), library_irt)
            end
        end
    end

    function finalize_irt_statistics!(
        prec_to_best_prob::Dictionary{UInt32, @NamedTuple{ best_prob::Float32,
                                                    best_ms_file_idx::UInt32,
                                                    best_scan_idx::UInt32,
                                                    best_library_irt::Float32,
                                                    consensus_library_irt::Float32,
                                                    median_library_irt::Float32,
                                                    mad_library_irt::Float32,
                                                    n::UInt16,
                                                    mz::Float32}},
        precursor_irt_values::Dictionary{UInt32, Vector{Float32}},
        mad_trim_multiplier::Float32)
        for (precursor_idx, val) in prec_to_best_prob
            irts = get(precursor_irt_values, precursor_idx, Float32[])
            if isempty(irts)
                prec_to_best_prob[precursor_idx] = (
                    best_prob = val.best_prob,
                    best_ms_file_idx = val.best_ms_file_idx,
                    best_scan_idx = val.best_scan_idx,
                    best_library_irt = val.best_library_irt,
                    consensus_library_irt = val.consensus_library_irt,
                    median_library_irt = val.median_library_irt,
                    mad_library_irt = val.mad_library_irt,
                    n = zero(UInt16),
                    mz = val.mz
                )
                continue
            end

            med = median(irts)
            abs_dev = abs.(irts .- med)
            mad_val = median(abs_dev)
            trimmed_irts = if mad_val == 0f0
                irts
            else
                filter(x -> abs(x - med) <= mad_trim_multiplier * mad_val, irts)
            end
            trimmed_median = median(trimmed_irts)
            trimmed_mean = mean(trimmed_irts)
            consensus_irt = length(trimmed_irts) == 1 ? trimmed_irts[1] : trimmed_mean

            prec_to_best_prob[precursor_idx] = (
                best_prob = val.best_prob,
                best_ms_file_idx = val.best_ms_file_idx,
                best_scan_idx = val.best_scan_idx,
                best_library_irt = val.best_library_irt,
                consensus_library_irt = Float32(consensus_irt),
                median_library_irt = Float32(trimmed_median),
                mad_library_irt = Float32(mad_val),
                n = UInt16(length(trimmed_irts)),
                mz = val.mz
            )
        end
    end
    # Initialize dictionary to store best precursor matches
    prec_to_best_prob = Dictionary{UInt32, @NamedTuple{ best_prob::Float32,
                                                        best_ms_file_idx::UInt32,
                                                        best_scan_idx::UInt32,
                                                        best_library_irt::Float32,
                                                        consensus_library_irt::Float32,
                                                        median_library_irt::Float32,
                                                        mad_library_irt::Float32,
                                                        n::UInt16,
                                                        mz::Float32}}()

    # First pass: collect best matches and mean library iRT

    n_precursors_vec = Vector{UInt64}()
    for psms_path in psms_paths #For each data frame 
        psms = Arrow.Table(psms_path)
        
        # Get the original file index from the PSM data
        if isempty(psms[:ms_file_idx])
            continue  # Skip empty files
        end
        file_idx = first(psms[:ms_file_idx])  # All PSMs in a file should have the same ms_file_idx

        # Check if RT model exists for this file
        if !haskey(rt_to_library_irt, file_idx)
            @warn "No RT model found for file index $file_idx, skipping"
            continue
        end

        push!(n_precursors_vec, length(psms[:precursor_idx]))
        #One row for each precursor
        readPSMs!(
            prec_to_best_prob,
            psms[:precursor_idx],
            psms[:q_value],
            psms[:prob],
            psms[:rt],
            psms[:scan_idx],
            psms[:ms_file_idx],
            rt_to_library_irt[file_idx],
            max_q_val
        )
    end
    
    # Handle case where no valid files were processed
    if isempty(n_precursors_vec)
        @warn "No valid PSM files found for cross-run analysis"
        return prec_to_best_prob
    end
    
    max_precursors = maximum(n_precursors_vec)
    # Filter to top N precursors by probability
    sort!(prec_to_best_prob, by = x->x[:best_prob], alg=PartialQuickSort(1:max_precursors), rev = true);
    N = 0
    for key in collect(keys(prec_to_best_prob))
        N += 1
        if N > max_precursors
            delete!(prec_to_best_prob, key)
        end
    end

    # Calculate trimmed consensus statistics
    finalize_irt_statistics!(prec_to_best_prob, precursor_irt_values, mad_trim_multiplier)

    return prec_to_best_prob #[(prob, idx) for (idx, prob) in sort(collect(top_probs), rev=true)]
end
