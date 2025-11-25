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
                                    rt_to_library_irt::Dict{Int64, RtConversionModel},
                                    rt_uncertainty::Dict{Int64, Float32};
                                    max_q_val::Float32=0.01f0,
                                    prior_variance::Float32=0.01f0)
    -> Dictionary{UInt32, NamedTuple}

Identify and collect best precursor matches across multiple runs for retention time calibration.

# Arguments
- `psms_paths`: Paths to PSM files from first pass search
- `prec_mzs`: Vector of precursor m/z values
- `rt_to_library_irt`: Dictionary mapping file indices to RT→library_iRT conversion models
- `rt_uncertainty`: Mapping from ms_file_idx to RT variance estimates (e.g., spline residuals)
- `max_q_val`: Maximum q-value threshold for considering PSMs
- `prior_variance`: Small variance term added to stabilize weights when limited runs are available

# Returns
Dictionary mapping precursor indices to NamedTuple containing:
- `best_prob`: Highest probability score
- `best_ms_file_idx`: File index with best match
- `best_scan_idx`: Scan index of best match
- `best_library_irt`: Library iRT value of best match
- `mean_library_irt`: Weighted mean library iRT across qualifying matches
- `var_library_irt`: Weighted variance accumulator across qualifying matches
- `total_rt_weight`: Sum of inverse-variance RT weights across runs
- `n`: Number of qualifying matches
- `mz`: Precursor m/z value

# Process
1. First pass: Collects best matches and calculates mean library iRT for each precursor across the runs
2. Filters to top N precursors by probability
3. Second pass: Calculates library iRT variance for remaining precursors
"""
function get_best_precursors_accross_runs(
                         psms_paths::Vector{String},
                         prec_mzs::AbstractVector{Float32},
                         rt_to_library_irt::Dict{Int64, RtConversionModel},
                         rt_uncertainty::Dict{Int64, Float32};
                         max_q_val::Float32 = 0.01f0,
                         prior_variance::Float32 = 0.01f0
                         )

    function readPSMs!(
        prec_to_best_prob::Dictionary{UInt32, @NamedTuple{ best_prob::Float32,
                                                    best_ms_file_idx::UInt32,
                                                    best_scan_idx::UInt32,
                                                    best_library_irt::Float32,
                                                    mean_library_irt::Union{Missing, Float32},
                                                    var_library_irt::Union{Missing, Float32},
                                                    n::Union{Missing, UInt16},
                                                    mz::Float32,
                                                    total_rt_weight::Float32}},
        precursor_idxs::AbstractVector{UInt32},
        q_values::AbstractVector{Float16},
        probs::AbstractVector{Float32},
        rts::AbstractVector{Float32},
        scan_idxs::AbstractVector{UInt32},
        ms_file_idxs::AbstractVector{UInt32},
        rt_to_library_irt::RtConversionModel,
        rt_uncertainty::Dict{Int64, Float32},
        max_q_val::Float32,
        prior_variance::Float32)

        function rt_weight(ms_file_idx::UInt32)
            variance = get(rt_uncertainty, Int(ms_file_idx), 0.0f0) + prior_variance
            variance = max(variance, eps(Float32))
            return Float32(inv(variance))
        end

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
            mean_library_irt = passed_q_val ? library_irt : missing
            var_library_irt = zero(Float32)
            total_rt_weight = passed_q_val ? rt_weight(ms_file_idx) : 0.0f0
            mz = prec_mzs[precursor_idx]

            #Has the precursor been encountered in a previous raw file?
            #Keep a running mean library_irt for instances below q-val threshold
            if haskey(prec_to_best_prob, precursor_idx)
                # Update existing precursor entry
                best_prob, best_ms_file_idx, best_scan_idx, best_library_irt, old_mean_library_irt, var_library_irt, old_n, mz, old_weight = prec_to_best_prob[precursor_idx]

                # Update best match if current is better
                if (best_prob < prob)
                    best_prob = prob
                    best_library_irt = library_irt
                    best_scan_idx = scan_idx
                    best_ms_file_idx = ms_file_idx
                end

                if passed_q_val
                    new_weight = rt_weight(ms_file_idx)
                    mean_library_irt = coalesce(old_mean_library_irt, library_irt)
                    total_weight = old_weight + new_weight
                    mean_library_irt += (library_irt - mean_library_irt) * (new_weight / total_weight)
                    total_rt_weight = total_weight
                    n = old_n + one(UInt16)
                else
                    mean_library_irt = old_mean_library_irt
                    total_rt_weight = old_weight
                    n = old_n
                end
                prec_to_best_prob[precursor_idx] = (
                                                best_prob = best_prob,
                                                best_ms_file_idx = best_ms_file_idx,
                                                best_scan_idx = best_scan_idx,
                                                best_library_irt = best_library_irt,
                                                mean_library_irt = mean_library_irt,
                                                var_library_irt = var_library_irt,
                                                n = n,
                                                mz = mz,
                                                total_rt_weight = total_rt_weight)
            else
                # Create new precursor entry
                val = (best_prob = prob,
                        best_ms_file_idx = ms_file_idx,
                        best_scan_idx = scan_idx,
                        best_library_irt = library_irt,
                        mean_library_irt = mean_library_irt,
                        var_library_irt = var_library_irt,
                        n = n,
                        mz = mz,
                        total_rt_weight = total_rt_weight)
                insert!(prec_to_best_prob, precursor_idx, val)
            end
        end
    end
    function getVariance!(
        prec_to_best_prob::Dictionary{UInt32, @NamedTuple{ best_prob::Float32,
                                                    best_ms_file_idx::UInt32,
                                                    best_scan_idx::UInt32,
                                                    best_library_irt::Float32,
                                                    mean_library_irt::Union{Missing, Float32},
                                                    var_library_irt::Union{Missing, Float32},
                                                    n::Union{Missing, UInt16},
                                                    mz::Float32,
                                                    total_rt_weight::Float32}},
        precursor_idxs::AbstractVector{UInt32},
        q_values::AbstractVector{Float16},
        rts::AbstractVector{Float32},
        ms_file_idxs::AbstractVector{UInt32},
        rt_to_library_irt::RtConversionModel,
        rt_uncertainty::Dict{Int64, Float32},
        max_q_val::Float32,
        prior_variance::Float32)
        function rt_weight(ms_file_idx::UInt32)
            variance = get(rt_uncertainty, Int(ms_file_idx), 0.0f0) + prior_variance
            variance = max(variance, eps(Float32))
            return Float32(inv(variance))
        end
        for row in eachindex(precursor_idxs)
            # Skip PSMs that don't pass q-value threshold
            q_value = q_values[row]

            # Get precursor info
            precursor_idx = precursor_idxs[row]
            library_irt = rt_to_library_irt(rts[row])
            ms_file_idx = UInt32(ms_file_idxs[row])

            if q_value > max_q_val
                continue
            end
            if haskey(prec_to_best_prob, precursor_idx)
                # Update variance calculation
                best_prob, best_ms_file_idx, best_scan_idx, best_library_irt, mean_library_irt, var_library_irt, n, mz, total_rt_weight = prec_to_best_prob[precursor_idx]
                if total_rt_weight > 0
                    weight = rt_weight(ms_file_idx)
                    delta = library_irt - coalesce(mean_library_irt, library_irt)
                    var_library_irt += weight * delta^2
                end
                prec_to_best_prob[precursor_idx] = (
                    best_prob = best_prob,
                    best_ms_file_idx= best_ms_file_idx,
                    best_scan_idx = best_scan_idx,
                    best_library_irt = best_library_irt,
                    mean_library_irt = mean_library_irt,
                    var_library_irt = var_library_irt,
                    n = n,
                    mz = mz,
                    total_rt_weight = total_rt_weight),

            end
        end
    end
    # Initialize dictionary to store best precursor matches
    prec_to_best_prob = Dictionary{UInt32, @NamedTuple{ best_prob::Float32,
                                                        best_ms_file_idx::UInt32,
                                                        best_scan_idx::UInt32,
                                                        best_library_irt::Float32,
                                                        mean_library_irt::Union{Missing, Float32},
                                                        var_library_irt::Union{Missing, Float32},
                                                        n::Union{Missing, UInt16},
                                                        mz::Float32,
                                                        total_rt_weight::Float32}}()

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
            rt_uncertainty,
            max_q_val,
            prior_variance
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

    # Second pass: calculate library iRT variance for remaining precursors
    for psms_path in psms_paths #For each data frame
        psms = Arrow.Table(psms_path)

        # Get the original file index from the PSM data
        if isempty(psms[:ms_file_idx])
            continue  # Skip empty files
        end
        file_idx = first(psms[:ms_file_idx])  # All PSMs in a file should have the same ms_file_idx

        # Check if RT model exists for this file
        if !haskey(rt_to_library_irt, file_idx)
            continue  # Skip files without RT models (already warned in first pass)
        end

        #One row for each precursor
        getVariance!(
            prec_to_best_prob,
            psms[:precursor_idx],
            psms[:q_value],
            psms[:rt],
            psms[:ms_file_idx],
            rt_to_library_irt[file_idx],
            rt_uncertainty,
            max_q_val,
            prior_variance
        )
    end 

    
    return prec_to_best_prob #[(prob, idx) for (idx, prob) in sort(collect(top_probs), rev=true)]
end
