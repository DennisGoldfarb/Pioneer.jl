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
                                    prec_mzs::AbstractVector{<:AbstractFloat},
                                    rt_irt::Dict{Int64, RtConversionModel},
                                    is_decoy::AbstractVector{Bool},
                                    sequences::AbstractVector,
                                    structural_mods::AbstractVector;
                                    max_q_val::Float32=0.01f0,
                                    pep_threshold::Float32=0.9f0,
                                    fdr_scale_factor::Float32=1.0f0)
    -> Dictionary{UInt32, NamedTuple}

Identify and collect high-confidence precursors across multiple runs for retention time calibration.

# Arguments
- `psms_paths`: Paths to PSM files from first pass search
- `prec_mzs`: Vector of precursor m/z values
- `rt_irt`: Dictionary mapping file indices to RT-iRT conversion models
- `is_decoy`: Indicator vector for decoy precursors
- `sequences`: Precursor peptide sequences
- `structural_mods`: Structural modification annotations for each precursor
- `max_q_val`: Maximum q-value threshold for contributing to RT statistics
- `pep_threshold`: Posterior error probability threshold applied per run
- `fdr_scale_factor`: Scale factor correcting target/decoy imbalance when estimating PEP

# Returns
Dictionary mapping precursor indices to NamedTuple containing:
- `best_prob`: Run-level probability after peptide backpropagation
- `best_ms_file_idx`: File index with best match
- `best_scan_idx`: Scan index of best match
- `best_irt`: iRT value of best match
- `mean_irt`: Sum of iRT observations below `max_q_val`
- `var_irt`: Sum of squared deviations for iRT observations
- `n`: Count of iRT observations contributing to statistics
- `mz`: Precursor m/z value

# Process
1. For each run, keep the best score per precursor/isotope combination.
2. Combine isotope scores into precursor probabilities and group precursors by modified peptide (sequence + structural mods).
3. Backpropagate peptide confidence to precursors and estimate precursor-level PEP within the run.
4. Retain precursors passing the PEP threshold, accumulating RT statistics across runs.
5. Limit retained precursors to the top-N by probability and compute cross-run iRT variance for qualifying precursors.
"""
function get_best_precursors_accross_runs(
                         psms_paths::Vector{String},
                         prec_mzs::AbstractVector{<:AbstractFloat},
                         rt_irt::Dict{Int64, RtConversionModel},
                         is_decoy::AbstractVector{Bool},
                         sequences::AbstractVector,
                         structural_mods::AbstractVector;
                         max_q_val::Float32 = 0.01f0,
                         pep_threshold::Float32 = 0.9f0,
                         fdr_scale_factor::Float32 = 1.0f0
                         )

    combine_probability_evidence(probs::AbstractVector{Float32}) = begin
        log_prod = 0.0f0
        has_positive = false
        @inbounds for prob in probs
            if prob <= 0f0
                continue
            end
            has_positive = true
            clamped = min(prob, 1.0f0 - 1f-6)
            log_prod += log1p(-clamped)
        end
        if !has_positive
            return 0.0f0
        end
        combined = 1.0f0 - exp(log_prod)
        return clamp(combined, 0.0f0, 1.0f0 - 1f-6)
    end

    normalize_isotopes(value) = begin
        if value === missing
            return (Int8(-1), Int8(-1))
        elseif value isa Tuple
            return (Int8(first(value)), Int8(last(value)))
        else
            return (Int8(-1), Int8(-1))
        end
    end

    function peptide_key_for_precursor(prec_idx::UInt32)
        idx = Int(prec_idx)
        seq_val = sequences[idx]
        struct_val = structural_mods[idx]
        seq_key = seq_val === missing ? "" : String(seq_val)
        struct_key = struct_val === missing ? "" : String(struct_val)
        return (is_decoy[idx], seq_key, struct_key)
    end

    prec_to_best_prob = Dictionary{UInt32, @NamedTuple{
        best_prob::Float32,
        best_ms_file_idx::UInt32,
        best_scan_idx::UInt32,
        best_irt::Float32,
        mean_irt::Union{Missing, Float32},
        var_irt::Union{Missing, Float32},
        n::Union{Missing, UInt16},
        mz::Float32
    }}()

    n_precursors_vec = Vector{UInt64}()

    for psms_path in psms_paths
        psms = Arrow.Table(psms_path)
        if isempty(psms[:ms_file_idx])
            continue
        end

        file_idx = Int(first(psms[:ms_file_idx]))
        if !haskey(rt_irt, file_idx)
            @warn "No RT model found for file index $file_idx, skipping"
            continue
        end

        precursor_idxs = psms[:precursor_idx]
        prob_col = psms[:prob]
        q_values = psms[:q_value]
        rt_vals = psms[:rt]
        scan_idxs = psms[:scan_idx]
        iso_available = hasproperty(psms, :isotopes_captured)
        iso_vals = iso_available ? psms[:isotopes_captured] : nothing

        push!(n_precursors_vec, UInt64(length(precursor_idxs)))

        trace_best = Dict{Tuple{UInt32, Tuple{Int8, Int8}}, Float32}()
        best_scan_idx = Dict{UInt32, UInt32}()
        best_psm_prob = Dict{UInt32, Float32}()
        best_irt = Dict{UInt32, Float32}()
        irt_sum = Dict{UInt32, Float32}()
        irt_count = Dict{UInt32, UInt16}()

        rt_model = rt_irt[file_idx]

        @inbounds for row in eachindex(precursor_idxs)
            prec_idx = precursor_idxs[row]
            prob = Float32(prob_col[row])
            iso_key = iso_available ? normalize_isotopes(iso_vals[row]) : (Int8(-1), Int8(-1))
            combo_key = (prec_idx, iso_key)
            if prob > get(trace_best, combo_key, 0f0)
                trace_best[combo_key] = prob
            end

            if prob > get(best_psm_prob, prec_idx, -Inf32)
                best_psm_prob[prec_idx] = prob
                best_scan_idx[prec_idx] = scan_idxs[row]
                best_irt[prec_idx] = Float32(rt_model(rt_vals[row]))
            end

            if q_values[row] <= max_q_val
                irt_val = Float32(rt_model(rt_vals[row]))
                irt_sum[prec_idx] = get(irt_sum, prec_idx, 0f0) + irt_val
                irt_count[prec_idx] = get(irt_count, prec_idx, UInt16(0)) + UInt16(1)
            end
        end

        if isempty(trace_best)
            continue
        end

        precursor_to_trace_probs = Dict{UInt32, Vector{Float32}}()
        for ((prec_idx, _), prob) in trace_best
            push!(get!(precursor_to_trace_probs, prec_idx, Float32[]), prob)
        end

        precursor_combined = Dict{UInt32, Float32}()
        for (prec_idx, prob_vec) in precursor_to_trace_probs
            precursor_combined[prec_idx] = combine_probability_evidence(prob_vec)
        end

        precursor_to_peptide_key = Dict{UInt32, Tuple{Bool, String, String}}()
        peptide_to_precursor_probs = Dict{Tuple{Bool, String, String}, Vector{Float32}}()
        for (prec_idx, prec_prob) in precursor_combined
            key = peptide_key_for_precursor(prec_idx)
            precursor_to_peptide_key[prec_idx] = key
            push!(get!(peptide_to_precursor_probs, key, Float32[]), prec_prob)
        end

        peptide_probs = Dict{Tuple{Bool, String, String}, Float32}()
        for (key, prob_vec) in peptide_to_precursor_probs
            peptide_probs[key] = combine_probability_evidence(prob_vec)
        end

        final_precursor_probs = Dict{UInt32, Float32}()
        for (prec_idx, prec_prob) in precursor_combined
            key = precursor_to_peptide_key[prec_idx]
            peptide_prob = peptide_probs[key]
            final_precursor_probs[prec_idx] = min(prec_prob, peptide_prob)
        end

        prec_ids = collect(keys(final_precursor_probs))
        if isempty(prec_ids)
            continue
        end

        scores = Float32[final_precursor_probs[pid] for pid in prec_ids]
        labels = Bool[!is_decoy[Int(pid)] for pid in prec_ids]
        peps = Vector{Float32}(undef, length(prec_ids))
        get_PEP!(scores, labels, peps; fdr_scale_factor=fdr_scale_factor)

        passing_precursors = Set{UInt32}()
        for (i, prec_idx) in enumerate(prec_ids)
            if peps[i] <= pep_threshold
                push!(passing_precursors, prec_idx)
            end
        end

        if isempty(passing_precursors)
            continue
        end

        for prec_idx in passing_precursors
            final_prob = final_precursor_probs[prec_idx]
            best_scan = get(best_scan_idx, prec_idx, UInt32(0))
            best_irt_val = get(best_irt, prec_idx, 0f0)
            sum_irt = get(irt_sum, prec_idx, 0f0)
            count = get(irt_count, prec_idx, UInt16(0))
            mz_value = Float32(coalesce(prec_mzs[Int(prec_idx)], 0f0))

            if haskey(prec_to_best_prob, prec_idx)
                old = prec_to_best_prob[prec_idx]
                best_prob_old = old.best_prob
                best_scan_old = old.best_scan_idx
                best_file_old = old.best_ms_file_idx
                best_irt_old = old.best_irt
                if final_prob > best_prob_old
                    best_prob_old = final_prob
                    best_scan_old = best_scan
                    best_file_old = UInt32(file_idx)
                    best_irt_old = best_irt_val
                end
                mean_sum_old = Float32(something(old.mean_irt, 0f0))
                n_old = UInt16(something(old.n, UInt16(0)))
                var_old = Float32(something(old.var_irt, 0f0))
                mean_sum_new = mean_sum_old + sum_irt
                n_new = n_old + count
                new_entry = (
                    best_prob = best_prob_old,
                    best_ms_file_idx = best_file_old,
                    best_scan_idx = best_scan_old,
                    best_irt = best_irt_old,
                    mean_irt = mean_sum_new,
                    var_irt = var_old,
                    n = n_new,
                    mz = old.mz,
                )
                prec_to_best_prob[prec_idx] = new_entry
            else
                new_entry = (
                    best_prob = final_prob,
                    best_ms_file_idx = UInt32(file_idx),
                    best_scan_idx = best_scan,
                    best_irt = best_irt_val,
                    mean_irt = sum_irt,
                    var_irt = 0f0,
                    n = count,
                    mz = mz_value,
                )
                insert!(prec_to_best_prob, prec_idx, new_entry)
            end
        end
    end

    if isempty(n_precursors_vec)
        @warn "No valid PSM files found for cross-run analysis"
        return prec_to_best_prob
    end

    max_precursors = maximum(n_precursors_vec)
    if max_precursors > 0
        sort!(prec_to_best_prob, by = x -> x[:best_prob], alg = PartialQuickSort(1:max_precursors), rev = true)
        N = 0
        for key in collect(keys(prec_to_best_prob))
            N += 1
            if N > max_precursors
                delete!(prec_to_best_prob, key)
            end
        end
    end

    if isempty(prec_to_best_prob)
        return prec_to_best_prob
    end

    for psms_path in psms_paths
        psms = Arrow.Table(psms_path)
        if isempty(psms[:ms_file_idx])
            continue
        end

        file_idx = Int(first(psms[:ms_file_idx]))
        if !haskey(rt_irt, file_idx)
            continue
        end

        precursor_idxs = psms[:precursor_idx]
        q_values = psms[:q_value]
        rt_vals = psms[:rt]
        rt_model = rt_irt[file_idx]

        @inbounds for row in eachindex(precursor_idxs)
            if q_values[row] > max_q_val
                continue
            end

            prec_idx = precursor_idxs[row]
            if !haskey(prec_to_best_prob, prec_idx)
                continue
            end

            entry = prec_to_best_prob[prec_idx]
            n = UInt16(something(entry.n, UInt16(0)))
            if n <= 1
                continue
            end

            mean_sum = Float32(something(entry.mean_irt, 0f0))
            mean_val = mean_sum / Float32(n)
            current_var = Float32(something(entry.var_irt, 0f0))
            irt_val = Float32(rt_model(rt_vals[row]))

            updated_entry = (
                best_prob = entry.best_prob,
                best_ms_file_idx = entry.best_ms_file_idx,
                best_scan_idx = entry.best_scan_idx,
                best_irt = entry.best_irt,
                mean_irt = entry.mean_irt,
                var_irt = current_var + (irt_val - mean_val)^2,
                n = entry.n,
                mz = entry.mz,
            )
            prec_to_best_prob[prec_idx] = updated_entry
        end
    end

    return prec_to_best_prob
end
