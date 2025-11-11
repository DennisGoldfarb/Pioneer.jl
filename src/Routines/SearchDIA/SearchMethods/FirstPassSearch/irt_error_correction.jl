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

const CANONICAL_AMINO_ACIDS = collect("ACDEFGHIKLMNPQRSTVWY")
const CANONICAL_AA_TO_INDEX = Dict{Char, Int}(aa => idx for (idx, aa) in enumerate(CANONICAL_AMINO_ACIDS))
"""
    correct_first_pass_irt_errors!(search_context::SearchContext;
                                   q_value_threshold::Float32 = 0.01f0,
                                   min_samples::Int = length(CANONICAL_AMINO_ACIDS) + 1)

Estimate systematic iRT prediction errors from high-confidence PSMs and
update the observed iRT dictionary stored in the search context.

The routine loads the cached first-pass PSM Arrow tables, filters to target
PSMs with q-values less than or equal to `q_value_threshold`, and computes
signed iRT errors (`predicted - observed`). Sequence features are generated
from peptide sequences by counting the occurrences of canonical amino acids
and distinct modification-specific residues. A linear system is solved to
estimate regression
coefficients that map these counts to iRT errors. If insufficient training
data are available or the design matrix is rank-deficient, the correction is
skipped and library iRT values are retained.

Returns `true` when a regression model is successfully fit; otherwise returns
`false` and leaves `search_context.irt_obs` populated with the original
library iRT values.
"""
function correct_first_pass_irt_errors!(search_context::SearchContext;
                                        q_value_threshold::Float32 = 0.01f0,
                                        min_samples::Int = length(CANONICAL_AMINO_ACIDS) + 1)
    precursors = getPrecursors(getSpecLib(search_context))
    library_irt = Vector{Float32}(getIrt(precursors))
    n_precursors = length(library_irt)
    if n_precursors == 0
        empty!(search_context.irt_obs)
        return false
    end

    is_decoy = Vector{Bool}(getIsDecoy(precursors))
    aa_feature_matrix = compute_precursor_aa_features(precursors)

    psms_paths = getFirstPassPsms(getMSData(search_context))
    training_precursors = Int[]
    training_errors = Float64[]

    for (ms_file_idx, psms_path) in enumerate(psms_paths)
        if isempty(psms_path) || !isfile(psms_path)
            continue
        end

        rt_model = getRtIrtModel(search_context, ms_file_idx)
        psms = try
            Arrow.Table(psms_path)
        catch
            continue
        end

        q_values = psms[:q_value]
        precursor_indices = psms[:precursor_idx]
        predicted_irts = psms[:irt_predicted]
        retention_times = psms[:rt]

        for row_idx in 1:length(q_values)
            if q_values[row_idx] > q_value_threshold
                continue
            end

            precursor_idx = Int(precursor_indices[row_idx])
            if precursor_idx < 1 || precursor_idx > n_precursors
                continue
            end
            if is_decoy[precursor_idx]
                continue
            end

            observed_irt = Float64(rt_model(retention_times[row_idx]))
            predicted_irt = Float64(predicted_irts[row_idx])
            push!(training_precursors, precursor_idx)
            push!(training_errors, predicted_irt - observed_irt)
        end
    end

    success = fit_irt_error_model!(search_context, library_irt, aa_feature_matrix,
                                   training_precursors, training_errors,
                                   min_samples)
    return success
end

function compute_precursor_aa_features(precursors)
    sequences_raw = getSequence(precursors)
    n_precursors = length(sequences_raw)

    tokens_per_precursor = Vector{Vector{String}}(undef, n_precursors)
    extra_tokens = String[]
    seen_extra_tokens = Set{String}()

    for prec_idx in 1:n_precursors
        tokens = parse_sequence_tokens(sequences_raw[prec_idx])
        tokens_per_precursor[prec_idx] = tokens
        for token in tokens
            if is_canonical_token(token)
                continue
            end
            if !in(token, seen_extra_tokens)
                push!(extra_tokens, token)
                push!(seen_extra_tokens, token)
            end
        end
    end

    n_features = length(CANONICAL_AMINO_ACIDS) + length(extra_tokens)
    features = zeros(Float64, n_precursors, n_features)

    extra_token_to_index = Dict{String, Int}()
    for (offset, token) in enumerate(extra_tokens)
        extra_token_to_index[token] = length(CANONICAL_AMINO_ACIDS) + offset
    end

    for prec_idx in 1:n_precursors
        feature_row = @view features[prec_idx, :]
        for token in tokens_per_precursor[prec_idx]
            if length(token) == 1
                aa = token[1]
                feature_idx = get(CANONICAL_AA_TO_INDEX, aa, 0)
                if feature_idx != 0
                    feature_row[feature_idx] += 1.0
                    continue
                end
            end

            feature_idx = get(extra_token_to_index, token, 0)
            if feature_idx != 0
                feature_row[feature_idx] += 1.0
            end
        end
    end

    return features
end

function parse_sequence_tokens(seq::Missing)
    return String[]
end

function parse_sequence_tokens(seq::AbstractString)
    tokens = String[]
    seq_str = String(seq)
    i = firstindex(seq_str)
    last = lastindex(seq_str)

    while i <= last
        ch = seq_str[i]
        if ch == '(' || ch == ')' || ch == ' '
            i = nextind(seq_str, i)
            continue
        end

        base = uppercase(ch)
        next_i = nextind(seq_str, i)
        token = String(base)

        if next_i <= last && seq_str[next_i] == '('
            close_idx = findnext(c -> c == ')', seq_str, next_i)
            if close_idx !== nothing
                inner_start = nextind(seq_str, next_i)
                inner_end = prevind(seq_str, close_idx)
                if inner_start <= inner_end
                    mod_content = uppercase(seq_str[inner_start:inner_end])
                    token = string(base, "(", mod_content, ")")
                    next_i = nextind(seq_str, close_idx)
                end
            end
        end

        push!(tokens, token)
        i = next_i
    end

    return tokens
end

is_canonical_token(token::String) = length(token) == 1 && haskey(CANONICAL_AA_TO_INDEX, token[1])

function fit_irt_error_model!(search_context::SearchContext,
                              library_irt::Vector{Float32},
                              aa_feature_matrix::Matrix{Float64},
                              training_precursors::Vector{Int},
                              training_errors::Vector{Float64},
                              min_samples::Int)
    n_features = size(aa_feature_matrix, 2) + 1
    n_samples = length(training_precursors)

    success = false
    coefficients = zeros(Float64, n_features)

    if n_samples >= max(min_samples, n_features)
        design_matrix = Array{Float64}(undef, n_samples, n_features)
        for (row_idx, precursor_idx) in enumerate(training_precursors)
            design_matrix[row_idx, 1] = 1.0
            @inbounds design_matrix[row_idx, 2:end] .= aa_feature_matrix[precursor_idx, :]
        end

        if LinearAlgebra.rank(design_matrix) == n_features
            coefficients .= design_matrix \ training_errors
            success = true
        end
    end

    empty!(search_context.irt_obs)
    for precursor_idx in 1:length(library_irt)
        base_irt = Float64(library_irt[precursor_idx])
        corrected = base_irt
        if success
            features = aa_feature_matrix[precursor_idx, :]
            predicted_error = coefficients[1] + sum(coefficients[2:end] .* features)
            corrected -= predicted_error
        end
        setPredIrt!(search_context, UInt32(precursor_idx), Float32(corrected))
    end

    return success
end
