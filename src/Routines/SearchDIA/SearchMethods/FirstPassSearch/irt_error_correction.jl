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
const TERMINAL_MOD_BASES = Set(['n', 'c'])
const EMPTY_SEQUENCE_MODIFICATIONS = NamedTuple{(:position, :residue, :mod_name), Tuple{Int, Char, String}}[]

normalize_residue(ch::Char) = uppercase(String(ch))[1]

function normalize_mod_residue(ch::Char)
    return ch in TERMINAL_MOD_BASES ? ch : normalize_residue(ch)
end

const SequenceModification = NamedTuple{(:position, :residue, :mod_name), Tuple{Int, Char, String}}

function parse_structural_mods(::Missing)
    return SequenceModification[]
end

function parse_structural_mods(mods::AbstractString)
    isempty(mods) && return SequenceModification[]
    mod_regex = r"\((\d+),([^,]+),([^\)]+)\)"
    parsed = SequenceModification[]

    for match in eachmatch(mod_regex, mods)
        position = parse(Int, match.captures[1])
        residue_token = strip(match.captures[2])
        residue_char = isempty(residue_token) ? '?' : first(residue_token)
        mod_name = strip(match.captures[3])
        push!(parsed, (position = position,
                       residue = normalize_mod_residue(residue_char),
                       mod_name = mod_name))
    end

    return parsed
end

structural_mod_token(residue::Char, mod_name::String) = string(residue, "(", mod_name, ")")

function determine_static_mods(sequences::Vector{String},
                               positional_mods::Vector{Vector{SequenceModification}},
                               nonpositional_mods::Vector{Vector{SequenceModification}})
    residue_site_totals = Dict{Char, Int}()
    for seq in sequences
        for ch in seq
            residue = normalize_residue(ch)
            residue_site_totals[residue] = get(residue_site_totals, residue, 0) + 1
        end
    end
    residue_site_totals['n'] = max(get(residue_site_totals, 'n', 0), length(sequences))
    residue_site_totals['c'] = max(get(residue_site_totals, 'c', 0), length(sequences))

    mod_counts = Dict{Tuple{Char, String}, Int}()
    for mods in positional_mods
        for mod in mods
            key = (mod.residue, mod.mod_name)
            mod_counts[key] = get(mod_counts, key, 0) + 1
        end
    end
    for mods in nonpositional_mods
        for mod in mods
            key = (mod.residue, mod.mod_name)
            mod_counts[key] = get(mod_counts, key, 0) + 1
        end
    end

    static_mods = Set{Tuple{Char, String}}()
    for (key, count) in mod_counts
        total_sites = get(residue_site_totals, key[1], length(sequences))
        if total_sites != 0 && count == total_sites
            push!(static_mods, key)
        end
    end

    return static_mods
end
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
and distinct modification-specific residues derived from the structural
modification annotations in the spectral library. A linear system is solved to
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
    structural_mods_raw = getStructuralMods(precursors)
    n_precursors = length(sequences_raw)

    sequences = Vector{String}(undef, n_precursors)
    positional_mods = Vector{Vector{SequenceModification}}(undef, n_precursors)
    nonpositional_mods = Vector{Vector{SequenceModification}}(undef, n_precursors)

    for prec_idx in 1:n_precursors
        seq = String(sequences_raw[prec_idx])
        sequences[prec_idx] = seq

        mods = parse_structural_mods(structural_mods_raw[prec_idx])
        positional = SequenceModification[]
        nonpositional = SequenceModification[]
        for mod in mods
            if haskey(CANONICAL_AA_TO_INDEX, mod.residue)
                push!(positional, mod)
            else
                push!(nonpositional, mod)
            end
        end
        positional_mods[prec_idx] = positional
        nonpositional_mods[prec_idx] = nonpositional
    end

    static_mods = determine_static_mods(sequences, positional_mods, nonpositional_mods)

    extra_tokens = String[]
    seen_extra_tokens = Set{String}()
    for mods in positional_mods
        for mod in mods
            key = (mod.residue, mod.mod_name)
            if in(key, static_mods)
                continue
            end
            token = structural_mod_token(mod.residue, mod.mod_name)
            if !in(token, seen_extra_tokens)
                push!(extra_tokens, token)
                push!(seen_extra_tokens, token)
            end
        end
    end
    for mods in nonpositional_mods
        for mod in mods
            key = (mod.residue, mod.mod_name)
            if in(key, static_mods)
                continue
            end
            token = structural_mod_token(mod.residue, mod.mod_name)
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
        seq = sequences[prec_idx]

        mods_by_position = Dict{Int, Vector{SequenceModification}}()
        for mod in positional_mods[prec_idx]
            push!(get!(mods_by_position, mod.position, SequenceModification[]), mod)
        end

        for (pos, ch) in enumerate(seq)
            residue = normalize_residue(ch)
            canonical_idx = get(CANONICAL_AA_TO_INDEX, residue, 0)
            has_variable_mod = false
            mods_here = get(mods_by_position, pos, EMPTY_SEQUENCE_MODIFICATIONS)

            for mod in mods_here
                if in((mod.residue, mod.mod_name), static_mods)
                    continue
                end
                has_variable_mod = true
                token = structural_mod_token(mod.residue, mod.mod_name)
                feature_idx = get(extra_token_to_index, token, 0)
                if feature_idx != 0
                    feature_row[feature_idx] += 1.0
                end
            end

            if !has_variable_mod && canonical_idx != 0
                feature_row[canonical_idx] += 1.0
            end
        end

        for mod in nonpositional_mods[prec_idx]
            if in((mod.residue, mod.mod_name), static_mods)
                continue
            end
            token = structural_mod_token(mod.residue, mod.mod_name)
            feature_idx = get(extra_token_to_index, token, 0)
            if feature_idx != 0
                feature_row[feature_idx] += 1.0
            end
        end
    end

    return features
end

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
