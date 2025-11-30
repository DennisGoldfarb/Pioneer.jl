# Copyright (C) 2025 Nathan Wamsley
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
    parse_structural_modifications(structural_mods::Union{Missing, AbstractString})

Convert structural modification string into a mapping from 1-based positions to
modification names. Returns an empty dictionary when no modifications are
present. Terminal modifications (n/c) are ignored because the refinement model
is amino-acid-centric.
"""
function parse_structural_modifications(structural_mods::Union{Missing, AbstractString})
    if ismissing(structural_mods) || isempty(structural_mods)
        return Dict{Int, String}()
    end

    mod_regex = r"\((\d+),([A-Z]|[nc]),([^,\)]+)\)"
    mods_by_position = Dict{Int, String}()

    for m in eachmatch(mod_regex, structural_mods)
        position = parse(Int, m.captures[1])
        aa = first(m.captures[2])
        if aa == 'n' || aa == 'c'
            continue
        end
        mods_by_position[position] = m.captures[3]
    end

    return mods_by_position
end

"""
    sanitize_mod_name(mod_name::String)

Sanitize a modification name so it can be embedded in a symbol-friendly feature
name. Non-alphanumeric characters are replaced with underscores and trimmed.
"""
function sanitize_mod_name(mod_name::String)
    sanitized = replace(mod_name, r"[^A-Za-z0-9]+" => "_")
    sanitized = strip(sanitized, '_')
    return isempty(sanitized) ? "mod" : sanitized
end

"""
    collect_feature_keys(sequences, structural_mods)

Identify all unique (amino acid, modification) combinations present in the
input sequences. Unmodified residues are labeled "unmodified" when present.
"""
function collect_feature_keys(
    sequences::Vector{String},
    structural_mods::AbstractVector
)::Vector{Tuple{Char, String}}
    feature_keys = Set{Tuple{Char, String}}()

    for (seq, mods_str) in zip(sequences, structural_mods)
        mods = parse_structural_modifications(mods_str)
        for (pos, aa) in enumerate(seq)
            if aa ∉ STANDARD_AAS
                continue
            end
            mod_name = get(mods, pos, "unmodified")
            push!(feature_keys, (aa, mod_name))
        end
    end

    sorted_keys = collect(feature_keys)
    sort!(sorted_keys, by = x -> (x[1], x[2]))
    return sorted_keys
end

"""
    prepare_features_dataframe(sequences::Vector{String},
                                structural_mods::AbstractVector,
                                library_irt::Vector{Float32},
                                irt_errors::Vector{Float32})
        -> Tuple{DataFrame, Vector{Symbol}, Vector{Tuple{Char, String}}}

Create feature matrix for iRT refinement model training with separate counts
for every observed amino acid and modification pairing. The function returns the
feature DataFrame along with the feature term symbols and their corresponding
feature keys so coefficients can be mapped back to specific modified residues.

# Arguments
- `sequences`: Peptide sequences
- `structural_mods`: Structural modification annotations aligned with sequences
- `library_irt`: Library iRT predictions
- `irt_errors`: Observed errors (library - observed)

# Returns
`(DataFrame, feature_terms, feature_keys)`
"""
function prepare_features_dataframe(
    sequences::Vector{String},
    structural_mods::AbstractVector,
    library_irt::Vector{Float32},
    irt_errors::Vector{Float32}
)::Tuple{DataFrame, Vector{Symbol}, Vector{Tuple{Char, String}}}
    n = length(sequences)

    feature_keys = collect_feature_keys(sequences, structural_mods)
    feature_symbols = Dict{Tuple{Char, String}, Symbol}(
        key => Symbol("count_$(key[1])_$(sanitize_mod_name(key[2]))") for key in feature_keys
    )

    features = Dict{Symbol, Vector{Float64}}()
    for sym in values(feature_symbols)
        features[sym] = zeros(Float64, n)
    end

    for (i, (seq, mods_str)) in enumerate(zip(sequences, structural_mods))
        mods = parse_structural_modifications(mods_str)
        for (pos, aa) in enumerate(seq)
            if aa ∉ STANDARD_AAS
                continue
            end
            mod_name = get(mods, pos, "unmodified")
            key = (aa, mod_name)
            if haskey(feature_symbols, key)
                features[feature_symbols[key]][i] += 1.0
            end
        end
    end

    features[:error] = Float64.(irt_errors)

    feature_terms = [feature_symbols[key] for key in feature_keys]
    return DataFrame(features), feature_terms, feature_keys
end

"""
    fit_irt_refinement_model(sequences::Vector{String},
                              irt_predicted::Vector{Float32},
                              irt_observed::Vector{Float32};
                              ms_file_idx::Int=1,
                              min_psms::Int=20,
                              train_fraction::Float64=0.67)
        -> Union{IrtRefinementModel, Nothing}

Train linear regression model to predict library iRT prediction errors.

# Workflow
1. Filter to sequences with sufficient PSMs (min_psms)
2. Split into training (train_fraction) and validation sets
3. Train linear model: error ~ counts for each observed amino-acid/
   modification combination
4. Evaluate on validation set
5. If validation MAE improves, retrain on full dataset
6. Return model if refinement helps, nothing otherwise

# Arguments
- `sequences`: Peptide sequences
- `structural_mods`: Structural modification annotations aligned with sequences
- `irt_predicted`: Library iRT predictions
- `irt_observed`: Observed iRT from RT alignment
- `ms_file_idx`: File index for logging
- `min_psms`: Minimum PSMs required (default: 20)
- `train_fraction`: Training set fraction (default: 0.67)

# Returns
`IrtRefinementModel` if refinement improves MAE, `nothing` otherwise

# Model Details
- Features: Counts for every observed amino acid + modification pairing
- Response: error = irt_predicted - irt_observed
- Algorithm: Ordinary least squares (GLM.lm)
- Validation: MAE on held-out validation set
"""
function fit_irt_refinement_model(
    sequences::Vector{String},
    structural_mods::AbstractVector,
    irt_predicted::Vector{Float32},
    irt_observed::Vector{Float32};
    ms_file_idx::Int=1,
    min_psms::Int=20,
    train_fraction::Float64=0.67
)::Union{IrtRefinementModel, Nothing}

    n = length(sequences)

    # Check minimum data requirement
    if n < min_psms
        @debug_l1 "File $ms_file_idx: Insufficient PSMs ($n < $min_psms), skipping iRT refinement"
        return nothing
    end

    # Calculate errors (library - observed; positive when library overestimates)
    irt_errors = irt_predicted .- irt_observed

    # Prepare feature matrix
    features_df, feature_terms, feature_keys = prepare_features_dataframe(
        sequences,
        structural_mods,
        irt_predicted,
        irt_errors
    )

    # Train/validation split
    n_train = round(Int, n * train_fraction)
    n_val = n - n_train

    if n_val < 5
        @debug_l1 "File $ms_file_idx: Insufficient validation data ($n_val < 5), skipping iRT refinement"
        return nothing
    end

    # Random shuffle
    indices = randperm(n)
    train_idx = indices[1:n_train]
    val_idx = indices[(n_train+1):end]

    train_df = features_df[train_idx, :]
    val_df = features_df[val_idx, :]

    # Build formula (modified amino acid counts only)
    formula_str = "error ~ 1"
    if !isempty(feature_terms)
        formula_str *= " + " * join(string.(feature_terms), " + ")
    end
    formula = @eval @formula($(Meta.parse(formula_str)))

    # Train model
    model = lm(formula, train_df)
    coef_values = coef(model)

    intercept = Float32(coef_values[1])
    r2_train = Float32(r2(model))

    # Validate
    feature_matrix = isempty(feature_terms) ? zeros(Float64, n_val, 0) : Matrix(val_df[:, feature_terms])
    val_matrix = hcat(ones(n_val), feature_matrix)
    val_predictions = val_matrix * coef_values

    val_errors = val_df.error
    ss_res = sum((val_errors .- val_predictions).^2)
    ss_tot = sum((val_errors .- mean(val_errors)).^2)
    r2_val = Float32(1 - ss_res / ss_tot)

    # Calculate MAEs
    mae_original = Float32(mean(abs.(val_errors)))
    mae_refined = Float32(mean(abs.(val_errors .- val_predictions)))

    # Decision: use refinement if MAE improves
    use_refinement = mae_refined < mae_original
    mae_improvement = mae_original - mae_refined
    mae_improvement_pct = (mae_improvement / mae_original) * 100.0

    if use_refinement
        @user_info "File $ms_file_idx iRT Refinement ENABLED: " *
                  "Training R²=$(round(r2_train, digits=4)), " *
                  "Validation R²=$(round(r2_val, digits=4)), " *
                  "MAE: $(round(mae_original, digits=4)) → $(round(mae_refined, digits=4)) " *
                  "(Δ=$(round(mae_improvement, digits=4)), $(round(mae_improvement_pct, digits=2))% improvement)"
    else
        @user_info "File $ms_file_idx iRT Refinement DISABLED: " *
                  "MAE: $(round(mae_original, digits=4)) → $(round(mae_refined, digits=4)) " *
                  "(Δ=$(round(mae_improvement, digits=4)), $(round(mae_improvement_pct, digits=2))%) " *
                  "- No improvement, using library iRT"
    end

    # Retrain on full dataset if refinement helps
    if use_refinement
        @debug_l1 "File $ms_file_idx: Retraining on full dataset"

        final_model = lm(formula, features_df)
        final_coef_values = coef(final_model)

        final_intercept = Float32(final_coef_values[1])

        feature_weights = Dict{Tuple{Char, String}, Float32}()
        for (i, key) in enumerate(feature_keys)
            feature_weights[key] = Float32(final_coef_values[1 + i])
        end

        return IrtRefinementModel(
            true,
            feature_weights,
            final_intercept,
            mae_original,
            mae_refined,
            r2_train,
            r2_val
        )
    else
        return nothing
    end
end

"""
    add_refined_irt_column!(psms_path::String,
                            refinement_model::Union{IrtRefinementModel, Nothing},
                            search_context::SearchContext;
                            batch_size::Int=100_000)

Add :refined_irt column to PSMs Arrow file using refinement model.

Uses batch processing for memory efficiency.

# Arguments
- `psms_path`: Path to PSMs Arrow file
- `refinement_model`: iRT refinement model (if nothing, copies :irt_predicted)
- `search_context`: SearchContext for accessing spectral library
- `batch_size`: Rows per batch (default 100k)

# Details
- If model exists: applies refinement to each sequence using its structural
  modifications
- If model is nothing: copies :irt_predicted to :refined_irt
- Uses ColumnOperations.add_column_to_file! for streaming

# Print Statements
Adds @user_info statements to confirm operation
"""
function add_refined_irt_column!(
    psms_path::String,
    refinement_model::Union{IrtRefinementModel, Nothing},
    search_context::SearchContext;
    batch_size::Int=100_000
)
    # Create FileReference
    ref = create_reference(psms_path, PSMFileReference)

    # Get sequences and structural modifications from library
    precursors = getPrecursors(getSpecLib(search_context))
    sequences = getSequence(precursors)
    structural_mods = getStructuralMods(precursors)

    # Define compute function
    compute_fn = if !isnothing(refinement_model) && refinement_model.use_refinement
        @user_info "  Adding :refined_irt column using refinement model..."

        (df_batch::DataFrame) -> begin
            refined = Vector{Float32}(undef, nrow(df_batch))
            for i in 1:nrow(df_batch)
                row = df_batch[i, :]
                seq = sequences[row.precursor_idx]
                mods = structural_mods[row.precursor_idx]
                lib_irt = row.irt_predicted
                refined[i] = refinement_model(seq, mods, lib_irt)
            end
            return refined
        end
    else
        @user_info "  Adding :refined_irt column (no refinement, copying :irt_predicted)..."

        (df_batch::DataFrame) -> Float32.(df_batch.irt_predicted)
    end

    # Use ColumnOperations infrastructure
    add_column_to_file!(ref, :refined_irt, compute_fn; batch_size=batch_size)
    @user_info "  ✓ Successfully added :refined_irt column to PSMs file"

    return nothing
end
