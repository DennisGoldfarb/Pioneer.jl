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
    prepare_features_dataframe(sequences::Vector{String},
                                library_irt::Vector{Float32},
                                irt_errors::Vector{Float32}) -> DataFrame

Create feature matrix for iRT refinement model training.

# Features (21 total)
- 20 columns: count_A, count_R, ..., count_V (amino acid counts)
- 1 column: irt_predicted (library iRT value)
- 1 response: error (library_irt - observed_irt)

# Arguments
- `sequences`: Peptide sequences
- `library_irt`: Library iRT predictions
- `irt_errors`: Observed errors (library - observed)

# Returns
DataFrame with 21 feature columns + error response
"""
function prepare_features_dataframe(
    sequences::Vector{String},
    library_irt::Vector{Float32},
    irt_errors::Vector{Float32}
)::DataFrame
    n = length(sequences)

    # Initialize feature matrix: 20 AAs + library_irt (use Float64 for GLM compatibility)
    features = Dict{Symbol, Vector{Float64}}()

    # Create columns for each standard amino acid
    for aa in STANDARD_AAS
        features[Symbol("count_$(aa)")] = zeros(Float64, n)
    end

    # Count amino acids in each sequence
    for (i, seq) in enumerate(sequences)
        for aa in seq
            col_name = Symbol("count_$(aa)")
            if haskey(features, col_name)
                features[col_name][i] += 1.0
            end
        end
    end

    # Add library iRT and error columns (convert to Float64 for GLM compatibility)
    features[:irt_predicted] = Float64.(library_irt)
    features[:error] = Float64.(irt_errors)

    return DataFrame(features)
end

"""
    fit_irt_refinement_model(sequences::Vector{String},
                              irt_predicted::Vector{Float32},
                              irt_observed::Vector{Float32};
                              ms_file_idx::Int=1,
                              min_psms::Int=20)
        -> Union{IrtRefinementModel, Nothing}

Train linear regression model to predict library iRT prediction errors using all
high-confidence target precursors (no cross validation).

# Workflow
1. Filter to sequences with sufficient PSMs (min_psms)
2. Train linear model: error ~ irt_predicted + count_A + ... + count_V
3. Evaluate MAE on the full training set
4. Return model if refinement helps, nothing otherwise

# Arguments
- `sequences`: Peptide sequences
- `irt_predicted`: Library iRT predictions
- `irt_observed`: Observed iRT from RT alignment
- `ms_file_idx`: File index for logging
- `min_psms`: Minimum PSMs required (default: 20)

# Returns
`IrtRefinementModel` if refinement improves MAE, `nothing` otherwise

# Model Details
- Features: 20 AA counts + library_irt (21 total)
- Response: error = irt_predicted - irt_observed
- Algorithm: Ordinary least squares (GLM.lm)
- Evaluation: MAE computed on the full filtered dataset
"""
function fit_irt_refinement_model(
    sequences::Vector{String},
    irt_predicted::Vector{Float32},
    irt_observed::Vector{Float32};
    ms_file_idx::Int=1,
    min_psms::Int=20
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
    features_df = prepare_features_dataframe(sequences, irt_predicted, irt_errors)

    # Build formula (20 AAs + irt_predicted)
    aa_terms = [Symbol("count_$(aa)") for aa in STANDARD_AAS]
    formula_str = "error ~ irt_predicted + " * join(string.(aa_terms), " + ")
    formula = @eval @formula($(Meta.parse(formula_str)))

    # Train model
    model = lm(formula, features_df)
    coef_values = coef(model)

    intercept = Float32(coef_values[1])
    irt_coef = Float32(coef_values[2])
    aa_coefficients = Float32.(coef_values[3:end])
    r2_train = Float32(r2(model))

    # Evaluate on full dataset
    full_matrix = hcat(ones(n), features_df[:, :irt_predicted], Matrix(features_df[:, aa_terms]))
    full_predictions = full_matrix * coef_values

    mae_original = Float32(mean(abs.(features_df.error)))
    mae_refined = Float32(mean(abs.(features_df.error .- full_predictions)))

    # Decision: use refinement if MAE improves
    use_refinement = mae_refined < mae_original
    mae_improvement = mae_original - mae_refined
    mae_improvement_pct = (mae_improvement / mae_original) * 100.0

    if use_refinement
        @user_info "File $ms_file_idx iRT Refinement ENABLED: " *
                  "Training R²=$(round(r2_train, digits=4)), " *
                  "MAE: $(round(mae_original, digits=4)) → $(round(mae_refined, digits=4)) " *
                  "(Δ=$(round(mae_improvement, digits=4)), $(round(mae_improvement_pct, digits=2))% improvement)"
    else
        @user_info "File $ms_file_idx iRT Refinement DISABLED: " *
                  "MAE: $(round(mae_original, digits=4)) → $(round(mae_refined, digits=4)) " *
                  "(Δ=$(round(mae_improvement, digits=4)), $(round(mae_improvement_pct, digits=2))%) " *
                  "- No improvement, using library iRT"
    end

    if use_refinement
        # Create Dict for callable model
        aa_weights = Dict{Char, Float32}()
        for (i, aa) in enumerate(STANDARD_AAS)
            aa_weights[aa] = Float32(coef_values[2 + i])
        end

        return IrtRefinementModel(
            true,
            aa_weights,
            intercept,
            irt_coef,
            mae_original,
            mae_refined,
            r2_train,
            r2_train
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
- If model exists: applies refinement to each sequence
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

    # Get sequences from library
    precursors = getPrecursors(getSpecLib(search_context))
    sequences = getSequence(precursors)

    # Define compute function
    compute_fn = if !isnothing(refinement_model) && refinement_model.use_refinement
        @user_info "  Adding :refined_irt column using refinement model..."

        (df_batch::DataFrame) -> begin
            refined = Vector{Float32}(undef, nrow(df_batch))
            for i in 1:nrow(df_batch)
                row = df_batch[i, :]
                seq = sequences[row.precursor_idx]
                lib_irt = row.irt_predicted
                refined[i] = refinement_model(seq, lib_irt)
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
