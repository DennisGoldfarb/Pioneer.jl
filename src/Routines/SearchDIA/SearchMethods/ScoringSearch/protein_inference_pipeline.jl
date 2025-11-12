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

#==========================================================
Pipeline Operations for Protein Inference
==========================================================#

"""
    add_peptide_metadata(precursors::LibraryPrecursors)

Add peptide sequence and protein information from precursor library.
Enriches PSM data with metadata needed for inference.
"""
function add_peptide_metadata(precursors::LibraryPrecursors)
    desc = "add_peptide_metadata"
    
    op = function(df)
        # Get library data
        all_sequences = getSequence(precursors)::AbstractVector{String}
        all_accessions = getAccessionNumbers(precursors)::AbstractVector{String}
        all_is_decoys = getIsDecoy(precursors)::AbstractVector{Bool}
        all_entrap_ids = getEntrapmentGroupId(precursors)::AbstractVector{UInt8}
        
        # Extract precursor_idx column with type assertion for performance
        precursor_idx = df.precursor_idx::AbstractVector{UInt32}
        n_rows = length(precursor_idx)
        
        # Add columns if they don't exist
        if !hasproperty(df, :sequence)
            sequences = Vector{String}(undef, n_rows)
            for i in 1:n_rows
                sequences[i] = all_sequences[precursor_idx[i]]
            end
            df.sequence = sequences
        end
        
        if !hasproperty(df, :accession_numbers)
            accessions = Vector{String}(undef, n_rows)
            for i in 1:n_rows
                accessions[i] = all_accessions[precursor_idx[i]]
            end
            df.accession_numbers = accessions
        end
        
        # Add target/decoy status if needed
        if !hasproperty(df, :is_decoy)
            is_decoy_vec = Vector{Bool}(undef, n_rows)
            for i in 1:n_rows
                is_decoy_vec[i] = all_is_decoys[precursor_idx[i]]
            end
            df.is_decoy = is_decoy_vec
        end
        
        # Add entrapment group if needed
        if !hasproperty(df, :entrap_id)
            entrap_vec = Vector{UInt8}(undef, n_rows)
            for i in 1:n_rows
                entrap_vec[i] = all_entrap_ids[precursor_idx[i]]
            end
            df.entrap_id = entrap_vec
        end
        
        return df
    end
    
    return desc => op
end

"""
    validate_peptide_data()

Ensure all required columns exist for protein inference.
"""
function validate_peptide_data()
    desc = "validate_peptide_data"
    
    op = function(df)
        required_cols = [:sequence, :accession_numbers, :is_decoy, 
                        :entrap_id, :precursor_idx]
        
        missing_cols = setdiff(required_cols, Symbol.(names(df)))
        if !isempty(missing_cols)
            error("Missing required columns for protein inference: $missing_cols")
        end
        
        return df
    end
    
    return desc => op
end

#==========================================================
Inference Wrapper Functions
==========================================================#

"""
    apply_inference_to_dataframe(df::DataFrame, precursors::LibraryPrecursors)

Apply the core infer_proteins algorithm to a prepared DataFrame.
Returns structured inference results using type-safe protein inference.
"""
function apply_inference_to_dataframe(df::DataFrame, precursors::LibraryPrecursors)
    if nrow(df) == 0
        # Return empty result for empty input
        return InferenceResult(
            Dictionary{PeptideKey, ProteinKey}()
        )
    end
    
    # Extract unique peptide-protein pairs
    unique_pairs = unique(df, [:sequence, :accession_numbers, :is_decoy, :entrap_id])
    
    # Convert to format expected by infer_proteins
    proteins_vec = Vector{ProteinKey}(undef, nrow(unique_pairs))
    peptides_vec = Vector{PeptideKey}(undef, nrow(unique_pairs))
    
    for (i, row) in enumerate(eachrow(unique_pairs))
        proteins_vec[i] = ProteinKey(
            row.accession_numbers,
            !row.is_decoy,  # Convert is_decoy to is_target
            row.entrap_id
        )
        
        peptides_vec[i] = PeptideKey(
            row.sequence,
            !row.is_decoy,  # Convert is_decoy to is_target
            row.entrap_id
        )
    end
    
    # Call inference algorithm
    return infer_proteins(proteins_vec, peptides_vec)
end

#==========================================================
Post-Inference Pipeline Operations
==========================================================#

"""
    add_inferred_protein_column(inference_result::InferenceResult)

Add inferred protein group assignments to PSMs based on inference results.
"""
function add_inferred_protein_column(inference_result::InferenceResult)
    desc = "add_inferred_protein_column"
    
    op = function(df)
        # Extract columns with type assertions for performance
        sequences = df.sequence::AbstractVector{String}
        is_decoy = df.is_decoy::AbstractVector{Bool}
        entrap_ids = df.entrap_id::AbstractVector{UInt8}
        accession_numbers = df.accession_numbers::AbstractVector{String}
        
        n_rows = length(sequences)
        inferred_proteins = Vector{Union{Missing, String}}(undef, n_rows)
        
        for i in 1:n_rows
            pep_key = PeptideKey(
                sequences[i],
                !is_decoy[i],
                entrap_ids[i]
            )
            
            if haskey(inference_result.peptide_to_protein, pep_key)
                protein_key = inference_result.peptide_to_protein[pep_key]
                inferred_proteins[i] = protein_key.name
            else
                # Fallback to original protein
                inferred_proteins[i] = missing#accession_numbers[i]
            end
        end
        
        df.inferred_protein_group = inferred_proteins
        return df
    end
    
    return desc => op
end

"""
    add_quantification_flag(inference_result::InferenceResult)

Add flag indicating whether peptide should be used for quantification.

Peptides present in the inference result (unique peptides assigned to minimal protein set)
are marked as usable for quantification. Peptides not in the result (shared peptides or
filtered out) are marked as NOT usable for quantification.
"""
function add_quantification_flag(inference_result::InferenceResult)
    desc = "add_quantification_flag"

    op = function(df)
        # Extract columns with type assertions for performance
        sequences = df.sequence::AbstractVector{String}
        is_decoy = df.is_decoy::AbstractVector{Bool}
        entrap_ids = df.entrap_id::AbstractVector{UInt8}

        n_rows = length(sequences)
        use_for_quant = Vector{Bool}(undef, n_rows)

        for i in 1:n_rows
            pep_key = PeptideKey(
                sequences[i],
                !is_decoy[i],
                entrap_ids[i]
            )

            # Peptide is usable for quantification if it's in the inference results
            # (i.e., it's a unique peptide assigned to a protein)
            use_for_quant[i] = haskey(inference_result.peptide_to_protein, pep_key)
        end

        df.use_for_protein_quant = use_for_quant
        return df
    end

    return desc => op
end

"""
    group_psms_by_protein(df::DataFrame)

Transform PSMs into protein groups by aggregating peptides.
Returns a DataFrame with one row per protein group.
"""
const PROTEIN_SCORE_EPS = 1f-6

"""
    combine_probability_evidence(probs::AbstractVector{Float32})

Combine probabilities assuming independent evidence, returning a probability in (0, 1).
Values are clamped slightly below 1.0 to avoid singularities in downstream logarithms.
"""
function combine_probability_evidence(probs::AbstractVector{Float32})
    log_prod = 0.0f0
    has_positive = false
    for prob in probs
        if prob <= 0f0
            continue
        end
        has_positive = true
        clamped = min(prob, 1.0f0 - PROTEIN_SCORE_EPS)
        log_prod += log1p(-clamped)
    end

    if !has_positive
        return 0.0f0
    end

    combined = 1.0f0 - exp(log_prod)
    return clamp(combined, 0.0f0, 1.0f0 - PROTEIN_SCORE_EPS)
end

"""
    accumulate_log_evidence(probs::AbstractVector{Float32})

Accumulate log-evidence style scores from probabilities using `-log(1 - p)`.
Probabilities at or below zero contribute no evidence, while values near one are
clamped to maintain numerical stability.
"""
function accumulate_log_evidence(probs::AbstractVector{Float32})
    total = 0.0f0
    for prob in probs
        if prob <= 0f0
            continue
        end
        clamped = min(prob, 1.0f0 - PROTEIN_SCORE_EPS)
        total += -log1p(-clamped)
    end
    return total
end

function group_psms_by_protein(df::DataFrame)
    if nrow(df) == 0
        # Return empty protein groups DataFrame with expected schema
        return DataFrame(
            protein_name = String[],
            target = Bool[],
            entrap_id = UInt8[],
            n_peptides = Int64[],
            peptide_list = String[],
            pg_score = Float32[],
            any_common_peps = Bool[]
        )
    end

    # Determine which probability column to use for protein scoring
    if hasproperty(df, :MBR_boosted_prec_prob)
        prob_col = :MBR_boosted_prec_prob
    else
        prob_col = :prec_prob
    end

    # Group by protein
    grouped = groupby(df, [:inferred_protein_group, :target, :entrap_id])
    
    # Aggregate to protein groups
    protein_groups = combine(grouped) do gdf
       
        quant_mask = gdf.use_for_protein_quant .== true
        quant_indices = findall(quant_mask)

        if isempty(quant_indices)
            has_common = any((gdf.missed_cleavage .== 0) .& (gdf.Mox .== 0))
            return DataFrame(
                n_peptides = 0,
                peptide_list = "",
                pg_score = 0.0f0,
                any_common_peps = has_common
            )
        end

        sequences = gdf.sequence
        precursor_idx = gdf.precursor_idx
        probabilities = gdf[!, prob_col]
        has_structural_mods = hasproperty(gdf, :structural_mods)
        has_isotopic_mods = hasproperty(gdf, :isotopic_mods)

        structural_mods = has_structural_mods ? gdf.structural_mods : nothing
        isotopic_mods = has_isotopic_mods ? gdf.isotopic_mods : nothing

        seen_peptides = Set{String}()
        peptide_order = String[]
        precursor_best = Dict{UInt32, Float32}()
        precursor_to_mod = Dict{UInt32, Tuple{String, String, String}}()

        for idx in quant_indices
            seq = sequences[idx]
            if !(seq in seen_peptides)
                push!(peptide_order, seq)
                push!(seen_peptides, seq)
            end

            struct_mod = ""
            if has_structural_mods
                val = structural_mods[idx]
                struct_mod = String(something(val, ""))
            end

            iso_mod = ""
            if has_isotopic_mods
                val = isotopic_mods[idx]
                iso_mod = String(something(val, ""))
            end

            key = (seq, struct_mod, iso_mod)
            prec_idx = UInt32(precursor_idx[idx])
            prob = Float32(probabilities[idx])

            if !haskey(precursor_best, prec_idx) || prob > precursor_best[prec_idx]
                precursor_best[prec_idx] = prob
            end
            precursor_to_mod[prec_idx] = key
        end

        mod_to_precursor_probs = Dict{Tuple{String, String, String}, Vector{Float32}}()
        for (prec_idx, prob) in precursor_best
            key = precursor_to_mod[prec_idx]
            push!(get!(mod_to_precursor_probs, key, Float32[]), prob)
        end

        mod_scores = Dict{Tuple{String, String, String}, Float32}()
        for (key, prob_vec) in mod_to_precursor_probs
            mod_scores[key] = combine_probability_evidence(prob_vec)
        end

        peptide_to_mod_probs = Dict{String, Vector{Float32}}()
        for ((seq, _, _), mod_prob) in mod_scores
            push!(get!(peptide_to_mod_probs, seq, Float32[]), mod_prob)
        end

        peptide_scores = Dict{String, Float32}()
        for (seq, mod_probs) in peptide_to_mod_probs
            peptide_scores[seq] = combine_probability_evidence(mod_probs)
        end

        peptide_probs_ordered = Float32[]
        for seq in peptide_order
            if haskey(peptide_scores, seq)
                push!(peptide_probs_ordered, peptide_scores[seq])
            end
        end

        pg_score = accumulate_log_evidence(peptide_probs_ordered)

        has_common = any(
            quant_mask .&
            (gdf.missed_cleavage .== 0) .&
            (gdf.Mox .== 0)
        )
        
        DataFrame(
            n_peptides = length(peptide_scores),
            peptide_list = join(peptide_order, ";"),
            pg_score = pg_score,
            any_common_peps = has_common
        )
    end
    
    # Rename the grouping column
    rename!(protein_groups, :inferred_protein_group => :protein_name)
    
    return protein_groups
end

"""
    filter_by_min_peptides(min_peptides::Int)

Filter protein groups that don't meet minimum peptide requirement.
"""
function filter_by_min_peptides(min_peptides::Int)
    desc = "filter_by_min_peptides(min=$min_peptides)"
    
    op = function(df)
        if hasproperty(df, :n_peptides)
            # For protein groups
            filter!(row -> row.n_peptides >= min_peptides, df)
        end
        return df
    end
    
    return desc => op
end

"""
    add_protein_features(protein_catalog::Dict)

Add protein-level features like peptide coverage.
"""
function add_protein_features(protein_catalog::Dict)
    desc = "add_protein_features"
    
    op = function(df)
        if !hasproperty(df, :protein_name)
            return df
        end
        
        n_rows = nrow(df)
        n_possible = Vector{Int64}(undef, n_rows)
        peptide_coverage = Vector{Float32}(undef, n_rows)
        
        for i in 1:n_rows
            key = (
                protein_name = df.protein_name[i],
                target = df.target[i],
                entrap_id = df.entrap_id[i]
            )
            
            if haskey(protein_catalog, key)
                possible_peptides = protein_catalog[key]
                n_possible[i] = length(possible_peptides)
                
                # Calculate coverage
                if n_possible[i] > 0
                    peptide_coverage[i] = Float32(df.n_peptides[i]) / Float32(n_possible[i])
                else
                    peptide_coverage[i] = 0.0f0
                end
            else
                n_possible[i] = 0
                peptide_coverage[i] = 0.0f0
            end
        end
        
        df.n_possible_peptides = n_possible
        df.peptide_coverage = peptide_coverage
        
        return df
    end
    
    return desc => op
end

#==========================================================
High-Level Interface
==========================================================#

"""
    perform_protein_inference_pipeline(psm_refs, output_folder, precursors, protein_catalog; kwargs...)

Perform protein inference using a composable pipeline approach.

# Arguments
- `psm_refs`: Vector of PSM file references
- `output_folder`: Directory for protein group output
- `precursors`: Library precursor information
- `protein_catalog`: Pre-computed protein-to-peptide mappings
- `min_peptides`: Minimum peptides per protein group (default: 2)

# Returns
- `pg_refs`: Vector of protein group file references
- `psm_to_pg_mapping`: Dictionary mapping PSM paths to protein group paths
"""
function perform_protein_inference_pipeline(
    psm_refs::Vector{PSMFileReference},
    output_folder::String,
    precursors::LibraryPrecursors,
    protein_catalog::Dict;
    min_peptides::Int = 2
)
    # Ensure output folder exists
    !isdir(output_folder) && mkpath(output_folder)
    
    # Build pre-inference pipeline
    pre_inference_pipeline = TransformPipeline() |>
        add_peptide_metadata(precursors) |>
        validate_peptide_data()
    
    # Process each file
    pg_refs = ProteinGroupFileReference[]
    psm_to_pg_mapping = Dict{String, String}()
    
    for (idx, psm_ref) in ProgressBar(collect(enumerate(psm_refs)))
        if !exists(psm_ref)
            continue
        end
        
        # Step 1: Apply pre-inference pipeline to add necessary columns
        apply_pipeline!(psm_ref, pre_inference_pipeline)
        
        # Step 2: Run inference on the prepared PSMs
        prepared_df = load_dataframe(psm_ref)
        inference_result = apply_inference_to_dataframe(prepared_df, precursors)
        
        # Step 3: Update PSMs with inference results
        psm_update_pipeline = TransformPipeline() |>
            add_inferred_protein_column(inference_result) |> #Protein group assigned to each peptide
            add_quantification_flag(inference_result) #Whether or not the peptide should be used for protein quant and inference (non ambiguous)
        
        # Apply updates to the same PSM file (which now has necessary columns)
        apply_pipeline!(psm_ref, psm_update_pipeline)
        
        # Step 4: Create protein groups
        # Reload updated PSMs
        updated_psms = load_dataframe(psm_ref)
        
        # Group by protein
        protein_groups_df = group_psms_by_protein(updated_psms)
        
        # Build post-inference pipeline
        post_inference_pipeline = TransformPipeline() |>
            filter_by_min_peptides(min_peptides) |>
            add_protein_features(protein_catalog)
        
        # Apply post-processing
        initial_rows = nrow(protein_groups_df)
        for (desc, op) in post_inference_pipeline.operations
            protein_groups_df = op(protein_groups_df)
            #@debug_l1 "Pipeline operation on protein groups" operation=desc rows_before=initial_rows rows_after=nrow(protein_groups_df)
            initial_rows = nrow(protein_groups_df)
        end
        
        # Write protein groups
        pg_filename = "protein_groups_$(lpad(idx, 3, '0')).arrow"
        pg_path = joinpath(output_folder, pg_filename)

        if nrow(protein_groups_df) > 0
            # Add file_idx for later merge/split operations
            protein_groups_df[!, :file_idx] = fill(Int64(idx), nrow(protein_groups_df))
            writeArrow(pg_path, protein_groups_df)
            pg_ref = ProteinGroupFileReference(pg_path)
            push!(pg_refs, pg_ref)
            psm_to_pg_mapping[file_path(psm_ref)] = pg_path
        else
            @user_warn "No protein groups to write after filtering" file_idx=idx pg_path=pg_path
        end
    end
    
    return pg_refs, psm_to_pg_mapping
end