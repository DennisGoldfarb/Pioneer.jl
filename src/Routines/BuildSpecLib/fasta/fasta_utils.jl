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

const TERMINAL_MUTATION_MAP = Dict{Char,Char}(
    'G' => 'L',
    'A' => 'L',
    'V' => 'L',
    'L' => 'V',
    'I' => 'V',
    'F' => 'L',
    'M' => 'L',
    'P' => 'L',
    'W' => 'L',
    'S' => 'T',
    'C' => 'S',
    'T' => 'S',
    'Y' => 'S',
    'H' => 'S',
    'K' => 'L',
    'R' => 'L',
    'Q' => 'N',
    'E' => 'D',
    'N' => 'Q',
    'D' => 'E',
)

"""
    PeptideSequenceSet

A data structure for efficiently storing and comparing peptide sequences with I/L equivalence.

The structure replaces all isoleucine (I) with leucine (L) in stored sequences to treat these
amino acids as equivalent for the purpose of comparisons.

# Fields
- `sequences::Set{String}`: Set of stored sequences with I replaced by L

# Methods
- `PeptideSequenceSet()`: Constructor to create an empty set
- `PeptideSequenceSet(fasta_entries::Vector{FastaEntry})`: Constructor to initialize from FASTA entries
- `push!(pss::PeptideSequenceSet, seq::AbstractString)`: Add a sequence to the set
- `in(seq::AbstractString, pss::PeptideSequenceSet)`: Check if a sequence exists in the set

# Examples
```julia
# Create an empty set
pss = PeptideSequenceSet()

# Add sequences
push!(pss, "PEPTIDE")
push!(pss, "PEPTLDE")  # Contains L instead of I

# Check membership (both return true due to I/L equivalence)
"PEPTIDE" in pss  # true
"PEPTLDE" in pss  # true

# Initialize from FASTA entries
pss = PeptideSequenceSet(fasta_entries)
```
"""
struct PeptideSequenceSet
    sequences::Set{Tuple{String, UInt8}}
    # Constructor
    function PeptideSequenceSet()
        new(Set{Tuple{String, UInt8}}())
    end
    function PeptideSequenceSet(fasta_entries::Vector{FastaEntry})
        pss = PeptideSequenceSet()
        # Add each sequence to the set
        for entry in fasta_entries
            push!(pss, get_sequence(entry), get_charge(entry))
        end
        return pss
    end
end

getSeqSet(s::PeptideSequenceSet) = s.sequences

import Base: push!
function push!(pss::PeptideSequenceSet, seq::AbstractString, charge::UInt8)
    # Replace 'I' with 'L' in the sequence and add to the set
    push!(pss.sequences, (replace(seq, 'I' => 'L'), charge))
    return pss
end

import Base: in
function in(seq_charge::Tuple{String, UInt8}, pss::PeptideSequenceSet)
    # Check if the modified sequence is in the set
    return (replace(first(seq_charge), 'I' => 'L'), last(seq_charge)) ∈ getSeqSet(pss)
end

function collect_modification_flags(
    seq_length::Int,
    structural_mods::Union{Missing, Vector{PeptideMod}},
    isotopic_mods::Union{Missing, Vector{PeptideMod}},
)
    modified_positions = falses(seq_length)
    n_term_modded = Ref(false)
    c_term_modded = Ref(false)

    function mark_mods!(mods)
        if ismissing(mods)
            return
        end
        for mod in mods
            if mod.aa == 'n'
                n_term_modded[] = true
            elseif mod.aa == 'c'
                c_term_modded[] = true
            else
                pos = Int(mod.position)
                if 1 <= pos <= seq_length
                    modified_positions[pos] = true
                end
            end
        end
    end

    mark_mods!(structural_mods)
    mark_mods!(isotopic_mods)

    return modified_positions, n_term_modded[], c_term_modded[]
end

function merge_modification_flags!(
    modified_positions::Vector{Bool},
    n_term_modded::Base.RefValue{Bool},
    c_term_modded::Base.RefValue{Bool},
    structural_mods::Union{Missing, Vector{PeptideMod}},
    isotopic_mods::Union{Missing, Vector{PeptideMod}},
)
    function mark_mods!(mods)
        if ismissing(mods)
            return
        end
        for mod in mods
            if mod.aa == 'n'
                n_term_modded[] = true
            elseif mod.aa == 'c'
                c_term_modded[] = true
            else
                pos = Int(mod.position)
                if 1 <= pos <= length(modified_positions)
                    modified_positions[pos] = true
                end
            end
        end
    end

    mark_mods!(structural_mods)
    mark_mods!(isotopic_mods)
end

function n_term_candidate_positions(
    seq_length::Int,
    modified_positions::Vector{Bool},
    n_term_modded::Bool,
)
    candidates = Int[]
    if seq_length <= 0
        return candidates
    elseif seq_length == 1
        push!(candidates, 1)
        return candidates
    end

    pos2_modded = modified_positions[2]
    pos3_exists = seq_length >= 3
    pos3_modded = pos3_exists ? modified_positions[3] : false
    terminal_modded = modified_positions[1] || n_term_modded

    if !pos2_modded
        push!(candidates, 2)
        if pos3_exists && !pos3_modded
            push!(candidates, 3)
        end
        if !terminal_modded
            push!(candidates, 1)
        end
    else
        if pos3_exists
            push!(candidates, 3)
            if pos3_modded
                if !terminal_modded
                    push!(candidates, 1)
                end
            else
                if !terminal_modded
                    push!(candidates, 1)
                end
            end
        else
            if !terminal_modded
                push!(candidates, 1)
            end
        end
        if terminal_modded
            push!(candidates, 2)
        end
    end

    if !pos2_modded && !terminal_modded
        # allow fallback to position 1 if duplicate occurs
        if !(1 in candidates)
            push!(candidates, 1)
        end
    elseif pos2_modded && !(2 in candidates)
        push!(candidates, 2)
    end

    # Remove duplicates while preserving order and ensure indices in bounds
    unique_candidates = Int[]
    for pos in candidates
        if 1 <= pos <= seq_length && !(pos in unique_candidates)
            push!(unique_candidates, pos)
        end
    end
    return unique_candidates
end

function c_term_candidate_positions(
    seq_length::Int,
    modified_positions::Vector{Bool},
    c_term_modded::Bool,
)
    candidates = Int[]
    if seq_length <= 0
        return candidates
    elseif seq_length == 1
        push!(candidates, 1)
        return candidates
    end

    primary = seq_length - 1
    primary = max(primary, 1)
    before_primary = primary - 1
    before_exists = before_primary >= 1
    before_modded = before_exists ? modified_positions[before_primary] : false
    primary_modded = modified_positions[primary]
    terminal_modded = modified_positions[seq_length] || c_term_modded

    if !primary_modded
        push!(candidates, primary)
        if before_exists && !before_modded
            push!(candidates, before_primary)
        end
        if !terminal_modded
            push!(candidates, seq_length)
        end
    else
        if before_exists
            push!(candidates, before_primary)
            if !before_modded && !terminal_modded
                push!(candidates, seq_length)
            elseif before_modded && !terminal_modded
                push!(candidates, seq_length)
            end
        else
            if !terminal_modded
                push!(candidates, seq_length)
            end
        end
        if terminal_modded
            push!(candidates, primary)
        end
    end

    if !primary_modded && !terminal_modded
        if !(seq_length in candidates)
            push!(candidates, seq_length)
        end
    elseif primary_modded && !(primary in candidates)
        push!(candidates, primary)
    end

    unique_candidates = Int[]
    for pos in candidates
        if 1 <= pos <= seq_length && !(pos in unique_candidates)
            push!(unique_candidates, pos)
        end
    end
    return unique_candidates
end

function attempt_mutation!(
    seq_chars::Vector{Char},
    candidate_pos::Int,
    charges::AbstractVector{UInt8},
    sequences_set::PeptideSequenceSet,
    mutated_positions::Vector{Int},
)
    if any(==(candidate_pos), mutated_positions)
        return false
    end
    if candidate_pos < 1 || candidate_pos > length(seq_chars)
        return false
    end

    original_aa = seq_chars[candidate_pos]
    new_aa = get(TERMINAL_MUTATION_MAP, original_aa, nothing)
    if new_aa === nothing || new_aa == original_aa
        return false
    end

    new_chars = copy(seq_chars)
    new_chars[candidate_pos] = new_aa
    mutated_sequence = String(new_chars)
    for charge in charges
        if (mutated_sequence, charge) ∈ sequences_set
            return false
        end
    end

    seq_chars[candidate_pos] = new_aa
    push!(mutated_positions, candidate_pos)
    return true
end

function mutate_termini!(
    seq_chars::Vector{Char},
    charges::AbstractVector{UInt8},
    sequences_set::PeptideSequenceSet,
    modified_positions::Vector{Bool},
    n_term_modded::Bool,
    c_term_modded::Bool,
)
    mutated_positions = Int[]
    seq_length = length(seq_chars)

    n_candidates = n_term_candidate_positions(seq_length, modified_positions, n_term_modded)
    for pos in n_candidates
        if attempt_mutation!(seq_chars, pos, charges, sequences_set, mutated_positions)
            break
        end
    end

    c_candidates = c_term_candidate_positions(seq_length, modified_positions, c_term_modded)
    for pos in c_candidates
        if attempt_mutation!(seq_chars, pos, charges, sequences_set, mutated_positions)
            break
        end
    end

    return mutated_positions
end

function generate_decoy_sequence(
    sequence::AbstractString,
    charges::AbstractVector{UInt8},
    sequences_set::PeptideSequenceSet,
    modified_positions::Vector{Bool},
    n_term_modded::Bool,
    c_term_modded::Bool,
)
    seq_chars = collect(sequence)
    mutate_termini!(seq_chars, charges, sequences_set, modified_positions, n_term_modded, c_term_modded)
    mutated_sequence = String(seq_chars)
    if mutated_sequence == sequence
        return nothing
    end
    return mutated_sequence
end

"""
    add_entrapment_sequences(
        target_fasta_entries::Vector{FastaEntry}, 
        entrapment_r::UInt8;
        max_shuffle_attempts::Int64 = 20
    )::Vector{FastaEntry}

Add entrapment sequences to a set of target peptides with modification handling.
Creates shuffled sequences while properly adjusting modification positions.

# Parameters
- `target_fasta_entries::Vector{FastaEntry}`: Vector of target peptide entries (with modifications) to generate entrapment sequences for
- `entrapment_r::UInt8`: Number of entrapment sequences to generate per target
- `max_shuffle_attempts::Int64`: Maximum attempts to generate unique shuffled sequence (default: 20)

# Returns
- `Vector{FastaEntry}`: Combined vector of original entries and their entrapment sequences

# Details
For each target peptide:
1. Creates `entrapment_r` shuffled versions using `shuffle_sequence!()`
2. Adjusts modification positions to match shuffled sequence using `adjust_mod_positions`
3. Preserves C-terminal amino acid to maintain enzymatic cleavage properties  
4. Ensures each shuffled sequence is unique (I/L equivalence considered)
5. Sets entrapment_group_id to indicate the entrapment group
6. Maintains original metadata (base_target_id, base_pep_id, etc.) for tracking
7. Properly handles both structural and isotopic modifications

# Examples
```julia
# Create 2 entrapment sequences per target
entries_with_entrapment = add_entrapment_sequences(target_entries, UInt8(2))

# Create 3 entrapment sequences with more shuffle attempts for difficult sequences
entries_with_entrapment = add_entrapment_sequences(
    target_entries, 
    UInt8(3), 
    max_shuffle_attempts=50
)
```

# Notes
Entrapment sequences help assess false discovery rates in peptide identification.
The function uses I/L equivalence when checking for sequence uniqueness.
"""
function add_entrapment_sequences(
    target_fasta_entries::Vector{FastaEntry}, 
    entrapment_r::UInt8;
    max_shuffle_attempts::Int64 = 20,
    fixed_chars::Vector{Char} = Vector{Char}(),
    entrapment_method::String = "shuffle"
)::Vector{FastaEntry}
    
    # Pre-allocate output vector
    entrapment_fasta_entries = Vector{FastaEntry}(
        undef, 
        length(target_fasta_entries) * entrapment_r
    )
    
    # Counters for tracking fallback to shuffle
    total_sequences = length(target_fasta_entries) * entrapment_r
    fallback_to_shuffle_count = 0
    
    # Track unique sequences
    sequences_set = PeptideSequenceSet(target_fasta_entries)#Set{String}()
    #sizehint!(sequences_set, length(entrapment_fasta_entries) + length(target_fasta_entries))
    #union!(sequences_set, target_sequences)
    
    n = 1

    shuffle_seq = ShuffleSeq(
        "",
        Vector{Char}(undef, 255),
        Vector{UInt8}(undef, 255),
        Vector{UInt8}(undef, 255),
        zero(UInt8),
        zero(UInt8),
       fixed_chars#['R','K']#Vector{Char}()
    )
    for target_entry in target_fasta_entries
        for entrapment_group_id in 1:entrapment_r
            n_shuffle_attempts = 0
            # Try the specified method first
            new_sequence = shuffle_sequence!(shuffle_seq, get_sequence(target_entry); method=entrapment_method)
            
            # If it creates a duplicate and we're using reverse, fall back to shuffle
            if (new_sequence, get_charge(target_entry)) ∈ sequences_set && entrapment_method == "reverse"
                @debug_l2 "Reverse created duplicate for entrapment of $(get_sequence(target_entry)), falling back to shuffle"
                fallback_to_shuffle_count += 1
                # Fall back to shuffle since reverse is deterministic
                while n_shuffle_attempts < max_shuffle_attempts
                    new_sequence = shuffle_sequence!(shuffle_seq, get_sequence(target_entry); method="shuffle")
                    
                    if (new_sequence, get_charge(target_entry)) ∉ sequences_set
                        break
                    end
                    n_shuffle_attempts += 1
                end
            elseif (new_sequence, get_charge(target_entry)) ∈ sequences_set
                # For shuffle method, keep trying with shuffle
                while n_shuffle_attempts < max_shuffle_attempts
                    new_sequence = shuffle_sequence!(shuffle_seq, get_sequence(target_entry); method="shuffle")
                    
                    if (new_sequence, get_charge(target_entry)) ∉ sequences_set
                        break
                    end
                    n_shuffle_attempts += 1
                end
            end
            
            #Make sure the entrapment sequence is unique (I and L are equivalent)
            if (new_sequence, get_charge(target_entry)) ∉ sequences_set
                    # Get sequence length for modification adjustment
                    seq_length = UInt8(length(get_sequence(target_entry)))
                    
                    # Adjust modification positions based on sequence shuffling
                    adjusted_structural_mods = adjust_mod_positions(
                        get_structural_mods(target_entry),
                        shuffle_seq.new_positions,
                        seq_length
                    )
                    
                    adjusted_isotopic_mods = adjust_mod_positions(
                        get_isotopic_mods(target_entry),
                        shuffle_seq.new_positions,
                        seq_length
                    )
                    
                    entrapment_fasta_entries[n] = FastaEntry(
                        get_id(target_entry),
                        get_description(target_entry),
                        get_gene(target_entry),
                        get_protein(target_entry),
                        get_organism(target_entry),
                        get_proteome(target_entry),
                        new_sequence,
                        get_start_idx(target_entry),
                        adjusted_structural_mods, #structural_mods - now properly adjusted
                        adjusted_isotopic_mods,   #isotopic_mods - now properly adjusted
                        get_charge(target_entry),
                        get_base_target_id(target_entry), # inherit base_target_id for tracking
                        get_base_pep_id(target_entry),
                        entrapment_group_id,
                        false
                    )
                    n += 1
                    push!(sequences_set, new_sequence, get_charge(target_entry))
            elseif n_shuffle_attempts >= max_shuffle_attempts
                @user_warn "Max shuffle attempts exceeded for $(get_sequence(target_entry))"
            end
        end
    end
    
    # Report statistics if using reverse method for entrapment
    #=
    if entrapment_method == "reverse" && total_sequences > 0
        if fallback_to_shuffle_count > 0
            @user_warn "Entrapment generation statistics for REVERSE method:"
            @user_warn "  Total entrapment sequences attempted: $total_sequences"
            @user_warn "  Sequences where reverse created duplicates: $fallback_to_shuffle_count"
            @user_warn "  Sequences successfully reversed: $(total_sequences - fallback_to_shuffle_count)"
            @user_warn "  Fallback rate: $(round(100.0 * fallback_to_shuffle_count / total_sequences, digits=1))%"
        end
    end
    =#
    
    return vcat(target_fasta_entries, entrapment_fasta_entries[1:n-1])
end

"""
    add_entrapment_sequences_grouped(
        target_fasta_entries::Vector{FastaEntry},
        entrapment_r::UInt8;
        max_shuffle_attempts::Int64 = 20,
        fixed_chars::Vector{Char} = Vector{Char}(),
        entrapment_method::String = "shuffle"
    )::Vector{FastaEntry}

Group-aware entrapment generation that ensures all modification variants of the
same base peptide sequence receive the same set of entrapment sequences.

# Parameters
- `target_fasta_entries::Vector{FastaEntry}`: Vector of peptide entries (after modifications)
- `entrapment_r::UInt8`: Number of entrapment sequences to generate per base sequence
- `max_shuffle_attempts::Int64`: Max attempts to find a unique shuffled sequence
- `fixed_chars::Vector{Char}`: Optional characters to keep fixed when shuffling
- `entrapment_method::String`: "shuffle" or "reverse" (reverse may fall back to shuffle)

# Returns
- `Vector{FastaEntry}`: Combined vector of original entries and grouped entrapment entries

# Details
Algorithm:
1. Group entries by base sequence (ignoring modifications)
2. For each base sequence, generate `entrapment_r` unique entrapment sequences once
3. Reuse those sequences for all modification variants in the group, adjusting mod positions
4. Maintain I/L equivalence when checking uniqueness (via PeptideSequenceSet)
"""
function add_entrapment_sequences_grouped(
    target_fasta_entries::Vector{FastaEntry},
    entrapment_r::UInt8;
    max_shuffle_attempts::Int64 = 20,
    fixed_chars::Vector{Char} = Vector{Char}(),
    entrapment_method::String = "shuffle"
)::Vector{FastaEntry}

    # Track existing sequences (I/L equivalence) including charges
    sequences_set = PeptideSequenceSet(target_fasta_entries)

    # Prepare shuffler
    shuffle_seq = ShuffleSeq(
        "",
        Vector{Char}(undef, 255),
        Vector{UInt8}(undef, 255),
        Vector{UInt8}(undef, 255),
        zero(UInt8),
        zero(UInt8),
        fixed_chars
    )

    # Group entries by base sequence (ignore mods)
    groups = Dict{String, Vector{Int}}()
    for (idx, entry) in enumerate(target_fasta_entries)
        base_seq = get_sequence(entry)
        if !haskey(groups, base_seq)
            groups[base_seq] = Vector{Int}()
        end
        push!(groups[base_seq], idx)
    end

    # High-level diagnostics
    n_entries = length(target_fasta_entries)
    n_groups = length(groups)
    avg_variants = n_groups == 0 ? 0.0 : round(n_entries / n_groups, digits=2)

    entrapments_out = Vector{FastaEntry}()
    fallback_to_shuffle_count = 0
    total_attempted = length(groups) * Int(entrapment_r)

    exhausted_groups = 0
    sample_logged = 0
    for (base_seq, idxs) in groups
        # Collect unique charges observed among variants (usually 0 at this stage)
        charges = unique([get_charge(target_fasta_entries[i]) for i in idxs])

        # Generate unique entrapment sequences for this base sequence
        entrap_seqs = Vector{String}()
        entrap_positions = Vector{Vector{UInt8}}()
        for entrapment_group_id in 1:entrapment_r
            n_shuffle_attempts = 0

            # Start with requested method
            new_sequence = shuffle_sequence!(shuffle_seq, base_seq; method=entrapment_method)

            # If duplicate: reverse may fall back to shuffle; shuffle keeps trying
            needs_retry = any(((new_sequence, c) ∈ sequences_set) for c in charges)
            if needs_retry && entrapment_method == "reverse"
                @user_warn "Reverse duplicate for entrapment of $base_seq; fallback to shuffle"
                fallback_to_shuffle_count += 1
            end

            while needs_retry && n_shuffle_attempts < max_shuffle_attempts
                new_sequence = shuffle_sequence!(shuffle_seq, base_seq; method="shuffle")
                needs_retry = any(((new_sequence, c) ∈ sequences_set) for c in charges)
                n_shuffle_attempts += 1
            end

            if needs_retry
                exhausted_groups += 1
                continue
            end

            # Snapshot positions mapping for consistent mod adjustments across variants
            positions_copy = Vector{UInt8}(shuffle_seq.new_positions)
            push!(entrap_seqs, new_sequence)
            push!(entrap_positions, positions_copy)

            # Reserve the sequence globally for all observed charges
            for c in charges
                push!(sequences_set, new_sequence, c)
            end
        end

        # Optional per-group diagnostics (debug level)
        if sample_logged < 5
            trunc_seq = length(base_seq) > 18 ? (base_seq[1:18] * "…") : base_seq
            sample_logged += 1
        end

        # Create entrapment entries for all variants using the same entrapment specs
        for idx in idxs
            target_entry = target_fasta_entries[idx]
            seq_length = UInt8(length(base_seq))
            for i in 1:length(entrap_seqs)
                adjusted_structural_mods = adjust_mod_positions(
                    get_structural_mods(target_entry),
                    entrap_positions[i],
                    seq_length
                )
                adjusted_isotopic_mods = adjust_mod_positions(
                    get_isotopic_mods(target_entry),
                    entrap_positions[i],
                    seq_length
                )

                push!(entrapments_out, FastaEntry(
                    get_id(target_entry),
                    get_description(target_entry),
                    get_gene(target_entry),
                    get_protein(target_entry),
                    get_organism(target_entry),
                    get_proteome(target_entry),
                    entrap_seqs[i],
                    get_start_idx(target_entry),
                    adjusted_structural_mods,
                    adjusted_isotopic_mods,
                    get_charge(target_entry),
                    get_base_target_id(target_entry),
                    get_base_pep_id(target_entry),
                    UInt8(i),
                    false
                ))
            end
        end
    end

    # Report statistics if using reverse method for entrapment
    if entrapment_method == "reverse" && total_attempted > 0
        if fallback_to_shuffle_count > 0
            @user_warn "Entrapment generation (GROUPED) stats for REVERSE:"
            @user_warn "  Total entrapments attempted: $total_attempted"
            @user_warn "  Reverse duplicates: $fallback_to_shuffle_count"
            @user_warn "  Fallback rate: $(round(100.0 * fallback_to_shuffle_count / total_attempted, digits=1))%"
        end
    end

    # General summary

    return vcat(target_fasta_entries, entrapments_out)
end

mutable struct ShuffleSeq
    old_sequence::String 
    new_sequence::Vector{Char}
    new_positions::Vector{UInt8}
    movable_positions::Vector{UInt8}
    n_movable::UInt8
    sequence_length::UInt8
    fixed_chars::Vector{Char}
end

"""
    resetSequence!(shuffle_sequence::ShuffleSeq, sequence::String)

Resets the 'old_sequence' and 'new_sequence' attributes of the `ShuffleSeq` object.
and resets the 'sequence_length' attribute of the 'shuffle_sequence'
"""
function resetSequence!(shuffle_sequence::ShuffleSeq, sequence::String)

    #Fills the character string for the 'new_sequence' to match 'sequence'
    #and resets the 'sequence_length' attribute of the 'shuffle_sequence'
    shuffle_sequence.old_sequence = sequence 
    shuffle_sequence.sequence_length = length(sequence)
    for i in range(one(UInt8), UInt8(shuffle_sequence.sequence_length))
        shuffle_sequence.new_sequence[i] = sequence[i]
        shuffle_sequence.new_positions[i] = i
    end
    shuffle_sequence.n_movable = zero(UInt8) 
    return nothing
end

"""
    fillMovablePositions!(ss::ShuffleSeq, sequence) 
    
Fills the character string for the 'new_sequence' to match 'sequence'
and resets the 'sequence_length' attribute of the 'shuffle_sequence'
"""
function fillMovablePositions!(shuffle_sequence::ShuffleSeq)

    shuffle_sequence.n_movable = zero(UInt8)
    #shuffle_sequence.sequence_length-1 because the last amino-acid is fixed 
    for i in range(one(UInt8), shuffle_sequence.sequence_length-1)
        if shuffle_sequence.old_sequence[i] ∉ shuffle_sequence.fixed_chars
            shuffle_sequence.n_movable += one(UInt8)
            # If the character is fixed, keep it in the same position
            shuffle_sequence.movable_positions[shuffle_sequence.n_movable] = i
        end
    end
    return nothing
end

function permuteNewPositions!(shuffle_sequence::ShuffleSeq)
    perm = randperm(shuffle_sequence.n_movable)
    # Update new_positions based on the permutation
    for (new_idx, old_idx) in enumerate(perm)
        # Update sequence 
        shuffle_sequence.new_sequence[shuffle_sequence.movable_positions[new_idx]] = 
            shuffle_sequence.old_sequence[shuffle_sequence.movable_positions[old_idx]]
        # Update positions 
        shuffle_sequence.new_positions[shuffle_sequence.movable_positions[new_idx]] = 
            shuffle_sequence.movable_positions[old_idx]
    end
    return nothing
end

function reverseMovablePositions!(shuffle_sequence::ShuffleSeq)
    # Reverse the movable positions (all but the last amino acid)
    n = shuffle_sequence.n_movable
    
    # Create reversed mapping
    for i in 1:n
        old_pos = shuffle_sequence.movable_positions[i]
        new_pos = shuffle_sequence.movable_positions[n - i + 1]
        
        # Update sequence - place character from position i at position (n-i+1)
        shuffle_sequence.new_sequence[old_pos] = 
            shuffle_sequence.old_sequence[new_pos]
        
        # Update position mapping
        shuffle_sequence.new_positions[old_pos] = new_pos
    end
    return nothing
end

function shuffle_sequence!(
    shuffle_sequence::ShuffleSeq,
    sequence::String;
    method::String = "shuffle"
)
    # Reset the sequence and positions
    resetSequence!(shuffle_sequence, sequence)
    
    # Fill movable positions
    fillMovablePositions!(shuffle_sequence)
    
    # Apply the selected decoy generation method
    if method == "shuffle"
        permuteNewPositions!(shuffle_sequence)
    elseif method == "reverse"
        reverseMovablePositions!(shuffle_sequence)
    else
        error("Unknown decoy method: $method. Must be 'shuffle' or 'reverse'")
    end
    
    return String(shuffle_sequence.new_sequence[1:shuffle_sequence.sequence_length])
end

"""
Enhanced version of shuffle_fast_with_positions that keeps specified characters fixed.
Pre-allocated vectors are passed in to avoid allocations.
"""
function shuffle_fast_with_positions_and_fixed_chars!(
    s::String, 
    positions::Vector{UInt8}, 
    fixed_chars::Set{Char},
    fixed_positions::Vector{Int},  # Pre-allocated
    movable_positions::Vector{Int},  # Pre-allocated
    temp_positions::Vector{UInt8}  # Pre-allocated for temporary storage
)
    ss = sizeof(s)
    l = length(s)
    
    # Count fixed and movable positions
    n_fixed = 0
    n_movable = 0
    
    # Create indices vector for byte positions
    v = Vector{Int}(undef, l)
    i = 1
    for j in 1:l
        v[j] = i
        i = nextind(s, i)
    end
    
    # Identify fixed and movable positions
    p = pointer(s)
    for j in 1:l
        c = Char(unsafe_load(p, v[j]))
        if j == l || c in fixed_chars  # Last position or fixed character
            n_fixed += 1
            fixed_positions[n_fixed] = j
        else
            n_movable += 1
            movable_positions[n_movable] = j
        end
    end
    
    # Generate permutation only for movable positions
    if n_movable > 0
        perm = randperm(n_movable)
        
        # Copy positions for movable characters to temp storage
        for idx in 1:n_movable
            temp_positions[idx] = positions[movable_positions[idx]]
        end
        
        # Apply permutation to positions
        for (new_idx, old_idx) in enumerate(perm)
            positions[movable_positions[new_idx]] = temp_positions[old_idx]
        end
        
        # Build the output string
        u = Vector{UInt8}(undef, ss)
        
        # Fill in the shuffled string
        for j in 1:l
            # Check if j is in fixed_positions (up to n_fixed)
            is_fixed = false
            for k in 1:n_fixed
                if fixed_positions[k] == j
                    is_fixed = true
                    break
                end
            end
            
            if is_fixed
                # Keep fixed characters in place
                u[v[j]] = unsafe_load(p, v[j])
            else
                # Find position in movable_positions
                idx = 0
                for k in 1:n_movable
                    if movable_positions[k] == j
                        idx = k
                        break
                    end
                end
                source_pos = movable_positions[perm[idx]]
                u[v[j]] = unsafe_load(p, v[source_pos])
            end
        end
        
        return String(u)
    else
        # All positions are fixed, return original string
        return s
    end
end

function adjust_mod_positions(
    mods::Union{Missing, Vector{PeptideMod}}, 
    positions::Vector{UInt8},
    seq_length::UInt8
)::Union{Missing, Vector{PeptideMod}}
    # If no modifications or missing, return as is
    if ismissing(mods) || isempty(mods)
        return mods
    end
    
    # Create new vector for adjusted modifications
    adjusted_mods = Vector{PeptideMod}(undef, length(mods))
    
    # Create a reverse mapping for efficiency
    # This maps original positions to new positions
    reverse_mapping = Vector{UInt8}(undef, seq_length)
    for new_pos in 1:seq_length
        orig_pos = positions[new_pos]
        if orig_pos <= seq_length
            reverse_mapping[orig_pos] = UInt8(new_pos)
        end
    end
    
    for (i, mod) in enumerate(mods)
        position = mod.position
        aa = mod.aa
        mod_name = mod.mod_name
        
        # Special case for N-terminal and C-terminal modifications
        if aa == 'n'
            # N-terminal modifications stay at position 1
            adjusted_mods[i] = PeptideMod(UInt8(1), 'n', mod_name)
        elseif aa == 'c'
            # C-terminal modifications stay at the end
            adjusted_mods[i] = PeptideMod(seq_length, 'c', mod_name)
        else
            # For normal residue modifications, use the reverse mapping
            if position <= seq_length
                adjusted_mods[i] = PeptideMod(reverse_mapping[position], aa, mod_name)
            else
                # Edge case handling if position is somehow out of bounds
                adjusted_mods[i] = mod
            end
        end
    end
    sort!(adjusted_mods)
    return adjusted_mods
end


"""
    add_decoy_sequences(target_fasta_entries::Vector{FastaEntry}; kwargs...)

Creates decoy sequences by mutating residues adjacent to each peptide terminus
according to a fixed substitution map.

# Parameters
- `target_fasta_entries::Vector{FastaEntry}`: Vector of target peptide entries to generate decoys for
- Additional keyword arguments (`max_shuffle_attempts`, `fixed_chars`,
  `decoy_method` and their `_`-prefixed variants) are accepted for backward
  compatibility but ignored by the termini-mutation workflow.

# Returns
- `Vector{FastaEntry}`: Sorted vector containing both original entries and their decoys

# Details
For each target peptide:
1. Determine mutation candidates near both termini using modification-aware rules.
2. Apply the substitution map (``GAVLIFMPWSCTYHKRQEND → LLLVVLLLLTSSSSLLNDQE``) to the
   first candidate that yields a unique sequence.
3. Mutate each terminus independently while preventing duplicates.
4. Preserve existing modification metadata and label the entry as a decoy.
5. Return the combined target and decoy peptides sorted by sequence.

# Notes
- Honors N- and C-terminal modifications when deciding mutation fallback order.
- Uses I/L equivalence when checking for sequence uniqueness.
- Keeps modification positions unchanged because the sequence length is constant.
"""
function add_decoy_sequences(
    target_fasta_entries::Vector{FastaEntry};
    kwargs...
    )
    kw = Dict{Symbol, Any}(kwargs)
    decoy_method = get(kw, :decoy_method, get(kw, :_decoy_method, "termini_mutation"))
    get(kw, :max_shuffle_attempts, get(kw, :_max_shuffle_attempts, 20))
    get(kw, :fixed_chars, get(kw, :_fixed_chars, Vector{Char}()))

    if decoy_method != "termini_mutation"
        @user_warn "Decoy method $(decoy_method) is no longer supported; defaulting to termini_mutation."
    end
    decoy_fasta_entries = Vector{FastaEntry}(undef, length(target_fasta_entries))
    sequences_set = PeptideSequenceSet(target_fasta_entries)

    n = 1
    for target_entry in target_fasta_entries
        target_sequence = get_sequence(target_entry)
        charge = get_charge(target_entry)
        seq_length = length(target_sequence)
        modified_positions, n_term_modded, c_term_modded = collect_modification_flags(
            seq_length,
            get_structural_mods(target_entry),
            get_isotopic_mods(target_entry),
        )

        decoy_sequence = generate_decoy_sequence(
            target_sequence,
            [charge],
            sequences_set,
            modified_positions,
            n_term_modded,
            c_term_modded,
        )

        if isnothing(decoy_sequence)
            @user_warn "Unable to generate unique decoy for $(target_sequence). Skipping."
            continue
        end

        decoy_fasta_entries[n] = FastaEntry(
            get_id(target_entry),
            get_description(target_entry),
            get_gene(target_entry),
            get_protein(target_entry),
            get_organism(target_entry),
            get_proteome(target_entry),
            decoy_sequence,
            get_start_idx(target_entry),
            get_structural_mods(target_entry),
            get_isotopic_mods(target_entry),
            get_charge(target_entry),
            get_base_target_id(target_entry),
            get_base_pep_id(target_entry),
            get_entrapment_pair_id(target_entry),
            true,
        )

        n += 1
        push!(sequences_set, decoy_sequence, get_charge(target_entry))
    end

    return sort(vcat(target_fasta_entries, decoy_fasta_entries[1:n-1]), by = x -> get_sequence(x))
end

"""
    add_decoy_sequences_grouped(target_fasta_entries::Vector{FastaEntry}; kwargs...)::Vector{FastaEntry}

Group-aware decoy generation using the termini mutation strategy so that all
modification variants of the same base peptide share a single decoy sequence.

# Parameters
- `target_fasta_entries::Vector{FastaEntry}`: Peptide entries to generate decoys for (typically includes targets and entrapments)
- Additional keyword arguments (`max_shuffle_attempts`, `fixed_chars`,
  `decoy_method` and their `_`-prefixed variants) are accepted for backward
  compatibility but ignored by the termini-mutation workflow.

# Returns
- `Vector{FastaEntry}`: Sorted vector with both original entries and their decoys

# Details
Algorithm:
1. Group entries by base sequence, ignoring modifications.
2. Aggregate modification information across variants to select mutation candidates.
3. Generate a single decoy sequence per base peptide using the termini mutation rules.
4. Reuse this sequence for all variants while preserving their modification metadata and setting `is_decoy = true`.
"""
function add_decoy_sequences_grouped(
    target_fasta_entries::Vector{FastaEntry};
    kwargs...
)::Vector{FastaEntry}
    kw = Dict{Symbol, Any}(kwargs)
    decoy_method = get(kw, :decoy_method, get(kw, :_decoy_method, "termini_mutation"))
    get(kw, :max_shuffle_attempts, get(kw, :_max_shuffle_attempts, 20))
    get(kw, :fixed_chars, get(kw, :_fixed_chars, Vector{Char}()))

    if decoy_method != "termini_mutation"
        @user_warn "Decoy method $(decoy_method) is no longer supported; defaulting to termini_mutation."
    end

    sequences_set = PeptideSequenceSet(target_fasta_entries)

    groups = Dict{String, Vector{Int}}()
    for (idx, entry) in enumerate(target_fasta_entries)
        base_seq = get_sequence(entry)
        if !haskey(groups, base_seq)
            groups[base_seq] = Vector{Int}()
        end
        push!(groups[base_seq], idx)
    end

    decoy_entries = Vector{FastaEntry}()
    for (base_seq, idxs) in groups
        charges = unique([get_charge(target_fasta_entries[i]) for i in idxs])
        seq_length = length(base_seq)
        modified_positions = falses(seq_length)
        n_term_modded = Ref(false)
        c_term_modded = Ref(false)

        for idx in idxs
            entry = target_fasta_entries[idx]
            merge_modification_flags!(
                modified_positions,
                n_term_modded,
                c_term_modded,
                get_structural_mods(entry),
                get_isotopic_mods(entry),
            )
        end

        decoy_sequence = generate_decoy_sequence(
            base_seq,
            charges,
            sequences_set,
            modified_positions,
            n_term_modded[],
            c_term_modded[],
        )

        if isnothing(decoy_sequence)
            @user_warn "Unable to generate unique decoy for $base_seq in grouped mode. Skipping."
            continue
        end

        for c in charges
            push!(sequences_set, decoy_sequence, c)
        end

        for idx in idxs
            target_entry = target_fasta_entries[idx]
            push!(decoy_entries, FastaEntry(
                get_id(target_entry),
                get_description(target_entry),
                get_gene(target_entry),
                get_protein(target_entry),
                get_organism(target_entry),
                get_proteome(target_entry),
                decoy_sequence,
                get_start_idx(target_entry),
                get_structural_mods(target_entry),
                get_isotopic_mods(target_entry),
                get_charge(target_entry),
                get_base_target_id(target_entry),
                get_base_pep_id(target_entry),
                get_entrapment_pair_id(target_entry),
                true,
            ))
        end
    end

    return sort(vcat(target_fasta_entries, decoy_entries), by = x -> get_sequence(x))
end

"""
    combine_shared_peptides(peptides::Vector{FastaEntry})::Vector{FastaEntry}

Combines entries that share identical peptide sequences by concatenating their protein accessions.

# Parameters
- `peptides::Vector{FastaEntry}`: Vector of peptide entries that may contain duplicates

# Returns
- `Vector{FastaEntry}`: Vector of unique peptide entries with concatenated metadata

# Details
This function:
1. Identifies peptides with identical sequences (considering I/L as equivalent)
2. For shared peptides, combines their protein accessions with semicolon separators
3. Also combines proteome identifiers and descriptions if multiple exist
4. Preserves other metadata from the first encountered instance of each sequence
5. Returns a vector containing only unique peptide sequences

# Examples
```julia
# Original entries with shared sequence
entries = [
    FastaEntry("P1", "desc1", "geneA", "protA", "human", "human", "PEPTIDE", 1, missing, missing, 0, 0, 0, 0, 0, 0, false),
    FastaEntry("P2", "desc2", "geneB", "protB", "human", "human", "PEPTIDE", 2, missing, missing, 0, 0, 0, 0, 0, 0, false),
    FastaEntry("P3", "desc3", "geneC", "protC", "human", "human", "UNIQUE", 3, missing, missing, 0, 0, 0, 0, 0, 0, false)
]

# Combine shared peptides
combined = combine_shared_peptides(entries)
# Results in 2 entries:
# 1. FastaEntry("P1;P2", "desc1;desc2", "geneA;geneB", "protA;protB", "human;human", "human;human", "PEPTIDE", 1, missing, missing, 0, 0, 0, 0, 0, 0, false)
# 2. FastaEntry("P3", "desc3", "geneC", "protC", "human", "human", "UNIQUE", 3, missing, missing, 0, 0, 0, 0, 0, 0, false)
```

# Notes
This function helps handle peptides that map to multiple proteins while maintaining
unique sequences in the library. It's particularly useful in bottom-up proteomics
where shared peptides are common.
"""
function combine_shared_peptides(peptides::Vector{FastaEntry})
    seq_to_fasta_entry = Dictionary{String, FastaEntry}()
    n = 0
    a = 0
    base_pep_id = one(UInt32)
    for peptide in peptides
        sequence = get_sequence(peptide)
        sequence_il_equiv = replace(sequence, 'I' => 'L')
        if haskey(seq_to_fasta_entry, sequence_il_equiv)
            a += 1
            fasta_entry = seq_to_fasta_entry[sequence_il_equiv]
            accession = get_id(peptide)*";"*get_id(fasta_entry)
            proteome = get_proteome(peptide)*";"*get_proteome(fasta_entry)
            description = get_description(peptide)*";"*get_description(fasta_entry)
            gene = get_gene(peptide)*";"*get_gene(fasta_entry)
            protein = get_protein(peptide)*";"*get_protein(fasta_entry)
            organism = get_organism(peptide)*";"*get_organism(fasta_entry)
            seq_to_fasta_entry[sequence_il_equiv] = FastaEntry(
                                                        accession,
                                                        description,
                                                        gene,
                                                        protein,
                                                        organism,
                                                        proteome,
                                                        get_sequence(fasta_entry),
                                                        get_start_idx(fasta_entry),
                                                        get_structural_mods(fasta_entry),
                                                        get_isotopic_mods(fasta_entry),
                                                        get_charge(fasta_entry),
                                                        get_base_target_id(fasta_entry), # preserve base_target_id
                                                        base_pep_id,
                                                        get_entrapment_pair_id(fasta_entry), 
                                                        is_decoy(fasta_entry)
                                                        )
            base_pep_id += one(UInt32)
        else
            n += 1
            insert!(seq_to_fasta_entry, sequence_il_equiv, peptide)
        end
    end
    fasta_entries = Vector{FastaEntry}(undef, length(seq_to_fasta_entry))
    i = 1
    for (key, value) in pairs(seq_to_fasta_entry)
        fasta_entries[i] = value
        i += 1
    end

    return fasta_entries
end

function assign_base_pep_ids!(fasta_entries::Vector{FastaEntry})
    """
    Assign sequential base_pep_id values starting from 1.
    Called after add_mods to identify unique peptides (sequence + modifications).
    Each peptide variant gets a unique base_pep_id for tracking through charge variants.
    
    Returns:
    - Int: Number of entries processed
    """
    
    for i in 1:length(fasta_entries)
        entry = fasta_entries[i]
        
        # Create new FastaEntry with sequential base_pep_id
        fasta_entries[i] = FastaEntry(
            get_id(entry),
            get_description(entry),
            get_gene(entry),
            get_protein(entry),
            get_organism(entry),
            get_proteome(entry),
            get_sequence(entry),
            get_start_idx(entry),
            get_structural_mods(entry),
            get_isotopic_mods(entry),
            get_charge(entry),
            get_base_target_id(entry), # preserve base_target_id
            UInt32(i),               # base_pep_id - sequential assignment
            get_entrapment_pair_id(entry),
            is_decoy(entry)
        )
    end
    
    return length(fasta_entries)
end


function assign_base_target_ids!(fasta_entries::Vector{FastaEntry})
    for i in 1:length(fasta_entries)
        entry = fasta_entries[i]
        # Create new FastaEntry with assigned base_target_id
        fasta_entries[i] = FastaEntry(
            get_id(entry),
            get_description(entry),
            get_gene(entry),
            get_protein(entry),
            get_organism(entry),
            get_proteome(entry),
            get_sequence(entry),
            get_start_idx(entry),
            get_structural_mods(entry),
            get_isotopic_mods(entry),
            get_charge(entry),
            UInt32(i),   # assign grouped base_target_id
            get_base_pep_id(entry),    # preserve existing base_pep_id
            get_entrapment_pair_id(entry),
            is_decoy(entry)
        )
    end
    
    return length(fasta_entries)
end
