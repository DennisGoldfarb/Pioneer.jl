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

const TERMINAL_MUTATION_SOURCE = collect("GAVLIFMPWSCTYHKRQEND")
const TERMINAL_MUTATION_TARGET = collect("LLLVVLLLLTSSSSLLNDQE")
const TERMINAL_MUTATION_MAP = Dict(zip(TERMINAL_MUTATION_SOURCE, TERMINAL_MUTATION_TARGET))

struct MutationContext
    modified_positions::Set{Int}
    n_term_modified::Bool
    c_term_modified::Bool
end

MutationContext() = MutationContext(Set{Int}(), false, false)

function _add_mods!(ctx::MutationContext, mods::Union{Missing,Vector{PeptideMod}})
    if mods === missing
        return ctx
    end
    for mod in mods
        aa = getAA(mod)
        if aa == 'n'
            ctx.n_term_modified = true
        elseif aa == 'c'
            ctx.c_term_modified = true
        else
            push!(ctx.modified_positions, Int(getPosition(mod)))
        end
    end
    return ctx
end

function MutationContext(structural_mods::Union{Missing,Vector{PeptideMod}}, isotopic_mods::Union{Missing,Vector{PeptideMod}})
    ctx = MutationContext()
    _add_mods!(ctx, structural_mods)
    _add_mods!(ctx, isotopic_mods)
    return ctx
end

struct MutationCandidate
    idx::Int
    force::Bool
end

function _candidate_order(side::Symbol, n::Int, ctx::MutationContext)
    base = side === :N ? 2 : n - 1
    fallback = side === :N ? 3 : n - 2
    terminal = side === :N ? 1 : n

    candidates = MutationCandidate[
        MutationCandidate(base, false),
        MutationCandidate(fallback, false),
        MutationCandidate(terminal, false),
    ]

    terminal_modified = if side === :N
        ctx.n_term_modified || (terminal >= 1 && terminal <= n && terminal in ctx.modified_positions)
    else
        ctx.c_term_modified || (terminal >= 1 && terminal <= n && terminal in ctx.modified_positions)
    end

    if terminal_modified && base >= 1 && base <= n
        push!(candidates, MutationCandidate(base, true))
    end

    filtered = MutationCandidate[]
    for cand in candidates
        if cand.idx < 1 || cand.idx > n
            continue
        end
        push!(filtered, cand)
    end
    return filtered
end

function _apply_candidate(seq_chars::Vector{Char}, cand::MutationCandidate, ctx::MutationContext)
    idx = cand.idx
    if idx < 1 || idx > length(seq_chars)
        return nothing
    end
    if idx in ctx.modified_positions && !cand.force
        return nothing
    end
    aa = seq_chars[idx]
    new_aa = get(TERMINAL_MUTATION_MAP, aa, nothing)
    if new_aa === nothing || new_aa == aa
        return nothing
    end
    mutated = copy(seq_chars)
    mutated[idx] = new_aa
    return mutated, Dict(idx => new_aa)
end

function _merge_mutations(m1::Dict{Int,Char}, m2::Dict{Int,Char})
    if isempty(m1)
        return copy(m2)
    end
    merged = copy(m1)
    for (k, v) in m2
        merged[k] = v
    end
    return merged
end

function generate_terminal_decoy(
    sequence::AbstractString,
    ctx::MutationContext,
    exists_fn::Function;
    allow_skip_c_term::Bool = true,
)
    n = length(sequence)
    if n == 0
        return nothing
    end
    seq_chars = collect(sequence)
    n_candidates = _candidate_order(:N, n, ctx)
    c_candidates = _candidate_order(:C, n, ctx)
    if allow_skip_c_term && isempty(c_candidates)
        c_candidates = [MutationCandidate(0, true)]
    elseif allow_skip_c_term
        push!(c_candidates, MutationCandidate(0, true))
    end

    for n_cand in n_candidates
        n_result = _apply_candidate(seq_chars, n_cand, ctx)
        if n_result === nothing
            continue
        end
        n_chars, n_mut = n_result
        for c_cand in c_candidates
            if c_cand.idx != 0 && c_cand.idx == n_cand.idx
                continue
            end
            if c_cand.idx == 0
                final_chars = n_chars
                final_mut = copy(n_mut)
            else
                c_result = _apply_candidate(n_chars, c_cand, ctx)
                if c_result === nothing
                    continue
                end
                final_chars, c_mut = c_result
                final_mut = _merge_mutations(n_mut, c_mut)
            end
            if isempty(final_mut)
                continue
            end
            final_seq = String(final_chars)
            if exists_fn(final_seq)
                continue
            end
            return final_seq, final_mut
        end
    end
    return nothing
end

function update_mods(mods::Union{Missing,Vector{PeptideMod}}, mutated_positions::Dict{Int,Char})
    if mods === missing || isempty(mods)
        return mods
    end
    updated = Vector{PeptideMod}(undef, length(mods))
    changed = false
    for (idx, mod) in enumerate(mods)
        aa = getAA(mod)
        pos = Int(getPosition(mod))
        if (aa == 'n') || (aa == 'c') || !(pos in mutated_positions)
            updated[idx] = mod
        else
            updated[idx] = PeptideMod(mod.position, mutated_positions[pos], mod.mod_name)
            changed = true
        end
    end
    return changed ? sort!(updated) : mods
end

const MOD_STRING_REGEX = r"\((\d+),([A-Za-z]),([^,\)]+)\)"

function mutation_context_from_string(mod_string::AbstractString)
    ctx = MutationContext()
    if isempty(mod_string)
        return ctx
    end
    for m in eachmatch(MOD_STRING_REGEX, mod_string)
        pos = parse(Int, m.captures[1])
        aa = first(m.captures[2])
        if aa == 'n'
            ctx.n_term_modified = true
        elseif aa == 'c'
            ctx.c_term_modified = true
        else
            push!(ctx.modified_positions, pos)
        end
    end
    return ctx
end

function update_mod_string(mod_string::AbstractString, mutated_positions::Dict{Int,Char})
    if isempty(mod_string) || isempty(mutated_positions)
        return mod_string
    end
    parts = String[]
    last_index = 1
    for m in eachmatch(MOD_STRING_REGEX, mod_string)
        push!(parts, mod_string[last_index:m.offset-1])
        pos = parse(Int, m.captures[1])
        aa = first(m.captures[2])
        name = m.captures[3]
        new_aa = (aa == 'n' || aa == 'c') ? aa : get(mutated_positions, pos, aa)
        push!(parts, "($(pos),$(new_aa),$(name))")
        last_index = m.offset + length(m.match)
    end
    push!(parts, mod_string[last_index:end])
    return join(parts)
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
    add_decoy_sequences(target_fasta_entries::Vector{FastaEntry}; max_shuffle_attempts::Int64 = 20)

Create terminal-mutation decoys by mutating residues adjacent to the peptide termini.

# Parameters
- `target_fasta_entries::Vector{FastaEntry}`: Target peptide entries to generate decoys for
- `max_shuffle_attempts::Int64`: (Unused; maintained for API compatibility)

# Returns
- `Vector{FastaEntry}`: Sorted vector containing both original entries and their decoys

# Details
For each target peptide:
1. Mutates the residue adjacent to the N-terminus using the predefined substitution map
2. Mutates the residue adjacent to the C-terminus using the same map (with fallbacks when modified)
3. Avoids residues with modifications unless required by the fallback rules
4. Skips candidate sequences that already exist in the target set
5. Updates modification annotations when a mutated residue carries a modification
6. Preserves peptide metadata while marking decoys with `is_decoy = true`

Decoy sequences are considered unique using I/L equivalence as in `PeptideSequenceSet`.
"""
function add_decoy_sequences(
    target_fasta_entries::Vector{FastaEntry};
    _max_shuffle_attempts::Int64 = 20,
    _fixed_chars::Vector{Char} = Vector{Char}(),
    _decoy_method::String = "shuffle"
    )
    decoy_fasta_entries = Vector{FastaEntry}(undef, length(target_fasta_entries))
    sequences_set = PeptideSequenceSet(target_fasta_entries)

    n = 1
    for target_entry in target_fasta_entries
        target_sequence = get_sequence(target_entry)
        charge = get_charge(target_entry)
        ctx = MutationContext(get_structural_mods(target_entry), get_isotopic_mods(target_entry))
        exists_fn = seq -> (seq, charge) ∈ sequences_set

        result = generate_terminal_decoy(target_sequence, ctx, exists_fn)

        if result === nothing
            @user_warn "Unable to generate decoy via terminal mutations for $(target_sequence)"
            continue
        end

        decoy_sequence, mutated_positions = result
        adjusted_structural_mods = update_mods(get_structural_mods(target_entry), mutated_positions)
        adjusted_isotopic_mods = update_mods(get_isotopic_mods(target_entry), mutated_positions)

        decoy_fasta_entries[n] = FastaEntry(
            get_id(target_entry),
            get_description(target_entry),
            get_gene(target_entry),
            get_protein(target_entry),
            get_organism(target_entry),
            get_proteome(target_entry),
            decoy_sequence,
            get_start_idx(target_entry),
            adjusted_structural_mods,
            adjusted_isotopic_mods,
            get_charge(target_entry),
            get_base_target_id(target_entry),
            get_base_pep_id(target_entry),
            get_entrapment_pair_id(target_entry),
            true
        )

        n += 1
        push!(sequences_set, decoy_sequence, charge)
    end

    return sort(vcat(target_fasta_entries, decoy_fasta_entries[1:n-1]), by = x -> get_sequence(x))
end

"""
    add_decoy_sequences_grouped(
        target_fasta_entries::Vector{FastaEntry};
        max_shuffle_attempts::Int64 = 20,
        fixed_chars::Vector{Char} = Vector{Char}(),
        decoy_method::String = "shuffle"
    )::Vector{FastaEntry}

Generate terminal-mutation decoys for groups of peptides sharing the same base
sequence while respecting modification constraints across all variants.

# Parameters
- `target_fasta_entries::Vector{FastaEntry}`: Peptide entries (targets and/or entrapments)
- `max_shuffle_attempts::Int64`: (Unused; preserved for compatibility)
- `fixed_chars::Vector{Char}`: (Unused with terminal mutations)
- `decoy_method::String`: (Unused; maintained for compatibility)

# Returns
- `Vector{FastaEntry}`: Sorted vector with both original entries and their decoys

# Details
1. Groups peptides by unmodified base sequence
2. Aggregates modification context across all variants in the group
3. Generates a single decoy sequence per group using terminal mutations
4. Applies the resulting amino-acid substitutions to all modification variants
5. Updates modification annotations when mutated residues carry modifications
6. Ensures decoy uniqueness using `PeptideSequenceSet`
"""
function add_decoy_sequences_grouped(
    target_fasta_entries::Vector{FastaEntry};
    _max_shuffle_attempts::Int64 = 20,
    _fixed_chars::Vector{Char} = Vector{Char}(),
    _decoy_method::String = "shuffle"
)::Vector{FastaEntry}

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

        ctx = MutationContext()
        for idx in idxs
            entry = target_fasta_entries[idx]
            _add_mods!(ctx, get_structural_mods(entry))
            _add_mods!(ctx, get_isotopic_mods(entry))
        end

        exists_fn = seq -> any(((seq, c) ∈ sequences_set) for c in charges)
        result = generate_terminal_decoy(base_seq, ctx, exists_fn)

        if result === nothing
            @user_warn "Unable to generate grouped decoy via terminal mutations for $(base_seq)"
            continue
        end

        decoy_sequence, mutated_positions = result

        for c in charges
            push!(sequences_set, decoy_sequence, c)
        end

        for idx in idxs
            target_entry = target_fasta_entries[idx]
            adjusted_structural_mods = update_mods(get_structural_mods(target_entry), mutated_positions)
            adjusted_isotopic_mods = update_mods(get_isotopic_mods(target_entry), mutated_positions)

            push!(decoy_entries, FastaEntry(
                get_id(target_entry),
                get_description(target_entry),
                get_gene(target_entry),
                get_protein(target_entry),
                get_organism(target_entry),
                get_proteome(target_entry),
                decoy_sequence,
                get_start_idx(target_entry),
                adjusted_structural_mods,
                adjusted_isotopic_mods,
                get_charge(target_entry),
                get_base_target_id(target_entry),
                get_base_pep_id(target_entry),
                get_entrapment_pair_id(target_entry),
                true
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
