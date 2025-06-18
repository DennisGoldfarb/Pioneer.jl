function parseMods(mods_string::AbstractString)::Base.RegexMatchIterator
    #Example: "1(7,M,Oxidation)(13,K,AnExampleMod)"
    mods_regex = r"(?<=\().*?(?=\))"
    return eachmatch(mods_regex, mods_string)
end

function insert_at_indices(original::S, insertions::Vector{Tuple{String, UInt8}}) where {S<:AbstractString}
    # Convert the original string into an array of characters for easier manipulation
    char_array = collect(original)

    # Sort the insertions by index in ascending order
    sorted_insertions = sort(insertions, by = x -> x[2])

    # Adjust the index for each insertion
    offset = 0
    for (substr, idx) in sorted_insertions
        # Adjust the index with the current offset
        insertion_index = idx + offset
        # Insert each character of the substring at the specified index
        for (i, char) in enumerate(substr)
            insert!(char_array, insertion_index + i, char)
        end
        # Update the offset by the length of the inserted substring
        offset += length(substr)
    end

    # Join the array of characters back into a single string
    return join(char_array)
end
function getModifiedSequence(
    sequence::S,
    isotope_mods::String,
    structural_mods::String,
    charge::UInt8) where {S<:AbstractString}

    mods = structural_mods*isotope_mods
    mods = [("("*getModName(mod.match)*")", getModIndex(mod.match)) for mod in parseMods(mods)]
    return "_"*insert_at_indices(sequence, mods)*"_."*string(charge)
end
function getModifiedSequence(
    sequence::S,
    isotope_mods::String,
    structural_mods::Missing,
    charge::UInt8) where {S<:AbstractString}

    mods = isotope_mods
    mods = [("("*getModName(mod.match)*")", getModIndex(mod.match)) for mod in parseMods(mods)]
    return "_"*insert_at_indices(sequence, mods)*"_."*string(charge)
end

function getModifiedSequence(
    sequence::S,
    isotope_mods::Missing,
    structural_mods::String,
    charge::UInt8) where {S<:AbstractString}

    mods = structural_mods
    mods = [("("*getModName(mod.match)*")", getModIndex(mod.match)) for mod in parseMods(mods)]
    return "_"*insert_at_indices(sequence, mods)*"_."*string(charge)
end
function getModifiedSequence(
    sequence::S,
    isotope_mods::Missing,
    structural_mods::Missing,
    charge::UInt8) where {S<:AbstractString}

    mods = ""
    mods = [("("*getModName(mod.match)*")", getModIndex(mod.match)) for mod in parseMods(mods)]
    return "_"*insert_at_indices(sequence, mods)*"_."*string(charge)
end

"""
    join_dict_values(dict::Dict{String,String}, key_string::Union{Missing,String})

Return a semicolon separated string of dictionary values corresponding to the
accession numbers contained in `key_string`. Missing entries are ignored.
"""
function join_dict_values(dict::Dict{String,String}, key_string::Union{Missing,String})
    if ismissing(key_string)
        return ""
    end
    accessions = split(String(key_string), ';')
    mapped = [get(dict, acc, "") for acc in accessions]
    nonempty = filter(!isempty, mapped)
    return join(nonempty, ';')
end

#Assume sorted by protein,peptide keys. Do this in batches and write a long and wide form .csv without
#loading the entire table into memory. 
function writePrecursorCSV(
    long_precursors_path::String,
    file_names::Vector{String},
    normalized::Bool;
    gene_dict::Dict{String,String} = Dict{String,String}(),
    protein_name_dict::Dict{String,String} = Dict{String,String}(),
    write_csv::Bool = true,
    batch_size::Int64 = 2000000)

    function makeWideFormat(
        longdf::DataFrame,
        cols::AbstractVector{Symbol},
        normalized::Bool)

        value_col = normalized ? :peak_area_normalized : :peak_area

        return unstack(longdf,
                    cols,
                    :file_name,
                    value_col;
                    combine = sum)      # or combine = maximum, first, etc.

    end
    precursors_long = DataFrame(Arrow.Table(long_precursors_path))
    n_rows = size(precursors_long, 1)

    # Map accession numbers to gene and protein names
    precursors_long[!,:gene_names] = join_dict_values.(Ref(gene_dict), precursors_long.accession_numbers)
    precursors_long[!,:protein_names] = join_dict_values.(Ref(protein_name_dict), precursors_long.accession_numbers)


    out_dir, arrow_path = splitdir(long_precursors_path)
    long_precursors_path = joinpath(out_dir,"precursors_long.tsv")
    wide_precursors_path = joinpath(out_dir,"precursors_wide.tsv")
    wide_precursors_arrow_path = joinpath(out_dir,"precursors_wide.arrow")
    if isfile(wide_precursors_arrow_path)
        rm(wide_precursors_arrow_path)
    end
    wide_columns = ["species"
    "gene_names"
    "protein_names"
    "accession_numbers"
    "inferred_protein_group"
    "sequence"
    "charge"
    "structural_mods"
    "isotopic_mods"
    "prec_mz"
    "global_qval"
    "use_for_protein_quant"
    "precursor_idx"
    "target"
    ]

    long_columns_exclude = [:isotopes_captured, :scan_idx, :protein_idx, :weight, :ms_file_idx]
    select!(precursors_long, Not(long_columns_exclude))
    rename!(precursors_long, [:new_best_scan => :apex_scan,
                              :prec_prob => :run_specific_prob,
                              :pg_score => :run_specific_pg_score,
                              :isotopes_captured_traces => :isotopes_captured,
                              :precursor_fraction_transmitted_traces => :precursor_fraction_transmitted])
    # order columns
    select!(precursors_long, [:file_name,
                             :species,
                             :gene_names,
                             :protein_names,
                             :accession_numbers,
                             :inferred_protein_group,
                             :sequence,
                             :charge,
                             :structural_mods,
                             :isotopic_mods,
                             :prec_mz,
                             :missed_cleavage,
                             :global_prob,
                             :run_specific_prob,
                             :global_qval,
                             :run_specific_qval,
                             :peak_area,
                             :peak_area_normalized,
                             :points_integrated,
                             :precursor_fraction_transmitted,
                             :isotopes_captured,
                             #:irt_obs,
                             :rt,
                             :apex_scan,
                             :global_pg_score,
                             :run_specific_pg_score,
                             :use_for_protein_quant,
                             :precursor_idx,
                             :target])

    sorted_columns = vcat(wide_columns, file_names)
    open(long_precursors_path,"w") do io1
        open(wide_precursors_path, "w") do io2
            #Make file headers 
            if write_csv
                println(io1,join(names(precursors_long),"\t"))
                println(io2,join(sorted_columns,"\t"))
            end
            batch_start_idx, batch_end_idx = 1, min(batch_size+1,n_rows)
            n_writes = 0
            while batch_start_idx <= n_rows
                #For the wide format, can't split a precursor between two batches.
                last_pid = precursors_long[batch_end_idx,:precursor_idx]
                while batch_end_idx < n_rows
                    if precursors_long[batch_end_idx + 1,:precursor_idx] != last_pid
                        break
                    end
                    batch_end_idx += 1
                end

                subdf =  precursors_long[range(batch_start_idx, batch_end_idx),:]
                batch_start_idx = batch_end_idx + 1
                batch_end_idx = min(batch_start_idx + batch_size, n_rows)
                if write_csv 
                    CSV.write(io1, subdf, append=true, header=false, delim="\t")
                end
                subunstack = makeWideFormat(subdf, Symbol.(wide_columns), normalized)
                col_names = names(subunstack)
                for fname in file_names
                    if fname ∉ col_names
                        subunstack[!,fname] .= Vector{Union{Missing,Float32}}(missing, nrow(subunstack))
                    end
                end
                if write_csv
                    CSV.write(io2, subunstack[!,sorted_columns], append=true,header=false,delim="\t")
                end
                if iszero(n_writes)
                    if isfile(wide_precursors_arrow_path)
                        rm(wide_precursors_arrow_path)
                    end
                    open(wide_precursors_arrow_path, "w") do io
                        Arrow.write(io, subunstack[!,sorted_columns]; file=false)  # file=false creates stream format
                    end
                else
                    Arrow.append(wide_precursors_arrow_path,subunstack[!,sorted_columns])
                end
                n_writes += 1
            end
        end
    end
    if write_csv == false
        rm(long_precursors_path, force = true)
        rm(wide_precursors_path, force = true)
    end
    return wide_precursors_arrow_path
end
function writeProteinGroupsCSV(
    long_pg_path::String,
    sequences::AbstractVector{<:AbstractString},
    isotope_mods::AbstractVector{Union{Missing, String}},
    structural_mods::AbstractVector{Union{String,Missing}},
    precursor_charge::AbstractVector{UInt8},
    file_names::Vector{String};
    gene_dict::Dict{String,String} = Dict{String,String}(),
    protein_name_dict::Dict{String,String} = Dict{String,String}(),
    write_csv::Bool = true,
    batch_size::Int64 = 2000000)

    function makeWideFormat(
        longdf::DataFrame)
        # First create a DataFrame with the non-abundance columns we want to keep
        metadata_df = unique(longdf[:, [:species, :gene_names, :protein_names, :protein, :target]])#, :n_peptides]])
        
        # Create the abundance wide format
        abundance_df = unstack(longdf,
            [:species, :protein, :target],
            :file_name, :abundance)
            
        # Join the metadata with the abundance data
        return leftjoin(metadata_df, abundance_df, on=[:species, :protein, :target])
    end

    protein_groups_long = DataFrame(Arrow.Table(long_pg_path))
    n_rows = size(protein_groups_long, 1)

    # Map accession numbers to gene names
    protein_groups_long[!,:gene_names] = join_dict_values.(Ref(gene_dict), protein_groups_long.protein)
    protein_groups_long[!,:protein_names] = join_dict_values.(Ref(protein_name_dict), protein_groups_long.protein)

    out_dir, arrow_path = splitdir(long_pg_path)
    long_protein_groups_path = joinpath(out_dir,"protein_groups_long.tsv")
    wide_protein_groups_path = joinpath(out_dir,"protein_groups_wide.tsv")
    wide_protein_groups_arrow = joinpath(out_dir,"protein_groups_wide.arrow")
    if isfile(wide_protein_groups_arrow)
        rm(wide_protein_groups_arrow)
    end
    # Update wide columns to include n_peptides
    wide_columns = ["species", "gene_names", "protein_names", "protein", "target"]

    select!(protein_groups_long, [:file_name,
                             :species,
                             :gene_names,
                             :protein_names,
                             :protein,
                             :target,
                             :entrap_id,
                             :peptides,
                             :n_peptides,
                             :global_qval,
                             :run_specific_qval,
                             :abundance])

    sorted_columns = vcat(wide_columns, file_names)
    open(long_protein_groups_path,"w") do io1
        open(wide_protein_groups_path, "w") do io2
            #Make file headers 
            if write_csv 
                println(io1,join(names(protein_groups_long),"\t"))
                println(io2,join(sorted_columns,"\t"))
            end
            batch_start_idx, batch_end_idx = 1, min(batch_size+1,n_rows)
            n_writes = 0
            while batch_start_idx <= n_rows
                #For the wide format, can't split a precursor between two batches.
                last_protein_group = protein_groups_long[batch_end_idx,:protein]
                while batch_end_idx < n_rows
                    if  protein_groups_long[batch_end_idx + 1,:protein] != last_protein_group
                        break
                    end
                    batch_end_idx += 1
                end
                subdf = protein_groups_long[range(batch_start_idx, batch_end_idx),:]
                batch_start_idx = batch_end_idx + 1
                batch_end_idx = min(batch_start_idx + batch_size, n_rows)
                subdf[!,:modified_sequence] = Vector{String}(undef, size(subdf, 1))
                for i in range(1, size(subdf, 1))
                    peptides = subdf[i,:peptides]
                    modified_sequences = Vector{String}(undef, length(peptides))
                    for (j, pid) in enumerate(peptides)
                        if ismissing(pid)
                            modified_sequences[j] = ""
                            continue
                        end
                        modified_sequences[j] = getModifiedSequence(
                            sequences[pid],
                            isotope_mods[pid],
                            structural_mods[pid],
                            precursor_charge[pid])
                    end
                    subdf[i,:modified_sequence] = join(filter(!isempty, modified_sequences),';')
                end
                subdf[!,:peptides] = subdf[!,:modified_sequence]
                select!(subdf, Not([:modified_sequence]))
                if write_csv 
                    CSV.write(io1, subdf, append=true, header=false, delim="\t")
                end
                subunstack = makeWideFormat(subdf)
                col_names = names(subunstack)
                for fname in file_names
                    if fname ∉ col_names
                        subunstack[!,fname] .= missing
                    end
                end
                if write_csv
                    CSV.write(io2, subunstack[!,sorted_columns], append=true,header=false,delim="\t")
                end

                if iszero(n_writes)
                    if isfile(wide_protein_groups_arrow)
                        rm(wide_protein_groups_arrow)
                    end
                    open(wide_protein_groups_arrow, "w") do io
                        Arrow.write(io, subunstack[!,sorted_columns]; file=false)  # file=false creates stream format
                    end
                else
                    Arrow.append(wide_protein_groups_arrow,subunstack[!,sorted_columns])
                end
                n_writes += 1
            end
        end
        if write_csv == false
            rm(long_protein_groups_path, force = true)
            rm(wide_protein_groups_path, force = true)
        end
    end
    return wide_protein_groups_arrow
end