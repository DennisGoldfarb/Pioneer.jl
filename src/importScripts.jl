function importScripts()
    #package_root = dirname(dirname(dirname(@__DIR__)))
    package_root = dirname(@__DIR__)
    function get_julia_files(dir::String)
        julia_files = String[]
        for (root, _, files) in walkdir(dir)
            for file in files
                if endswith(file, ".jl")
                    push!(julia_files, joinpath(root, file))
                end
            end
        end
        return julia_files
    end
    function include_files!(files_loded::Set{String}, file_dir::String, file_names::Vector{String})
        file_paths = [joinpath(file_dir, fname) for fname in file_names]
        [include(fpath) for fpath in file_paths if fpath ∉ files_loaded]
        push!(files_loaded, file_paths...)
        return nothing 
    end
    files_loaded = Set{String}()
    include_files!(
        files_loaded, 
        joinpath(package_root, "src","utils", "quadTransmissionModeling"),
        [
            "quadTransmissionModel.jl",
            "generalGaussModel.jl",
            "noQuadModel.jl",
            "RazoQuadModel.jl",
            "SplineQuadModel.jl",
            "binIsotopeRatioData.jl",
            "SquareQuadModel.jl"
        ]
    )

    include(joinpath(package_root, "src", "Routines","SearchDIA", "ParseInputs", "parseParams.jl"))
    push!(files_loaded, joinpath(package_root, "src", "Routines","SearchDIA", "ParseInputs", "parseParams.jl"))

    include_files!(
        files_loaded, 
        joinpath(package_root, "src","structs"),
        [
            "MassSpecData.jl",
            "ChromObject.jl",
            "ArrayDict.jl",
            "Counter.jl",
            "Ion.jl",
            "LibraryIon.jl",
            "LibraryFragmentIndex.jl",                                                                    
            "IsotopeTraceType.jl",
            "MatchIon.jl",
            "SparseArray.jl",
            "FragBoundModel.jl",
            "RetentionTimeIndex.jl",
            "MassErrorModel.jl",                                                                  
            "RetentionTimeConversionModel.jl"
            ]
    )

    
    # Utilities/ML
    include_files!(
        files_loaded,
        joinpath(package_root, "src", "utils", "ML"),
        [
            "percolatorSortOf.jl",
            "piecewiseLinearFunction.jl",
            "probitRegression.jl",
            "spectralLinearRegression.jl",
            "uniformBasisCubicSpline.jl",
            "wittakerHendersonSmoothing.jl",
            "libraryBSpline.jl"
        ]
    )


    # Utils
    include_files!(
        files_loaded,
        joinpath(package_root, "src", "utils"),
        [
            "isotopes.jl",
            "isotopeSplines.jl",
            "maxLFQ.jl",
                        "writeArrow.jl",
            "proteinInference.jl"
        ]
    )

    # PSMs
    include_files!(
        files_loaded,
        joinpath(package_root, "src", "Routines", "SearchDIA", "PSMs"),
        [
            "PSM.jl",
            "spectralDistanceMetrics.jl",
            "UnscoredPSMs.jl",
            "ScoredPSMs.jl"
        ]
    )

        
    #Search Method 
    include(joinpath(package_root, "src", "Routines", "SearchDIA", "SearchMethods", "SearchTypes.jl"))
    push!(files_loaded, joinpath(package_root, "src", "Routines", "SearchDIA", "SearchMethods", "SearchTypes.jl"))

    include(joinpath(package_root, "src", "Routines", "SearchDIA", "CommonSearchUtils", "selectTransitions", "selectTransitions.jl"))
    push!(files_loaded, joinpath(package_root, "src", "Routines", "SearchDIA", "CommonSearchUtils", "selectTransitions", "selectTransitions.jl"))

    #[println(fpath) for fpath in collect(files_loaded)]

    [include(jfile) for jfile in get_julia_files(joinpath(package_root, "src", "Routines", "SearchDIA", "CommonSearchUtils")) if jfile ∉ files_loaded]

    [include(jfile) for jfile in get_julia_files(joinpath(package_root, "src", "Routines", "SearchDIA", "ParseInputs")) if jfile ∉ files_loaded]

    [include(jfile) for jfile in get_julia_files(joinpath(package_root, "src", "Routines", "SearchDIA", "SearchMethods")) if jfile ∉ files_loaded]

    [include(jfile) for jfile in get_julia_files(joinpath(package_root, "src", "Routines", "SearchDIA", "WriteOutputs")) if jfile ∉ files_loaded]

    include(joinpath(package_root, "src", "Routines", "SearchDIA", "LibrarySearch.jl"))


    
    # BuildSpecLib
    root_path = joinpath(package_root, "src", "Routines", "BuildSpecLib")
    [include(fname) for fname in get_julia_files(joinpath(package_root, "src", "structs", "KoinaStructs"))]
    # FASTA processing
    include(joinpath(root_path, "fasta", "fasta_parser.jl"))
    include(joinpath(root_path, "fasta", "fasta_digest.jl"))
    include(joinpath(root_path, "fasta", "fasta_utils.jl"))
    # Fragment handling
    include(joinpath(root_path, "fragments", "get_frag_bounds.jl"))
    include(joinpath(root_path, "fragments", "fragment_parse.jl"))
    include(joinpath(root_path, "fragments", "fragment_index.jl"))
    include(joinpath(root_path, "fragments", "fragment_annotation.jl"))
    include(joinpath(root_path, "fragments", "fragment_predict.jl"))
    # Koina integration
    include(joinpath(root_path, "koina", "koina_api.jl"))
    include(joinpath(root_path, "koina", "koina_batch_prep.jl"))
    include(joinpath(root_path, "koina", "koina_batch_parse.jl"))
    # Utilities
    include(joinpath(root_path, "utils", "io.jl"))
    include(joinpath(root_path, "utils", "estimate_collision_ev.jl"))
    include(joinpath(root_path, "utils", "math.jl"))
    include(joinpath(root_path, "utils", "get_mz.jl"))
    include(joinpath(root_path, "utils", "parse_isotope_mods.jl"))
    include(joinpath(root_path, "utils", "check_params.jl"))
    #structs 
    include(joinpath(root_path, "structs", "EmpiricalLibrary.jl"))
    include(joinpath(root_path, "utils", "parse_mods.jl"))
    # Library building
    include(joinpath(root_path, "build", "build_poin_lib.jl"))
    #Chronologer Methods 
    include(joinpath(root_path, "chronologer", "chronologer_prep.jl"))
    include(joinpath(root_path, "chronologer", "chronologer_predict.jl"))
    include(joinpath(root_path, "chronologer", "chronologer_parse.jl"))
end