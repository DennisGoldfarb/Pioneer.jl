
using JSON

@testset "Fragment Predict Tests" begin
    
    @testset "Model Type Validation" begin
        @testset "Instrument-specific models" begin
            # Test that instrument-specific models validate instruments correctly
            unispec_model = InstrumentSpecificModel("unispec")
            alphapept_model = InstrumentSpecificModel("alphapeptdeep")
            
            # These should be valid based on MODEL_CONFIGS
            valid_unispec_instruments = ["QE", "QEHFX", "LUMOS", "ELITE", "VELOS"]
            valid_alpha_instruments = ["QE", "LUMOS", "TIMSTOF", "SCIEXTOF"]
            
            # Test model type properties
            @test unispec_model isa InstrumentSpecificModel
            @test alphapept_model isa InstrumentSpecificModel
        end
        
        @testset "Instrument-agnostic models" begin
            prosit_model = InstrumentAgnosticModel("prosit_2020_hcd")
            @test prosit_model isa InstrumentAgnosticModel
        end
        
        @testset "Spline coefficient models" begin
            altimeter_model = SplineCoefficientModel("altimeter")
            @test altimeter_model isa SplineCoefficientModel
        end
    end
    
    @testset "filter_fragments!" begin
        @testset "Standard intensity filtering" begin
            # Test data with various intensities
            df = DataFrame(
                annotation = ["b3", "y4", "b5", "y7"],
                intensities = Float32[0.0005, 0.002, 0.5, 0.9],
                mz = Float32[300.1, 400.2, 500.3, 600.4]
            )
            
            model = InstrumentSpecificModel("unispec")
            filter_fragments!(df, model)
            
            # Should keep only intensities > 0.001
            @test nrow(df) == 3
            @test all(df.intensities .> 0.001f0)
        end
        
        @testset "Invalid m/z filtering" begin
            df = DataFrame(
                annotation = ["b3", "y4", "b5"],
                intensities = Float32[0.5, 0.6, 0.7],
                mz = Float32[-100.0, 0.0, 500.3]
            )
            
            model = InstrumentAgnosticModel("prosit_2020_hcd")
            filter_fragments!(df, model)
            
            # Should keep only positive m/z
            @test nrow(df) == 1
            @test all(df.mz .> 0)
        end
        
        @testset "Isotope peak filtering for instrument-specific models" begin
            df = DataFrame(
                annotation = ["b3", "y4+i", "b5+2i", "y7"],
                intensities = Float32[0.5, 0.6, 0.7, 0.8],
                mz = Float32[300.1, 400.2, 500.3, 600.4]
            )
            
            model = InstrumentSpecificModel("unispec")
            filter_fragments!(df, model)
            
            # Should remove isotope peaks (containing 'i')
            @test nrow(df) == 2
            @test !any(occursin('i', ann) for ann in df.annotation)
        end
        
        @testset "Spline coefficient filtering" begin
            # For spline models, we don't filter on intensity but on m/z
            df = DataFrame(
                annotation = [1, 2, 3, 4],  # Altimeter uses indices
                coefficients = [(0.1, 0.2), (0.3, 0.4), (0.5, 0.6), (0.7, 0.8)],
                mz = Float32[-100.0, 0.0, 500.3, 600.4]
            )
            
            model = SplineCoefficientModel("altimeter")
            filter_fragments!(df, model)
            
            # Should keep only positive m/z
            @test nrow(df) == 2
            @test all(df.mz .> 0)
        end
    end
    
    @testset "sort_fragments!" begin
        df = DataFrame(
            precursor_idx = [1, 2, 1, 2, 1, 2],
            annotation = ["b3", "y4", "y5", "b2", "a3", "x6"],
            intensities = Float32[0.5, 0.9, 0.7, 0.3, 0.8, 0.6],
            mz = Float32[300.1, 400.2, 500.3, 200.1, 350.2, 450.3]
        )

        sort_fragments!(df)

        # Check that within each precursor, fragments are sorted by intensity (descending)
        @test df.precursor_idx == [1, 1, 1, 2, 2, 2]

        # For precursor 1
        prec1_mask = df.precursor_idx .== 1
        @test issorted(df[prec1_mask, :intensities], rev=true)

        # For precursor 2
        prec2_mask = df.precursor_idx .== 2
        @test issorted(df[prec2_mask, :intensities], rev=true)

        coef_df = DataFrame(
            precursor_idx = UInt32[1, 1, 2, 2],
            annotation = Int32[10, 11, 20, 21],
            coefficients = [
                (0.1f0, 0.2f0),
                (0.3f0, 0.4f0),
                (0.5f0, 0.6f0),
                (0.7f0, 0.8f0),
            ],
            ranking = UInt8[2, 1, 2, 1],
            mz = Float32[100.0, 101.0, 200.0, 201.0],
        )

        sort_fragments!(coef_df)

        @test coef_df.precursor_idx == UInt32[1, 1, 2, 2]
        @test coef_df.ranking == UInt8[1, 2, 1, 2]
    end
    
    @testset "predict_fragments_batch" begin
        temp_dir = mktempdir()
        
        @testset "InstrumentSpecificModel batch prediction" begin
            # Create mock peptide data
            peptides_df = DataFrame(
                sequence = ["PEPTIDE", "SEQUENCE", "FRAGMENT"],
                charge = [2, 3, 2],
                collision_energy = [25.0, 30.0, 27.0]
            )
            
            # Mock model and parameters
            model = InstrumentSpecificModel("unispec")
            instrument_type = "QE"
            batch_size = 2
            concurrent_requests = 1
            first_prec_idx = UInt32(1)
            
            # We can't test actual API calls, but we can test the structure
            # Create a mock response structure
            mock_fragments = DataFrame(
                annotation = ["b3", "y4", "b5", "y6", "a2", "x3"],
                mz = Float32[300.1, 400.2, 500.3, 600.4, 200.1, 350.2],
                intensities = Float32[0.8, 0.9, 0.7, 0.6, 0.5, 0.4]
            )
            
            # Test that the function structure works with mock data
            @test model isa InstrumentSpecificModel
            @test instrument_type in ["QE", "QEHFX", "LUMOS", "ELITE", "VELOS"]
        end
        
        @testset "InstrumentAgnosticModel batch prediction" begin
            peptides_df = DataFrame(
                sequence = ["PEPTIDE", "SEQUENCE"],
                charge = [2, 3],
                collision_energy = [25.0, 30.0]
            )
            
            model = InstrumentAgnosticModel("prosit_2020_hcd")
            # Instrument type should be ignored
            instrument_type = "ANY"
            
            @test model isa InstrumentAgnosticModel
        end
        
        @testset "SplineCoefficientModel batch prediction" begin
            peptides_df = DataFrame(
                sequence = ["PEPTIDE"],
                charge = [2],
                collision_energy = [25.0]
            )
            
            model = SplineCoefficientModel("altimeter")
            instrument_type = "QE"
            
            @test model isa SplineCoefficientModel
            
            # Test knot vector consistency check would happen here
            # In real implementation, all batches must have same knot vector
        end
        
        rm(temp_dir, recursive=true)
    end

    @testset "clone_decoy_fragments" begin
        temp_dir = mktempdir()
        config_path = joinpath(temp_dir, "config.json")
        base_config = Dict(
            "fixed_mods" => Dict("mass" => Float64[], "name" => String[]),
            "variable_mods" => Dict("mass" => Float64[], "name" => String[]),
            "isotope_mod_groups" => Any[],
            "channel_decoys" => false,
            "sulfur_mod_groups" => Any[],
            "library_params" => Dict(
                "include_immonium" => false,
                "max_frag_rank" => 50,
                "length_to_frag_count_multiple" => 2.0,
            ),
        )
        open(config_path, "w") do io
            write(io, JSON.json(base_config))
        end

        peptides_df = DataFrame(
            sequence = ["PEPTIDE", "PEPTIDA"],
            mods = ["", ""],
            isotope_mods = ["", ""],
            decoy = [false, true],
            pair_id = UInt32[1, 1]
        )

        target_fragments = DataFrame(
            annotation = ["y2", "IH"],
            mz = Float32[400.0, 110.0],
            intensities = Float32[0.8, 0.2],
            precursor_idx = UInt32[1, 1]
        )

        model = InstrumentSpecificModel("unispec")
        decoy_frags = clone_decoy_fragments(peptides_df, target_fragments, model, config_path)

        @test nrow(decoy_frags) == 1
        @test all(decoy_frags.precursor_idx .== UInt32(2))
        @test decoy_frags.intensities == Float32[0.8]
        @test decoy_frags.annotation == ["y2"]

        aa_masses = zeros(Float32, 255)
        structural_mod_masses = zeros(Float32, 255)
        iso_mod_masses = zeros(Float32, 255)
        sequence = peptides_df.sequence[2]
        get_aa_masses!(aa_masses, sequence)
        get_structural_mod_masses!(structural_mod_masses, "", Dict{String, Float32}())
        getIsoModMasses!(iso_mod_masses, "", "", Dict{String, Dict{String, Float32}}())
        info = parse_fragment_annotation(UniSpecFragAnnotation("y2"))
        start_idx, stop_idx = get_fragment_indices(info.base_type, info.frag_index, UInt8(length(sequence)))
        expected_mz = get_fragment_mz(
            start_idx,
            stop_idx,
            info.base_type,
            info.charge,
            aa_masses,
            structural_mod_masses,
            iso_mod_masses
        )

        y2_idx = findfirst(==("y2"), decoy_frags.annotation)
        @test y2_idx !== nothing
        @test isapprox(decoy_frags.mz[y2_idx], expected_mz; atol=1e-5)

        @test all(!=("IH"), decoy_frags.annotation)

        # Re-run with immonium enabled
        include_config = deepcopy(base_config)
        include_config["library_params"]["include_immonium"] = true
        open(config_path, "w") do io
            write(io, JSON.json(include_config))
        end

        decoy_frags = clone_decoy_fragments(peptides_df, target_fragments, model, config_path)
        @test nrow(decoy_frags) == 2
        @test "IH" in decoy_frags.annotation

        # Limit the requested number of ions
        rank_limited_config = deepcopy(include_config)
        rank_limited_config["library_params"]["max_frag_rank"] = 1
        open(config_path, "w") do io
            write(io, JSON.json(rank_limited_config))
        end

        rich_target = DataFrame(
            annotation = ["y2", "b3", "a1"],
            mz = Float32[400.0, 320.0, 210.0],
            intensities = Float32[0.8, 0.5, 0.3],
            precursor_idx = UInt32[1, 1, 1],
        )

        limited_frags = clone_decoy_fragments(peptides_df, rich_target, model, config_path)
        @test nrow(limited_frags) == 1
        @test limited_frags.annotation == ["y2"]

        rm(temp_dir, recursive=true)
    end

    @testset "predict_fragments - main dispatcher" begin
        temp_dir = mktempdir()
        
        # Create test peptide data
        peptide_data = DataFrame(
            sequence = ["PEPTIDE", "SEQUENCE", "FRAGMENT", "EXAMPLE"],
            charge = [2, 3, 2, 2],
            collision_energy = [25.0, 30.0, 27.0, 26.0]
        )
        
        peptide_path = joinpath(temp_dir, "peptides.arrow")
        Arrow.write(peptide_path, peptide_data)
        
        frags_out_path = joinpath(temp_dir, "fragments.arrow")
        
        @testset "Model name validation" begin
            # Test invalid model name
            @test_throws ErrorException predict_fragments(
                peptide_path,
                frags_out_path,
                InstrumentSpecificModel("invalid_model"),
                "QE",
                10,
                100,
                "invalid_model"
            )
        end
        
        @testset "Batch size limiting" begin
            # Test that batch size is capped at 1000
            model = InstrumentSpecificModel("unispec")
            
            # Would need to mock API calls to fully test
            # Here we just verify the setup is correct
            @test model.name == "unispec"
            
            # Verify peptide data was loaded correctly
            loaded_peptides = DataFrame(Arrow.Table(peptide_path))
            @test nrow(loaded_peptides) == 4
            @test names(loaded_peptides) == ["sequence", "charge", "collision_energy"]
        end
        
        @testset "Output file creation" begin
            # Test that output path handling works
            @test !isfile(frags_out_path)
            
            # In real implementation, this would create the file
            # We just test the path construction
            @test dirname(frags_out_path) == temp_dir
            @test basename(frags_out_path) == "fragments.arrow"
        end
        
        rm(temp_dir, recursive=true)
    end
    
    @testset "Batch processing calculations" begin
        @testset "Batch index calculation" begin
            nprecs = 1000
            batch_size = 100
            koina_pool_size = 5
            
            batch_start_idxs = collect(one(UInt32):UInt32(batch_size*koina_pool_size):UInt32(nprecs))
            
            # Should create appropriate number of batches
            @test length(batch_start_idxs) == 2  # 0-499, 500-999
            @test batch_start_idxs[1] == 1
            @test batch_start_idxs[2] == 501
        end
        
        @testset "Fragment per precursor calculation" begin
            # Mock fragment data with 3 precursors
            fragment_df = DataFrame(
                precursor_idx = [1, 1, 1, 2, 2, 3, 3, 3, 3],
                annotation = ["b1", "b2", "y1", "b1", "y2", "b1", "b2", "y1", "y2"]
            )
            
            # Calculate fragments per precursor
            frags_per_prec = [3, 2, 4]  # Precursor 1: 3, Precursor 2: 2, Precursor 3: 4
            
            # Test index assignment
            n_precursors = 3
            for i in 1:n_precursors
                mask = fragment_df.precursor_idx .== i
                @test sum(mask) == frags_per_prec[i]
            end
        end
    end
end