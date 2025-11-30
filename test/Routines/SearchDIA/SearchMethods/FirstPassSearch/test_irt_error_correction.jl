using Test
using Pioneer
using Arrow
using DataFrames
using Plots
using Random

struct MockSearchStructures <: Pioneer.SearchDataStructures end

struct MockPrecursors
    sequences::Vector{String}
    is_decoy::Vector{Bool}
    irt::Vector{Float32}
end

struct MockSpectralLibrary <: Pioneer.SpectralLibrary
    precursors::MockPrecursors
end

Pioneer.getPrecursors(lib::MockSpectralLibrary) = lib.precursors
Pioneer.getSequence(precursors::MockPrecursors) = precursors.sequences
Pioneer.getIsDecoy(precursors::MockPrecursors) = precursors.is_decoy
Pioneer.getIrt(precursors::MockPrecursors) = precursors.irt

function build_mock_context(psms_df::DataFrame, sequences::Vector{String},
                            is_decoy::Vector{Bool}, library_irt::Vector{Float32})
    tmp_dir = mktempdir()
    psms_path = joinpath(tmp_dir, "psms.arrow")
    Arrow.write(psms_path, psms_df)

    ms_reference = Pioneer.ArrowTableReference([psms_path])
    Pioneer.setFirstPassPsms!(ms_reference, 1, psms_path)

    library = MockSpectralLibrary(MockPrecursors(sequences, is_decoy, library_irt))
    temp_structures = [MockSearchStructures()]

    context = Pioneer.SearchContext(library, temp_structures, ms_reference,
                                    1, length(sequences), 1)
    Pioneer.setRtIrtMap!(context, Pioneer.IdentityModel(), 1)
    return context
end

function compute_expected_irt(sequence::AbstractString, intercept::Float64,
                              aa_weights::Dict{Char, Float64}, base_irt::Float32)
    err = intercept
    for aa in uppercase(sequence)
        err += get(aa_weights, aa, 0.0)
    end
    return Float32(Float64(base_irt) - err)
end

@testset "First pass iRT error correction" begin
    canonical = collect("ACDEFGHIKLMNPQRSTVWY")
    intercept = 0.25
    aa_weights = Dict{Char, Float64}(aa => 0.05 * idx for (idx, aa) in enumerate(canonical))

    target_sequences = [string(aa) for aa in canonical]
    push!(target_sequences, join(canonical))
    decoy_sequence = "WYAC"
    sequences = vcat(target_sequences, decoy_sequence)

    library_irt = Float32.(100.0:1.0:100.0 + length(sequences) - 1)
    is_decoy = vcat(fill(false, length(target_sequences)), true)

    errors = [intercept + sum(get(aa_weights, aa, 0.0) for aa in uppercase(seq))
              for seq in sequences]

    ms_file_idx = UInt32(1)
    psms_df = DataFrame(
        ms_file_idx = fill(ms_file_idx, length(sequences)),
        scan_idx = UInt32.(1:length(sequences)),
        precursor_idx = UInt32.(1:length(sequences)),
        rt = Float32.(library_irt .- errors),
        irt_predicted = library_irt,
        q_value = vcat(fill(Float32(0.001), length(target_sequences)), Float32(0.1)),
        score = zeros(Float32, length(sequences)),
        prob = zeros(Float32, length(sequences)),
        scan_count = fill(UInt32(1), length(sequences))
    )

    ctx = build_mock_context(psms_df, sequences, is_decoy, library_irt)
    Pioneer.correct_first_pass_irt_errors!(ctx)

    for (idx, seq) in enumerate(sequences)
        expected = compute_expected_irt(seq, intercept, aa_weights, library_irt[idx])
        @test isapprox(Pioneer.getPredIrt(ctx, UInt32(idx)), expected; atol=1f-5)
    end

    # Fallback when insufficient training samples
    limited_psms = psms_df[1:2, :]
    limited_psms.q_value .= Float32(0.2)
    ctx_fail = build_mock_context(limited_psms, sequences, is_decoy, library_irt)
    Pioneer.correct_first_pass_irt_errors!(ctx_fail)

    for (idx, _) in enumerate(sequences)
        @test Pioneer.getPredIrt(ctx_fail, UInt32(idx)) == library_irt[idx]
    end
end

@testset "First pass rescoring after iRT correction" begin
    canonical = collect("ACDEFGHIKLMNPQRSTVWY")
    sequences = vcat([string(aa) for aa in canonical], "WYAC")
    is_decoy = vcat(fill(false, length(canonical)), true)
    library_irt = Float32.(100.0:1.0:100.0 + length(sequences) - 1)

    params_path = Pioneer.asset_path("example_config", "defaultSearchParams.json")
    pioneer_params = Pioneer.parse_pioneer_parameters(params_path)
    first_params = Pioneer.FirstPassSearchParameters(pioneer_params)

    initial_psms = DataFrame(
        ms_file_idx = UInt32[1],
        scan_idx = UInt32[1],
        precursor_idx = UInt32[1],
        rt = Float32[100.0],
        irt_predicted = Float32[100.0],
        q_value = Float32[0.01],
        score = Float32[0.0],
        prob = Float32[0.5],
        scan_count = UInt32[1]
    )

    ctx = build_mock_context(initial_psms, sequences, is_decoy, library_irt)
    results = Pioneer.init_search_results(first_params, ctx)

    n_precursors = length(sequences)
    per_precursor = 4
    n_rows = n_precursors * per_precursor
    precursor_idxs = repeat(UInt32.(1:n_precursors), inner = per_precursor)
    scan_idxs = UInt32.(1:n_rows)
    ms_file_col = fill(UInt32(1), n_rows)
    rt_vals = Float32.(collect(1:n_rows))
    base_pred = library_irt[Int.(precursor_idxs)]

    rescoring_psms = DataFrame(
        spectral_contrast = Float32.(0.1:0.1:0.1 * n_rows),
        city_block = Float32.(0.2:0.1:0.2 + 0.1 * (n_rows - 1)),
        entropy_score = Float32.(0.3:0.1:0.3 + 0.1 * (n_rows - 1)),
        scribe = Float32.(0.4:0.1:0.4 + 0.1 * (n_rows - 1)),
        percent_theoretical_ignored = Float32.(range(0.0f0, length = n_rows, stop = 0.5f0)),
        charge2 = Float32.(rand(Float32, n_rows)),
        poisson = Float32.(rand(Float32, n_rows)),
        irt_error = Float16.(abs.(rt_vals .- base_pred)),
        missed_cleavage = Float32.(rand(Float32, n_rows)),
        Mox = Float32.(rand(Float32, n_rows)),
        TIC = Float32.(range(10f0, length = n_rows, stop = 20f0)),
        y_count = Float32.(rand(Float32, n_rows)),
        err_norm = Float32.(rand(Float32, n_rows)),
        spectrum_peak_count = Float32.(rand(Float32, n_rows) .* 10f0),
        intercept = Float32.(ones(n_rows)),
        ms_file_idx = ms_file_col,
        score = zeros(Float32, n_rows),
        precursor_idx = precursor_idxs,
        scan_idx = scan_idxs,
        q_value = fill(Float16(0.5), n_rows),
        log2_summed_intensity = Float32.(rand(Float32, n_rows) .* 20f0),
        irt = rt_vals,
        rt = rt_vals,
        irt_predicted = Float32.(base_pred),
        target = repeat(vcat(fill(true, per_precursor - 1), false), n_precursors),
        prob = zeros(Float32, n_rows)
    )

    Pioneer.cache_psms_for_rescoring!(results, ctx, rescoring_psms, 1)

    new_predictions = library_irt .+ 5f0
    for idx in 1:n_precursors
        Pioneer.setPredIrt!(ctx, UInt32(idx), new_predictions[idx])
    end

    Pioneer.rescore_first_pass_psms!(ctx, results, first_params)

    rescore_path = results.rescoring_psm_paths[1]
    updated_rescore = DataFrame(Arrow.Table(rescore_path))
    expected_preds = Float32.(new_predictions[Int.(updated_rescore.precursor_idx)])
    @test all(updated_rescore.irt_predicted .== expected_preds)
    expected_errors = Float16.(abs.(Float32.(updated_rescore.rt) .- expected_preds))
    @test all(updated_rescore.irt_error .== expected_errors)

    final_path = Pioneer.getFirstPassPsms(Pioneer.getMSData(ctx), 1)
    final_psms = DataFrame(Arrow.Table(final_path))
    @test !isempty(final_psms)
    final_expected = Float32.(new_predictions[Int.(final_psms.precursor_idx)])
    @test all(final_psms.irt_predicted .== final_expected)
end
