using Test
using Pioneer
using Arrow
using DataFrames

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
