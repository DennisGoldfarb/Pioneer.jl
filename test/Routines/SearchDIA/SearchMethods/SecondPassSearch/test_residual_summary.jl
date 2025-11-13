using Test
using DataFrames
using LinearAlgebra
using Statistics
using Pioneer: ComplexScoredPSM, init_summary_columns!, get_summary_scores!, RtConversionModel,
                RESIDUAL_PROFILE_LENGTH

@testset "Residual profile summaries" begin
    K = RESIDUAL_PROFILE_LENGTH
    profiles = [
        Float32[0.45, 0.12, 0.18, 0.04, 0.06, 0.03],
        Float32[0.22, 0.35, 0.08, 0.02, 0.05, 0.07],
        Float32[0.38, 0.07, 0.16, 0.05, 0.03, 0.02]
    ]

    weights = Float32[0.4, 0.9, 0.6]
    gof = Float16[1.1, 1.5, 1.3]
    matched_ratio = Float16[0.2, 0.4, 0.3]
    fitted_manhattan = Float16[1.2, 1.6, 1.4]
    fitted_spectral = Float16[0.9, 1.1, 1.0]
    spectral_contrast = Float16[0.95, 1.05, 0.98]
    scribe = Float16[0.8, 1.0, 0.95]
    y_counts = UInt8[5, 6, 5]

    psms_vec = ComplexScoredPSM{Float32, Float16}[
        ComplexScoredPSM(
            UInt8(1), UInt8(1), UInt8(3), UInt8(3),
            UInt8(2), UInt8(2), UInt8(4), y_counts[i], UInt8(1), UInt8(0), UInt8(0),
            Float16(0.0),
            Float16(-0.5),
            Float16(0.1),
            spectral_contrast[i],
            fitted_spectral[i],
            gof[i],
            Float16(0.2),
            Float16(0.15),
            fitted_manhattan[i],
            matched_ratio[i],
            Float16(0.05),
            scribe[i],
            copy(profiles[i]),
            weights[i],
            UInt32(42),
            UInt32(7),
            UInt32(i)
        ) for i in 1:length(profiles)
    ]

    df = DataFrame(psms_vec)
    df.rt = Float32[101.0, 101.4, 100.9]
    df.best_scan = falses(nrow(df))

    init_summary_columns!(df)

    grouped = groupby(df, :precursor_idx)
    rt_model = RtConversionModel()

    @test length(grouped) == 1
    group = first(grouped)

    get_summary_scores!(
        group,
        group[!, :weight],
        group[!, :gof],
        group[!, :matched_ratio],
        group[!, :fitted_manhattan_distance],
        group[!, :fitted_spectral_contrast],
        group[!, :scribe],
        group[!, :y_count],
        rt_model
    )

    @test df.residual_profile isa Vector{Vector{Float32}}
    @test df.residual_profile[2] == profiles[2]

    apex_index = argmax(group[!, :weight])
    apex_row = group[apex_index, :]

    matrix = Matrix{Float32}(undef, length(profiles), K)
    for (row, profile) in enumerate(profiles)
        matrix[row, :] = profile
    end
    corr_matrix = cor(matrix; dims=1)
    corr_matrix = Float32.(corr_matrix)
    corr_matrix .= isfinite.(corr_matrix) .* corr_matrix

    total_pairs = K * (K - 1)
    off_diag_sum = 0.0f0
    negative_count = 0
    for i in 1:K, j in 1:K
        i == j && continue
        val = corr_matrix[i, j]
        off_diag_sum += val
        negative_count += val < 0.0f0
    end

    expected_mean = off_diag_sum / Float32(total_pairs)
    expected_neg_fraction = Float32(negative_count) / Float32(total_pairs)
    eigenvalues = eigvals(Symmetric(corr_matrix))
    expected_ratio = Float32(maximum(eigenvalues) / sum(eigenvalues))

    @test isapprox(apex_row.residual_corr_mean, expected_mean; atol=1f-6)
    @test isapprox(apex_row.residual_corr_negative_fraction, expected_neg_fraction; atol=1f-6)
    @test isapprox(apex_row.residual_corr_dom_eig_ratio, expected_ratio; atol=1f-6)
end
