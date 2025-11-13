using Test
using DataFrames
using Pioneer
using Pioneer: SparseArray, computeFittedMetricsFor, init_summary_columns!, get_summary_scores!, RtConversionModel

@testset "Fragment compensation log-likelihood" begin
    H = SparseArray{Int64, Float32}(2)
    H.n_vals = 2
    H.m = 2
    H.n = 1
    H.rowval = Int64[1, 2]
    H.colval = UInt16[1, 1]
    H.nzval = Float32[0.6, 0.4]
    H.matched = Bool[true, false]
    H.isotope = UInt8[0, 0]
    H.x = Float32[0.5, 0.3]
    H.colptr = Int64[1, 3]

    w = Float32[1.0]
    r = Float32[0.1, 0.1]
    included = [1, 2]

    metrics = computeFittedMetricsFor(w, H, r, 1, included)
    (_, _, _, _, _, _, _, _, frag_loglik, _, _, _) = metrics

    fitted_peaks = Float32.(w[1] .* H.nzval[included])
    shadow_peaks = fitted_peaks .- r[H.rowval[included]]
    total_counts = sum(shadow_peaks)
    fitted_norm = fitted_peaks ./ sum(fitted_peaks)
    shadow_norm = shadow_peaks ./ total_counts
    prior_strength = Float32(Pioneer.FRAG_COMP_PRIOR_STRENGTH)
    prior_offset = Float32(Pioneer.FRAG_COMP_PRIOR_OFFSET)
    alpha = prior_strength .* fitted_norm .+ prior_offset
    alpha0 = sum(alpha)
    log_comb = lgamma(total_counts + 1f0) - sum(lgamma.(shadow_peaks .+ 1f0))
    log_prior = lgamma(alpha0) - lgamma(total_counts + alpha0)
    log_post = sum(lgamma.(shadow_peaks .+ alpha) .- lgamma.(alpha))
    expected_loglik = log_comb + log_prior + log_post

    @test isapprox(frag_loglik, expected_loglik; atol=1f-3)

    psms = DataFrame(
        weight = Float32[0.2, 0.5, 0.3],
        gof = Float16[1, 2, 3],
        matched_ratio = Float16[-1, 0, 1],
        frag_compensation_loglik = Float16.([-0.5, -0.2, -0.3]),
        fitted_manhattan_distance = Float16[-2, -1, -0.5],
        fitted_spectral_contrast = Float16[0.1, 0.2, 0.3],
        scribe = Float16[0.2, 0.4, 0.3],
        y_count = UInt8[2, 3, 1],
        rt = Float32[10.0, 10.1, 10.2]
    )
    psms.best_scan = falses(3)

    init_summary_columns!(psms)

    get_summary_scores!(
        SubDataFrame(psms, 1:3),
        psms.weight,
        psms.gof,
        psms.matched_ratio,
        psms.frag_compensation_loglik,
        psms.fitted_manhattan_distance,
        psms.fitted_spectral_contrast,
        psms.scribe,
        psms.y_count,
        RtConversionModel()
    )

    apex = argmax(psms.weight)
    vals = Float32.(psms.frag_compensation_loglik)
    expected_mean = sum(vals) / length(vals)
    expected_std = sqrt(max(sum(vals .^ 2) / length(vals) - expected_mean^2, 0f0))

    @test psms.best_scan[apex]
    @test psms.mean_frag_compensation_loglik[apex] ≈ expected_mean atol=1f-6
    @test psms.std_frag_compensation_loglik[apex] ≈ expected_std atol=1f-6
    @test all(psms.mean_frag_compensation_loglik[setdiff(1:3, apex)] .== 0f0)
    @test all(psms.std_frag_compensation_loglik[setdiff(1:3, apex)] .== 0f0)
end
