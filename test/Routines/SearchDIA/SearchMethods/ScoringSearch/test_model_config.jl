using Test
using Pioneer

@testset "ScoringSearch model configurations" begin
    configs = Pioneer.create_model_configurations()

    # Expect five models after addition of simplified probit
    @test length(configs) == 5

    names = getfield.(configs, :name)
    @test "ProbitRegression" in names
    @test "ProbitRegressionSimple" in names

    # Intercept should appear only in probit model feature sets
    model_by_name = Dict(name => cfg for (name, cfg) in zip(names, configs))
    @test :intercept in model_by_name["ProbitRegression"].features
    @test :intercept in model_by_name["ProbitRegressionSimple"].features
    @test :intercept ∉ model_by_name["SimpleLightGBM"].features
    @test :intercept ∉ model_by_name["AdvancedLightGBM"].features
    @test :intercept ∉ model_by_name["SuperSimplified"].features

    weight_trend_features = [
        :weight_scribe_slope,
        :weight_scribe_rho,
        :weight_fitted_spectral_contrast_slope,
        :weight_fitted_spectral_contrast_rho,
        :weight_matched_ratio_slope,
        :weight_matched_ratio_rho,
    ]

    for feature in weight_trend_features
        @test feature in model_by_name["SimpleLightGBM"].features
        @test feature in model_by_name["AdvancedLightGBM"].features
    end
end
