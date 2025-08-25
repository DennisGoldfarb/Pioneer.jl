using Test, DataFrames

@testset "Sundial parse" begin
    response = Dict(
        "outputs" => [
            Dict(
                "name" => "coefficients",
                "datatype" => "FP32",
                "shape" => [2, 5],
                "data" => Float32[1,2,3,4,5,6,7,8,9,10]
            ),
        ]
    )
    result = parse_koina_batch(RetentionTimeModel("sundial"), response)
    @test size(result.fragments, 1) == 2
    @test result.fragments.bias == Float32[5,10]
    @test result.fragments.coefficients[1] == (1f0,2f0,3f0,4f0)
end
