using Test
using DataFrames
using Pioneer: get_training_data_for_iteration!

@testset "filter paired precursors" begin
    psms = DataFrame(
        prob = Float32[0.9, 0.8, 0.1],
        target = Bool[true, true, false],
        q_value = Float32[0.0, 0.0, 1.0],
        MBR_is_best_decoy = Bool[false, false, false],
        MBR_max_pair_prob = Float32[0.9, 0.8, 0.1],
        MBR_is_missing = Bool[false, false, false],
        MBR_is_paired = Bool[false, true, false]
    )

    res1 = get_training_data_for_iteration!(psms, 1, true, 1.0f0, 0.2f0, 2.0f0, false)
    @test nrow(res1) == 2

    res2 = get_training_data_for_iteration!(psms, 2, true, 1.0f0, 0.2f0, 2.0f0, false)
    @test nrow(res2) == 2

    res_last = get_training_data_for_iteration!(psms, 2, true, 1.0f0, 0.2f0, 2.0f0, true)
    @test nrow(res_last) == 3
end
