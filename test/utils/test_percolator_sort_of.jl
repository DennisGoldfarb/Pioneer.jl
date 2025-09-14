using Test
using DataFrames
using CSV
using Pioneer: get_training_data_for_iteration!, write_cv_debug_files

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

@testset "write cv fold debug files" begin
    prob2_train = Float32[0.1, 0.2]
    prob2_test = Float32[0.3, 0.4]
    prob3_train = Float32[0.5, 0.6]
    prob3_test = Float32[0.7, 0.8]
    fold_indices = Dict(UInt8(1) => [1], UInt8(2) => [2])
    train_indices = Dict(UInt8(1) => [2], UInt8(2) => [1])

    psms1 = DataFrame(
        precursor_idx = UInt32[1, 2],
        sequence = ["PEPTIDE", "PEPTIDER"],
        mods = ["", ""],
        run = ["run1", "run2"],
        charge = [2, 3],
        target = Bool[true, false],
        MBR_transfer_candidate = Bool[false, true],
    )
    mktempdir() do dir
        write_cv_debug_files(psms1, prob2_train, prob2_test,
                             prob3_train, prob3_test,
                             fold_indices, train_indices, dir)
        df_train = DataFrame(CSV.File(joinpath(dir, "cv_fold_1_train.tsv"); delim='\t'))
        df_test = DataFrame(CSV.File(joinpath(dir, "cv_fold_1_test.tsv"); delim='\t'))
        @test :sequence in names(df_train)
        @test :mods in names(df_train)
        @test df_train.prob_2nd_iteration[1] == prob2_train[2]
        @test df_test.prob_2nd_iteration[1] == prob2_test[1]
        @test df_train.prob_2nd_iteration[1] != df_test.prob_2nd_iteration[1]
    end

    psms2 = DataFrame(
        precursor_idx = UInt32[1, 2],
        structural_mods = ["", "(Oxidation[M])"],
        isotopic_mods = ["", ""],
        run = ["run1", "run2"],
        charge = [2, 3],
        target = Bool[true, false],
        MBR_transfer_candidate = Bool[false, true],
    )
    mktempdir() do dir
        write_cv_debug_files(psms2, prob2_train, prob2_test,
                             prob3_train, prob3_test,
                             fold_indices, train_indices, dir)
        df = DataFrame(CSV.File(joinpath(dir, "cv_fold_1_train.tsv"); delim='\t'))
        @test :sequence ∉ names(df)
        @test :mods in names(df)
    end
end
