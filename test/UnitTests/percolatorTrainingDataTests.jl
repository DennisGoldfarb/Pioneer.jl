# Copyright (C) 2024 Nathan Wamsley
#
# This file is part of Pioneer.jl
#
# Pioneer.jl is free software: you can redistribute it and/or modify
# it under the terms of the GNU Affero General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU Affero General Public License for more details.
#
# You should have received a copy of the GNU Affero General Public License
# along with this program. If not, see <https://www.gnu.org/licenses/>.

using Test
using DataFrames
using Pioneer: get_training_data_for_iteration!

@testset "training data filtering" begin
    psms = DataFrame(
        prob = Float32[0.9, 0.8, 0.7, 0.6],
        target = Bool[true, false, true, false],
        q_value = zeros(Float32, 4),
        MBR_is_best_decoy = falses(4),
        MBR_max_pair_prob = Float32[0.9, 0.9, 0.9, 0.9],
        MBR_is_missing = falses(4),
        passed_first_search = Bool[true, true, false, false]
    )

    iter1 = get_training_data_for_iteration!(psms, 1, true, 0.01f0, 0.2f0, 0.9f0, false)
    @test all(iter1.passed_first_search)
    @test nrow(iter1) == 2

    iter2 = get_training_data_for_iteration!(psms, 2, true, 0.01f0, 0.2f0, 0.9f0, false)
    @test all(iter2.passed_first_search)
    @test nrow(iter2) == 2

    iter3 = get_training_data_for_iteration!(psms, 3, true, 0.01f0, 0.2f0, 0.9f0, true)
    @test nrow(iter3) == 4
end
