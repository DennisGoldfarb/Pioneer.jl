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

"""
    initialize_prob_group_features!(psms, match_between_runs)

Ensure the probability- and match-between-runs-related columns exist on the
provided PSM table.

The columns are created if they are missing and reset to zero/false so that
subsequent rescoring passes start from a known state.
"""
function initialize_prob_group_features!(
    psms::AbstractDataFrame,
    match_between_runs::Bool,
)
    n = nrow(psms)
    psms[!, :prob] = zeros(Float32, n)
    psms[!, :q_value] = zeros(Float64, n)

    if match_between_runs
        psms[!, :MBR_max_pair_prob] = zeros(Float32, n)
        psms[!, :MBR_best_irt_diff] = zeros(Float32, n)
        psms[!, :MBR_log2_weight_ratio] = zeros(Float32, n)
        psms[!, :MBR_log2_explained_ratio] = zeros(Float32, n)
        psms[!, :MBR_rv_coefficient] = zeros(Float32, n)
        psms[!, :MBR_is_best_decoy] = trues(n)
        psms[!, :MBR_num_runs] = zeros(Int32, n)
        psms[!, :MBR_transfer_candidate] = falses(n)
        psms[!, :MBR_is_missing] = falses(n)
    end

    return psms
end
