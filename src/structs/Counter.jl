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

mutable struct Counter{I,C<:Unsigned}
    ids::Vector{I}
    counts::Vector{C}
    size::Int64
    matches::Int64
    function Counter(I::DataType, C::DataType, size::Int) #where {I,C<:Unsigned}
        new{I, C}(zeros(I, size), zeros(C, size), 1, 0)
    end
end

getCount(c::Counter{I,C}, id::I) where {I,C<:Unsigned} = c.counts[id]
getSize(c::Counter{I,C}) where {I,C<:Unsigned} = c.size
incSize!(c::Counter{I,C}) where {I,C<:Unsigned} = c.size += 1
getID(c::Counter{I,C}, idx::Int) where {I,C<:Unsigned} = c.ids[idx]

function update!(c::Counter{I,C}, id::I, pred_intensity::C) where {I,C<:Unsigned}
    @inbounds c.counts[id] += pred_intensity;
    return nothing
end

function reset!(c::Counter{I,C}, id::I) where {I,C<:Unsigned}
    c.counts[id] = zero(T);
    return nothing
end
#=
function inc!(c::Counter{I,C}, id::I, pred_intensity::C) where {I,C<:Unsigned} 
    @inbounds @fastmath begin 
        if c.counts[id]===zero(C)
            c.ids[c.size] = id
            c.size += 1#no_previous_encounter
        end
        c.counts[id] += pred_intensity;
    end
    return nothing
end
=#
function inc!(c::Counter{I,C}, id::I, pred_intensity::C) where {I,C<:Unsigned}
    @inbounds @fastmath begin
        no_previous_encounter = c.counts[id]===zero(C)
        c.ids[c.size] = id
        c.size += no_previous_encounter
        c.counts[id] += pred_intensity;
    end
    return nothing
end

function or_inc!(c::Counter{I,C}, id::I, frag_score::C) where {I,C<:Unsigned}
    @inbounds @fastmath begin
        prev_score = c.counts[id]
        no_previous_encounter = prev_score === zero(C)
        c.ids[c.size] = id
        c.size += no_previous_encounter
        c.counts[id] = prev_score | frag_score
    end
    return nothing
end

function or_merge!(dest::Counter{I,C}, source::Counter{I,C}) where {I,C<:Unsigned}
    @inbounds @fastmath for i in 1:(getSize(source) - 1)
        id = source.ids[i]
        or_inc!(dest, id, source.counts[id])
    end
    return nothing
end

const FRAG_SCORE_WEIGHTS = UInt8[1, 1, 2, 2, 4, 4, 8]

@inline function convert_frag_score(score::C) where {C<:Unsigned}
    weighted_score = zero(C)
    @inbounds @fastmath for (idx, weight) in pairs(FRAG_SCORE_WEIGHTS)
        weighted_score += ((score >> (idx - 1)) & one(C)) * weight
    end
    return weighted_score
end

import Base.sort!
function sortCounter!(counter::Counter{I,C}) where {I,C<:Unsigned}
    sort!(
                @view(counter.ids[1:counter.matches]), 
                by = id -> counter.counts[id],
                rev = true,
                alg=PartialQuickSort(counter.matches)
        )
    return nothing
end

function reset!(c::Counter{I,C}) where {I,C<:Unsigned}
    #for i in 1:(getSize(c) - 1)
    @turbo for i in 1:(getSize(c) - 1)
        c.counts[c.ids[i]] = zero(T);
        c.ids[i] = zero(I)
    end
    c.size, c.matches = 1, 0
    return nothing
end

function countFragMatches(c::Counter{I,C}, min_count::C) where {I,C<:Unsigned}
    c.matches = 0
    @inbounds for i in 1:(getSize(c) - 1)
        id = c.ids[i]
        weighted_score = convert_frag_score(getCount(c, id))
        if weighted_score >= min_count
            c.ids[c.matches + 1] = c.ids[i]
            c.matches += 1
        end
    end
    return 0#c.matches
end

function countFragMatches(
    smoothed::Counter{I,C},
    smoothed_min::C,
    current::Counter{I,C},
    current_min::C,
)
    smoothed.matches = 0
    @inbounds for i in 1:(getSize(smoothed) - 1)
        id = smoothed.ids[i]
        smoothed_score = convert_frag_score(getCount(smoothed, id))
        if smoothed_score >= smoothed_min
            current_score = convert_frag_score(getCount(current, id))
            if current_score >= current_min
                smoothed.ids[smoothed.matches + 1] = id
                smoothed.matches += 1
            end
        end
    end
    return 0
end
