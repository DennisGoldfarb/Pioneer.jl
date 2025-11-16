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

const DEBUG_PRECURSOR_CACHE = Base.RefValue{Union{Nothing, UInt32}}(nothing)
const DEBUG_PRECURSOR_CACHE_INITIALIZED = Base.RefValue(false)

"""
    _parse_debug_precursor_idx()

Parse the `PIONEER_DEBUG_PRECURSOR_IDX` environment variable (if present) into a
`UInt32`. Returns `nothing` when unset or invalid.
"""
function _parse_debug_precursor_idx()
    raw = get(ENV, "PIONEER_DEBUG_PRECURSOR_IDX", nothing)
    if raw === nothing
        return nothing
    end
    stripped = strip(raw)
    isempty(stripped) && return nothing
    parsed = tryparse(UInt64, stripped)
    parsed === nothing && return nothing
    return UInt32(parsed)
end

"""
    get_debug_precursor_idx()

Return the cached debug precursor index (if configured).
"""
function get_debug_precursor_idx()
    if !DEBUG_PRECURSOR_CACHE_INITIALIZED[]
        DEBUG_PRECURSOR_CACHE[] = _parse_debug_precursor_idx()
        DEBUG_PRECURSOR_CACHE_INITIALIZED[] = true
    end
    return DEBUG_PRECURSOR_CACHE[]
end

"""
    should_trace_precursor(precursor_idx)

Return `true` if `precursor_idx` matches the configured debug precursor.
"""
should_trace_precursor(precursor_idx::Integer) = begin
    debug_idx = get_debug_precursor_idx()
    debug_idx === nothing && return false
    return UInt32(precursor_idx) == debug_idx
end

"""
    log_precursor_trace(stage; kwargs...)
    log_precursor_trace(stage, precursor_idx; kwargs...)

Emit a formatted log message when tracing is enabled for the configured
precursor. Additional keyword arguments are appended as key/value pairs for
context.
"""
function log_precursor_trace(stage::AbstractString; kwargs...)
    debug_idx = get_debug_precursor_idx()
    debug_idx === nothing && return
    log_precursor_trace(stage, debug_idx; kwargs...)
end

function log_precursor_trace(stage::AbstractString, precursor_idx::Integer; kwargs...)
    should_trace_precursor(precursor_idx) || return
    io = IOBuffer()
    print(io, "[precursor-trace] stage=", stage, " precursor_idx=", UInt32(precursor_idx))
    for (name, value) in kwargs
        value === nothing && continue
        print(io, " ", name, "=")
        if value isa AbstractFloat
            print(io, round(Float64(value); digits=6))
        else
            print(io, value)
        end
    end
    @user_info String(take!(io))
end
