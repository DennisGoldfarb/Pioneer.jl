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
    enableExponentialSearchLogging(path::AbstractString; append::Bool=false)

Enable logging of exponential search diagnostics to a TSV file located at `path`.
When `append=true`, appends to an existing file (writing the header only if the
file is empty). Returns the canonicalized path to the log file.
"""
function enableExponentialSearchLogging(path::AbstractString; append::Bool=false)
    mkpath(dirname(path))
    mode = append ? "a" : "w"
    io = open(path, mode)
    file_nonempty = append && isfile(path) && filesize(path) > 0
    logger = ExponentialSearchLogger(io, abspath(path), file_nonempty)
    lock(EXP_SEARCH_LOG_LOCK) do
        if EXP_SEARCH_LOGGER[] !== nothing
            close(EXP_SEARCH_LOGGER[].io)
        end
        EXP_SEARCH_LOGGER[] = logger
        if !file_nonempty
            println(io, join(EXP_SEARCH_LOG_HEADER, '\t'))
            flush(io)
            logger.header_written = true
        end
    end
    return logger.path
end

"""
    disableExponentialSearchLogging()

Close the exponential search diagnostics log file if it is open and disable logging.
"""
function disableExponentialSearchLogging()
    lock(EXP_SEARCH_LOG_LOCK) do
        if EXP_SEARCH_LOGGER[] !== nothing
            close(EXP_SEARCH_LOGGER[].io)
            EXP_SEARCH_LOGGER[] = nothing
        end
    end
    return nothing
end

"""
    isExponentialSearchLoggingEnabled() -> Bool

Return `true` if exponential search diagnostics logging is currently enabled.
"""
isExponentialSearchLoggingEnabled() = EXP_SEARCH_LOGGER[] !== nothing

const EXP_SEARCH_LOG_HEADER = (
    "initial_step_size",
    "prev_lower_index",
    "prev_upper_index",
    "final_lower_index",
    "final_upper_index",
    "distance_from_prev_lower",
    "distance_from_prev_upper",
    "exp_search_steps",
    "total_binary_steps",
    "window_start_index",
    "window_stop_index",
    "fragments_examined",
    "bins_examined",
    "target_mz",
    "peak_density",
    "peaks_remaining",
    "total_peaks"
)

mutable struct ExponentialSearchLogger
    io::IO
    path::String
    header_written::Bool
end

mutable struct ExponentialSearchTrace
    prev_lower::UInt32
    prev_upper::UInt32
    step_size::UInt32
    exp_steps::UInt16
    final_lower::UInt32
    final_upper::UInt32
    window_start_idx::UInt32
    window_stop_idx::UInt32
    fragments_examined::UInt32
    bins_examined::UInt32
    target_mz::Float32
    peak_density::Float32
    peaks_remaining::Int32
    total_peaks::Int32
    binary_steps::UInt32
end

function ExponentialSearchTrace(prev_lower::UInt32, prev_upper::UInt32, step_size::UInt32,
                                target_mz::Float32,
                                peak_density::Float32, peaks_remaining::Int32,
                                total_peaks::Int32)
    ExponentialSearchTrace(
        prev_lower,
        prev_upper,
        step_size,
        UInt16(0),
        prev_lower,
        prev_upper,
        UInt32(0),
        UInt32(0),
        UInt32(0),
        UInt32(0),
        target_mz,
        peak_density,
        peaks_remaining,
        total_peaks,
        UInt32(0)
    )
end

const EXP_SEARCH_LOG_LOCK = ReentrantLock()
const EXP_SEARCH_LOGGER = Ref{Union{Nothing, ExponentialSearchLogger}}(nothing)

function log_exponential_search_metrics(trace::ExponentialSearchTrace)
    logger = EXP_SEARCH_LOGGER[]
    logger === nothing && return nothing
    lock(EXP_SEARCH_LOG_LOCK) do
        io = logger.io
        if !logger.header_written
            println(io, join(EXP_SEARCH_LOG_HEADER, '\t'))
            logger.header_written = true
        end
        row = (
            trace.step_size,
            trace.prev_lower,
            trace.prev_upper,
            trace.final_lower,
            trace.final_upper,
            Int64(trace.final_lower) - Int64(trace.prev_lower),
            Int64(trace.final_upper) - Int64(trace.prev_upper),
            trace.exp_steps,
            trace.binary_steps,
            trace.window_start_idx,
            trace.window_stop_idx,
            trace.fragments_examined,
            trace.bins_examined,
            trace.target_mz,
            trace.peak_density,
            trace.peaks_remaining,
            trace.total_peaks
        )
        println(io, join(row, '\t'))
        flush(io)
    end
    return nothing
end
