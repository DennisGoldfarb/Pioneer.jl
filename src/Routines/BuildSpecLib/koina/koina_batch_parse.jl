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
Parse results for standard intensity prediction models.
"""
function parse_koina_batch(model::InstrumentSpecificModel,
                          response::Dict{String,Any})::KoinaBatchResult{Nothing}
    df = DataFrame()
    n_precs, n_frags = first(response["outputs"])["shape"]
    
    for col in response["outputs"]
        col_name = Symbol(col["name"])
        if col_name ∈ [:intensities, :mz]
            df[!, col_name] = Float32.(col["data"])
        else
            df[!, :annotation] = string.(col["data"])
        end
    end
    
    return KoinaBatchResult(df, Int64(n_frags), nothing)
end

"""
Parse results for instrument-agnostic models.
"""
function parse_koina_batch(model::InstrumentAgnosticModel,
                          response::Dict{String,Any})::KoinaBatchResult{Nothing}
    # Currently same as InstrumentSpecificModel
    parse_koina_batch(InstrumentSpecificModel(model.name), response)
end

"""
Parse results for spline coefficient models.
"""
function parse_koina_batch(model::SplineCoefficientModel,
                          response::Dict{String,Any})::KoinaBatchResult{Vector{Float32}}
    df = DataFrame()
    n_precs, n_coef_per_frag, n_frags = response["outputs"][4]["shape"]
    knot_vector = Float32[]
    
    for col in response["outputs"]
        col_name = Symbol(col["name"])
        if col_name == :coefficients
            flat_coeffs = Float32.(col["data"])
            coefs = Vector{NTuple{n_coef_per_frag, Float32}}(undef, n_precs * n_frags)
            
            idx = 1
            for i in 1:n_precs
                for j in 1:n_frags
                    coefs[idx] = ntuple(k -> 
                        flat_coeffs[(i-1)*n_coef_per_frag*n_frags + j + (k-1)*n_frags],
                        n_coef_per_frag)
                    idx += 1
                end
            end
            df[!, :coefficients] = coefs
            
        elseif col_name == :knots
            knot_vector = Float32.(col["data"])
        elseif col_name == :mz
            df[!, col_name] = Float32.(col["data"])
        elseif col_name == :annotations
            df[!, :annotation] = Int32.(col["data"])
        end
    end
    
    return KoinaBatchResult(df, Int64(n_frags), knot_vector)
end

"""
Parse results for retention time prediction models.
"""
function parse_koina_batch(model::RetentionTimeModel,
                          response::Dict{String,Any})::KoinaBatchResult{Nothing}
    df = DataFrame()

    for col in response["outputs"]
        col_name = Symbol(col["name"])
        if col_name == :rt
            df[!, col_name] = Float32.(col["data"])
        elseif col_name == :coefficients
            n_precs, n_coef_plus_bias = col["shape"]
            flat = Float32.(col["data"])
            coef_count = n_coef_plus_bias - 1
            coefs = Vector{NTuple{coef_count, Float32}}(undef, n_precs)
            bias = Vector{Float32}(undef, n_precs)
            for i in 1:n_precs
                start = (i-1)*n_coef_plus_bias + 1
                coefs[i] = ntuple(j -> flat[start + j - 1], coef_count)
                bias[i] = flat[start + coef_count]
            end
            df[!, :coefficients] = coefs
            df[!, :bias] = bias
        end
    end

    return KoinaBatchResult(df, 1, nothing)
end