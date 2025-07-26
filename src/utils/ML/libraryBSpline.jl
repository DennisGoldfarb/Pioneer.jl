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
    B(x, k, i, t)

Evaluate B-spline basis function of degree k for the i-th basis function.
"""
function B(x::T, k::Int, i::Int, t::NTuple{N,T}) where {N,T<:AbstractFloat}
    # Base case for k=0 (constant basis function)
    if k == 0
        return T(t[i] ≤ x < t[i+1])
    end
    
    # First term: uses nodes t[i] through t[i+k]
    c1 = if t[i+k] == t[i]
        zero(T)
    else
        ((x - t[i]) / (t[i+k] - t[i])) * B(x, k-1, i, t)
    end
    
    # Second term: uses nodes t[i+1] through t[i+k+1]
    c2 = if t[i+k+1] == t[i+1]
        zero(T)
    else
        ((t[i+k+1] - x) / (t[i+k+1] - t[i+1])) * B(x, k-1, i+1, t)
    end
    
    return c1 + c2
end

"""
    splevl(x, t, c, k)

Evaluate B-spline with knots t, coefficients c, and degree k at point x.
"""
function splevl(x::T, knots::NTuple{N,T}, c::NTuple{M,T}, k::Int) where {M,N,T<:AbstractFloat}
    # Number of basis functions is (length(t) - k - 1)
    n = length(knots) - k - 1
    
    # Input validation
    @assert n ≥ k+1 && length(c) ≥ n "Invalid input sizes"
    
    # Sum up the contributions from each basis function
    v = zero(T)
    for i in 1:n
        v += c[i] * B(x, k, i, knots)
    end
    return v
end

"""
    getSplineQuadrature(dtype::Type)

Get Quadrature nodes and weights as static vector to integrate B-Splines
"""
function getSplineQuadrature(dtype::Type)
    x, w = gausslegendre(20)
    return SVector{20}(dtype.(x)), SVector{20}(dtype.(w))
end

"""
    splevl(knotw, c, d, gqw, gqx)

Use Fast Gaussian Quadrature to compute the definite integral of a B-spline
over its domain
"""
function splint(knots::NTuple{N, T}, 
                c::NTuple{M, T}, 
                d::Int,
                gqx::SVector{20, T},
                gqw::SVector{20, T}) where {M,N,T<:AbstractFloat}
    i_eval = zero(T)
    for i in range(1, length(gqx))
        i_eval += splevl(gqx[i], knots, c, d)*gqw[i]
    end
    return i_eval
end

"""
    getSplineQuadrature(dtype::Type, x0::AbstractFloat, x1::AbstractFloat)

Get quadrature weights and nodes as static arrays. Defaults to 20. 
"""
function getSplineQuadrature(dtype::Type, x0::AbstractFloat, x1::AbstractFloat)
    x, w = gausslegendre(20)
    ws = (x1 - x0)/2
    for i in range(1, length(x))
        x[i] = x0 + (x[i] + 1)*(x1- x0)/2
        w[i] = ws*w[i]
    end
    return SVector{20}(dtype.(x)), SVector{20}(dtype.(w))
end

function getSplineAreas(knots::NTuple{N, T}, 
                        coefficients::AbstractVector{NTuple{M, T}},
                        d::Int,
                        gqx::SVector{20, T},
                        gqw::SVector{20, T}) where {M,N,T<:AbstractFloat}
    intensities = Vector{T}(undef, length(coefficients))
    for i in range(1, length(coefficients))
        intensities[i] = (
            splint(knots, coefficients[i], d, gqx, gqw)
        )
    end
    return intensities
end
#=
knots = (6.0f0, 13.0f0, 20.0f0, 27.0f0, 34.0f0, 41.0f0, 48.0f0, 55.0f0)
c = (1.6181915f-6, 7.382022f-6, 7.887343f-5, 0.00023642876f0)
t = 37.0f0
degree = 3
@btime splevl(t, knots, c, degree)
gqx, gqw = getSplineQuadrature(Float32, 20.0, 40.0);
@btime splint(knots, c, 3, gqx, gqw)
=#




"""
    BSplineEvaluator
Struct to hold precomputed values for efficient B-spline evaluation.
"""
struct BSplineEvaluator{N,T<:AbstractFloat}
    knots::NTuple{N,T}
    k::Int
    n::Int
    denoms1::Vector{T}
    denoms2::Vector{T}
    workspace::Vector{T}
    workspace2::Vector{T}
end

function BSplineEvaluator(knots::NTuple{N,T}, k::Int) where {N,T<:AbstractFloat}
    n = length(knots) - k - 1
    denoms1, denoms2 = precompute_denominators(knots, k)
    workspace = zeros(T, n + k)
    workspace2 = zeros(T, n + k)
    return BSplineEvaluator(knots, k, n, denoms1, denoms2, workspace, workspace2)
end


"""
    precompute_denominators(knots, k)
Precompute all denominators used in B-spline recursion since knots are fixed.
"""
function precompute_denominators(knots::NTuple{N,T}, k::Int) where {N,T<:AbstractFloat}
    n = length(knots) - k - 1
    denoms1 = Vector{T}(undef, n * k)
    denoms2 = Vector{T}(undef, n * k)
    
    idx = 1
    for deg in 1:k
        for i in 1:n
            if i + deg - 1 ≤ length(knots) - 1
                denoms1[idx] = knots[i+deg] - knots[i]
                denoms2[idx] = knots[i+deg+1] - knots[i+1]
            else
                denoms1[idx] = zero(T)
                denoms2[idx] = zero(T)
            end
            idx += 1
        end
    end
    
    return denoms1, denoms2
end


"""
    evaluate!(evaluator, x, c)
Evaluate B-spline with precomputed evaluator. Most efficient method.
"""
function evaluate!(evaluator::BSplineEvaluator{N,T}, x::T, c::NTuple{M,T}) where {N,M,T<:AbstractFloat}
    (; knots, k, n, denoms1, denoms2, workspace, workspace2) = evaluator
    
    # Initialize degree 0 basis functions
    @inbounds for i in 1:length(workspace)-1
        workspace[i] = T(knots[i] ≤ x < knots[i+1])
    end
    workspace[end] = zero(T)
    
    # Build up higher degree basis functions
    for deg in 1:k
        fill!(workspace2, zero(T))
        base_idx = (deg - 1) * n
        
        @inbounds for i in 1:n
            idx = base_idx + i
            
            if denoms1[idx] != zero(T)
                workspace2[i] += ((x - knots[i]) / denoms1[idx]) * workspace[i]
            end
            
            if denoms2[idx] != zero(T)
                workspace2[i] += ((knots[i+deg+1] - x) / denoms2[idx]) * workspace[i+1]
            end
        end
        
        workspace, workspace2 = workspace2, workspace
    end
    
    # Sum contributions
    result = zero(T)
    @inbounds for i in 1:n
        result += c[i] * workspace[i]
    end
    
    return result
end
