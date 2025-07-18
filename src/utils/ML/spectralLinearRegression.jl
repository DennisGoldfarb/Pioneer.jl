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

abstract type RegularizationType end
struct L1Norm <: RegularizationType end
struct L2Norm <: RegularizationType end
struct NoNorm <: RegularizationType end

function getRegL1(λ::T, xk::T, ::NoNorm) where T<:AbstractFloat
    return zero(Float32)
end

function getRegL1(λ::T, xk::T, ::L1Norm) where T<:AbstractFloat
    return λ
end

function getRegL1(λ::T, xk::T, ::L2Norm) where T<:AbstractFloat
    return λ*Float32(2)*xk
end


function getRegL2(λ::T, xk::T, ::NoNorm) where T<:AbstractFloat
    return zero(Float32)
end

function getRegL2(λ::T, xk::T, ::L1Norm) where T<:AbstractFloat
    return zero(Float32)
end

function getRegL2(λ::T, xk::T, ::L2Norm) where T<:AbstractFloat
    return Float32(2)*λ
end


function updateResiduals!(Hs::SparseArray{Ti, T}, r::Vector{T}, col::Int64, X1, X0) where {Ti<:Integer, T<:AbstractFloat}
    @inbounds @fastmath for i in Hs.colptr[col]:(Hs.colptr[col + 1] - 1)
        row_val = Hs.rowval[i]
        nz_val = Hs.nzval[i]
        r[row_val] += nz_val*(X1 - X0)
    end
end

function getDerivatives!(Hs::SparseArray{Ti, T}, 
                            r::Vector{Float32}, 
                            col::Int64, 
                            δ::Float32, 
                            λ::Float32,
                            xk::Float32,
                            regularization_type::RegularizationType
                            ) where {Ti<:Integer,T<:AbstractFloat}
    L1 = zero(Float32)
    L2 = zero(Float32)
    @inbounds @fastmath for i in Hs.colptr[col]:(Hs.colptr[col + 1] - 1)
        rval = r[Hs.rowval[i]]
        hsval = Hs.nzval[i]
        RS = (1 + (rval/δ)^2)
        # Quake's Fast Inverse Square Root Algorithm 
        R = RS
        int32 = reinterpret(UInt32, R)
        int32 = 0x5f3759df - int32 >> 1 
        R = reinterpret(Float32, int32)
        R *= 1.5f0 - RS * 0.5f0 * R^2
        HSVAL_R = hsval * R

        L1 += HSVAL_R * rval
        L2 += (hsval) * HSVAL_R * ((R)^(2))
    end

    return Float32(L1) + getRegL1(λ, xk, regularization_type), Float32(L2) + getRegL2(λ, xk, regularization_type)
end

function getL1(Hs::SparseArray{Ti, Float32}, 
                r::Vector{Float32}, 
                col::Int64, 
                δ::Float32, 
                λ::Float32,
                xk::Float32,
                regularization_type::RegularizationType) where {Ti<:Integer}
    L1 = zero(Float32)
    @inbounds @fastmath for i in Hs.colptr[col]:(Hs.colptr[col + 1] - 1)
        rval = r[Hs.rowval[i]]
        hsval = Hs.nzval[i]
        RS = (1 + (rval/δ)^2)
        # Quake's Fast Inverse Square Root Algorithm 
        R = RS
        int32 = reinterpret(UInt32, R)
        int32 = 0x5f3759df - int32 >> 1 #Magic 
        R = reinterpret(Float32, int32)
        R *= 1.5f0 - RS * 0.5f0 * R^2
        L1 += rval * hsval * R
    end
    return Float32(L1) + getRegL1(λ, xk, regularization_type)
end

function newton_bisection!(Hs::SparseArray{Ti, T}, 
                            r::Vector{T}, 
                            X₁::Vector{T}, 
                            col::Int64, 
                            δ::T, 
                            λ::T,
                            max_iter_newton::Int64, 
                            max_iter_bisection::Int64,
                            accuracy_newton::T,
                            accuracy_bisection::T,
                            regularization_type::RegularizationType,
                            rel_tol::T = T(0.01)) where {Ti<:Integer,T<:AbstractFloat}
    n = 0
    X_init = X₁[col]  # Initial estimate before optimization
    X0 = X₁[col]      # Previous estimate at each iteration
    _ranbisection_ = false  # Flag to indicate if bisection was used
    # Track maximum estimates for bisection bounds if Newton fails
    max_l1, max_x1 = typemax(T), typemax(T)

    @inbounds begin
        # Newton-Raphson iterations
        while (n < max_iter_newton)
            # Compute first and second derivatives
            L1, L2 = getDerivatives!(Hs, r, col, δ, λ, X₁[col], regularization_type)
            update_rule = (L1)/L2
            
            # Check for numerical issues
            if isnan(update_rule)
                n = max_iter_newton
                break
            end

            # Track positive L1 with smallest absolute value for bisection bounds
            if (sign(L1) == 1) & (L1 < max_l1)
                max_x1, max_l1 = X₁[col], L1
            end

            # Newton-Raphson update (constrained to be non-negative)
            X0 = X₁[col] 
            X₁[col] = max(X₁[col] - update_rule, zero(T))
            n += 1

            # Update residuals
            updateResiduals!(Hs, r, col, X₁[col], X0)

            # Check convergence
            abs_change = abs(X₁[col] - X0)
            
            if !iszero(X0)  # Can check relative change
                rel_change = abs_change / abs(X0)
                if rel_change < rel_tol
                    break
                end
            else
                # For zero weights, use absolute tolerance only
                if abs_change < accuracy_newton
                    break
                end
            end
        end
        
        # If Newton's method failed to converge, use bisection
        if n == max_iter_newton
            _ranbisection_ = true
            # Reset to zero (lower bound)
            X0 = X₁[col]
            X₁[col] = zero(T)
            updateResiduals!(Hs, r, col, X₁[col], X0)
            L1 = getL1(Hs, r, col, δ, λ, X₁[col], regularization_type)
            
            # Only use bisection if minimum is at X₁[col] > 0
            if sign(L1) != 1
                _ = bisection!(Hs, r, X₁, col, δ, λ, zero(T), 
                            min(max(max_x1, zero(Float32)), Float32(1e11)),
                            L1,  
                            max_iter_bisection,
                            accuracy_bisection,
                            regularization_type)
            end
            return X₁[col] - X_init
        else
            # Convergence reached
            return X₁[col] - X_init
        end
    end
end

function bisection!(Hs::SparseArray{Ti, T}, 
                    r::Vector{T}, 
                    X₁::Vector{T},
                    col::Int64, 
                    δ::T, 
                    λ::T, 
                    a::T, 
                    b::T, 
                    fa::Float32,
                    max_iter::Int64, 
                    accuracy_bisection::T,
                    regularization_type::RegularizationType) where {Ti<:Integer,T<:AbstractFloat}
    n = 0
    c = (a + b)/2
    # Update residuals for new X₁[col] value 
    updateResiduals!(Hs, r, col, c, X₁[col])
    X0 = X₁[col]
    X₁[col] = c
    X_init, X0 = X₁[col],  X₁[col]
    while (n < max_iter)

        # Evaluate first partial derivative
        fc = getL1(Hs, r, col, δ, λ, X₁[col], regularization_type)
        # Bisection Rule
        if (sign(fc) != sign(fa))
            b, fb = c, fc
        else
            a, fa = c, fc
        end

        c, X0 = (a + b)/2, X₁[col]
        X₁[col] = c
        # Update residuals
        updateResiduals!(Hs, r, col, X₁[col], X0)

        # Stopping Criterion
        abs(X₁[col] - X0) < accuracy_bisection ? break : nothing
        n += 1
    end
    return X₁[col] - X_init
end

function solveHuber!(Hs::SparseArray{Ti, T}, 
                        r::Vector{T}, 
                        X₁::Vector{T}, 
                        δ::T,
                        λ::T,
                        max_iter_newton::Int64, 
                        max_iter_bisection::Int64, 
                        max_iter_outer::Int64,
                        accuracy_newton::T,
                        accuracy_bisection::T,
                        relative_convergence_threshold::T,  # max_diff from config: max relative change in weights
                        regularization_type::RegularizationType) where {Ti<:Integer,T<:AbstractFloat}
    
    # Recommended values for solver parameters:
    # max_iter_newton = 25
    # max_iter_bisection = 100
    # max_iter_outer = max(1000, Hs.n*5)
    
    # Use user-provided convergence threshold as relative tolerance for Newton's method
    newton_rel_tol = relative_convergence_threshold
    
    # Initialize iteration counter
    i = 0
    while i < max_iter_outer
        _diff = T(0)
        for col in range(1, Hs.n)
            
            # Update coefficient
            δx = abs(newton_bisection!(Hs, r, X₁, col, δ, λ,
                                        max_iter_newton, 
                                        max_iter_bisection,
                                        accuracy_newton,
                                        accuracy_bisection,
                                        regularization_type,
                                        newton_rel_tol))

            # Track maximum relative change
            if !iszero(X₁[col])
                rel_change = δx / abs(X₁[col])
                if rel_change > _diff
                    _diff = rel_change
                end
            end
        end
        
        # Check convergence
        if _diff < relative_convergence_threshold
            break
        end  
        i += 1
    end
    
    return nothing
end