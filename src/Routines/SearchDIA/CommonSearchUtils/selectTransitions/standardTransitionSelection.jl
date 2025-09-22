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

function _select_transitions_impl!(
    transitions::Vector{DetailedFrag{Float32}},
    ::StandardTransitionSelection,
    prec_estimation_type::PrecEstimation,
    transition_idx::Int64,
    lookup::LibraryFragmentLookup,
    scan_to_prec_idx::UnitRange{Int64},
    precursors_passed_scoring::Vector{UInt32},
    prec_mzs::AbstractArray{Float32},
    prec_charges::AbstractArray{UInt8},
    prec_sulfur_counts::AbstractArray{UInt8},
    prec_irts::AbstractArray{Float32},
    iso_splines::IsotopeSplineModel,
    quad_transmission_func::QuadTransmissionFunction,
    precursor_transmission::Vector{Float32},
    isotopes::Vector{Float32},
    n_frag_isotopes::Int64,
    max_frag_rank::UInt8,
    iRT::Float32,
    iRT_tol::Float32,
    frag_mz_bounds::Tuple{Float32, Float32};
    isotope_err_bounds::Tuple{I, I} = (3, 1),
    block_size::Int64 = 10000,
    min_fraction_transmitted::Float32 = 0.0f0,
) where {I<:Integer}

    mz_low_offset = Float32(NEUTRON * first(isotope_err_bounds))
    mz_high_offset = Float32(NEUTRON * last(isotope_err_bounds))
    processed_precursors = 0

    for i in scan_to_prec_idx
        # Get precursor properties
        prec_idx = precursors_passed_scoring[i]
        prec_charge = prec_charges[prec_idx]
        prec_mz = prec_mzs[prec_idx]

        # Check retention time tolerance
        if abs(prec_irts[prec_idx] - iRT) > iRT_tol
            continue
        end

        # Handle isotope errors
        mz_low = getPrecMinBound(quad_transmission_func) - mz_low_offset / prec_charge
        mz_high = getPrecMaxBound(quad_transmission_func) + mz_high_offset / prec_charge

        if (prec_mz < mz_low) | (prec_mz > mz_high)
            continue
        end

        prec_sulfur_count = prec_sulfur_counts[prec_idx]
        prec_isotope_set = getPrecursorIsotopeSet(prec_mz, prec_charge, quad_transmission_func)
        last(prec_isotope_set) < 0 && continue

        fraction_transmitted = getPrecursorFractionTransmitted!(
            precursor_transmission,
            iso_splines,
            prec_isotope_set,
            quad_transmission_func,
            prec_mz,
            prec_charge,
            prec_sulfur_count,
        )

        fraction_transmitted < min_fraction_transmitted && continue

        processed_precursors += 1

        transition_idx = @inline fillTransitionListPrecomputed!(
            transitions,
            prec_estimation_type,
            getPrecFragRange(lookup, prec_idx),
            getFragments(lookup),
            getSplineData(lookup, prec_charge, prec_mz),
            prec_mz,
            prec_charge,
            prec_sulfur_count,
            transition_idx,
            precursor_transmission,
            prec_isotope_set,
            isotopes,
            n_frag_isotopes,
            max_frag_rank,
            iso_splines,
            frag_mz_bounds,
            block_size,
        )
    end

    return transition_idx, processed_precursors
end
