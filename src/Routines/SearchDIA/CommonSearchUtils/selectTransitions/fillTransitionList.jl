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
    fillTransitionList!(transitions, prec_estimation_type, args...) -> Int64

Fill a transition list with fragment ions and their isotopes based on the specified
precursor estimation strategy.

# Arguments
- `transitions::Vector{DetailedFrag{Float32}}`: Vector to store generated transitions
- `prec_estimation_type::PrecEstimation`: Strategy for isotope pattern calculation
- `precursor_fragment_range::UnitRange{UInt64}`: Range of fragments to process
- `fragment_ions::Vector{AltimeterFragment}`: The library fragment ions
- `nce::Union{Missing, Float32}`: Normalized collision energy (if applicable)
- `knots::Union{Missing, NTuple{M, Float32}}`: Spline knots for intensity prediction
- `prec_mz::Float32`: Precursor m/z
- `prec_charge::UInt8`: Precursor charge
- `prec_sulfur_count::UInt8`: Number of sulfur atoms in precursor
- `transition_idx::Int64`: Current transition index
- `quad_transmission_func::QuadTransmissionFunction`: Quadrupole transmission function
- `precursor_transmission::Vector{Float32}`: Buffer for precursor transmission values
- `isotopes::Vector{Float32}`: Buffer for isotope calculations
- `n_frag_isotopes::Int64`: Number of fragment isotopes to consider
- `max_frag_rank::UInt8`: Maximum fragment rank to include
- `frag_iso_cutoff::Float32`: Fraction of the most intense isotope across the
  precursor below which fragment isotopes are discarded
- `iso_splines::IsotopeSplineModel`: Model for isotope pattern prediction
- `frag_mz_bounds::Tuple{Float32, Float32}`: m/z bounds for fragments
- `block_size::Int64`: Size for array growth when needed

# Returns
- `Int64`: Updated transition index
"""
function fillTransitionList!(transitions::Vector{DetailedFrag{Float32}}, 
                            prec_estimation_type::PrecEstimation,
                            precursor_fragment_range::UnitRange{UInt64},
                            fragment_ions::Vector{F},
                            spline_data::G,
                            prec_mz::Float32,
                            prec_charge::UInt8,
                            prec_sulfur_count::UInt8,
                            transition_idx::Int64, 
                            quad_transmission_func::QuadTransmissionFunction,
                            precursor_transmission::Vector{Float32},
                            isotopes::Vector{Float32},
                            n_frag_isotopes::Int64,
                            max_frag_rank::UInt8,
                            frag_iso_cutoff::Float32,
                            iso_splines::IsotopeSplineModel,
                            frag_mz_bounds::Tuple{Float32, Float32},
                            block_size::Int64)::Int64 where {G<:IntensityDataType, F <: AltimeterFragment}#where {T,U,V,W<:AbstractFloat,I<:Integer}

    # Calculate precursor isotope transmission and range
    getPrecursorIsotopeTransmission!(precursor_transmission, prec_mz, prec_charge, quad_transmission_func)
    prec_isotope_set = getPrecursorIsotopeSet(prec_mz, prec_charge, quad_transmission_func)
    frag_iso_idx_range = range(0, min(n_frag_isotopes - 1, last(prec_isotope_set)))

    # First pass: collect isotope intensities for each fragment while tracking
    # the current global maximum. Stop generating isotopes early when they fall
    # below the cutoff relative to this running maximum.
    global_max_intensity = zero(Float32)
    frag_iso_data = Vector{Tuple{UInt64, Vector{Float32}}}()
    for frag_idx in precursor_fragment_range
        frag = fragment_ions[frag_idx]
        frag.rank > max_frag_rank && continue

        frag_isos = Float32[]
        global_max_intensity = getFragIsotopes!(
            prec_estimation_type,
            isotopes,
            precursor_transmission,
            prec_isotope_set,
            frag_iso_idx_range,
            iso_splines,
            prec_mz,
            prec_charge,
            prec_sulfur_count,
            frag,
            spline_data,
            global_max_intensity,
            frag_iso_cutoff,
            frag_isos,
        )
        !isempty(frag_isos) && push!(frag_iso_data, (frag_idx, frag_isos))
    end

    # Second pass: emit transitions using a cutoff relative to the global maximum
    for (frag_idx, frag_isos) in frag_iso_data
        frag = fragment_ions[frag_idx]
        iso_range = 0:(length(frag_isos) - 1)
        transition_idx = addTransitionIsotopes!(transitions, transition_idx,
                                                frag, frag_isos, iso_range,
                                                frag_mz_bounds, block_size,
                                                frag_iso_cutoff, global_max_intensity)
    end

    return transition_idx
end


"""
Helper function to add transitions for each isotope of a fragment.
"""
function addTransitionIsotopes!(transitions::Vector{DetailedFrag{Float32}},
                                transition_idx::Int64,
                                frag::AltimeterFragment,
                                isotopes::Vector{Float32},
                                frag_iso_idx_range::UnitRange{Int64},
                                frag_mz_bounds::Tuple{Float32, Float32},
                                block_size::Int64,
                                frag_iso_cutoff::Float32,
                                global_max_intensity::Float32)::Int64
    for iso_idx in frag_iso_idx_range
        intensity = isotopes[iso_idx + 1]
        intensity < frag_iso_cutoff * global_max_intensity && break

        frag_mz = Float32(frag.mz + iso_idx * NEUTRON/frag.frag_charge)

        # Skip if outside m/z bounds
        (frag_mz < first(frag_mz_bounds) || frag_mz > last(frag_mz_bounds)) && continue

        transition_idx += 1
        transitions[transition_idx] = DetailedFrag(
            frag.prec_id, frag_mz, Float16(intensity),
            frag.ion_type, frag.is_y, frag.is_b, frag.is_p, iso_idx > 0,
            frag.frag_charge, frag.ion_position, frag.prec_charge,
            frag.rank, frag.sulfur_count
        )

        ensureTransitionCapacity!(transitions, transition_idx, block_size)
    end
    return transition_idx
end

function getFragIsotopes!(
                            ::PartialPrecCapture,
                            tmp_isotopes::Vector{Float32},
                            precursor_transmition::Vector{Float32},
                            prec_isotope_set::Tuple{Int64, Int64},
                            frag_iso_idx_range::UnitRange{Int64},
                            iso_splines::IsotopeSplineModel,
                            prec_mz::Float32,
                            prec_charge::UInt8,
                            prec_sulfur_count::UInt8,
                            frag::F,
                            spline_data::G,
                            current_max::Float32,
                            frag_iso_cutoff::Float32,
                            out::Vector{Float32}) where {F<:AltimeterFragment, G<:IntensityDataType}
    #Reset relative abundances of isotopes to zero
    fill!(tmp_isotopes, zero(eltype(tmp_isotopes)))
    #Predicted total fragment ion intensity (sum of fragment isotopes)
    total_fragment_intensity =  getIntensity(frag, spline_data)

    getFragAbundance!(
                    tmp_isotopes,
                    precursor_transmition,
                    iso_splines,
                    prec_mz,
                    prec_charge,
                    prec_sulfur_count,
                    frag
                    )

    threshold = frag_iso_cutoff * current_max
    for iso_idx in frag_iso_idx_range
        intensity = total_fragment_intensity * tmp_isotopes[iso_idx + 1]
        intensity < threshold && break
        push!(out, intensity)
        current_max = max(current_max, intensity)
        threshold = frag_iso_cutoff * current_max
    end
    return current_max
end

function getFragIsotopes!(
                            ::FullPrecCapture,
                            tmp_isotopes::Vector{Float32},
                            precursor_transmition::Vector{Float32},
                            prec_isotope_set::Tuple{Int64, Int64},
                            frag_iso_idx_range::UnitRange{Int64},
                            iso_splines::IsotopeSplineModel,
                            prec_mz::Float32,
                            prec_charge::UInt8,
                            prec_sulfur_count::UInt8,
                            frag::F,
                            spline_data::G,
                            current_max::Float32,
                            frag_iso_cutoff::Float32,
                            out::Vector{Float32}) where {F<:AltimeterFragment, G<:IntensityDataType}

    #Predicted total fragment ion intensity (sum of fragment isotopes)
    total_fragment_intensity = getIntensity(frag, spline_data)
    frag_mz = getMz(frag)
    frag_charge = getPrecCharge(frag)
    frag_nsulfur = Int64(getSulfurCount(frag))
    threshold = frag_iso_cutoff * current_max
    for iso_idx in frag_iso_idx_range
        intensity = iso_splines(
                                min(frag_nsulfur, 5),
                                iso_idx,
                                frag_mz*frag_charge
                                )*total_fragment_intensity
        intensity < threshold && break
        push!(out, intensity)
        current_max = max(current_max, intensity)
        threshold = frag_iso_cutoff * current_max
    end
    return current_max
end

# Compatibility methods without early stopping
function getFragIsotopes!(
                            capture::PartialPrecCapture,
                            frag_isotopes::Vector{Float32},
                            precursor_transmition::Vector{Float32},
                            prec_isotope_set::Tuple{Int64, Int64},
                            frag_iso_idx_range::UnitRange{Int64},
                            iso_splines::IsotopeSplineModel,
                            prec_mz::Float32,
                            prec_charge::UInt8,
                            prec_sulfur_count::UInt8,
                            frag::F,
                            spline_data::G) where {F<:AltimeterFragment, G<:IntensityDataType}
    tmp = similar(frag_isotopes, length(frag_iso_idx_range))
    empty!(frag_isotopes)
    getFragIsotopes!(capture, tmp, precursor_transmition, prec_isotope_set,
                     frag_iso_idx_range, iso_splines, prec_mz, prec_charge,
                     prec_sulfur_count, frag, spline_data, 0f0, 0f0,
                     frag_isotopes)
    return nothing
end

function getFragIsotopes!(
                            capture::FullPrecCapture,
                            frag_isotopes::Vector{Float32},
                            precursor_transmition::Vector{Float32},
                            prec_isotope_set::Tuple{Int64, Int64},
                            frag_iso_idx_range::UnitRange{Int64},
                            iso_splines::IsotopeSplineModel,
                            prec_mz::Float32,
                            prec_charge::UInt8,
                            prec_sulfur_count::UInt8,
                            frag::F,
                            spline_data::G) where {F<:AltimeterFragment, G<:IntensityDataType}
    tmp = similar(frag_isotopes, length(frag_iso_idx_range))
    empty!(frag_isotopes)
    getFragIsotopes!(capture, tmp, precursor_transmition, prec_isotope_set,
                     frag_iso_idx_range, iso_splines, prec_mz, prec_charge,
                     prec_sulfur_count, frag, spline_data, 0f0, 0f0,
                     frag_isotopes)
    return nothing
end

