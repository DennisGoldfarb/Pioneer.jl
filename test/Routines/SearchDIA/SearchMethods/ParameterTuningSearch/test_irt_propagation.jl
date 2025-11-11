using Test
using Pioneer
using DataFrames

struct MockSearchStructures <: Pioneer.SearchDataStructures end
struct MockMassSpecReference <: Pioneer.MassSpecDataReference end

struct MockPrecursors
    is_decoy::Vector{Bool}
    irt::Vector{Float32}
    charge::Vector{UInt8}
end

struct MockSpectralLibrary <: Pioneer.SpectralLibrary
    precursors::MockPrecursors
end

struct MockMassSpecData <: Pioneer.MassSpecData
    retention_times::Vector{Float32}
    tics::Vector{Float32}
end

Pioneer.getPrecursors(lib::MockSpectralLibrary) = lib.precursors
Pioneer.getIsDecoy(p::MockPrecursors) = p.is_decoy
Pioneer.getIrt(p::MockPrecursors) = p.irt
Pioneer.getCharge(p::MockPrecursors) = p.charge
Pioneer.getRetentionTimes(ms::MockMassSpecData) = ms.retention_times
Pioneer.getTICs(ms::MockMassSpecData) = ms.tics

function build_mock_context(precursors::MockPrecursors)
    library = MockSpectralLibrary(precursors)
    temp_structures = [MockSearchStructures()]
    ms_reference = MockMassSpecReference()
    Pioneer.SearchContext(library, temp_structures, ms_reference, 1, length(precursors.irt), 1)
end

function make_psm_table()
    DataFrame(
        precursor_idx = UInt32[1, 2, 3],
        scan_idx = UInt32[1, 2, 3],
        matched_ratio = Float16[10, 20, 30]
    )
end

function make_mock_spectra()
    MockMassSpecData(
        Float32[100.0, 200.0, 300.0],
        Float32[1.0e5, 2.0e5, 3.0e5]
    )
end

@testset "Parameter tuning uses corrected iRT predictions" begin
    precursors = MockPrecursors(
        [false, true, false],
        Float32[50.0, 60.0, 70.0],
        UInt8[2, 3, 2]
    )
    spectra = make_mock_spectra()
    psms = make_psm_table()
    context = build_mock_context(precursors)

    Pioneer.setPredIrt!(context, UInt32(1), 45.5f0)
    Pioneer.setPredIrt!(context, UInt32(2), 65.25f0)
    Pioneer.setPredIrt!(context, UInt32(3), 80.0f0)

    Pioneer.add_tuning_search_columns!(
        psms,
        context,
        spectra,
        Pioneer.getIsDecoy(precursors),
        Pioneer.getIrt(precursors),
        Pioneer.getCharge(precursors),
        Pioneer.getRetentionTimes(spectra),
        Pioneer.getTICs(spectra)
    )

    @test psms[!, :irt_predicted] ≈ Float32[45.5, 65.25, 80.0]
    @test psms[!, :target] == .!precursors.is_decoy
end

@testset "Parameter tuning falls back to library iRT when corrections missing" begin
    precursors = MockPrecursors(
        [false, false, true],
        Float32[10.0, 20.0, 30.0],
        UInt8[2, 2, 2]
    )
    spectra = make_mock_spectra()
    psms = make_psm_table()
    context = build_mock_context(precursors)

    Pioneer.setPredIrt!(context, UInt32(2), 18.0f0)

    Pioneer.add_tuning_search_columns!(
        psms,
        context,
        spectra,
        Pioneer.getIsDecoy(precursors),
        Pioneer.getIrt(precursors),
        Pioneer.getCharge(precursors),
        Pioneer.getRetentionTimes(spectra),
        Pioneer.getTICs(spectra)
    )

    expected = Float32[10.0, 18.0, 30.0]
    @test psms[!, :irt_predicted] ≈ expected
end
