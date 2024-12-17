"""
    process_psms!(psms::DataFrame, spectra::Arrow.Table, 
                 search_context::SearchContext, params::NceTuningSearchParameters)

Process and filter peptide-spectrum matches (PSMs) from search results to identify high-confidence matches.

# Arguments
- `psms`: DataFrame containing raw PSM results
- `spectra`: MS/MS spectral data
- `search_context`: Contains spectral library and search context
- `params`: NCE tuning search parameters including filtering thresholds

# Process
1. Adds analysis columns including target/decoy status, retention times, and charge states
2. Scores PSMs using presearch scoring model
3. Calculates q-values for false discovery rate control
4. Selects best PSM per precursor/scan combination based on SCRIBE score
5. Filters to retain only:
   - Target PSMs (non-decoy)
   - PSMs below q-value threshold
   - Best scoring PSM per precursor
6. Adds precursor m/z information

# Returns
Filtered DataFrame containing only high-confidence PSMs that passed all criteria.

# Side Effects
Modifies input `psms` DataFrame by adding columns and filtering rows.
"""
function process_psms!(
    psms::DataFrame,
    spectra::Arrow.Table,
    search_context::SearchContext,
    params::NceTuningSearchParameters
)
    # Add columns
    precursors = getPrecursors(getSpecLib(search_context))
    add_tuning_search_columns!(
        psms,
        spectra,
        getIsDecoy(precursors),
        getIrt(precursors),
        getCharge(precursors),
        spectra[:retentionTime],
        spectra[:TIC]
    )

    # Score PSMs
    score_presearch!(psms)
    get_qvalues!(psms[!,:prob], psms[!,:target], psms[!,:q_value])

    # Get best PSMs per precursor/scan
    spsms = combine(groupby(psms, [:precursor_idx, :scan_idx])) do group
        max_idx = argmax(group[!, :scribe])
        return group[max_idx:max_idx, :]
    end

    # Apply filters
    filter!(row -> row.target && row.q_value <= params.max_q_val, spsms)
    passing_precs = Set(spsms[!, :precursor_idx])
    filter!(row -> row.precursor_idx ∈ passing_precs, psms)

    # Add precursor info
    psms[!, :prec_mz] = [
        getMz(precursors)[pid]
        for pid in psms[!, :precursor_idx]
    ]
    
    # Select best PSMs
    psms[!, :best_psms] .= false
    for group in groupby(psms, :precursor_idx)
        best_idx = argmax(group[!, :scribe])
        group[best_idx, :best_psms] = true
    end
    filter!(row -> row.best_psms, psms)

    return psms
end

