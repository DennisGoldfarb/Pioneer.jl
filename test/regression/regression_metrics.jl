#!/usr/bin/env julia

using Arrow
using DataFrames
using Dates
using JSON
using Printf
using Statistics
using Downloads

include("metrics_helpers.jl")
include("entrapment_metrics.jl")
include("three_proteome_metrics.jl")
include("ftr_metrics.jl")
include("metrics_pipeline.jl")
using .RegressionMetricsHelpers: condition_columns, count_nonyeast_ids, count_species_ids, count_total_ids, count_yeast_ids, gene_names_column, mean_for_columns, quant_column_names_from_proteins, select_quant_columns, species_column, unique_species_value
using .EntrapmentMetrics: compute_entrapment_metrics
using .ThreeProteomeMetrics: experimental_design_entry, experimental_design_for_dataset, fold_change_metrics_for_table, gene_counts_metrics_by_run, load_experimental_design, load_three_proteome_designs, normalize_metric_label, run_groups_for_dataset, three_proteome_design_entry

const DEFAULT_METRIC_GROUPS = ["identification", "CV", "eFDR", "runtime"]
const RELEASE_METRICS_ROOT = "/storage1/fs1/d.goldfarb/Active/Automation/Pioneer/metrics/release"
const DEVELOP_METRICS_ROOT = "/storage1/fs1/d.goldfarb/Active/Automation/Pioneer/metrics/develop"

const REPORT_FILENAME = "metrics_report.md"

function read_required_table(path::AbstractString)
    isfile(path) || error("Required file not found: $path")
    DataFrame(Arrow.Table(path))
end

function runtime_minutes_from_report(path::AbstractString)
    if !isfile(path)
        @warn "Runtime report not found; skipping runtime metric" path=path
        return nothing
    end

    for line in eachline(path)
        if startswith(line, "Total Runtime:")
            m = match(r"^Total Runtime:\s*([0-9]+(?:\.[0-9]+)?)\s+minutes", line)
            if m !== nothing
                return parse(Float64, m.captures[1])
            end
        end
    end

    @warn "Total runtime line not found in report" path=path
    nothing
end

function compute_wide_metrics(
    df::DataFrame,
    quant_col_names::AbstractVector{<:Union{Symbol, String}};
    table_label::AbstractString = "wide_table",
    dataset_name::AbstractString = "dataset",
)
    existing_quant_cols = select_quant_columns(df, quant_col_names)
    runs = length(existing_quant_cols)
    available_cols = Symbol.(names(df))

    if runs == 0
        expected = Symbol.(quant_col_names)
        @warn "No quantification columns found in dataset" dataset=dataset_name table_label=table_label expected_quant_cols=expected available_cols=available_cols expected_quant_cols_list=join(expected, ", ") available_cols_list=join(available_cols, ", ")
        return (; runs = 0, complete_rows = 0, data_completeness = 0.0)
    end

    if length(existing_quant_cols) < length(quant_col_names)
        missing_cols = setdiff(Symbol.(quant_col_names), Symbol.(existing_quant_cols))
        expected = Symbol.(quant_col_names)
        available_quant_syms = Symbol.(existing_quant_cols)
        @warn "Missing quantification columns in dataset" dataset=dataset_name table_label=table_label missing_cols=missing_cols expected_quant_cols=expected available_quant_cols=available_quant_syms available_cols=available_cols missing_cols_list=join(missing_cols, ", ") expected_quant_cols_list=join(expected, ", ") available_quant_cols_list=join(available_quant_syms, ", ") available_cols_list=join(available_cols, ", ")
    end

    quant_data = df[:, existing_quant_cols]
    quant_matrix = Matrix(quant_data)

    complete_rows = sum(all(!ismissing, row) for row in eachrow(quant_matrix))
    non_missing_values = count(!ismissing, quant_matrix)
    total_cells = nrow(df) * runs

    (; runs, complete_rows, data_completeness = total_cells > 0 ? non_missing_values / total_cells : 0.0)
end

function compute_cv_metrics(
    df::DataFrame,
    quant_col_names::AbstractVector{<:Union{Symbol, String}};
    table_label::AbstractString = "wide_table",
    groups::Dict{String, Vector{String}} = Dict{String, Vector{String}}(),
)
    function cvs_for_columns(
        columns::AbstractVector{<:Union{Symbol, String}},
        label::AbstractString,
    )
        column_syms = Symbol.(columns)
        runs = length(columns)
        if runs == 0
            return (; runs = 0, rows_evaluated = 0, cvs = Float64[])
        end

        quant_data = df[:, column_syms]
        complete_data = dropmissing(quant_data)
        rows_evaluated = nrow(complete_data)

        first_row_values = rows_evaluated == 0 ? NamedTuple() : NamedTuple(complete_data[1, :])
        @info "Computing CVs for group" table_label=table_label group_label=label quant_columns=column_syms rows_evaluated=rows_evaluated first_row_values=first_row_values
        if rows_evaluated == 0
            return (; runs, rows_evaluated, cvs = Float64[])
        end

        quant_complete = Matrix(complete_data)
        computed_cvs = Float64[]
        for row in eachrow(quant_complete)
            mean_val = mean(row)
            if mean_val != 0
                push!(computed_cvs, std(row) / mean_val)
            end
        end

        (; runs, rows_evaluated, cvs = computed_cvs)
    end

    if isempty(groups)
        existing_quant_cols = select_quant_columns(df, quant_col_names)
        stats = cvs_for_columns(existing_quant_cols, "all_runs")
        median_cv = isempty(stats.cvs) ? 0.0 : median(stats.cvs)
        return (; runs = stats.runs, rows_evaluated = stats.rows_evaluated, median_cv)
    end

    all_runs = Set{Symbol}()
    all_cvs = Float64[]
    total_rows = 0
    for (label, runs) in pairs(groups)
        columns = select_quant_columns(df, runs)
        union!(all_runs, Symbol.(columns))
        stats = cvs_for_columns(columns, label)
        total_rows += stats.rows_evaluated
        append!(all_cvs, stats.cvs)
    end

    median_cv = isempty(all_cvs) ? 0.0 : median(all_cvs)
    (; runs = length(all_runs), rows_evaluated = total_rows, median_cv)
end

function normalize_metric_groups(groups)
    if groups isa AbstractVector
        return [String(g) for g in groups]
    end

    @warn "Invalid metric groups entry; falling back to empty list" groups_type=typeof(groups)
    String[]
end

function metric_preferences(config::Dict, dataset_name::AbstractString)
    entry = get(config, dataset_name, DEFAULT_METRIC_GROUPS)

    if entry isa AbstractDict
        groups = normalize_metric_groups(get(entry, "groups", DEFAULT_METRIC_GROUPS))
        return (; groups)
    end

    groups = normalize_metric_groups(entry)
    (; groups)
end

function metric_groups_for_search(config::Dict, search_name::AbstractString)
    entry = get(config, search_name, DEFAULT_METRIC_GROUPS)

    if entry isa AbstractDict
        groups = get(entry, "groups", DEFAULT_METRIC_GROUPS)
        return normalize_metric_groups(groups)
    elseif entry isa AbstractVector
        return normalize_metric_groups(entry)
    else
        @warn "Invalid metrics entry for search; using defaults" search=search_name entry_type=typeof(entry)
        return DEFAULT_METRIC_GROUPS
    end
end

function flatten_metrics(data; prefix::AbstractString = "")
    flattened = Dict{String, Any}()

    if data isa AbstractDict
        for (key, value) in data
            key_str = String(key)
            path = isempty(prefix) ? key_str : string(prefix, "/", key_str)
            if value isa AbstractDict
                merge!(flattened, flatten_metrics(value; prefix = path))
            else
                flattened[path] = value
            end
        end
    else
        flattened[prefix] = data
    end

    flattened
end

function parse_metrics_filename(filename::AbstractString; dataset_hint::AbstractString = "", search_hint::AbstractString = "")
    if startswith(filename, "metrics_") && endswith(filename, ".json")
        base = replace(filename, r"^metrics_" => "")
        base = replace(base, r"\.json$" => "")
        if !isempty(search_hint) && endswith(base, string("_", search_hint))
            dataset = replace(base, string("_", search_hint) => "")
            return dataset, search_hint
        end
        if !isempty(dataset_hint) && startswith(base, string(dataset_hint, "_"))
            search = replace(base, string(dataset_hint, "_") => "")
            return dataset_hint, search
        end
    end

    return nothing, nothing
end

function metrics_value_is_numeric(value)
    value isa Number && isfinite(float(value))
end

function format_metric_value(value)
    if value === nothing
        return "NA"
    elseif value isa Missing
        return "NA"
    elseif value isa Integer
        return string(value)
    elseif value isa AbstractFloat
        return @sprintf("%.4f", value)
    elseif value isa Number
        return string(value)
    else
        return string(value)
    end
end

function format_metric_delta(current, previous)
    if metrics_value_is_numeric(current) && metrics_value_is_numeric(previous)
        delta = float(current) - float(previous)
        return @sprintf("%+.4f", delta)
    end
    return "NA"
end

function collect_metrics_from_versioned_root(root::AbstractString)
    metrics = Dict{String, Dict{String, Dict{String, Any}}}()
    !isdir(root) && return metrics

    for version_dir in readdir(root; join = true)
        isdir(version_dir) || continue
        version = basename(version_dir)
        version_metrics = Dict{String, Dict{String, Any}}()

        for search_dir in readdir(version_dir; join = true)
            isdir(search_dir) || continue
            search_name = basename(search_dir)

            for file in readdir(search_dir; join = true)
                isfile(file) || continue
                filename = basename(file)
                !startswith(filename, "metrics_") && continue
                !endswith(filename, ".json") && continue

                dataset, search = parse_metrics_filename(filename; search_hint = search_name)
                if dataset === nothing || search === nothing
                    continue
                end

                raw = JSON.parsefile(file)
                flattened = flatten_metrics(raw)
                dataset_metrics = get!(version_metrics, dataset, Dict{String, Any}())
                dataset_metrics[search] = flattened
            end
        end

        metrics[version] = version_metrics
    end

    metrics
end

function collect_metrics_from_search_root(root::AbstractString; version_label::AbstractString = "develop")
    metrics = Dict{String, Dict{String, Dict{String, Any}}}()
    !isdir(root) && return metrics

    version_metrics = Dict{String, Dict{String, Any}}()
    for search_dir in readdir(root; join = true)
        isdir(search_dir) || continue
        search_name = basename(search_dir)

        for file in readdir(search_dir; join = true)
            isfile(file) || continue
            filename = basename(file)
            !startswith(filename, "metrics_") && continue
            !endswith(filename, ".json") && continue

            dataset, search = parse_metrics_filename(filename; search_hint = search_name)
            if dataset === nothing || search === nothing
                continue
            end

            raw = JSON.parsefile(file)
            flattened = flatten_metrics(raw)
            dataset_metrics = get!(version_metrics, dataset, Dict{String, Any}())
            dataset_metrics[search] = flattened
        end
    end

    metrics[version_label] = version_metrics
    metrics
end

function collect_metrics_from_results_root(root::AbstractString; version_label::AbstractString)
    metrics = Dict{String, Dict{String, Dict{String, Any}}}()
    !isdir(root) && return metrics

    version_metrics = Dict{String, Dict{String, Any}}()
    for dataset_dir in readdir(root; join = true)
        isdir(dataset_dir) || continue
        dataset_name = basename(dataset_dir)

        for file in readdir(dataset_dir; join = true)
            isfile(file) || continue
            filename = basename(file)
            !startswith(filename, "metrics_") && continue
            !endswith(filename, ".json") && continue

            dataset, search = parse_metrics_filename(filename; dataset_hint = dataset_name)
            if dataset === nothing || search === nothing
                continue
            end

            raw = JSON.parsefile(file)
            flattened = flatten_metrics(raw)
            dataset_metrics = get!(version_metrics, dataset, Dict{String, Any}())
            dataset_metrics[search] = flattened
        end
    end

    metrics[version_label] = version_metrics
    metrics
end

function collect_metrics_from_results_dirs(results_dirs::AbstractVector{<:AbstractString}; version_label::AbstractString)
    metrics = Dict{String, Dict{String, Dict{String, Any}}}()
    version_metrics = Dict{String, Dict{String, Any}}()

    for results_dir in results_dirs
        isdir(results_dir) || continue
        dataset_name = dataset_name_from_results(results_dir)
        for file in readdir(results_dir; join = true)
            isfile(file) || continue
            filename = basename(file)
            !startswith(filename, "metrics_") && continue
            !endswith(filename, ".json") && continue

            dataset, search = parse_metrics_filename(filename; dataset_hint = dataset_name)
            if dataset === nothing || search === nothing
                continue
            end

            raw = JSON.parsefile(file)
            flattened = flatten_metrics(raw)
            dataset_metrics = get!(version_metrics, dataset, Dict{String, Any}())
            dataset_metrics[search] = flattened
        end
    end

    metrics[version_label] = version_metrics
    metrics
end

function parse_version_tuple(version::AbstractString)
    cleaned = replace(version, r"^v" => "")
    parts = split(cleaned, ".")
    parsed = map(part -> tryparse(Int, part), parts)
    if any(x -> x === nothing, parsed)
        return (typemax(Int),)
    end
    return Tuple(parsed)
end

function ordered_release_versions(metrics::Dict{String, Dict{String, Dict{String, Any}}})
    versions = collect(keys(metrics))
    sort(versions; by = parse_version_tuple)
end

function build_metrics_report(metrics_by_version::Dict{String, Dict{String, Dict{String, Any}}}, version_order::Vector{String})
    datasets = Set{String}()
    searches = Dict{String, Set{String}}()

    for (version, dataset_metrics) in metrics_by_version
        for (dataset, search_metrics) in dataset_metrics
            push!(datasets, dataset)
            dataset_searches = get!(searches, dataset, Set{String}())
            for search in keys(search_metrics)
                push!(dataset_searches, search)
            end
        end
    end

    report = IOBuffer()
    println(report, "# Regression Metrics Report")
    println(report, "")
    println(report, "- Generated: $(Dates.format(now(), dateformat\"yyyy-mm-dd HH:MM:SS\"))")
    println(report, "- Versions: $(join(version_order, ", "))")
    println(report, "")

    for dataset in sort(collect(datasets))
        dataset_searches = get(searches, dataset, Set{String}())
        for search in sort(collect(dataset_searches))
            println(report, "## Dataset: $(dataset) — Search: $(search)")
            println(report, "")

            metric_keys = Set{String}()
            for version in version_order
                dataset_metrics = get(metrics_by_version, version, Dict{String, Dict{String, Any}}())
                search_metrics = get(get(dataset_metrics, dataset, Dict{String, Any}()), search, Dict{String, Any}())
                for key in keys(search_metrics)
                    push!(metric_keys, key)
                end
            end

            header_cells = ["Metric"]
            for (idx, version) in enumerate(version_order)
                push!(header_cells, version)
                idx > 1 && push!(header_cells, "Δ")
            end
            println(report, "| ", join(header_cells, " | "), " |")
            println(report, "|", join(fill("---", length(header_cells)), "|"), "|")

            for metric in sort(collect(metric_keys))
                row_cells = String[metric]
                previous_value = nothing
                for (idx, version) in enumerate(version_order)
                    dataset_metrics = get(metrics_by_version, version, Dict{String, Dict{String, Any}}())
                    search_metrics = get(get(dataset_metrics, dataset, Dict{String, Any}()), search, Dict{String, Any}())
                    value = get(search_metrics, metric, nothing)
                    push!(row_cells, format_metric_value(value))
                    if idx > 1
                        push!(row_cells, format_metric_delta(value, previous_value))
                    end
                    previous_value = value
                end
                println(report, "| ", join(row_cells, " | "), " |")
            end

            println(report, "")
        end
    end

    String(take!(report))
end

function current_commit_label()
    sha = get(ENV, "GITHUB_SHA", "")
    if isempty(sha)
        try
            sha = readchomp(`git rev-parse --short HEAD`)
        catch err
            @warn "Unable to resolve git commit SHA for report label" error=err
            sha = "current"
        end
    else
        sha = first(sha, min(length(sha), 7))
    end

    "commit $(sha)"
end

function write_report(report_text::AbstractString, results_root::AbstractString)
    mkpath(results_root)
    report_path = joinpath(results_root, REPORT_FILENAME)
    open(report_path, "w") do io
        write(io, report_text)
    end
    report_path
end

function post_report_to_pr(report_text::AbstractString)
    event_path = get(ENV, "GITHUB_EVENT_PATH", "")
    isempty(event_path) && return false
    !isfile(event_path) && return false

    token = get(ENV, "GITHUB_TOKEN", "")
    isempty(token) && return false

    repo = get(ENV, "GITHUB_REPOSITORY", "")
    isempty(repo) && return false

    event = JSON.parsefile(event_path)
    pull_request = get(event, "pull_request", nothing)
    pull_request === nothing && return false
    pr_number = get(pull_request, "number", nothing)
    pr_number === nothing && return false

    body = """
    ## Regression Metrics Report

    <details>
    <summary>View report</summary>

    $(report_text)
    </details>
    """

    url = "https://api.github.com/repos/$(repo)/issues/$(pr_number)/comments"
    headers = Dict(
        "Authorization" => "token $(token)",
        "Accept" => "application/vnd.github+json",
        "User-Agent" => "PioneerMetricsReporter",
    )

    response = Downloads.request(
        url,
        "POST",
        headers,
        JSON.json(Dict("body" => body)),
    )

    response.status in 200:299
end

function generate_regression_report(results_root::AbstractString; current_metrics::Dict{String, Dict{String, Dict{String, Any}}} = Dict{String, Dict{String, Dict{String, Any}}}())
    release_metrics = collect_metrics_from_versioned_root(RELEASE_METRICS_ROOT)
    develop_metrics = collect_metrics_from_search_root(DEVELOP_METRICS_ROOT; version_label = "develop")

    commit_label = current_commit_label()
    version_order = ordered_release_versions(release_metrics)
    push!(version_order, "develop")
    push!(version_order, commit_label)

    metrics_by_version = Dict{String, Dict{String, Dict{String, Any}}}()
    merge!(metrics_by_version, release_metrics)
    metrics_by_version["develop"] = get(develop_metrics, "develop", Dict{String, Dict{String, Any}}())

    if isempty(current_metrics)
        current_metrics = collect_metrics_from_results_root(results_root; version_label = commit_label)
    end
    metrics_by_version[commit_label] = get(current_metrics, commit_label, Dict{String, Dict{String, Any}}())

    report_text = build_metrics_report(metrics_by_version, version_order)
    report_path = write_report(report_text, results_root)
    posted = post_report_to_pr(report_text)
    @info "Generated regression metrics report" report_path=report_path posted_to_pr=posted
end

function dataset_dirs_from_root(root::AbstractString)
    !isdir(root) && return String[]
    filter(readdir(root; join=true)) do path
        isdir(path)
    end
end

function resolve_results_root()
    env_root = get(ENV, "PIONEER_RESULTS_DIR", "")
    arg_root = length(ARGS) >= 1 ? ARGS[1] : ""

    if !isempty(env_root)
        return env_root
    elseif !isempty(arg_root)
        return arg_root
    end

    error("Results directory is not specified; set PIONEER_RESULTS_DIR or pass a path argument")
end

function load_metrics_config(path::AbstractString)
    if isempty(path)
        @info "No metrics config provided; using defaults for all searches"
        return Dict{String, Any}()
    end

    if !isfile(path)
        @warn "Metrics config not found; using defaults for all searches" metrics_config_path=path
        return Dict{String, Any}()
    end

    try
        raw = JSON.parsefile(path)
        if raw isa AbstractDict
            return Dict{String, Any}(String(k) => v for (k, v) in raw)
        end
        @warn "Metrics config is not a dictionary; using defaults" metrics_config_path=path
        return Dict{String, Any}()
    catch err
        @warn "Failed to parse metrics config; using defaults" metrics_config_path=path error=err
        return Dict{String, Any}()
    end
end

function search_param_files(params_dir::AbstractString)
    !isdir(params_dir) && return String[]
    filter(readdir(params_dir; join=true)) do path
        isfile(path) && occursin(r"search.*\.json$", basename(path))
    end
end

function results_dir_from_param(path::AbstractString)
    try
        parsed = JSON.parsefile(path)
        results = get(parsed, "results", nothing)
        if results isa AbstractString
            return results
        end

        paths_block = get(parsed, "paths", nothing)
        if paths_block isa AbstractDict
            nested_results = get(paths_block, "results", nothing)
            if nested_results isa AbstractString
                return nested_results
            end
        end
    catch err
        @warn "Failed to parse search param file" param_path=path error=err
        return nothing
    end

    @warn "Search param file missing results entry" param_path=path
    nothing
end

function dataset_name_from_results(results_dir::AbstractString)
    parent_dir = dirname(results_dir)
    if isempty(parent_dir) || parent_dir == results_dir
        return basename(results_dir)
    end

    basename(parent_dir)
end

function output_metrics_path(results_dir::AbstractString, dataset_name::AbstractString, search_name::AbstractString)
    filename = string("metrics_", dataset_name, "_", search_name, ".json")
    joinpath(results_dir, filename)
end

function is_log_file(path::AbstractString)
    lowercasepath = lowercase(path)
    endswith(lowercasepath, ".log") || endswith(lowercasepath, ".log.gz")
end

function cleanup_entrapment_dir(entrapment_dir::AbstractString)
    isdir(entrapment_dir) || return

    for entry in readdir(entrapment_dir; join = true)
        if isfile(entry)
            endswith(entry, ".png") && continue
            rm(entry; force = true, recursive = true)
        else
            rm(entry; force = true, recursive = true)
        end
    end
end

function cleanup_results_dir(results_dir::AbstractString, metrics_path::AbstractString)
    if !isdir(results_dir)
        @warn "Results directory missing during cleanup" results_dir = results_dir
        return
    end

    metrics_path_abs = abspath(metrics_path)

    for entry in readdir(results_dir; join = true)
        entry_abs = abspath(entry)
        if entry_abs == metrics_path_abs
            continue
        end

        if isfile(entry)
            if is_log_file(entry)
                continue
            end
        elseif isdir(entry) && basename(entry) == "entrapment_analysis"
            cleanup_entrapment_dir(entry)
            continue
        end

        rm(entry; force = true, recursive = true)
    end
end

function archive_results(
    results_dir::AbstractString,
    metrics_path::AbstractString;
    archive_root::AbstractString = "",
    dataset_name::AbstractString = "dataset",
    search_name::AbstractString = "search",
)
    if isempty(archive_root)
        @info "Archive root not provided; skipping move of regression outputs" results_dir=results_dir
        return
    end

    target_dir = joinpath(archive_root, "results", dataset_name)
    mkpath(target_dir)

    if isfile(metrics_path)
        target_metrics = joinpath(target_dir, basename(metrics_path))
        mv(metrics_path, target_metrics; force = true)
        metrics_path = target_metrics
    end

    entrapment_dir = joinpath(results_dir, "entrapment_analysis")
    if isdir(entrapment_dir)
        for entry in readdir(entrapment_dir; join = true)
            if isfile(entry) && endswith(lowercase(entry), ".png")
                mv(entry, joinpath(target_dir, basename(entry)); force = true)
            end
        end
    end

    for entry in readdir(results_dir; join = true)
        if isfile(entry) && is_log_file(entry)
            mv(entry, joinpath(target_dir, basename(entry)); force = true)
        end
    end

    try
        rm(results_dir; force = true, recursive = true)
    catch err
        @warn "Failed to remove original results directory after archiving" results_dir=results_dir error=err
    end

    @info "Archived regression outputs" results_dir=results_dir metrics_path=metrics_path target_dir=target_dir
end

function compute_metrics_for_params_dir(
    params_dir::AbstractString;
    metrics_config_path::AbstractString = "",
    experimental_design_path::AbstractString = "",
    three_proteome_designs_path::AbstractString = "",
    archive_root::AbstractString = "",
)
    isdir(params_dir) || error("Params directory does not exist: $params_dir")

    param_files = search_param_files(params_dir)
    isempty(param_files) && error("No search*.json files found in params directory $params_dir")

    metrics_config = load_metrics_config(metrics_config_path)
    experimental_design = isempty(experimental_design_path) ? Dict{String, Any}() : load_experimental_design(experimental_design_path)

    dataset_entries = NamedTuple{(:search_name, :results_dir, :dataset_name)}[]
    dataset_paths = Dict{String, String}()

    for param_file in param_files
        results_dir = results_dir_from_param(param_file)
        results_dir === nothing && error("Missing results entry in search param file $param_file")

        dataset_name = dataset_name_from_results(results_dir)
        search_name = replace(basename(param_file), ".json" => "")
        push!(dataset_entries, (; search_name, results_dir, dataset_name))
        haskey(dataset_paths, dataset_name) || (dataset_paths[dataset_name] = results_dir)
    end

    isempty(dataset_entries) && error("No valid parameter files remain after parsing results paths")

    three_proteome_designs = nothing

    for entry in dataset_entries
        metric_groups = metric_groups_for_search(metrics_config, entry.search_name)
        normalized_groups = Set(replace.(lowercase.(metric_groups), "-" => "_"))
        need_three_proteome = ("fold_change" in normalized_groups) || ("three_proteome" in normalized_groups)

        if need_three_proteome && three_proteome_designs === nothing && !isempty(three_proteome_designs_path)
            three_proteome_designs = load_three_proteome_designs(three_proteome_designs_path)
        end

        metrics = compute_dataset_metrics(
            entry.results_dir,
            entry.dataset_name;
            metric_groups = metric_groups,
            experimental_design = experimental_design,
            three_proteome_designs = three_proteome_designs,
            dataset_paths = dataset_paths,
        )

        output_path = output_metrics_path(entry.results_dir, entry.dataset_name, entry.search_name)

        metrics === nothing && begin
            @warn "No metrics produced for search; archiving logs only" search=entry.search_name dataset=entry.dataset_name results_dir=entry.results_dir
            cleanup_results_dir(entry.results_dir, output_path)
            archive_results(
                entry.results_dir,
                output_path;
                archive_root = archive_root,
                dataset_name = entry.dataset_name,
                search_name = entry.search_name,
            )
            continue
        end

        if metrics isa AbstractDict && isempty(metrics)
            @warn "No metrics produced for search; archiving logs only" search=entry.search_name dataset=entry.dataset_name results_dir=entry.results_dir
            cleanup_results_dir(entry.results_dir, output_path)
            archive_results(
                entry.results_dir,
                output_path;
                archive_root = archive_root,
                dataset_name = entry.dataset_name,
                search_name = entry.search_name,
            )
            continue
        end

        open(output_path, "w") do io
            JSON.print(io, metrics)
        end

        cleanup_results_dir(entry.results_dir, output_path)
        archive_results(
            entry.results_dir,
            output_path;
            archive_root = archive_root,
            dataset_name = entry.dataset_name,
            search_name = entry.search_name,
        )
    end

    report_root = if !isempty(archive_root)
        joinpath(archive_root, "results")
    else
        first_results_dir = first(dataset_entries).results_dir
        dataset_dir = dirname(first_results_dir)
        dirname(dataset_dir)
    end

    current_metrics = if !isempty(archive_root)
        collect_metrics_from_results_root(report_root; version_label = current_commit_label())
    else
        unique_results_dirs = unique([entry.results_dir for entry in dataset_entries])
        collect_metrics_from_results_dirs(unique_results_dirs; version_label = current_commit_label())
    end

    generate_regression_report(report_root; current_metrics = current_metrics)
end

function main()
    params_dir_override = get(ENV, "PIONEER_PARAMS_DIR", "")
    if !isempty(params_dir_override)
        metrics_config_path = get(ENV, "PIONEER_METRICS_FILE", "")
        experimental_design_path = get(ENV, "PIONEER_EXPERIMENTAL_DESIGN", "")
        three_proteome_designs_path = get(ENV, "PIONEER_THREE_PROTEOME_DESIGNS", "")
        archive_root = get(ENV, "PIONEER_ARCHIVE_ROOT", "")

        compute_metrics_for_params_dir(
            params_dir_override;
            metrics_config_path = metrics_config_path,
            experimental_design_path = experimental_design_path,
            three_proteome_designs_path = three_proteome_designs_path,
            archive_root = archive_root,
        )
        return
    end

    results_dir = resolve_results_root()
    isdir(results_dir) || error("Results directory does not exist: $results_dir")
    dataset_dirs = dataset_dirs_from_root(results_dir)
    isempty(dataset_dirs) && error("No dataset directories found in $results_dir")

    dataset_paths = Dict{String, String}(basename(path) => path for path in dataset_dirs)

    @info "Using regression results directory" results_dir=results_dir

    metrics_config_path = get(
        ENV,
        "PIONEER_METRICS_CONFIG",
        joinpath(@__DIR__, "..", "..", "pioneer-regression-configs", "metrics_config.json"),
    )
    metric_group_config = if isfile(metrics_config_path)
        JSON.parsefile(metrics_config_path)
    else
        @info "No metrics_config.json found; using default metric groups" metrics_config_path=metrics_config_path
        Dict{String, Any}()
    end

    experimental_design_path = get(
        ENV,
        "PIONEER_EXPERIMENTAL_DESIGN",
        joinpath(@__DIR__, "..", "..", "pioneer-regression-configs", "experimental_designs"),
    )
    experimental_design = load_experimental_design(experimental_design_path)

    three_proteome_designs_path = get(
        ENV,
        "PIONEER_THREE_PROTEOME_DESIGNS",
        joinpath(@__DIR__, "..", "..", "pioneer-regression-configs", "experimental_designs"),
    )
    three_proteome_designs = nothing

    dataset_dirs = filter(dataset_dirs) do path
        dataset_name = basename(path)
        preferences = metric_preferences(metric_group_config, dataset_name)
        if isempty(preferences.groups)
            @info "Skipping dataset without requested metrics" dataset=dataset_name
            return false
        end
        return true
    end

    isempty(dataset_dirs) && error("No dataset directories remain after filtering")

    for dataset_dir in dataset_dirs
        dataset_name = basename(dataset_dir)

        preferences = metric_preferences(metric_group_config, dataset_name)

        metric_groups = preferences.groups
        normalized_groups = Set(replace.(lowercase.(metric_groups), "-" => "_"))
        need_three_proteome = ("fold_change" in normalized_groups) || ("three_proteome" in normalized_groups)
        if need_three_proteome && three_proteome_designs === nothing
            three_proteome_designs = load_three_proteome_designs(three_proteome_designs_path)
        end

        metrics = compute_dataset_metrics(
            dataset_dir,
            dataset_name;
            metric_groups = metric_groups,
            experimental_design = experimental_design,
            three_proteome_designs = three_proteome_designs,
            dataset_paths = dataset_paths,
        )
        metrics === nothing && continue

        output_path = joinpath(dataset_dir, "metrics_$(dataset_name).json")
        open(output_path, "w") do io
            JSON.print(io, metrics)
        end
    end

    generate_regression_report(results_dir)
end

main()
