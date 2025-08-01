function main_pioneer(argv=ARGS)::Cint
    s = ArgParseSettings()
    s.prog = "pioneer"
    s.description = "Pioneer - Mass Spectrometry Data Analysis"
    s.epilog = "Subcommands:\n" *
        "  search\n" *
        "  predict\n" *
        "  empirical\n" *
        "  search-config\n" *
        "  predict-config\n" *
        "  empirical-config\n" *
        "  convert-mzml\n\n" *
        "For subcommand-specific help: pioneer <subcommand> --help"
    @add_arg_table s begin
        "--threads"
            help = "Set number of Julia threads"
            arg_type = Int
            required = false
    end

    parsed, extra = parse_args(argv, s; as_symbols=true, skip_extra_args=true)

    if haskey(parsed, :threads) && parsed[:threads] !== nothing
        ENV["JULIA_NUM_THREADS"] = string(parsed[:threads])
    elseif get(ENV, "JULIA_NUM_THREADS", "") == ""
        ENV["JULIA_NUM_THREADS"] = "auto"
    end

    if isempty(extra)
        println("Subcommand required\n")
        print_usage(s)
        return 1
    end

    cmd = lowercase(first(extra))
    args = extra[2:end]
    try
        if cmd == "search"
            return main_SearchDIA(args)
        elseif cmd == "predict"
            return main_BuildSpecLib(args)
        elseif cmd == "empirical"
            return main_ParseSpecLib(args)
        elseif cmd == "search-config"
            return main_GetSearchParams(args)
        elseif cmd == "predict-config"
            return main_GetBuildLibParams(args)
        elseif cmd == "empirical-config"
            return main_GetParseSpecLibParams(args)
        elseif cmd == "convert-mzml"
            return main_convertMzML(args)
        else
            println("Unknown subcommand: $cmd")
            return 1
        end
    catch
        Base.invokelatest(Base.display_error, Base.catch_stack())
        return 1
    end
    return 0
end
