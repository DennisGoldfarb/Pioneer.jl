function main_pioneer()::Cint
    if isempty(ARGS)
        println("Usage: pioneer <subcommand> [args...]")
        println("Available subcommands: searchdia, buildspeclib, parsespeclib, getsearchparams, getbuildlibparams, convertmzml")
        return 1
    end
    cmd = lowercase(ARGS[1])
    args = ARGS[2:end]
    try
        if cmd == "searchdia"
            return main_SearchDIA(args)
        elseif cmd == "buildspeclib"
            return main_BuildSpecLib(args)
        elseif cmd == "parsespeclib"
            return main_ParseSpecLib(args)
        elseif cmd == "getsearchparams"
            return main_GetSearchParams(args)
        elseif cmd == "getbuildlibparams"
            return main_GetBuildLibParams(args)
        elseif cmd == "convertmzml"
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
