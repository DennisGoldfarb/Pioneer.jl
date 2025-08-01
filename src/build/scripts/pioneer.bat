@echo off
setlocal enabledelayedexpansion

if "%JULIA_NUM_THREADS%"=="" set JULIA_NUM_THREADS=auto

set SCRIPT_DIR=%~dp0
set CLI=%SCRIPT_DIR%bin\PioneerCLI.exe
set ARGS=

:parse
if "%~1"=="" goto run
if "%~1"=="--threads" (
    if "%~2"=="" (
        echo Error: --threads requires a value
        exit /b 1
    )
    set JULIA_NUM_THREADS=%~2
    shift
    shift
    goto parse
)
for /f "tokens=1,2 delims==" %%A in ("%~1") do (
    if "%%A"=="--threads" (
        set JULIA_NUM_THREADS=%%B
        shift
        goto parse
    )
)
if defined ARGS (
    set ARGS=%ARGS% %~1
) else (
    set ARGS=%~1
)
shift
goto parse

:run
if defined ARGS (
    "%CLI%" %ARGS%
) else (
    "%CLI%"
)
