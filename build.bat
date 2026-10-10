@echo off
REM ===========================================================================
REM  build.bat - build the "nightsky" observing planner with MSVC.
REM
REM  Steps:
REM    1. Generate the binary catalogue nightsky.dat from the CSV sources
REM       (catalog_raw.csv + catalog_ngc.csv) using gencat.py. Python is a
REM       BUILD-TIME tool only; nightsky.exe itself never needs it.
REM    2. Compile nightsky.cpp optimized for speed + size into nightsky.exe.
REM
REM  Ship BOTH nightsky.exe and nightsky.dat together: the program loads its
REM  catalogue only from that file. To change the object list, edit a CSV and
REM  re-run this script (or just "python gencat.py") -- no C++ recompile needed.
REM
REM  Optimization flags:
REM    /O2   maximize speed        /GL   whole-program optimization
REM    /Os   favor small code      /MT   static CRT (standalone exe, no VC DLLs)
REM    /link /OPT:REF /OPT:ICF /LTCG  strip unused code, fold identical code
REM ===========================================================================
setlocal
call "C:\Program Files\Microsoft Visual Studio\18\Community\VC\Auxiliary\Build\vcvars64.bat" >nul

echo [1/2] Generating binary catalogue nightsky.dat...
python gencat.py
if %errorlevel% neq 0 (
    echo    gencat.py failed - is Python on PATH?
    if not exist nightsky.dat (
        echo    ERROR: no nightsky.dat and generator failed. Aborting.
        exit /b 1
    )
    echo    Reusing existing nightsky.dat.
)

echo [2/2] Compiling nightsky.exe...
cl /nologo /std:c++20 /EHsc /O2 /Os /GL /MT /DNDEBUG ^
   /Fe:nightsky.exe nightsky.cpp windowsapp.lib ^
   /link /OPT:REF /OPT:ICF /LTCG
if %errorlevel% neq 0 (
    echo    Compilation FAILED.
    exit /b 1
)

REM Clean up intermediate object file.
if exist nightsky.obj del nightsky.obj

echo.
echo Build complete: nightsky.exe (+ nightsky.dat)
echo Try:  nightsky.exe 10
endlocal
