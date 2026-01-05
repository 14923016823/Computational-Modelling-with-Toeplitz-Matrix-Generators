@echo off
:: --- STEP 1: AUTO-ELEVATE TO ADMIN ---
net session >nul 2>&1
if %errorLevel% == 0 (
    goto :admin_ok
) else (
    echo Requesting Administrator privileges...
    powershell -Command "Start-Process -FilePath '%0' -Verb RunAs"
    exit /b
)

:admin_ok
setlocal enabledelayedexpansion

:: --- STEP 2: USER INPUT FOR PATH ---
echo ===================================================
echo   Intel oneAPI Roofline Automator (Full Path Mode)
echo ===================================================
echo Example: C:\Users\Name\Documents\Project\MVM.cpp
set /p FULL_CPP_PATH="Enter the FULL path to your .cpp file: "

:: Extract the Directory, Filename, and Extension from the path
for %%i in ("%FULL_CPP_PATH%") do (
    set "SOURCE_DIR=%%~dpi"
    set "FILE_NAME=%%~ni"
)

:: Change drive and directory to where the source file is
cd /d "%SOURCE_DIR%"

set EXE_NAME=%FILE_NAME%.exe
set PROJ_DIR=.\Advisor_Results_%FILE_NAME%

:: --- STEP 3: INITIALIZE VISUAL STUDIO & ONEAPI ---
echo [1/3] Initializing Environments...
set "VS2022INSTALLDIR=C:\Program Files (x86)\Microsoft Visual Studio\2022\BuildTools"

:: Run the Intel setvars script
call "C:\Program Files (x86)\Intel\oneAPI\setvars.bat" intel64 vs2022

:: --- STEP 4: COMPILE WITH ICPX ---
echo.
echo [2/3] Compiling %FILE_NAME%.cpp from %SOURCE_DIR%...
icpx -O3 -xHost -g "%FULL_CPP_PATH%" -o "%EXE_NAME%"

if %errorLevel% neq 0 (
    echo.
    echo [ERROR] Compilation failed.
    pause
    exit /b
)

:: --- STEP 5: RUN ADVISOR ROOFLINE ---
echo.
echo [3/3] Running Intel Advisor Roofline Analysis...
:: We run advisor from the source directory so the results folder is created there
advisor --collect=roofline --project-dir=%PROJ_DIR% -- .\%EXE_NAME%

echo.
echo ===================================================
echo ANALYSIS COMPLETE!
echo.
echo Source Dir: %SOURCE_DIR%
echo Results in: %SOURCE_DIR%%PROJ_DIR%
echo Binary:     %SOURCE_DIR%%EXE_NAME%
echo ===================================================
pause