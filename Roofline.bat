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
echo   Intel oneAPI Roofline Automator
echo ===================================================
set /p FULL_CPP_PATH="Enter the FULL path to your .cpp file: "

for %%i in ("%FULL_CPP_PATH%") do (
    set "SOURCE_DIR=%%~dpi"
    set "FILE_NAME=%%~ni"
)

cd /d "%SOURCE_DIR%"
set EXE_NAME=%FILE_NAME%.exe
set PROJ_DIR=%FILE_NAME%_Advisor

:: --- STEP 3: INITIALIZE ENVIRONMENTS ---
echo.
echo [1/3] Initializing Environments...

:: Set defaults
set "DEFAULT_VS=C:\Program Files (x86)\Microsoft Visual Studio\2022\BuildTools"
set "DEFAULT_ONEAPI=C:\Program Files (x86)\Intel\oneAPI\setvars.bat"

echo Default VS Path: %DEFAULT_VS%
set /p USER_VS="Press Enter to use default or specify new path: "
if "%USER_VS%"=="" (set "VS2022INSTALLDIR=%DEFAULT_VS%") else (set "VS2022INSTALLDIR=%USER_VS%")

echo.
echo Default oneAPI Path: %DEFAULT_ONEAPI%
set /p USER_ONEAPI="Press Enter to use default or specify new path: "
if "%USER_ONEAPI%"=="" (set "ONEAPI_PATH=%DEFAULT_ONEAPI%") else (set "ONEAPI_PATH=%USER_ONEAPI%")

call "%ONEAPI_PATH%" intel64 vs2022

:: --- STEP 4: COMPILE ---
echo.
echo [2/3] Compiling %FILE_NAME%.cpp...
icpx -O3 -xHost -g "%FULL_CPP_PATH%" -o "%EXE_NAME%"

if %errorLevel% neq 0 (
    echo [ERROR] Compilation failed.
    pause
    exit /b
)

:: --- STEP 5: RUN ADVISOR ---
echo.
echo [3/3] Running Intel Advisor Roofline Analysis...
advisor --collect=roofline --project-dir=.\%PROJ_DIR% -- .\%EXE_NAME%

echo.
echo ===================================================
echo ANALYSIS COMPLETE!
echo ===================================================

:: --- STEP 6: PLOTTING ---
set /p PLOT="Generate roofline plot? (y/n): "
if /i "%PLOT%"=="y" (
    echo Searching for latest result folder...
    
    REM Use REM instead of :: inside IF blocks to avoid syntax crashes
    set "E000_DIR=%SOURCE_DIR%%PROJ_DIR%\e000"
    set "LATEST_HS="

    REM Loop through all directories starting with 'hs'
    for /d %%D in ("!E000_DIR!\hs*") do (
        set "LATEST_HS=%%D"
    )

    if "!LATEST_HS!"=="" (
        echo [ERROR] No result folder (hsXXX) found in !E000_DIR!
        pause
    ) else (
        echo Found latest result: !LATEST_HS!
        
        REM Check for Python script using the directory of this batch file
        if exist "%~dp0Roofline.py" (
            set "PYTHON_SCRIPT=%~dp0Roofline.py"
        ) else (
            echo [WARNING] Roofline.py not found in %~dp0
            set /p PYTHON_SCRIPT="Please enter the FULL path to Roofline.py: "
        )

        echo Launching Python...
        REM Using "python" and quoting both the script and the data path
        python "!PYTHON_SCRIPT!" "!LATEST_HS!"
        
        if %errorlevel% neq 0 (
            echo [ERROR] Python script failed to execute.
            pause
        )
    )
)

pause