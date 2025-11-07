@echo off
REM ============================================================
REM MoleculeViewer - Complete Master Startup
REM Starts the Flask backend server and opens the web interface
REM ============================================================

title MoleculeViewer - Master Startup

echo.
echo ========================================================
echo        MOLECULEVIEWER - COMPLETE SYSTEM START
echo ========================================================
echo.

cd /d "%~dp0"

echo [1/4] Checking Python installation...
python --version >nul 2>&1
if errorlevel 1 (
    echo ERROR: Python is not installed or not in PATH
    echo Please install Python 3.8+ from python.org
    pause
    exit /b 1
)
echo      Python: OK

echo.
echo [2/4] Stopping any existing Flask servers...
taskkill /F /IM python.exe >nul 2>&1
timeout /t 1 /nobreak >nul

echo.
echo [3/4] Starting Flask backend server...
echo      URL: http://localhost:5000
echo.

REM Start Python in background and capture the PID
start "MoleculeViewer Backend" python -m flask --app app.api run --host=0.0.0.0 --port=5000

REM Wait for Flask to start
timeout /t 3 /nobreak >nul

echo.
echo [4/4] Opening MoleculeViewer in your browser...
timeout /t 1 /nobreak >nul
start http://localhost:5000

echo.
echo ========================================================
echo                    READY TO USE!
echo ========================================================
echo.
echo USAGE:
echo   Type in your browser:  chem:benzene:
echo   Or use SMILES:         chem:CCO:
echo   Or complex names:      chem:1-chloro-benzene:
echo.
echo STOP:
echo   Close the black "MoleculeViewer Backend" window
echo.
echo ========================================================
echo.

pause
