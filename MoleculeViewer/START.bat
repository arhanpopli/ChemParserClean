@echo off
REM ============================================================
REM MoleculeViewer - Master Startup Script
REM Starts the Flask backend server on http://localhost:5000
REM ============================================================

title MoleculeViewer Server

echo.
echo ========================================================
echo          MOLECULEVIEWER - Starting Server
echo ========================================================
echo.

cd /d "%~dp0"

echo [1/3] Checking Python installation...
python --version >nul 2>&1
if errorlevel 1 (
    echo ERROR: Python is not installed or not in PATH
    echo Please install Python 3.8+ from python.org
    pause
    exit /b 1
)
echo      Python: OK

echo.
echo [2/3] Starting Flask server...
echo      URL: http://localhost:5000
echo      Press Ctrl+C to stop the server
echo.
echo ========================================================
echo.

REM Start the Flask server using the exact command that worked
python -m flask --app app.api run --host=0.0.0.0 --port=5000

REM If we get here, the server stopped
echo.
echo ========================================================
echo Server stopped.
echo ========================================================
pause
