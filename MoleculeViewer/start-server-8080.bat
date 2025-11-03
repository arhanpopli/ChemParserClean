@echo off
REM MoleculeViewer Server Launcher - Port 8080
REM Starts the Flask server on port 8080

echo.
echo ====================================
echo MoleculeViewer Server - Port 8080
echo ====================================
echo.

REM Change to the MoleculeViewer directory
cd /d "%~dp0"

REM Check if Python is installed
python --version >nul 2>&1
if errorlevel 1 (
    echo ERROR: Python is not installed or not in PATH
    pause
    exit /b 1
)

echo [*] Starting MoleculeViewer Flask Server on port 8080...
echo [*] Available at: http://localhost:8080
echo [*] Press CTRL+C to stop
echo.

REM Start the server
python server_8080.py

echo.
echo Server stopped.
pause
