@echo off
REM MoleculeViewer Flask Server Launcher
REM This batch file starts the MoleculeViewer Flask server

echo.
echo ====================================
echo MoleculeViewer Server Launcher
echo ====================================
echo.

REM Change to the MoleculeViewer directory
cd /d "%~dp0"

REM Check if Python is installed
python --version >nul 2>&1
if errorlevel 1 (
    echo ERROR: Python is not installed or not in PATH
    echo Please install Python 3.9 or higher
    pause
    exit /b 1
)

echo [*] Starting MoleculeViewer Flask Server...
echo [*] Server will run on: http://localhost:5000
echo.

REM Start Flask server
python start.py

REM If the server exits, show a message
echo.
echo Server stopped.
pause
