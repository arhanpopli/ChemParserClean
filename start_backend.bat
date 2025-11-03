@echo off
REM ==============================================================
REM MoleculeViewer - Start Backend Server
REM ==============================================================
REM This script starts the Flask backend for MoleculeViewer
REM Location: http://localhost:5000
REM ==============================================================

echo.
echo  ==========================================
echo  MoleculeViewer Backend - Starting
echo  ==========================================
echo.

REM Change to MoleculeViewer directory
cd /d "%~dp0MoleculeViewer"

REM Check if Python is available
python --version >nul 2>&1
if errorlevel 1 (
    echo  ERROR: Python not found!
    echo.
    echo  Please install Python from https://www.python.org
    echo  Make sure to check "Add Python to PATH" during installation
    echo.
    pause
    exit /b 1
)

echo  Python found: 
python --version
echo.

REM Check if required packages are installed
echo  Checking dependencies...
python -m pip show flask > nul 2>&1
if errorlevel 1 (
    echo.
    echo  WARNING: Missing dependencies. Installing...
    echo.
    python -m pip install flask rdkit requests
    if errorlevel 1 (
        echo.
        echo  ERROR: Failed to install dependencies!
        pause
        exit /b 1
    )
)

echo  All dependencies OK!
echo.

REM Start the backend using Flask directly
echo  ==========================================
echo  Starting Flask backend on http://localhost:5000
echo  ==========================================
echo.
echo  Server is starting. Keep this window open!
echo  Press Ctrl+C to stop the server.
echo.

python run.py

if errorlevel 1 (
    echo.
    echo  ==========================================
    echo  ERROR: Backend failed to start!
    echo  ==========================================
    echo.
    pause
    exit /b 1
)

pause
