@echo off
REM ==============================================================
REM MoleculeViewer - Start Backend Server
REM ==============================================================
REM This script starts the Flask backend for MoleculeViewer
REM Location: http://localhost:5000
REM
REM Requirements:
REM   - Python 3.7+
REM   - pip packages: flask, rdkit, requests
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
python -c "import flask, rdkit, requests" >nul 2>&1
if errorlevel 1 (
    echo.
    echo  WARNING: Missing dependencies. Installing...
    echo.
    pip install flask rdkit requests
    if errorlevel 1 (
        echo.
        echo  ERROR: Failed to install dependencies!
        pause
        exit /b 1
    )
)

echo  All dependencies OK!
echo.

REM Start the backend
echo  ==========================================
echo  Starting Flask backend on http://localhost:5000
echo  ==========================================
echo.
echo  The server is starting. You should see output below:
echo.

REM Create a temporary Python script to run the server
setlocal enabledelayedexpansion
set "TEMP_SCRIPT=%TEMP%\mol_start_server.py"

(
echo import os
echo import sys
echo from app.api import app
echo.
echo try:
echo     print(' ✅ Backend module loaded successfully'^)
echo     print('^)
echo     print(' 🚀 Starting Flask server...'^)
echo     print(' 📍 Location: http://localhost:5000'^)
echo     print('^)
echo     print(' Available endpoints:'^)
echo     print('   - /img/smiles?smiles=CCO'^)
echo     print('   - /img/nomenclature?nomenclature=acetone'^)
echo     print('   - /health'^)
echo     print('^)
echo     print(' NOTE: Keep this window open while using the extension!'^)
echo     print('^)
echo     print(' =========================================='^)
echo     print('^)
echo     app.run(host='0.0.0.0', port=5000, debug=False, use_reloader=False^)
echo except Exception as e:
echo     print(' ❌ ERROR: Failed to start backend!'^)
echo     print(f' Details: {str(e)}'^)
echo     print('^)
echo     sys.exit(1^)
) > "%TEMP_SCRIPT%"

python "%TEMP_SCRIPT%"

if errorlevel 1 (
    echo.
    echo  ==========================================
    echo  ERROR: Backend failed to start!
    echo  ==========================================
    echo.
    del "%TEMP_SCRIPT%" >nul 2>&1
    pause
    exit /b 1
)

del "%TEMP_SCRIPT%" >nul 2>&1
pause
