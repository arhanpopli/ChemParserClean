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

python -c "
import os
import sys
from app.api import app

try:
    print(' ✅ Backend module loaded successfully')
    print()
    print(' 🚀 Starting Flask server...')
    print(' 📍 Location: http://localhost:5000')
    print()
    print(' Available endpoints:')
    print('   - /img/smiles?smiles=CCO')
    print('   - /img/nomenclature?nomenclature=acetone')
    print('   - /health')
    print()
    print(' NOTE: Keep this window open while using the extension!')
    print()
    print(' ==========================================')
    print()
    app.run(host='0.0.0.0', port=5000, debug=False, use_reloader=False)
except Exception as e:
    print(' ❌ ERROR: Failed to start backend!')
    print(f' Details: {str(e)}')
    print()
    sys.exit(1)
"

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
