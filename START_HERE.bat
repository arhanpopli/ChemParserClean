@echo off
REM ==============================================================
REM MoleculeViewer - Complete Startup (Backend + Instructions)
REM ==============================================================
REM This is the ONLY file you need to run!
REM Starts the Flask backend and shows everything you need
REM ==============================================================

color 0A
cls

echo.
echo  ╔════════════════════════════════════════════════════════╗
echo  ║          🚀 MoleculeViewer Complete Startup            ║
echo  ╚════════════════════════════════════════════════════════╝
echo.

REM Change to MoleculeViewer directory
cd /d "%~dp0MoleculeViewer"

REM ==============================================================
REM STEP 1: Check Python
REM ==============================================================
echo  [1/4] Checking Python installation...
python --version >nul 2>&1
if errorlevel 1 (
    color 0C
    echo.
    echo  ❌ ERROR: Python not found!
    echo.
    echo  Please install Python from https://www.python.org
    echo  ✓ Check "Add Python to PATH" during installation
    echo.
    pause
    exit /b 1
)

for /f "tokens=*" %%i in ('python --version 2^>^&1') do set PYTHON_VERSION=%%i
echo  ✅ Found: %PYTHON_VERSION%
echo.

REM ==============================================================
REM STEP 2: Check Dependencies
REM ==============================================================
echo  [2/4] Checking dependencies (flask, rdkit, requests)...
python -m pip show flask > nul 2>&1
if errorlevel 1 (
    echo  ⚠️  Installing missing dependencies...
    python -m pip install -q flask rdkit requests
    if errorlevel 1 (
        color 0C
        echo.
        echo  ❌ ERROR: Failed to install dependencies!
        echo.
        pause
        exit /b 1
    )
)
echo  ✅ All dependencies ready!
echo.

REM ==============================================================
REM STEP 3: Start Flask Backend
REM ==============================================================
echo  [3/4] Starting Flask backend...
echo.

python run.py >nul 2>&1 &
set BACKEND_PID=%ERRORLEVEL%

timeout /t 2 /nobreak >nul

REM Check if backend started successfully
python -c "import socket; s = socket.socket(); s.connect(('127.0.0.1', 5000)); s.close()" >nul 2>&1
if errorlevel 1 (
    color 0C
    echo.
    echo  ❌ ERROR: Backend failed to start!
    echo  Check that port 5000 is not already in use
    echo.
    pause
    exit /b 1
)

color 0A
echo  ✅ Backend is running on http://localhost:5000
echo.

REM ==============================================================
REM STEP 4: Show Ready Message
REM ==============================================================
echo  [4/4] Setup complete!
echo.
echo  ╔════════════════════════════════════════════════════════╗
echo  ║              🎉 READY TO USE! 🎉                      ║
echo  ╚════════════════════════════════════════════════════════╝
echo.
echo  ✅ Backend:     RUNNING on http://localhost:5000
echo  ✅ Extensions:  Ready to use
echo  ✅ ChatGPT:     Ready for molecules
echo.
echo  ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
echo.
echo  📋 NEXT STEPS:
echo.
echo  1. Reload Chrome Extension:
echo     • Go to: chrome://extensions/
echo     • Find: Chemistry Renderer
echo     • Click: 🔄 Reload button
echo.
echo  2. Open ChatGPT (new tab)
echo.
echo  3. Type a molecule with TRAILING COLON:
echo     • chem:benzene:
echo     • chem:aspirin:
echo     • chem:CCO:
echo.
echo  ⚠️  IMPORTANT: Keep this window open while using molecules!
echo  (The backend server runs in the background)
echo.
echo  📊 ENDPOINTS AVAILABLE:
echo     • /img/smiles?smiles=CCO
echo     • /img/nomenclature?nomenclature=acetone
echo     • /health
echo.
echo  ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
echo.
echo  🛑 To stop: Close this window or press Ctrl+C
echo.

REM Keep the backend running
python run.py
