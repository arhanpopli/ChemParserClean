@echo off
REM ============================================================
REM CHEMPARSER - START ALL SERVERS
REM Consolidated startup script for all required services
REM ============================================================

title CHEMPARSER - All Servers

setlocal enabledelayedexpansion

REM Get the root directory of this script
set ROOT_DIR=%~dp0
cd /d "%ROOT_DIR%"

echo.
echo ============================================================
echo   CHEMPARSER - COMPLETE SYSTEM STARTUP
echo ============================================================
echo   This script starts:
echo   1. MoleculeViewer Backend (Flask on port 5000)
echo   2. PubChem Server (Node.js on port 5002)
echo ============================================================
echo.

REM Check Python
echo [1/2] Checking Python installation...
python --version >nul 2>&1
if errorlevel 1 (
    echo ERROR: Python not found!
    echo Please install Python 3.8+ from python.org
    pause
    exit /b 1
)
for /f "tokens=*" %%i in ('python --version 2^>^&1') do set PYTHON_VERSION=%%i
echo      %PYTHON_VERSION% - OK

REM Check Node.js
echo [2/2] Checking Node.js installation...
node --version >nul 2>&1
if errorlevel 1 (
    echo ERROR: Node.js not found!
    echo Please install Node.js from nodejs.org
    pause
    exit /b 1
)
for /f "tokens=*" %%i in ('node --version 2^>^&1') do set NODE_VERSION=%%i
echo      Node.js !NODE_VERSION! - OK

echo.
echo Stopping any existing servers...
taskkill /F /IM python.exe >nul 2>&1
taskkill /F /IM node.exe >nul 2>&1
timeout /t 2 >nul

echo.
echo ============================================================
echo   STARTING BACKEND SERVERS
echo ============================================================
echo.

REM Terminal 1: PubChem Server (Node.js on port 5002)
echo Starting PubChem Server on port 5002...
cd /d "%ROOT_DIR%MoleculeViewer\pubchem"
start "PubChem Server (Node.js:5002)" /d "%ROOT_DIR%MoleculeViewer\pubchem" cmd /c "node server.js"
timeout /t 3 >nul

REM Terminal 2: MoleculeViewer Backend (Flask on port 5000)
echo Starting MoleculeViewer Backend on port 5000...
cd /d "%ROOT_DIR%"
start "MoleculeViewer Backend (Flask:5000)" /d "%ROOT_DIR%" cmd /c "python backend_server.py"
timeout /t 2 >nul

echo.
echo ============================================================
echo   SYSTEM STATUS
echo ============================================================
echo.
echo [✓] PubChem Server:      http://localhost:5002
echo [✓] MoleculeViewer:      http://localhost:5000
echo [✓] Extension:           Installed in Chrome
echo.
echo SERVICES RUNNING:
echo   - PubChem 3D Viewer (Node.js)
echo   - Chemical Structure Rendering (Flask)
echo.
echo NEXT STEPS:
echo   1. Reload the Chemistry Extension in Chrome
echo      Go to chrome://extensions/ and click reload
echo.
echo   2. Enable 3D Viewer in extension
echo      Click extension icon > Developer Options > Enable 3D Viewer
echo.
echo   3. Test on any webpage
echo      Type: chem:histamine:
echo      Should show 3D molecular viewer inline
echo.
echo DO NOT CLOSE THIS WINDOW - SERVERS WILL STOP
echo ============================================================
echo.

REM Keep this window open
pause
