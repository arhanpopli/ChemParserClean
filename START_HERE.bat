@echo off
REM ==============================================================
REM MoleculeViewer - Reliable Startup with Auto-Restart
REM ==============================================================
REM Backend stays running and auto-restarts if it crashes!
REM ==============================================================

setlocal enabledelayedexpansion

color 0A
cls

echo.
echo  ╔════════════════════════════════════════════════════════╗
echo  ║    🚀 MoleculeViewer - Starting (Auto-Restart ON)      ║
echo  ╚════════════════════════════════════════════════════════╝
echo.

REM Change to MoleculeViewer directory
cd /d "%~dp0MoleculeViewer" || (
    echo ERROR: Could not change to MoleculeViewer directory
    pause
    exit /b 1
)

REM ==============================================================
REM STEP 1: Check Python
REM ==============================================================
echo  [1/3] Checking Python...
python --version >nul 2>&1
if errorlevel 1 (
    color 0C
    echo  ❌ ERROR: Python not found!
    echo  Please install Python from https://www.python.org
    pause
    exit /b 1
)

for /f "tokens=*" %%i in ('python --version 2^>^&1') do set PYTHON_VERSION=%%i
echo  ✅ Found: %PYTHON_VERSION%
echo.

REM ==============================================================
REM STEP 2: Check Dependencies
REM ==============================================================
echo  [2/3] Checking dependencies...
python -m pip show flask > nul 2>&1
if errorlevel 1 (
    echo  ⚠️  Installing dependencies...
    python -m pip install -q flask rdkit requests >nul 2>&1
    if errorlevel 1 (
        color 0C
        echo  ❌ ERROR: Could not install dependencies!
        pause
        exit /b 1
    )
)
echo  ✅ All dependencies ready!
echo.

REM ==============================================================
REM STEP 3: Kill any existing Python process on port 5000
REM ==============================================================
echo  [3/3] Cleaning up old processes...
for /f "tokens=5" %%a in ('netstat -ano 2^>nul ^| findstr :5000') do (
    taskkill /PID %%a /F >nul 2>&1
)
timeout /t 1 /nobreak >nul

REM ==============================================================
REM STEP 4: Start Backend with Auto-Restart Loop
REM ==============================================================
set RESTART_COUNT=0

:RESTART_LOOP
set /a RESTART_COUNT+=1

if %RESTART_COUNT% gtr 1 (
    color 0C
    echo.
    echo  ⚠️  Backend crashed! Auto-restarting... (Attempt %RESTART_COUNT%)
    color 0A
    echo.
    timeout /t 2 /nobreak >nul
)

if %RESTART_COUNT% equ 1 (
    echo  Starting backend server...
    echo.
    echo  ✅ Backend RUNNING on http://localhost:5000
    echo.
    echo  📋 Keep this window OPEN while using molecules!
    echo.
    echo  🧪 To test:
    echo     1. Go to chrome://extensions/
    echo     2. Reload "Chemistry Renderer"
    echo     3. Type in ChatGPT: chem:benzene:
    echo.
    echo  ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
    echo.
)

REM Start the backend
python run.py

REM If we get here, the backend crashed - restart it
goto RESTART_LOOP
