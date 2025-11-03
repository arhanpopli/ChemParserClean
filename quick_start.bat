@echo off
REM ==============================================================
REM Quick Start - MoleculeViewer Full Setup
REM ==============================================================
REM This script opens both:
REM   1. Backend server (this window)
REM   2. Instructions in default browser
REM ==============================================================

echo.
echo  ==========================================
echo  MoleculeViewer - Quick Start
echo  ==========================================
echo.
echo  Starting backend server...
echo.

REM Get the directory where this script is located
set SCRIPT_DIR=%~dp0

REM Start backend in a new window and keep it running
start "MoleculeViewer Backend" cmd /k "%SCRIPT_DIR%start_backend.bat"

REM Wait a moment for server to start
timeout /t 3 /nobreak

echo.
echo  Backend started! Now reload your Chrome extension:
echo.
echo  1. Go to: chrome://extensions/
echo  2. Find: Chemistry Renderer
echo  3. Click: Reload button
echo.
echo  4. Go to ChatGPT
echo  5. Type: chem:benzene:
echo.
echo  (Note the trailing colon!)
echo.
pause
