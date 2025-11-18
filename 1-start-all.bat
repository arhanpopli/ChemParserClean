@echo off
REM ============================================================
REM CHEMPARSER - START ALL SERVERS (MAIN STARTUP)
REM ============================================================
REM This is the PRIMARY startup script for the entire system.
REM It starts all three servers: Mol2ChemFig, MoleculeViewer, PubChem
REM Then opens the unified interface.
REM ============================================================

title ChemParser - Complete System Startup
color 0A

setlocal enabledelayedexpansion

REM Get the root directory
set ROOT_DIR=%~dp0
cd /d "%ROOT_DIR%"

echo.
echo ============================================================
echo   CHEMPARSER - COMPLETE SYSTEM STARTUP
echo ============================================================
echo   Starting 3 servers in separate windows...
echo   1. Mol2ChemFig Server     (Flask, Port 5001)
echo   2. MoleculeViewer Backend (Node.js, Port 5000)
echo   3. PubChem Server        (Node.js, Port 5002)
echo   Then: Unified Interface in your browser
echo ============================================================
echo.

REM Kill existing processes to avoid port conflicts
echo [PRE] Cleaning up existing processes...
taskkill /F /IM node.exe >nul 2>&1
taskkill /F /IM python.exe >nul 2>&1
timeout /t 2 >nul

echo.
echo [1/4] Starting Mol2ChemFig Server (Port 5001)...
start "Mol2ChemFig (5001)" cmd /k python mol2chemfig_server.py
timeout /t 3 >nul

echo [2/4] Starting MoleculeViewer Backend (Port 5000)...
cd MoleculeViewer
start "MoleculeViewer (5000)" cmd /k node server.js
cd ..
timeout /t 3 >nul

echo [3/4] Starting PubChem Server (Port 5002)...
cd MoleculeViewer\pubchem
start "PubChem (5002)" cmd /k node server.js
cd ..\..
timeout /t 3 >nul

echo [4/4] Opening Unified Interface...
timeout /t 2 >nul
start http://localhost:5000/unified-interface.html

echo.
echo ============================================================
echo   ✓ ALL SERVERS STARTED SUCCESSFULLY
echo ============================================================
echo.
echo   Access the system at:
echo   → http://localhost:5000/unified-interface.html
echo.
echo   Individual servers:
echo   → Mol2ChemFig:      http://localhost:5001
echo   → MoleculeViewer:   http://localhost:5000
echo   → PubChem:          http://localhost:5002
echo.
echo   To stop servers: run util-stop-all.bat
echo   To check status: run util-status.bat
echo ============================================================
echo.
pause
