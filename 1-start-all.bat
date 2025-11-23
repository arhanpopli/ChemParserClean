@echo off
REM ============================================================
REM CHEMPARSER - START ALL SERVERS (MAIN STARTUP)
REM ============================================================
REM This is the PRIMARY startup script for the entire system.
REM It starts Docker, all servers, and opens the unified interface.
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
echo   Starting all services:
echo   0. Docker Backend         (Docker, Port 8000) - for mol2chemfig
echo   1. Mol2ChemFig Server     (Flask, Port 5001)
echo   2. MoleculeViewer Backend (Node.js, Port 5000)
echo   3. PubChem Server         (Node.js, Port 5002)
echo   Then: Unified Interface in your browser
echo ============================================================
echo.

REM Kill existing processes to avoid port conflicts
echo [PRE] Cleaning up existing processes...
taskkill /F /IM node.exe >nul 2>&1
taskkill /F /IM python.exe >nul 2>&1
timeout /t 2 >nul

echo.
echo [0/5] Starting Docker Backend (Port 8000)...
echo       NOTE: Docker is OPTIONAL - local mol2chemfig processing is available!
docker-compose up -d 2>nul
if errorlevel 1 (
    echo       [INFO] Docker not running - using LOCAL mol2chemfig processing
    echo       Local processing uses mol2chemfigPy3 + MiKTeX LaTeX
) else (
    echo       Docker backend starting (will be used as fallback)...
    timeout /t 5 >nul
)

echo.
echo [1/5] Starting Mol2ChemFig Server (Port 5001)...
start "Mol2ChemFig (5001)" /D "%ROOT_DIR%" cmd /k python mol2chemfig_server.py
timeout /t 3 >nul

echo [2/5] Starting MoleculeViewer Backend (Port 5000)...
start "MoleculeViewer (5000)" /D "%ROOT_DIR%MoleculeViewer" cmd /k node server.js
timeout /t 3 >nul

echo [3/5] Starting PubChem Server (Port 5002)...
start "PubChem (5002)" /D "%ROOT_DIR%MoleculeViewer\pubchem" cmd /k node server.js
timeout /t 3 >nul

echo [4/5] Opening Unified Interface...
timeout /t 2 >nul
start http://localhost:5000/unified-interface.html

echo.
echo ============================================================
echo   ALL SERVERS STARTED
echo ============================================================
echo.
echo   Access the system at:
echo   http://localhost:5000/unified-interface.html
echo.
echo   Individual servers:
echo   Docker Backend:    http://localhost:8000 (mol2chemfig core)
echo   Mol2ChemFig:       http://localhost:5001 (Flask wrapper)
echo   MoleculeViewer:    http://localhost:5000
echo   PubChem:           http://localhost:5002
echo.
echo   To stop servers: run util-stop-all.bat
echo   To check status: run util-status.bat
echo ============================================================
echo.
pause
