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
echo   0. Docker Backend         (Docker, Port 8000) - DEPRECATED
echo   1. Mol2ChemFig Server     (Flask, Port 1000)
echo   2. MoleculeViewer Backend (Node.js, Port 5000)
echo   3. PubChem Server         (Node.js, Port 5002)
echo   4. MolView PHP Server     (PHP, Port 8000) - Main viewer
echo   5. MolView Search API     (Node.js, Port 8001) - Unified search
echo   Then: Unified Interface in your browser
echo ============================================================
echo.

REM Kill existing processes to avoid port conflicts
echo [PRE] Cleaning up existing processes...
taskkill /F /IM node.exe >nul 2>&1
taskkill /F /IM python.exe >nul 2>&1
taskkill /F /IM php.exe >nul 2>&1
timeout /t 2 >nul

echo.
echo [SKIPPED] Docker Backend (Port 8000) - Using MolView PHP server instead
echo.

echo [1/5] Starting Mol2ChemFig Server (Port 1000)...
start "Mol2ChemFig (1000)" /D "%ROOT_DIR%" cmd /k python mol2chemfig_server.py
timeout /t 3 >nul

echo [2/5] Starting MoleculeViewer Backend (Port 5000)...
start "MoleculeViewer (5000)" /D "%ROOT_DIR%MoleculeViewer" cmd /k node server.js
timeout /t 3 >nul

echo [3/5] Starting PubChem Server (Port 5002)...
start "PubChem (5002)" /D "%ROOT_DIR%MoleculeViewer\pubchem" cmd /k node server.js
timeout /t 3 >nul

echo [4/5] Starting MolView PHP Server (Port 8000)...
start "MolView PHP (8000)" /D "%ROOT_DIR%Molview\molview" cmd /k php -S localhost:8000
timeout /t 3 >nul

echo [5/5] Starting MolView Search API (Port 8001)...
start "MolView Search (8001)" /D "%ROOT_DIR%Molview\molview" cmd /k node search-server.js
timeout /t 3 >nul

echo [6/6] Opening Unified Interface...
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
echo   MolView PHP:       http://localhost:8000 (main viewer)
echo   MolView Search:    http://localhost:8001 (unified search API)
echo   Mol2ChemFig:       http://localhost:1000 (Flask wrapper)
echo   MoleculeViewer:    http://localhost:5000
echo   PubChem:           http://localhost:5002
echo.
echo   Extension uses:
echo   - Port 8001 for ALL queries (with autocorrect)
echo   - Port 8000 for MolView embeds
echo.
echo   To stop servers: run util-stop-all.bat
echo   To check status: run util-status.bat
echo ============================================================
echo.
pause
