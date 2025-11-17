@echo off
REM ============================================================
REM ChemParser - Check Server Status
REM This script checks the status of all ChemParser servers
REM ============================================================

title ChemParser - Server Status

setlocal enabledelayedexpansion

echo.
echo ============================================================
echo   CHEMPARSER - SERVER STATUS CHECK
echo ============================================================
echo.

REM Get script directory
set SCRIPT_DIR=%~dp0
cd /d "%SCRIPT_DIR%"

echo Checking server availability...
echo.

REM Check MoleculeViewer (Port 5000)
echo [1/4] MoleculeViewer Server (Port 5000)
curl -s http://localhost:5000/health >nul 2>&1
if errorlevel 1 (
    echo   Status: [✗] OFFLINE
    echo   URL: http://localhost:5000
) else (
    echo   Status: [✓] ONLINE
    echo   URL: http://localhost:5000
)
echo.

REM Check Docker Backend (Port 8000)
echo [2/4] Mol2ChemFig Docker Backend (Port 8000)
curl -s http://localhost:8000 >nul 2>&1
if errorlevel 1 (
    echo   Status: [✗] OFFLINE
    echo   URL: http://localhost:8000
    echo   Note: Make sure Docker Desktop is running and containers are started
) else (
    echo   Status: [✓] ONLINE
    echo   URL: http://localhost:8000
)
echo.

REM Check Mol2ChemFig Server (Port 5001)
echo [3/4] Mol2ChemFig Server (Port 5001)
curl -s http://localhost:5001/health >nul 2>&1
if errorlevel 1 (
    echo   Status: [✗] OFFLINE
    echo   URL: http://localhost:5001
) else (
    echo   Status: [✓] ONLINE
    echo   URL: http://localhost:5001
)
echo.

REM Check PubChem Server (Port 5002)
echo [4/4] PubChem Server (Port 5002)
curl -s http://localhost:5002/health >nul 2>&1
if errorlevel 1 (
    echo   Status: [✗] OFFLINE
    echo   URL: http://localhost:5002
) else (
    echo   Status: [✓] ONLINE
    echo   URL: http://localhost:5002
)
echo.

echo ============================================================
echo   Docker Containers Status
echo ============================================================
echo.

docker ps --filter "name=m2cf" --format "table {{.Names}}\t{{.Status}}\t{{.Ports}}" 2>nul
if errorlevel 1 (
    echo   Docker not available or no containers running
)

echo.
echo ============================================================
echo   Process Information
echo ============================================================
echo.

echo Node.js processes:
tasklist /FI "IMAGENAME eq node.exe" /FI "WINDOWTITLE eq *ChemParser*" 2>nul | find "node.exe"
if errorlevel 1 (
    echo   No Node.js processes found
)

echo.
echo Python processes:
tasklist /FI "IMAGENAME eq python.exe" /FI "WINDOWTITLE eq *ChemParser*" 2>nul | find "python.exe"
if errorlevel 1 (
    echo   No Python processes found
)

echo.
echo ============================================================
echo   Quick Actions
echo ============================================================
echo.
echo   Start all servers:  start_all.bat
echo   Stop all servers:   stop_all.bat
echo   Open launcher:      launcher.html
echo ============================================================
echo.

pause
