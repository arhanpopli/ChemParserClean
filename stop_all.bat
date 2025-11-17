@echo off
REM ============================================================
REM ChemParser - Stop All Servers
REM This script stops all running ChemParser servers
REM ============================================================

title ChemParser - Stopping All Servers

setlocal enabledelayedexpansion

echo.
echo ============================================================
echo   CHEMPARSER - STOPPING ALL SERVERS
echo ============================================================
echo.

REM Get script directory
set SCRIPT_DIR=%~dp0
cd /d "%SCRIPT_DIR%"

echo [1/3] Stopping Docker containers...
echo.

REM Stop Docker containers
docker-compose down >nul 2>&1
if errorlevel 1 (
    echo   ⚠ No Docker containers found or Docker not running
) else (
    echo   ✓ Docker containers stopped
)

echo.
echo [2/3] Stopping Node.js servers...
echo.

REM Kill Node.js processes
taskkill /F /IM node.exe /FI "WINDOWTITLE eq *ChemParser*" >nul 2>&1
if errorlevel 1 (
    echo   ⚠ No Node.js servers found
) else (
    echo   ✓ Node.js servers stopped
)

echo.
echo [3/3] Stopping Python servers...
echo.

REM Kill Python processes
taskkill /F /IM python.exe /FI "WINDOWTITLE eq *ChemParser*" >nul 2>&1
if errorlevel 1 (
    echo   ⚠ No Python servers found
) else (
    echo   ✓ Python servers stopped
)

REM Wait for processes to fully terminate
timeout /t 2 >nul

echo.
echo ============================================================
echo   ALL SERVERS STOPPED
echo ============================================================
echo.
echo   All ChemParser services have been terminated.
echo   To restart, run: start_all.bat
echo ============================================================
echo.

pause
