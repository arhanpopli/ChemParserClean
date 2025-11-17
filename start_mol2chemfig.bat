@echo off
REM ============================================================
REM Start Mol2ChemFig Flask Server (Port 5001)
REM ============================================================

title ChemParser - Mol2ChemFig (Port 5001)

set SCRIPT_DIR=%~dp0
cd /d "%SCRIPT_DIR%"

echo.
echo ============================================================
echo   Starting Mol2ChemFig Server
echo ============================================================
echo   Port: 5001
echo   Type: Flask (Python)
echo ============================================================
echo.

python mol2chemfig_server.py

pause
