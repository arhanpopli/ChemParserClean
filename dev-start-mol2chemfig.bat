@echo off
REM ============================================================
REM CHEMPARSER - DEV: START MOL2CHEMFIG SERVER ONLY
REM ============================================================
REM Development use: Start only the Mol2ChemFig server (Port 5001)
REM Useful for testing/debugging without other servers
REM ============================================================

title ChemParser - Mol2ChemFig Server (Dev)

set ROOT_DIR=%~dp0
cd /d "%ROOT_DIR%"

echo.
echo ============================================================
echo   Mol2ChemFig Server (Development)
echo ============================================================
echo   Port: 5001
echo   Type: Flask Server (Python)
echo.
echo   Access at: http://localhost:5001
echo ============================================================
echo.

python mol2chemfig_server.py
