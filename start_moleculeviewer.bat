@echo off
REM ============================================================
REM Start MoleculeViewer Server (Node.js - Port 5000)
REM ============================================================

title ChemParser - MoleculeViewer (Port 5000)

set SCRIPT_DIR=%~dp0
cd /d "%SCRIPT_DIR%MoleculeViewer"

echo.
echo ============================================================
echo   Starting MoleculeViewer Server
echo ============================================================
echo   Port: 5000
echo   Type: Node.js
echo ============================================================
echo.

node server.js

pause
