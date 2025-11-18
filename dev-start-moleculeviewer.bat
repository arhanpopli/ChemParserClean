@echo off
REM ============================================================
REM CHEMPARSER - DEV: START MOLECULEVIEWER SERVER ONLY
REM ============================================================
REM Development use: Start only MoleculeViewer backend (Port 5000)
REM Useful for testing SVG generation without other servers
REM ============================================================

title ChemParser - MoleculeViewer Backend (Dev)

set ROOT_DIR=%~dp0
cd /d "%ROOT_DIR%"

echo.
echo ============================================================
echo   MoleculeViewer Backend Server (Development)
echo ============================================================
echo   Port: 5000
echo   Type: Node.js Server
echo.
echo   Access at: http://localhost:5000
echo ============================================================
echo.

cd MoleculeViewer
node server.js
