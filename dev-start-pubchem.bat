@echo off
REM ============================================================
REM CHEMPARSER - DEV: START PUBCHEM SERVER ONLY
REM ============================================================
REM Development use: Start only PubChem 3D viewer (Port 5002)
REM Useful for testing 3D molecule visualization without other servers
REM ============================================================

title ChemParser - PubChem 3D Server (Dev)

set ROOT_DIR=%~dp0
cd /d "%ROOT_DIR%"

echo.
echo ============================================================
echo   PubChem 3D Viewer Server (Development)
echo ============================================================
echo   Port: 5002
echo   Type: Node.js Server
echo.
echo   Access at: http://localhost:5002
echo ============================================================
echo.

cd MoleculeViewer\pubchem
node server.js
