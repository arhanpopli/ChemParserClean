@echo off
REM ============================================================
REM START PUBCHEM SERVER ONLY
REM Node.js server on port 5002 for 3D molecular viewing
REM ============================================================

title PubChem Server (3D Viewer) - Port 5002

cd /d "%~dp0"

echo.
echo ============================================================
echo   PUBCHEM 3D VIEWER SERVER
echo ============================================================
echo   Port: 5002
echo   URL:  http://localhost:5002
echo ============================================================
echo.

REM Check Node.js
node --version >nul 2>&1
if errorlevel 1 (
    echo ERROR: Node.js not found!
    echo Install from: https://nodejs.org/
    pause
    exit /b 1
)

REM Stop existing node processes
taskkill /F /IM node.exe >nul 2>&1
timeout /t 1 >nul

REM Navigate to PubChem directory and start
cd /d "%~dp0MoleculeViewer\pubchem"
echo Starting server from: %CD%
echo.

node .\server.js
