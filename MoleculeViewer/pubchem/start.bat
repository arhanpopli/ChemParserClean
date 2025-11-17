@echo off
echo ============================================
echo Starting PubChem Integration Server
echo ============================================
echo.
echo This server provides:
echo   - 2D molecular structure images from PubChem
echo   - 3D interactive molecular viewers
echo   - Direct PubChem API integration
echo.
echo Server will start on: http://localhost:5002
echo.
echo ============================================
echo.

cd /d "%~dp0"

REM Check if node_modules exists
if not exist "node_modules\" (
    echo Installing dependencies...
    call npm install
    echo.
)

echo Starting server...
node server.js

pause
