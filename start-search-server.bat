@echo off
REM Chemistry Renderer Search Server Starter
REM Starts the search server on port 8001

echo.
echo ===============================================
echo   Chemistry Renderer - Search Server Starter
echo ===============================================
echo.
echo Starting search server on port 8001...
echo.

cd /d "%~dp0"

REM Kill any existing node processes on port 8001
taskkill /F /IM node.exe >nul 2>&1

REM Wait a moment for the process to fully terminate
timeout /t 2 /nobreak >nul

REM Start the search server
echo.
echo ✓ Starting server...
echo.
node Molview\molview\search-server.js

REM If the script exits, display a message
echo.
echo Server stopped. Press any key to close this window.
pause >nul
