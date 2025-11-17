@echo off
REM ============================================================
REM START LAUNCHER CONTROL SERVER
REM This starts the Node.js backend that controls all servers
REM ============================================================

title ChemParser Launcher Control - Port 3000

cd /d "%~dp0"

echo.
echo ============================================================
echo   CHEMPARSER LAUNCHER CONTROL SERVER
echo ============================================================
echo   Starting on port 3000...
echo.

REM Check Node.js
node --version >nul 2>&1
if errorlevel 1 (
    echo ERROR: Node.js not found!
    echo Install from: https://nodejs.org/
    pause
    exit /b 1
)

REM Install dependencies if needed
if not exist "node_modules\express" (
    echo Installing dependencies...
    npm install express cors
    echo.
)

echo Starting control server...
echo.

REM Start launcher in new window with correct working directory
start "ChemParser Launcher Control" cmd /k "cd /d "%~dp0" && node launcher-server.js"

echo Waiting for control server to start...
timeout /t 2 /nobreak >nul

REM Check if server is responding (retry up to 10 times)
set /a attempts=0
:check_server
set /a attempts+=1
curl -s http://localhost:3000/health >nul 2>&1
if %ERRORLEVEL% EQU 0 (
    echo [SUCCESS] Launcher control server is running!
    goto server_ready
)
if %attempts% GEQ 10 (
    echo.
    echo [ERROR] Launcher control server failed to start!
    echo.
    echo Troubleshooting:
    echo   1. Check if port 3000 is already in use
    echo   2. Look at the "ChemParser Launcher Control" window for errors
    echo   3. Make sure Node.js is installed: node --version
    echo.
    pause
    exit /b 1
)
timeout /t 1 /nobreak >nul
goto check_server

:server_ready
echo.
echo Opening control panel in browser...
start http://localhost:3000/launcher.html

echo.
echo ============================================================
echo   SUCCESS! Control Panel is running
echo ============================================================
echo.
echo   Control Panel: http://localhost:3000/launcher.html
echo   Launcher Window: "ChemParser Launcher Control" (keep it open!)
echo.
echo   You can now click buttons in the browser to start/stop servers.
echo   DO NOT close the "ChemParser Launcher Control" window!
echo.
echo ============================================================
pause
