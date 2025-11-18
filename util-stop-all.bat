@echo off
REM ============================================================
REM CHEMPARSER - UTILITY: STOP ALL SERVERS
REM ============================================================
REM Kills all running ChemParser servers (Node.js and Python)
REM Safe to run anytime - won't crash if nothing is running
REM ============================================================

title ChemParser - Stop All Servers

echo.
echo ============================================================
echo   ChemParser - Stopping All Servers
echo ============================================================
echo.

echo [1/2] Stopping Node.js servers...
taskkill /F /IM node.exe >nul 2>&1
if errorlevel 1 (
    echo   (No Node.js processes running)
) else (
    echo   ✓ Node.js servers stopped
)

echo.
echo [2/2] Stopping Python servers...
taskkill /F /IM python.exe >nul 2>&1
if errorlevel 1 (
    echo   (No Python processes running)
) else (
    echo   ✓ Python servers stopped
)

echo.
echo ============================================================
echo   ✓ All servers stopped
echo ============================================================
echo.
pause
