@echo off
REM ============================================================
REM CHEMPARSER - UTILITY: CHECK SERVER STATUS
REM ============================================================
REM Quick health check of all servers
REM ============================================================

title ChemParser - Server Status Check

echo.
echo ============================================================
echo   ChemParser - Server Status Check
echo ============================================================
echo.

echo Checking servers...
echo.

REM Check Mol2ChemFig (Port 5001)
echo [1/3] Mol2ChemFig (Port 5001)...
timeout /t 1 >nul
netstat -ano | find ":5001" >nul 2>&1
if errorlevel 1 (
    echo   Status: ✗ OFFLINE
) else (
    echo   Status: ✓ ONLINE
)

REM Check MoleculeViewer (Port 5000)
echo.
echo [2/3] MoleculeViewer (Port 5000)...
timeout /t 1 >nul
netstat -ano | find ":5000" >nul 2>&1
if errorlevel 1 (
    echo   Status: ✗ OFFLINE
) else (
    echo   Status: ✓ ONLINE
)

REM Check PubChem (Port 5002)
echo.
echo [3/3] PubChem (Port 5002)...
timeout /t 1 >nul
netstat -ano | find ":5002" >nul 2>&1
if errorlevel 1 (
    echo   Status: ✗ OFFLINE
) else (
    echo   Status: ✓ ONLINE
)

echo.
echo ============================================================
echo   To start all servers, run: 1-start-all.bat
echo   To stop all servers, run: util-stop-all.bat
echo ============================================================
echo.
pause
