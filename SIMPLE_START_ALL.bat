@echo off
title ChemParser - Simple Starter
color 0A

echo.
echo ============================================================
echo   ChemParser - Simple Start All Servers
echo ============================================================
echo.
echo   This script starts all servers in separate windows.
echo   NO control panel needed - just direct server startup!
echo.
echo ============================================================
echo.

REM Check prerequisites
echo [1/4] Checking prerequisites...
echo.

REM Check Node.js
node --version >nul 2>&1
if %ERRORLEVEL% NEQ 0 (
    echo [ERROR] Node.js not installed!
    echo Install from: https://nodejs.org/
    pause
    exit /b 1
)
echo   [OK] Node.js found:
node --version

REM Check Python
python --version >nul 2>&1
if %ERRORLEVEL% NEQ 0 (
    echo [ERROR] Python not installed!
    echo Install from: https://www.python.org/
    pause
    exit /b 1
)
echo   [OK] Python found:
python --version

REM Check Docker
docker --version >nul 2>&1
if %ERRORLEVEL% NEQ 0 (
    echo [WARNING] Docker not installed or not running
    echo mol2chemfig will not start
)

echo.
echo ============================================================
echo [2/4] Creating .env file if needed...
echo ============================================================
echo.

if not exist ".env" (
    (
        echo # mol2chemfig ports
        echo BACKEND_PORT=8000
        echo FRONTEND_PORT=8080
    ) > .env
    echo [CREATED] .env file
) else (
    echo [OK] .env file exists
)

echo.
echo ============================================================
echo [3/4] Killing any existing servers...
echo ============================================================
echo.

taskkill /F /FI "WINDOWTITLE eq *MoleculeViewer*" 2>nul
taskkill /F /FI "WINDOWTITLE eq *PubChem*" 2>nul
taskkill /F /FI "WINDOWTITLE eq *mol2chemfig*" 2>nul
docker-compose down 2>nul

echo [OK] Cleaned up old processes

echo.
echo ============================================================
echo [4/4] Starting all servers...
echo ============================================================
echo.

REM Start MoleculeViewer (Port 5000)
echo [Starting] MoleculeViewer on port 5000...
start "MoleculeViewer Server - Port 5000" cmd /k "cd /d "%~dp0\MoleculeViewer" && node server.js"
timeout /t 2 /nobreak >nul

REM Start mol2chemfig Docker (Port 8000)
echo [Starting] mol2chemfig Docker on port 8000...
start "mol2chemfig Docker - Port 8000" cmd /k "cd /d "%~dp0" && docker-compose up"
timeout /t 3 /nobreak >nul

REM Start PubChem Server (Port 5002)
echo [Starting] PubChem Server on port 5002...
start "PubChem Server - Port 5002" cmd /k "cd /d "%~dp0" && python pubchem_server.py"
timeout /t 2 /nobreak >nul

echo.
echo ============================================================
echo   ALL SERVERS STARTED!
echo ============================================================
echo.
echo   You should see 3 new command windows:
echo     1. MoleculeViewer Server - Port 5000
echo     2. mol2chemfig Docker - Port 8000
echo     3. PubChem Server - Port 5002
echo.
echo   Testing servers in 5 seconds...
echo.
timeout /t 5 /nobreak >nul

echo ============================================================
echo   HEALTH CHECKS
echo ============================================================
echo.

REM Test each server
echo Testing MoleculeViewer (port 5000)...
curl -s http://localhost:5000/health >nul 2>&1
if %ERRORLEVEL% EQU 0 (
    echo   [OK] MoleculeViewer is responding
) else (
    echo   [FAIL] MoleculeViewer not responding
)

echo.
echo Testing mol2chemfig (port 8000)...
curl -s http://localhost:8000/ >nul 2>&1
if %ERRORLEVEL% EQU 0 (
    echo   [OK] mol2chemfig is responding
) else (
    echo   [FAIL] mol2chemfig not responding (Docker may still be starting...)
)

echo.
echo Testing PubChem (port 5002)...
curl -s http://localhost:5002/health >nul 2>&1
if %ERRORLEVEL% EQU 0 (
    echo   [OK] PubChem is responding
) else (
    echo   [FAIL] PubChem not responding
)

echo.
echo ============================================================
echo   QUICK TEST LINKS
echo ============================================================
echo.
echo   MoleculeViewer:
echo     http://localhost:5000/img/smiles?smiles=CCO
echo.
echo   mol2chemfig:
echo     http://localhost:8000/
echo.
echo   PubChem:
echo     http://localhost:5002/health
echo.
echo ============================================================
echo   TO STOP ALL SERVERS
echo ============================================================
echo.
echo   Run: stop_all.bat
echo   OR close each server window
echo.
echo ============================================================
echo.
echo   Press any key to close this window...
echo   (Servers will keep running in their own windows)
echo.
pause >nul
