@echo off
REM ============================================================
REM START ALL CHEMPARSER SERVERS
REM Starts all three servers: PubChem, MoleculeViewer, Mol2ChemFig
REM ============================================================

title ChemParser - Start All Servers

cd /d "%~dp0"

echo.
echo ============================================================
echo   CHEMPARSER - START ALL SERVERS
echo ============================================================
echo.

REM Kill any existing processes
echo Killing existing processes...
taskkill /F /IM node.exe >nul 2>&1
taskkill /F /IM python.exe >nul 2>&1
timeout /t 2 >nul

echo.
echo Starting servers...
echo.

REM Terminal 1: PubChem Server (Node.js on port 5002)
echo [1/3] Starting PubChem Server (Port 5002)...
start "PubChem Server (5002)" cmd /k cd MoleculeViewer\pubchem ^& node server.js
timeout /t 3 >nul

REM Terminal 2: MoleculeViewer Server (Node.js on port 5000)
echo [2/3] Starting MoleculeViewer Server (Port 5000)...
start "MoleculeViewer (5000)" cmd /k cd MoleculeViewer ^& node server.js
timeout /t 3 >nul

REM Terminal 3: Mol2ChemFig Server (Python on port 5001)
echo [3/3] Starting Mol2ChemFig Server (Port 5001)...
start "Mol2ChemFig (5001)" cmd /k python mol2chemfig_server.py
timeout /t 2 >nul

echo.
echo ============================================================
echo   ALL SERVERS STARTED!
echo ============================================================
echo.
echo Services now available:
echo   PubChem (3D):        http://localhost:5002
echo   MoleculeViewer:      http://localhost:5000
echo   Mol2ChemFig:         http://localhost:5001
echo.
echo DO NOT CLOSE THESE WINDOWS - Servers will stop!
echo.
pause

REM Get script directory
set SCRIPT_DIR=%~dp0
cd /d "%SCRIPT_DIR%"

echo [1/4] Checking prerequisites...
echo.

REM Check Python
python --version >nul 2>&1
if errorlevel 1 (
    echo ERROR: Python not found!
    echo Please install Python 3.8+ from python.org
    pause
    exit /b 1
)
for /f "tokens=*" %%i in ('python --version 2^>^&1') do set PYTHON_VERSION=%%i
echo   ✓ %PYTHON_VERSION%

REM Check Node.js
node --version >nul 2>&1
if errorlevel 1 (
    echo ERROR: Node.js not found!
    echo Please install Node.js from nodejs.org
    pause
    exit /b 1
)
for /f "tokens=*" %%i in ('node --version 2^>^&1') do set NODE_VERSION=%%i
echo   ✓ Node.js %NODE_VERSION%

REM Check Docker
docker --version >nul 2>&1
if errorlevel 1 (
    echo WARNING: Docker not found!
    echo Mol2ChemFig Docker backend will not start.
    echo Install Docker Desktop to use this feature.
    set DOCKER_AVAILABLE=0
) else (
    for /f "tokens=*" %%i in ('docker --version 2^>^&1') do set DOCKER_VERSION=%%i
    echo   ✓ !DOCKER_VERSION!
    set DOCKER_AVAILABLE=1
)

echo.
echo [2/4] Cleaning up existing processes...
echo.

REM Stop any existing servers
taskkill /F /IM python.exe /FI "WINDOWTITLE eq *ChemParser*" >nul 2>&1
taskkill /F /IM node.exe /FI "WINDOWTITLE eq *ChemParser*" >nul 2>&1
timeout /t 2 >nul

echo   ✓ Cleanup complete
echo.
echo [3/4] Starting Docker backend...
echo.

if %DOCKER_AVAILABLE%==1 (
    REM Check if .env file exists, create if not
    if not exist ".env" (
        echo Creating .env file...
        (
            echo BACKEND_CONTAINER_PORT=8000
            echo BACKEND_HOST_PORT=8000
            echo FRONTEND_CONTAINER_PORT=80
            echo FRONTEND_HOST_PORT=8080
        ) > .env
    )

    REM Start Docker containers
    docker-compose up -d
    if errorlevel 1 (
        echo WARNING: Failed to start Docker containers
        echo Please check Docker Desktop is running
    ) else (
        echo   ✓ Docker backend started on port 8000
    )
    timeout /t 3 >nul
) else (
    echo   ⚠ Skipping Docker backend (not available)
)

echo.
echo [4/4] Starting application servers...
echo.

REM Start MoleculeViewer (Node.js on port 5000)
echo Starting MoleculeViewer Server...
cd /d "%SCRIPT_DIR%MoleculeViewer"
start "ChemParser - MoleculeViewer (Port 5000)" cmd /k "node server.js"
cd /d "%SCRIPT_DIR%"
timeout /t 2 >nul
echo   ✓ MoleculeViewer started on port 5000

REM Start Mol2ChemFig Server (Flask on port 5001)
echo Starting Mol2ChemFig Server...
start "ChemParser - Mol2ChemFig (Port 5001)" cmd /k "python mol2chemfig_server.py"
timeout /t 2 >nul
echo   ✓ Mol2ChemFig Server started on port 5001

REM Start PubChem Server (Flask on port 5002)
echo Starting PubChem Server...
start "ChemParser - PubChem (Port 5002)" cmd /k "python pubchem_server.py"
timeout /t 2 >nul
echo   ✓ PubChem Server started on port 5002

echo.
echo ============================================================
echo   ALL SERVERS STARTED SUCCESSFULLY!
echo ============================================================
echo.
echo   Server Status:
if %DOCKER_AVAILABLE%==1 (
    echo   [✓] Docker Backend:      http://localhost:8000
)
echo   [✓] MoleculeViewer:      http://localhost:5000
echo   [✓] Mol2ChemFig Server:  http://localhost:5001
echo   [✓] PubChem Server:      http://localhost:5002
echo.
echo   Quick Links:
echo   - Master Launcher:       launcher.html
echo   - Test Mol2ChemFig:      http://localhost:5001
echo   - Test PubChem:          http://localhost:5002
echo   - Chrome Extension:      Load chem-extension folder
echo.
echo   To stop all servers, run: stop_all.bat
echo   To check status, run: status.bat
echo ============================================================
echo.

REM Open the launcher in default browser
echo Opening Master Launcher...
start launcher.html

pause
