@echo off
REM Docker Compose Launcher for Mol2chemfig
REM Starts the Mol2chemfig application using Docker Compose

echo.
echo ================================================
echo  Mol2chemfig Docker Compose Launcher
echo ================================================
echo.

REM Change to the project directory
cd /d "%~dp0"

REM Check if Docker is installed
docker --version >nul 2>&1
if errorlevel 1 (
    echo ERROR: Docker is not installed or not in PATH
    echo Please install Docker Desktop: https://www.docker.com/products/docker-desktop
    pause
    exit /b 1
)

REM Check if Docker daemon is running
docker ps >nul 2>&1
if errorlevel 1 (
    echo.
    echo ERROR: Docker daemon is not running!
    echo.
    echo Please start Docker Desktop and try again.
    echo.
    pause
    exit /b 1
)

echo [*] Docker is ready
echo.
echo [*] Starting Mol2chemfig services via Docker Compose...
echo.
echo Services will be available at:
echo   - Frontend: http://localhost:8080
echo   - Backend:  http://localhost:8000
echo.
echo [*] Press CTRL+C to stop all services
echo.

REM Start docker-compose
docker-compose up

echo.
echo Docker Compose stopped.
pause
