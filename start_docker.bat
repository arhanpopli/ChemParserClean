@echo off
REM ============================================================
REM Start Mol2ChemFig Docker Backend (Port 8000)
REM ============================================================

title ChemParser - Docker Backend (Port 8000)

set SCRIPT_DIR=%~dp0
cd /d "%SCRIPT_DIR%"

echo.
echo ============================================================
echo   Starting Mol2ChemFig Docker Backend
echo ============================================================
echo   Port: 8000
echo   Type: Docker Container
echo ============================================================
echo.

REM Check if .env exists, create if not
if not exist ".env" (
    echo Creating .env file...
    (
        echo BACKEND_CONTAINER_PORT=8000
        echo BACKEND_HOST_PORT=8000
        echo FRONTEND_CONTAINER_PORT=80
        echo FRONTEND_HOST_PORT=8080
    ) > .env
    echo .env file created
    echo.
)

echo Starting Docker containers...
docker-compose up -d

if errorlevel 1 (
    echo.
    echo ERROR: Failed to start Docker containers
    echo.
    echo Please check:
    echo 1. Docker Desktop is running
    echo 2. No other service is using port 8000
    echo 3. Docker Compose is installed
    pause
    exit /b 1
)

echo.
echo ============================================================
echo   Docker Backend Started Successfully
echo ============================================================
echo   Backend URL: http://localhost:8000
echo   Status: docker ps
echo ============================================================
echo.

docker ps --filter "name=m2cf"

echo.
echo Press any key to view logs (Ctrl+C to exit logs)...
pause >nul

docker-compose logs -f
