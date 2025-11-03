@echo off
REM MoleculeViewer - Start Script
REM This script sets up and runs the MoleculeViewer Flask application

setlocal enabledelayedexpansion

echo.
echo ================================
echo   MoleculeViewer - Startup
echo ================================
echo.

REM Change to project directory
cd /d "%~dp0"
echo Project Directory: %CD%
echo.

REM Check if Python is installed
python --version >nul 2>&1
if errorlevel 1 (
    echo [ERROR] Python is not installed or not in PATH
    echo Please install Python 3.9+ from https://www.python.org/
    pause
    exit /b 1
)

echo [OK] Python is installed
python --version
echo.

REM Create virtual environment if it doesn't exist
if not exist "venv" (
    echo [*] Creating virtual environment...
    python -m venv venv
    if errorlevel 1 (
        echo [ERROR] Failed to create virtual environment
        pause
        exit /b 1
    )
    echo [OK] Virtual environment created
    echo.
)

REM Activate virtual environment
echo [*] Activating virtual environment...
call venv\Scripts\activate.bat
if errorlevel 1 (
    echo [ERROR] Failed to activate virtual environment
    pause
    exit /b 1
)
echo [OK] Virtual environment activated
echo.

REM Install/upgrade requirements
echo [*] Installing/upgrading dependencies...
echo    (This may take several minutes on the first run)
pip install --upgrade pip setuptools wheel --quiet >nul 2>&1

pip install -r requirements.txt
if errorlevel 1 (
    echo [ERROR] Failed to install requirements
    pause
    exit /b 1
)
echo [OK] Dependencies installed
echo.

REM Start the Flask application
echo ================================
echo   Starting MoleculeViewer
echo ================================
echo.
echo Application will be available at:
echo   http://localhost:5000
echo.
echo Press Ctrl+C to stop the server
echo.

python run.py

REM If the app exits, pause so user can see any error messages
pause
