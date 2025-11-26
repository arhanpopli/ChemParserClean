@echo off
echo Starting MolView Local Server...
echo.

REM Get the directory where this batch file is located
set SCRIPT_DIR=%~dp0

REM Change to the MolView directory
cd /d "%SCRIPT_DIR%Molview\molview-2.4.6 - github repo, src"

echo Checking PHP installation...
php -v >nul 2>&1
if errorlevel 1 (
    echo Error: PHP is not installed or not in PATH
    echo Please install PHP and add it to your system PATH
    pause
    exit /b 1
)

echo Found PHP installation
echo.

REM Start PHP built-in server
echo Starting MolView server on port 8080...
echo Open your browser and go to http://localhost:8080
echo Press Ctrl+C to stop the server
echo.

php -S localhost:8080

pause