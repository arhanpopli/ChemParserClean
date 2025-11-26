@echo off
echo Starting MolView Python Server...
echo.

REM Get the directory where this batch file is located
set SCRIPT_DIR=%~dp0

REM Change to the script directory
cd /d "%SCRIPT_DIR%"

echo Launching Python-based MolView server...
echo Server will be available at http://localhost:5000
echo.

REM Start the Python server
python molview_python_server.py

pause