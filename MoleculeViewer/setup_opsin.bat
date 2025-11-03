@echo off
REM Setup OPSIN parser for MoleculeViewer
REM This script downloads and sets up OPSIN for IUPAC nomenclature parsing

echo.
echo =====================================================
echo OPSIN Parser Setup for MoleculeViewer
echo =====================================================
echo.

REM Check if Java is installed
java -version >nul 2>&1
if errorlevel 1 (
    echo.
    echo ERROR: Java is not installed or not in PATH
    echo.
    echo OPSIN requires Java to run. Please install Java from:
    echo https://www.java.com/en/download/
    echo.
    echo After installing Java, run this script again.
    echo.
    pause
    exit /b 1
)

echo ✓ Java found
java -version 2>&1 | findstr /r "version"
echo.

REM Check if opsin-cli.jar already exists
if exist "opsin-cli.jar" (
    echo ✓ OPSIN parser already present (opsin-cli.jar)
    echo.
    echo Setup complete!
    echo.
    pause
    exit /b 0
)

echo Downloading OPSIN parser...
echo.

REM Use Python to download the file (more reliable)
python -c "
import urllib.request
import sys

url = 'https://files.pythonhosted.org/packages/1e/ee/b44a1aa481bcbf7d3abffc51cf6de9c9fb6b6a5a1d5a5e5f5f5f5f5f5f5f/opsin-cli-2.8.0.jar'
filename = 'opsin-cli.jar'

try:
    print(f'Downloading from {url}...')
    urllib.request.urlretrieve(url, filename)
    print(f'✓ Downloaded successfully to {filename}')
except Exception as e:
    print(f'✗ Download failed: {e}')
    print('Please download manually from:')
    print('https://github.com/opsin/opsin/releases')
    sys.exit(1)
" 2>&1

if errorlevel 1 (
    echo.
    echo Failed to download. Trying alternative method...
    echo.
    powershell -Command "(New-Object Net.WebClient).DownloadFile('https://github.com/opsin/opsin/releases/download/v2.8.0/opsin-cli-2.8.0.jar', 'opsin-cli.jar')" 2>nul
    
    if not exist "opsin-cli.jar" (
        echo.
        echo DOWNLOAD FAILED
        echo.
        echo Please download OPSIN manually:
        echo 1. Visit: https://github.com/opsin/opsin/releases
        echo 2. Download: opsin-cli-*.jar
        echo 3. Copy to this directory
        echo 4. Rename to: opsin-cli.jar
        echo.
        pause
        exit /b 1
    )
)

echo.
echo Testing OPSIN...
echo.

java -jar opsin-cli.jar ethanol >nul 2>&1
if errorlevel 1 (
    echo Warning: OPSIN test failed
) else (
    echo ✓ OPSIN test successful
)

echo.
echo =====================================================
echo OPSIN Setup Complete!
echo =====================================================
echo.
echo You can now use IUPAC nomenclature parsing.
echo.
echo Configure in app/config.py:
echo   NOMENCLATURE_PARSER = 'auto'  (or 'opsin' for IUPAC only)
echo.
pause
