@echo off
echo ============================================================
echo   MOLVIEW SERVERS - Starting PHP + Search API
echo ============================================================

REM Kill existing MolView processes
echo [PRE] Cleaning up existing PHP and Node processes...
taskkill /F /IM php.exe >nul 2>&1
taskkill /F /FI "WINDOWTITLE eq MolView*" >nul 2>&1
timeout /t 2 >nul

REM Start MolView PHP Server (Port 8000)
echo [1/2] Starting MolView PHP Server (Port 8000)...
echo       This serves the embed viewer at /embed/v2/
start "MolView-PHP-8000" cmd /k "cd /d C:\Users\Kapil\Personal\PROJECTS\Chemparser\Molview\molview && php -S localhost:8000"
timeout /t 3 /nobreak >nul

REM Start MolView Search API (Port 8001)
echo [2/2] Starting MolView Search API (Port 8001)...
echo       This provides autocorrect and intelligent filtering
start "MolView-Search-8001" cmd /k "cd /d C:\Users\Kapil\Personal\PROJECTS\Chemparser\Molview\molview && node search-server.js"
timeout /t 3 /nobreak >nul

echo.
echo ============================================================
echo   MOLVIEW SERVERS STARTED!
echo ============================================================
echo.
echo   MolView PHP Server:    http://localhost:8000
echo   MolView Search API:    http://localhost:8001
echo.
echo   Test URLs:
echo   - Main viewer:         http://localhost:8000/
echo   - Embed (compound):    http://localhost:8000/embed/v2/?cid=2244
echo   - Embed (protein):     http://localhost:8000/embed/v2/?pdbid=1RHV
echo   - Embed (mineral):     http://localhost:8000/embed/v2/?codid=9004137
echo   - Search API:          http://localhost:8001/search?q=aspirin
echo.
echo   Extension uses:
echo   - Port 8001: ALL queries (autocorrect + intelligent filtering)
echo   - Port 8000: Embed URLs returned by search API
echo.
echo ============================================================
echo.
pause
