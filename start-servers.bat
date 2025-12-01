@echo off
echo ============================================================
echo   CHEMPARSER - Starting All Servers
echo ============================================================

REM Kill existing processes to avoid port conflicts
echo [PRE] Cleaning up existing processes...
taskkill /F /IM node.exe >nul 2>&1
taskkill /F /IM python.exe >nul 2>&1
taskkill /F /IM php.exe >nul 2>&1
timeout /t 2 >nul

REM Start Mol2ChemFig Server
echo [1/6] Starting Mol2ChemFig Server (Port 5001)...
start "Mol2ChemFig-5001" cmd /k "cd /d C:\Users\Kapil\Personal\PROJECTS\Chemparser && python mol2chemfig_server.py"
timeout /t 2 /nobreak >nul

REM Start MoleculeViewer Server
echo [2/6] Starting MoleculeViewer Server (Port 5000)...
start "MoleculeViewer-5000" cmd /k "cd /d C:\Users\Kapil\Personal\PROJECTS\Chemparser\MoleculeViewer && node server.js"
timeout /t 2 /nobreak >nul

REM Start PubChem Server
echo [3/6] Starting PubChem Server (Port 5002)...
start "PubChem-5002" cmd /k "cd /d C:\Users\Kapil\Personal\PROJECTS\Chemparser\MoleculeViewer\pubchem && node server.js"
timeout /t 2 /nobreak >nul

REM Start MolView PHP Server (Port 8000) - REQUIRED for embeds
echo [4/6] Starting MolView PHP Server (Port 8000)...
start "MolView-PHP-8000" cmd /k "cd /d C:\Users\Kapil\Personal\PROJECTS\Chemparser\Molview\molview && php -S localhost:8000"
timeout /t 3 /nobreak >nul

REM Start MolView Search API (Port 8001) - REQUIRED for autocorrect
echo [5/6] Starting MolView Search API (Port 8001)...
start "MolView-Search-8001" cmd /k "cd /d C:\Users\Kapil\Personal\PROJECTS\Chemparser\Molview\molview && node search-server.js"
timeout /t 3 /nobreak >nul

REM Open browser
echo [6/6] Opening Unified Interface...
timeout /t 2 /nobreak >nul
start http://localhost:5000/unified-interface.html

echo.
echo ============================================================
echo   ALL SERVERS STARTED!
echo ============================================================
echo.
echo   Mol2ChemFig:       http://localhost:5001
echo   MoleculeViewer:    http://localhost:5000
echo   PubChem:           http://localhost:5002
echo   MolView PHP:       http://localhost:8000  (embeds)
echo   MolView Search:    http://localhost:8001  (autocorrect API)
echo.
echo   Extension uses:
echo   - Port 8001 for ALL queries (autocorrect + intelligent filtering)
echo   - Port 8000 for MolView embed URLs
echo.
echo   Test Search API: http://localhost:8001/search?q=aspirin
echo.
echo ============================================================
echo.
pause
