@echo off
echo ============================================================
echo   CHEMPARSER - Starting All Servers
echo ============================================================

REM Start Docker (optional)
echo Starting Docker backend...
docker-compose up -d 2>nul

REM Start Mol2ChemFig Server
echo Starting Mol2ChemFig Server (Port 5001)...
start "Mol2ChemFig-5001" cmd /c "cd /d C:\Users\Kapil\Personal\PROJECTS\Chemparser && python mol2chemfig_server.py"

timeout /t 2 /nobreak >nul

REM Start MoleculeViewer Server
echo Starting MoleculeViewer Server (Port 5000)...
start "MoleculeViewer-5000" cmd /c "cd /d C:\Users\Kapil\Personal\PROJECTS\Chemparser\MoleculeViewer && node server.js"

timeout /t 2 /nobreak >nul

REM Start PubChem Server
echo Starting PubChem Server (Port 5002)...
start "PubChem-5002" cmd /c "cd /d C:\Users\Kapil\Personal\PROJECTS\Chemparser\MoleculeViewer\pubchem && node server.js"

timeout /t 3 /nobreak >nul

REM Open browser
echo Opening Unified Interface...
start http://localhost:5000/unified-interface.html

echo.
echo ============================================================
echo   All servers started!
echo ============================================================
echo.
pause
