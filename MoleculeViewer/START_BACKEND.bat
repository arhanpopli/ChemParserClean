@echo off
REM Start MoleculeViewer Flask Backend
REM This keeps the server running persistently

cd /d "C:\Users\Kapil\Personal\PROJECTS\Mol2chemfig\MoleculeViewer"

echo.
echo ========================================================================
echo  MOLECULEVIEWER BACKEND SERVER
echo ========================================================================
echo  Backend is starting on http://localhost:5000
echo  This window MUST STAY OPEN for the extension to work!
echo.
echo  DO NOT CLOSE THIS WINDOW!
echo ========================================================================
echo.

:loop
python start_server_simple.py
echo.
echo.
echo Restart in 3 seconds...
timeout /t 3
goto loop
