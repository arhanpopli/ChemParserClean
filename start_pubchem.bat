@echo off
echo ====================================
echo Starting PubChem Integration Server
echo ====================================
echo.
echo Installing/checking dependencies...
pip install -r requirements_pubchem.txt
echo.
echo Starting server on port 5002...
echo.
python pubchem_server.py
