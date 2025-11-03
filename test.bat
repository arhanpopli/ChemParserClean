@echo off
REM Test Script for Mol2ChemFig SVG Pipeline
REM Tests all main endpoints and displays results

echo ============================================
echo   MOL2CHEMFIG - COMPLETE SYSTEM TEST
echo ============================================
echo.

REM Check if Docker is running
echo [1/5] Checking Docker containers...
docker-compose ps >nul 2>&1
if %errorlevel% neq 0 (
    echo ERROR: Docker not running. Run: docker-compose up -d
    pause
    exit /b 1
)
echo OK - Containers running
echo.

REM Test SVG endpoint
echo [2/5] Testing SMILES to SVG endpoint...
powershell -Command "$body = @{textAreaData='C1=CC=CC=C1'} | ConvertTo-Json; $r = Invoke-WebRequest -Uri 'http://localhost:8000/m2cf/svg' -Method POST -Headers @{'Content-Type'='application/json'} -Body $body -ErrorAction SilentlyContinue; if ($r.StatusCode -eq 200) { Write-Host 'OK - SVG generated'; $j = $r.Content | ConvertFrom-Json; Write-Host ('Size: ' + $j.svg.Length + ' bytes') } else { Write-Host 'FAILED' }"
echo.

REM Test nomenclature endpoint
echo [3/5] Testing nomenclature to SVG endpoint...
powershell -Command "$body = @{nomenclature='water'} | ConvertTo-Json; $r = Invoke-WebRequest -Uri 'http://localhost:8000/m2cf/nomenclature-to-svg' -Method POST -Headers @{'Content-Type'='application/json'} -Body $body -ErrorAction SilentlyContinue; if ($r.StatusCode -eq 200) { $j = $r.Content | ConvertFrom-Json; if ($j.smiles) { Write-Host ('OK - Found: ' + $j.smiles) } else { Write-Host ('Status 200 but: ' + $j.error) } } else { Write-Host 'HTTP ' $r.StatusCode }"
echo.

REM Test submit endpoint
echo [4/5] Testing submit endpoint...
powershell -Command "$body = @{textAreaData='CC(=O)O'} | ConvertTo-Json; $r = Invoke-WebRequest -Uri 'http://localhost:8000/m2cf/submit' -Method POST -Headers @{'Content-Type'='application/json'} -Body $body -ErrorAction SilentlyContinue; if ($r.StatusCode -eq 200) { Write-Host 'OK - LaTeX + PDF generated' } else { Write-Host 'FAILED' }"
echo.

REM Test search endpoint
echo [5/5] Testing search endpoint...
powershell -Command "$body = @{searchTerm='water'} | ConvertTo-Json; $r = Invoke-WebRequest -Uri 'http://localhost:8000/m2cf/search' -Method POST -Headers @{'Content-Type'='application/json'} -Body $body -ErrorAction SilentlyContinue; if ($r.StatusCode -eq 200) { $j = $r.Content | ConvertFrom-Json; if ($j.smiles) { Write-Host ('OK - Found: ' + $j.smiles) } else { Write-Host ('Status 200 but: ' + $j.error) } } else { Write-Host 'HTTP ' $r.StatusCode }"
echo.

echo ============================================
echo   ALL TESTS COMPLETE
echo ============================================
echo.
echo Frontend available at: http://localhost:8080
echo Backend API at: http://localhost:8000
echo.
echo Next steps:
echo 1. Read INTEGRATION_GUIDE.md to integrate in your app
echo 2. Read API_ENDPOINTS.md for detailed endpoint docs
echo 3. Read FINAL_SUMMARY.md for complete overview
echo.
pause
