#!/usr/bin/env powershell
<#
.SYNOPSIS
    Quick test script for MoleculeViewer Node.js Server
.DESCRIPTION
    Tests SMILES and nomenclature endpoints
#>

Write-Host "🧪 Testing MoleculeViewer Node.js Server" -ForegroundColor Cyan
Write-Host "========================================`n" -ForegroundColor Cyan

$baseUrl = "http://localhost:5000"

# Test 1: Health check
Write-Host "Test 1: Health Check" -ForegroundColor Yellow
try {
    $response = Invoke-WebRequest -Uri "$baseUrl/health" -UseBasicParsing
    Write-Host "✅ Server is healthy: $($response.StatusCode)" -ForegroundColor Green
    Write-Host $response.Content | ConvertFrom-Json | ConvertTo-Json
} catch {
    Write-Host "❌ Health check failed: $_" -ForegroundColor Red
}

Write-Host "`n"

# Test 2: SMILES endpoint
Write-Host "Test 2: SMILES Endpoint (Ethanol: CCO)" -ForegroundColor Yellow
try {
    $response = Invoke-WebRequest -Uri "$baseUrl/img/smiles?smiles=CCO&width=300&height=200" -UseBasicParsing
    Write-Host "✅ SMILES endpoint works: $($response.StatusCode)" -ForegroundColor Green
    Write-Host "Content-Type: $($response.Headers['Content-Type'])" -ForegroundColor Cyan
    Write-Host "Response size: $($response.Content.Length) bytes" -ForegroundColor Cyan
} catch {
    Write-Host "❌ SMILES endpoint failed: $_" -ForegroundColor Red
}

Write-Host "`n"

# Test 3: Nomenclature endpoint
Write-Host "Test 3: Nomenclature Endpoint (Acetone)" -ForegroundColor Yellow
try {
    $response = Invoke-WebRequest -Uri "$baseUrl/img/nomenclature?nomenclature=acetone&width=300&height=200" -UseBasicParsing
    Write-Host "✅ Nomenclature endpoint works: $($response.StatusCode)" -ForegroundColor Green
    Write-Host "Content-Type: $($response.Headers['Content-Type'])" -ForegroundColor Cyan
    Write-Host "Response size: $($response.Content.Length) bytes" -ForegroundColor Cyan
} catch {
    Write-Host "❌ Nomenclature endpoint failed: $_" -ForegroundColor Red
}

Write-Host "`n"

# Test 4: Cache info
Write-Host "Test 4: Cache Information" -ForegroundColor Yellow
try {
    $response = Invoke-WebRequest -Uri "$baseUrl/cache-info" -UseBasicParsing
    Write-Host "✅ Cache info endpoint works: $($response.StatusCode)" -ForegroundColor Green
    $cacheInfo = $response.Content | ConvertFrom-Json
    Write-Host "Cached SVGs: $($cacheInfo.cachedSvgs)" -ForegroundColor Cyan
    Write-Host "Cache size: $($cacheInfo.totalCacheSize)" -ForegroundColor Cyan
} catch {
    Write-Host "❌ Cache info failed: $_" -ForegroundColor Red
}

Write-Host "`n========================================" -ForegroundColor Cyan
Write-Host "✅ All tests completed!" -ForegroundColor Green
Write-Host "`n📝 Next Steps:" -ForegroundColor Cyan
Write-Host "1. Reload the Chrome extension"
Write-Host "2. Test in ChatGPT with: chem:acetone or chem:CCO"
Write-Host "3. Images should render inline!" -ForegroundColor Green
