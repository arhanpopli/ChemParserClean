@echo off
title ChemistryLaTeX Server (Vercel)
echo.
echo ============================================
echo   ChemistryLaTeX Server - Vercel Dev Environment
echo ============================================
echo.
echo This will start the local Vercel development server.
echo Make sure you have the Vercel CLI installed (npx vercel).
echo.

cd /d "%~dp0server"
npx vercel dev

pause
