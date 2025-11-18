@echo off
REM Setup MCP Server for autonomous development

echo.
echo ============================================
echo  MCP Server Setup
echo ============================================
echo.

REM Check if Python is available
python --version >nul 2>&1
if errorlevel 1 (
    echo ERROR: Python not found. Please install Python.
    pause
    exit /b 1
)

echo [OK] Python found
echo.

REM Set paths
set PROJECT_ROOT=C:\Users\Kapil\Personal\PROJECTS\Mol2chemfig
set FLASK_PATH=%PROJECT_ROOT%\MoleculeViewer

echo Project root: %PROJECT_ROOT%
echo Flask path: %FLASK_PATH%
echo.

REM Test the MCP server
echo Testing MCP server...
python %PROJECT_ROOT%\mcp_server.py check_status
if errorlevel 1 (
    echo ERROR: MCP server test failed
    pause
    exit /b 1
)

echo.
echo [OK] MCP Server is ready!
echo.
echo Next steps:
echo 1. Run: python mcp_server.py start_flask
echo 2. Or run: python mcp_server.py run_tests
echo 3. Or run: python mcp_server.py check_status
echo.
pause
