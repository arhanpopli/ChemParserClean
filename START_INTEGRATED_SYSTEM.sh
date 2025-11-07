#!/usr/bin/env bash
# Quick Start: MoleculeViewer + Mol2ChemFig Integration

echo "🚀 Starting Integrated MoleculeViewer + Mol2ChemFig System"
echo "=========================================================="
echo ""

# Check if Docker is running
echo "✓ Checking Docker containers..."
if ! docker ps | grep -q "m2cf_backend"; then
    echo "⚠ Mol2ChemFig Docker backend not running!"
    echo "  Starting Docker..."
    cd "C:\Users\Kapil\Personal\PROJECTS\Mol2chemfig"
    docker-compose up -d
    sleep 5
fi

# Check if backend is responsive
echo "✓ Verifying mol2chemfig backend (localhost:8000)..."
if ! curl -s http://localhost:8000/m2cf/submit > /dev/null 2>&1; then
    echo "⚠ Backend not responding yet, waiting..."
    sleep 3
fi

# Start MoleculeViewer server
echo "✓ Starting MoleculeViewer on localhost:5000..."
cd "C:\Users\Kapil\Personal\PROJECTS\Mol2chemfig\MoleculeViewer"

# Kill any existing Python processes on port 5000
if lsof -Pi :5000 -sTCP:LISTEN -t >/dev/null 2>&1; then
    echo "  Stopping existing server..."
    lsof -ti:5000 | xargs kill -9
    sleep 1
fi

# Start the server in background
python run_server.py &
SERVER_PID=$!
sleep 3

echo ""
echo "🎉 System is ready!"
echo "=========================================================="
echo ""
echo "📊 MoleculeViewer with Mol2ChemFig Integration"
echo ""
echo "  URL: http://localhost:5000/"
echo ""
echo "  Features:"
echo "  - 📊 MoleculeViewer Tab: RDKit-based SMILES to SVG"
echo "  - 🧬 Mol2ChemFig Tab: LaTeX-quality ChemFig rendering"
echo ""
echo "  Backend Services:"
echo "  - MoleculeViewer: http://localhost:5000 ✓"
echo "  - Mol2ChemFig Docker: http://localhost:8000 ✓"
echo ""
echo "=========================================================="
echo ""
echo "Press CTRL+C to stop the server"
echo ""

wait $SERVER_PID
