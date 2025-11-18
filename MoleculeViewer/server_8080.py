#!/usr/bin/env python
"""
MoleculeViewer Server Launcher - Port 8080
Starts the Flask app on port 8080 instead of 5000
"""
import sys
import os

sys.path.insert(0, os.getcwd())
os.chdir(os.path.dirname(os.path.abspath(__file__)))

from app.api import app

print("\n" + "="*60)
print("  MoleculeViewer Flask Server - Port 8080")
print("="*60)
print("\nServer running on:")
print("  Local:   http://127.0.0.1:8080")
print("  Network: http://0.0.0.0:8080")
print("\nPress CTRL+C to stop\n")

try:
    app.run(
        host='0.0.0.0',
        port=8080,
        debug=False,
        use_reloader=False,
        threaded=True
    )
except KeyboardInterrupt:
    print("\n\nServer stopped.")
    sys.exit(0)
