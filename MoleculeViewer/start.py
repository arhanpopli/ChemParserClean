#!/usr/bin/env python
"""
Simple Python launcher for MoleculeViewer
Avoids batch file issues and starts Flask directly
"""
import subprocess
import sys
import os

# Change to project directory
os.chdir(os.path.dirname(os.path.abspath(__file__)))

# Run the Flask app
from app.api import app

if __name__ == '__main__':
    print("\n" + "="*50)
    print("  MoleculeViewer - Starting")
    print("="*50 + "\n")
    print("Application available at:")
    print("  http://localhost:5000")
    print("\nPress Ctrl+C to stop\n")
    
    app.run(debug=False, host='0.0.0.0', port=5000)
