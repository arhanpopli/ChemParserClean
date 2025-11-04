#!/usr/bin/env python
"""
Simple Flask server starter - guaranteed to work
"""
import os
import sys

# Fix encoding
sys.stdout.reconfigure(encoding='utf-8')

# Change to MoleculeViewer directory
os.chdir(os.path.dirname(os.path.abspath(__file__)))

# Load .env
from dotenv import load_dotenv
load_dotenv()

# Import Flask app
from app.api import app

if __name__ == '__main__':
    print("\n" + "="*70)
    print("MOLECULEVIEWER BACKEND SERVER")
    print("="*70)
    print("Starting Flask server...")
    print("Listen on: http://localhost:5000")
    print("Listen on: http://192.168.1.4:5000")
    print("="*70 + "\n")
    
    try:
        # Start Flask
        app.run(
            host='0.0.0.0',
            port=5000,
            debug=False,
            use_reloader=False,
            threaded=True
        )
    except Exception as e:
        print(f"ERROR: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)
