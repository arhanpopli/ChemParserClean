"""
Main entry point for MoleculeViewer application.
"""

import os
import sys
from dotenv import load_dotenv

# Fix encoding for emoji printing on Windows
sys.stdout.reconfigure(encoding='utf-8')

# Load environment variables from .env file
load_dotenv()

from app.api import app

if __name__ == '__main__':
    print("\n" + "="*60)
    print("MoleculeViewer Server Starting")
    print("="*60)
    print("Server URL: http://192.168.1.4:5000")
    print("Local: http://localhost:5000")
    print("Worldwide accessible links will be generated")
    print("="*60 + "\n")
    app.run(debug=False, host='0.0.0.0', port=5000, use_reloader=False)
