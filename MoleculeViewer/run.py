"""
Main entry point for MoleculeViewer application.
"""

import os
from dotenv import load_dotenv

# Load environment variables from .env file
load_dotenv()

from app.api import app

if __name__ == '__main__':
    print("\n" + "="*60)
    print("🚀 MoleculeViewer Server Starting")
    print("="*60)
    print("📍 Server: http://localhost:5000")
    print("🔗 Worldwide accessible cache links available")
    print("="*60 + "\n")
    
    try:
        app.run(debug=False, host='0.0.0.0', port=5000, use_reloader=False, threaded=True)
    except OSError as e:
        if "Address already in use" in str(e):
            print("\n❌ ERROR: Port 5000 is already in use!")
            print("Please close the other application using port 5000")
            import sys
            sys.exit(1)
        else:
            raise
    except Exception as e:
        print(f"\n❌ ERROR: {str(e)}")
        import sys
        sys.exit(1)
