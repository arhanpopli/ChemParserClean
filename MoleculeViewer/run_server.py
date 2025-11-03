#!/usr/bin/env python
"""Run MoleculeViewer server"""
import sys
import os

# Add the current directory to the path
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

try:
    from app.api import app
    
    if __name__ == '__main__':
        print("\n" + "="*50)
        print("  MoleculeViewer Flask Server")
        print("="*50)
        print("\nServer available at:")
        print("  http://localhost:5000")
        print("  http://127.0.0.1:5000")
        print("\nPress CTRL+C to stop\n")
        
        import socket
        sock = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
        result = sock.connect_ex(('127.0.0.1', 5000))
        sock.close()
        if result == 0:
            print("WARNING: Port 5000 is already in use!")
            sys.exit(1)
        
        app.run(debug=False, host='0.0.0.0', port=5000, use_reloader=False, threaded=True)
        
except Exception as e:
    print(f"\nERROR: {e}", file=sys.stderr)
    import traceback
    traceback.print_exc()
    sys.exit(1)
