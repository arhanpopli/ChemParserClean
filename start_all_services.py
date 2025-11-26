#!/usr/bin/env python3
"""
Unified Startup Script for ChemParser Suite
Starts all servers together: MolView, Mol2ChemFig, and PubChem
"""

import os
import sys
import subprocess
import threading
import time
import requests
import webbrowser
from pathlib import Path
import signal
import atexit

# Add project root to Python path
project_root = Path(__file__).parent
sys.path.insert(0, str(project_root))

# Global variables for server processes
SERVER_PROCESSES = []

def start_molview_server():
    """Start the MolView Python server"""
    try:
        # Check if Python MolView server exists
        molview_server_path = project_root / "molview_python_server.py"
        if molview_server_path.exists():
            print("[GAME_CONTROLLER] Starting MolView Python server on port 5000...")

            # Run the Python-based MolView server
            proc = subprocess.Popen([
                sys.executable, str(molview_server_path)
            ], stdout=subprocess.PIPE, stderr=subprocess.PIPE)

            SERVER_PROCESSES.append(proc)
            print("[CHECK] MolView server started on http://localhost:5000")
            return proc
        else:
            print("[WARNING] MolView Python server not found, attempting PHP version...")

            # Try to start PHP version if available
            molview_dir = project_root / "Molview" / "molview-2.4.6 - github repo, src"
            if molview_dir.exists():
                try:
                    # Check if PHP is available
                    result = subprocess.run(['php', '--version'],
                                          capture_output=True, text=True, timeout=5)
                    if result.returncode == 0:
                        print("[CHECK] PHP found, starting PHP MolView server...")
                        proc = subprocess.Popen([
                            'php', '-S', 'localhost:5000'
                        ], cwd=str(molview_dir), stdout=subprocess.PIPE, stderr=subprocess.PIPE)

                        SERVER_PROCESSES.append(proc)
                        print("[CHECK] PHP MolView server started on http://localhost:5000")
                        return proc
                    else:
                        print("[ERROR] PHP not found - MolView server unavailable")
                        return None
                except (subprocess.CalledProcessError, FileNotFoundError):
                    print("[ERROR] PHP not found - MolView server unavailable")
                    return None
            else:
                print("[ERROR] MolView directory not found - downloading required")
                return None
    except Exception as e:
        print(f"[ERROR] Error starting MolView server: {e}")
        return None

def start_mol2chemfig_server():
    """Start the Mol2ChemFig server"""
    try:
        mol2chemfig_server_path = project_root / "mol2chemfig_server.py"
        if mol2chemfig_server_path.exists():
            print("[TEST_TUBE] Starting Mol2ChemFig server on port 5001...")

            proc = subprocess.Popen([
                sys.executable, str(mol2chemfig_server_path)
            ], stdout=subprocess.PIPE, stderr=subprocess.PIPE)

            SERVER_PROCESSES.append(proc)
            print("[CHECK] Mol2ChemFig server started on http://localhost:5001")
            return proc
        else:
            print("[ERROR] Mol2ChemFig server not found")
            return None
    except Exception as e:
        print(f"[ERROR] Error starting Mol2ChemFig server: {e}")
        return None

def start_pubchem_server():
    """Start the PubChem server"""
    try:
        pubchem_server_path = project_root / "pubchem_server.py"
        if pubchem_server_path.exists():
            print("[WEB] Starting PubChem server on port 5002...")

            proc = subprocess.Popen([
                sys.executable, str(pubchem_server_path)
            ], stdout=subprocess.PIPE, stderr=subprocess.PIPE)

            SERVER_PROCESSES.append(proc)
            print("[CHECK] PubChem server started on http://localhost:5002")
            return proc
        else:
            print("[ERROR] PubChem server not found")
            return None
    except Exception as e:
        print(f"[ERROR] Error starting PubChem server: {e}")
        return None

def check_server_health(url, name, timeout=5):
    """Check if a server is responding"""
    try:
        response = requests.get(url, timeout=timeout)
        return response.status_code == 200
    except:
        return False

def start_all_servers():
    """Start all servers and monitor them"""
    print("[ROCKET] Starting ChemParser Suite...")
    print("="*60)

    # Start servers
    molview_proc = start_molview_server()
    mol2chemfig_proc = start_mol2chemfig_server()
    pubchem_proc = start_pubchem_server()

    print("\n[WAIT] Waiting for servers to initialize...")
    time.sleep(5)  # Give servers time to start

    # Check server health
    print("\n[CHECKLIST] Server Status:")
    checks = [
        ("MolView", "http://localhost:5000/", molview_proc),
        ("Mol2ChemFig", "http://localhost:5001/", mol2chemfig_proc),
        ("PubChem", "http://localhost:5002/", pubchem_proc)
    ]

    for name, url, proc in checks:
        if proc is not None:
            status = "[OK] Online" if check_server_health(url, name) else "[ERROR] Offline"
            print(f"  {name:12}: {status} ({url})")
        else:
            print(f"  {name:12}: [WARN] Not started (missing)")

    print("="*60)
    print("[LIGHTBULB] Access the unified interface at: http://localhost:5000/")
    print("[LIGHTBULB] Or individual servers at their respective ports")
    print("[LIGHTBULB] Press Ctrl+C to stop all servers")
    print("="*60)

    # Optionally open browser
    try:
        # Wait a bit more for full initialization
        time.sleep(3)
        webbrowser.open("http://localhost:5000/")
        print("[GLOBE] Browser opened to unified interface")
    except:
        print("[GLOBE] Could not open browser automatically")

    # Monitor servers
    try:
        # Keep main thread alive
        while True:
            time.sleep(1)
            # Check if any server died
            for i, proc in enumerate(SERVER_PROCESSES[:]):
                if proc.poll() is not None:
                    print(f"\n[WARNING] Server process {i} has stopped unexpectedly")
                    SERVER_PROCESSES.remove(proc)
                    if not SERVER_PROCESSES:
                        print("[ERROR] All servers have stopped")
                        break
    except KeyboardInterrupt:
        print("\n[STOP] Shutting down servers...")
        stop_all_servers()

def stop_all_servers():
    """Stop all running server processes"""
    for proc in SERVER_PROCESSES:
        try:
            proc.terminate()
            proc.wait(timeout=5)
        except:
            try:
                proc.kill()
            except:
                pass
    print("[CHECK] All servers stopped")

def signal_handler(sig, frame):
    """Handle Ctrl+C gracefully"""
    print('\n\n[HANDSHAKE] Gracefully shutting down...')
    stop_all_servers()
    sys.exit(0)

if __name__ == "__main__":
    # Register signal handler for graceful shutdown
    signal.signal(signal.SIGINT, signal_handler)
    signal.signal(signal.SIGTERM, signal_handler)
    
    # Ensure cleanup on exit
    atexit.register(stop_all_servers)
    
    start_all_servers()