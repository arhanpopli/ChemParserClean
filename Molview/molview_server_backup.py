#!/usr/bin/env python3
"""
MolView Local Server
This script sets up a local PHP server to run MolView locally
"""

import os
import sys
import subprocess
import threading
import time
import socket
from pathlib import Path
import webbrowser
import argparse

def check_port(port):
    """Check if a port is available"""
    with socket.socket(socket.AF_INET, socket.SOCK_STREAM) as s:
        return s.connect_ex(('localhost', port)) != 0

def find_available_port(start_port=8080):
    """Find an available port starting from start_port"""
    port = start_port
    while not check_port(port):
        port += 1
    return port

def start_php_server(directory, port):
    """Start PHP built-in server"""
    try:
        # Change to the directory
        os.chdir(directory)
        
        # Start PHP server
        cmd = ['php', '-S', f'localhost:{port}']
        process = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        
        return process
    except FileNotFoundError:
        print("Error: PHP is not installed or not in PATH")
        print("Please install PHP to run MolView locally")
        return None

def setup_molview_configs(molview_dir):
    """Set up any required configurations for MolView"""
    # Check if composer.json exists and install dependencies if needed
    composer_json = os.path.join(molview_dir, 'composer.json')
    if os.path.exists(composer_json):
        print("Installing PHP dependencies...")
        try:
            subprocess.run(['composer', 'install'], cwd=molview_dir, check=True)
        except (subprocess.CalledProcessError, FileNotFoundError):
            print("Warning: Composer not found or install failed. Some features might not work.")

    # Set up writable directories if needed
    cache_dirs = ['cache', 'tmp', 'log']
    for cache_dir in cache_dirs:
        full_path = os.path.join(molview_dir, cache_dir)
        if not os.path.exists(full_path):
            os.makedirs(full_path, exist_ok=True)
        # Make sure it's writable
        os.chmod(full_path, 0o755)

def main():
    parser = argparse.ArgumentParser(description='MolView Local Server')
    parser.add_argument('--port', type=int, default=8080, help='Port to run the server on (default: 8080)')
    parser.add_argument('--molview-dir', type=str, help='Path to MolView directory')
    parser.add_argument('--no-browser', action='store_true', help='Don\'t open browser automatically')
    parser.add_argument('--build', action='store_true', help='Run build process before starting')
    
    args = parser.parse_args()
    
    # Find MolView directory
    molview_dir = args.molview_dir
    
    if not molview_dir:
        # Try to find MolView directory automatically
        current_dir = Path.cwd()
        
        # Look in common locations
        possible_paths = [
            Path('./Molview/molview-2.4.6 - github repo, src'),
            Path('../Molview/molview-2.4.6 - github repo, src'),
            Path('./molview-2.4.6 - github repo, src'),
            Path('../molview-2.4.6 - github repo, src')
        ]
        
        for path in possible_paths:
            if path.exists():
                molview_dir = str(path.absolute())
                break
    
    if not molview_dir or not os.path.exists(molview_dir):
        print("Error: Could not find MolView directory")
        print("Please specify --molview-dir or run this script from the correct location")
        print("\nPossible locations:")
        for path in possible_paths:
            print(f"  - {path}")
        sys.exit(1)
    
    # Check if this is a valid MolView directory
    required_files = ['index.php', 'php/', 'src/']
    missing_files = []
    for req_file in required_files:
        if not os.path.exists(os.path.join(molview_dir, req_file)):
            missing_files.append(req_file)
    
    if missing_files:
        print(f"Error: Missing required MolView files: {missing_files}")
        sys.exit(1)
    
    print(f"MolView directory found: {molview_dir}")
    
    # Set up configurations
    setup_molview_configs(molview_dir)
    
    # If build flag is set, try to build MolView
    if args.build:
        print("Building MolView...")
        build_scripts = ['build.sh', 'Gruntfile.js', 'gulpfile.js', 'package.json']
        
        if os.path.exists(os.path.join(molview_dir, 'package.json')):
            try:
                # Install npm packages
                subprocess.run(['npm', 'install'], cwd=molview_dir, check=True)
                print("Installed npm packages")
            except subprocess.CalledProcessError:
                print("Warning: npm install failed")
        
        if os.path.exists(os.path.join(molview_dir, 'build.sh')):
            try:
                # Run build script
                subprocess.run(['bash', './build.sh'], cwd=molview_dir, check=True)
                print("Built MolView successfully")
            except (subprocess.CalledProcessError, FileNotFoundError):
                print("Warning: Build script failed or bash not available")
    
    # Find available port
    port = find_available_port(args.port)
    if port != args.port:
        print(f"Port {args.port} is busy, using {port} instead")
    
    # Start PHP server
    print(f"Starting MolView server on http://localhost:{port}")
    process = start_php_server(molview_dir, port)
    
    if not process:
        sys.exit(1)
    
    print(f"MolView server started successfully!")
    print(f"Access MolView at: http://localhost:{port}")
    print("Press Ctrl+C to stop the server")
    
    # Open browser if requested
    if not args.no_browser:
        time.sleep(2)  # Wait a bit for server to start
        webbrowser.open(f'http://localhost:{port}')
    
    try:
        # Wait for process to finish (Ctrl+C will interrupt)
        process.wait()
    except KeyboardInterrupt:
        print("\nShutting down server...")
        process.terminate()
        try:
            process.wait(timeout=5)
        except subprocess.TimeoutExpired:
            process.kill()
        print("Server stopped")

if __name__ == "__main__":
    main()