#!/usr/bin/env python3
"""
Comprehensive Test Runner for Chemparser Project
Tests all components: MoleculeViewer, mol2chemfig, and extension
"""

import requests
import subprocess
import time
import sys
import os
from pathlib import Path

# Configuration
PROJECT_ROOT = Path(__file__).parent
MOLECULE_VIEWER_PORT = 5000
MOL2CHEMFIG_PORT = 8000
MOL2CHEMFIG_WRAPPER_PORT = 5001

class Colors:
    GREEN = '\033[92m'
    RED = '\033[91m'
    YELLOW = '\033[93m'
    BLUE = '\033[94m'
    END = '\033[0m'
    BOLD = '\033[1m'

def print_header(text):
    print(f"\n{Colors.BOLD}{Colors.BLUE}{'='*70}{Colors.END}")
    print(f"{Colors.BOLD}{Colors.BLUE}{text}{Colors.END}")
    print(f"{Colors.BOLD}{Colors.BLUE}{'='*70}{Colors.END}\n")

def print_success(text):
    print(f"{Colors.GREEN}[OK] {text}{Colors.END}")

def print_error(text):
    print(f"{Colors.RED}[FAIL] {text}{Colors.END}")

def print_warning(text):
    print(f"{Colors.YELLOW}[WARN] {text}{Colors.END}")

def print_info(text):
    print(f"{Colors.BLUE}[INFO] {text}{Colors.END}")

def test_server(url, name, timeout=2):
    """Test if a server is running"""
    try:
        response = requests.get(url, timeout=timeout)
        if response.status_code in [200, 404]:  # 404 is ok for root
            print_success(f"{name} is running at {url}")
            return True
        else:
            print_warning(f"{name} returned status {response.status_code}")
            return False
    except requests.exceptions.ConnectionError:
        print_error(f"{name} is not running at {url}")
        return False
    except Exception as e:
        print_error(f"{name} error: {str(e)}")
        return False

def test_moleculeviewer():
    """Test MoleculeViewer server"""
    print_header("Testing MoleculeViewer Server (Port 5000)")

    base_url = f"http://localhost:{MOLECULE_VIEWER_PORT}"

    # Test health endpoint
    if not test_server(f"{base_url}/health", "MoleculeViewer health"):
        return False

    # Test SMILES rendering
    test_cases = [
        ("CCO", "Ethanol"),
        ("c1ccccc1", "Benzene"),
        ("CC(=O)C", "Acetone")
    ]

    passed = 0
    for smiles, name in test_cases:
        try:
            url = f"{base_url}/img/smiles?smiles={smiles}&width=300&height=200"
            response = requests.get(url, timeout=10)
            if response.status_code == 200 and 'svg' in response.headers.get('Content-Type', ''):
                print_success(f"{name} ({smiles}): SVG generated ({len(response.content)} bytes)")
                passed += 1
            else:
                print_error(f"{name} ({smiles}): Failed")
        except Exception as e:
            print_error(f"{name} ({smiles}): {str(e)}")

    print(f"\n{passed}/{len(test_cases)} tests passed\n")
    return passed == len(test_cases)

def test_mol2chemfig():
    """Test mol2chemfig Docker backend"""
    print_header("Testing mol2chemfig Backend (Port 8000)")

    base_url = f"http://localhost:{MOL2CHEMFIG_PORT}"

    # Test if backend is running
    if not test_server(base_url, "mol2chemfig backend"):
        print_info("mol2chemfig backend might not be running")
        print_info("Start it with: docker-compose up -d")
        return False

    return True

def test_extension_files():
    """Test if extension files exist"""
    print_header("Testing Extension Files")

    extension_dir = PROJECT_ROOT / "chem-extension"
    required_files = [
        "manifest.json",
        "content.js",
        "popup.html",
        "popup.js"
    ]

    all_exist = True
    for file in required_files:
        filepath = extension_dir / file
        if filepath.exists():
            print_success(f"{file} exists")
        else:
            print_error(f"{file} missing")
            all_exist = False

    return all_exist

def start_moleculeviewer():
    """Start MoleculeViewer server"""
    print_header("Starting MoleculeViewer Server")

    try:
        # Check if already running
        if test_server(f"http://localhost:{MOLECULE_VIEWER_PORT}/health", "MoleculeViewer", timeout=1):
            print_info("MoleculeViewer already running")
            return True

        # Start the server
        print_info("Starting MoleculeViewer...")
        moleculeviewer_dir = PROJECT_ROOT / "MoleculeViewer"

        subprocess.Popen(
            ["node", "server.js"],
            cwd=moleculeviewer_dir,
            stdout=subprocess.DEVNULL,
            stderr=subprocess.DEVNULL
        )

        # Wait for server to start
        for i in range(10):
            time.sleep(1)
            if test_server(f"http://localhost:{MOLECULE_VIEWER_PORT}/health", "MoleculeViewer", timeout=1):
                print_success("MoleculeViewer started successfully")
                return True

        print_error("MoleculeViewer failed to start within 10 seconds")
        return False

    except Exception as e:
        print_error(f"Failed to start MoleculeViewer: {str(e)}")
        return False

def check_dependencies():
    """Check if required dependencies are installed"""
    print_header("Checking Dependencies")

    deps = {
        "node": ["node", "--version"],
        "npm": ["npm", "--version"],
        "python": [sys.executable, "--version"],
        "docker": ["docker", "--version"]
    }

    all_installed = True
    for name, cmd in deps.items():
        try:
            result = subprocess.run(cmd, capture_output=True, text=True, timeout=5)
            if result.returncode == 0:
                version = result.stdout.strip() or result.stderr.strip()
                print_success(f"{name}: {version}")
            else:
                print_error(f"{name}: Not working properly")
                all_installed = False
        except FileNotFoundError:
            print_error(f"{name}: Not installed")
            all_installed = False
        except Exception as e:
            print_warning(f"{name}: Could not check ({str(e)})")

    return all_installed

def main():
    """Run all tests"""
    print_header("Chemparser Project Test Suite")

    results = {}

    # Check dependencies
    results['dependencies'] = check_dependencies()

    # Test extension files
    results['extension_files'] = test_extension_files()

    # Test or start MoleculeViewer
    if not test_server(f"http://localhost:{MOLECULE_VIEWER_PORT}/health", "MoleculeViewer", timeout=1):
        start_moleculeviewer()
        time.sleep(2)
    results['moleculeviewer'] = test_moleculeviewer()

    # Test mol2chemfig
    results['mol2chemfig'] = test_mol2chemfig()

    # Print summary
    print_header("Test Summary")

    for test_name, passed in results.items():
        if passed:
            print_success(f"{test_name}: PASSED")
        else:
            print_error(f"{test_name}: FAILED")

    total = len(results)
    passed = sum(1 for v in results.values() if v)

    print(f"\n{Colors.BOLD}Overall: {passed}/{total} test suites passed{Colors.END}\n")

    if passed == total:
        print_success("All tests passed! System is ready.")
        return 0
    else:
        print_warning("Some tests failed. Check the output above.")
        return 1

if __name__ == "__main__":
    sys.exit(main())
