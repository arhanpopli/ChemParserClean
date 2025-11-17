#!/usr/bin/env python3
"""
Test script to verify cache separation implementation
Tests both MoleculeViewer and mol2chemfig servers
"""

import os
import sys
import time
import requests
import json
from pathlib import Path

# ANSI color codes for terminal output
GREEN = '\033[92m'
RED = '\033[91m'
YELLOW = '\033[93m'
BLUE = '\033[94m'
RESET = '\033[0m'

def print_header(text):
    """Print a formatted header"""
    print(f"\n{BLUE}{'=' * 70}{RESET}")
    print(f"{BLUE}{text:^70}{RESET}")
    print(f"{BLUE}{'=' * 70}{RESET}\n")

def print_success(text):
    """Print success message"""
    print(f"{GREEN}✓ {text}{RESET}")

def print_error(text):
    """Print error message"""
    print(f"{RED}✗ {text}{RESET}")

def print_warning(text):
    """Print warning message"""
    print(f"{YELLOW}⚠ {text}{RESET}")

def print_info(text):
    """Print info message"""
    print(f"{BLUE}ℹ {text}{RESET}")

def check_server(url, name):
    """Check if a server is running"""
    try:
        response = requests.get(url, timeout=2)
        if response.status_code == 200:
            print_success(f"{name} is running on {url}")
            return True
        else:
            print_error(f"{name} returned status {response.status_code}")
            return False
    except requests.exceptions.RequestException as e:
        print_error(f"{name} is not running on {url}")
        print_info(f"Error: {str(e)}")
        return False

def check_cache_dir(path, name):
    """Check if cache directory exists"""
    if os.path.exists(path):
        print_success(f"{name} cache directory exists: {path}")
        return True
    else:
        print_warning(f"{name} cache directory does not exist yet: {path}")
        print_info("It will be created when the server generates its first file")
        return False

def count_cache_files(path):
    """Count files in cache directory"""
    if not os.path.exists(path):
        return 0
    return len([f for f in os.listdir(path) if os.path.isfile(os.path.join(path, f))])

def test_moleculeviewer():
    """Test MoleculeViewer server and cache"""
    print_header("MOLECULEVIEWER SERVER TEST")

    # Check server
    if not check_server("http://localhost:5000/health", "MoleculeViewer"):
        print_error("MoleculeViewer server is not running!")
        print_info("Start it with: cd MoleculeViewer && node server.js")
        return False

    # Check cache directory
    cache_path = Path("MoleculeViewer/cache/moleculeviewer")
    check_cache_dir(cache_path, "MoleculeViewer")

    # Test SVG generation
    print_info("\nTesting SVG generation...")
    test_smiles = "CCO"

    try:
        # Initial cache count
        initial_count = count_cache_files(cache_path)
        print_info(f"Initial cache files: {initial_count}")

        # Generate SVG
        response = requests.get(
            f"http://localhost:5000/img/smiles",
            params={"smiles": test_smiles},
            timeout=10
        )

        if response.status_code == 200:
            print_success(f"Generated SVG for {test_smiles}")

            # Wait a moment for file to be written
            time.sleep(0.5)

            # Check cache count
            final_count = count_cache_files(cache_path)
            print_info(f"Final cache files: {final_count}")

            if final_count > initial_count:
                print_success("Cache file was created successfully!")
            else:
                print_info("File may have been cached already (no new file created)")

            # Get cache info
            cache_info = requests.get("http://localhost:5000/cache-info").json()
            print_info(f"Total cached SVGs: {cache_info.get('cachedSvgs', 0)}")
            print_info(f"Cache size: {cache_info.get('totalCacheSize', '0 KB')}")

            return True
        else:
            print_error(f"Failed to generate SVG: Status {response.status_code}")
            return False

    except Exception as e:
        print_error(f"Error testing MoleculeViewer: {str(e)}")
        return False

def test_mol2chemfig():
    """Test mol2chemfig server and cache"""
    print_header("MOL2CHEMFIG SERVER TEST")

    # Check server
    if not check_server("http://localhost:5001/health", "mol2chemfig wrapper"):
        print_error("mol2chemfig wrapper server is not running!")
        print_info("Start it with: python mol2chemfig_server.py")
        return False

    # Check cache directory
    cache_path = Path("cache/mol2chemfig")
    check_cache_dir(cache_path, "mol2chemfig")

    # Test SVG generation
    print_info("\nTesting SVG generation...")
    test_smiles = "c1ccccc1"

    try:
        # Initial cache count
        initial_count = count_cache_files(cache_path)
        print_info(f"Initial cache files: {initial_count}")

        # Generate SVG
        response = requests.post(
            "http://localhost:5001/api/generate",
            json={"smiles": test_smiles, "return_format": "svg"},
            timeout=30
        )

        if response.status_code == 200:
            result = response.json()
            if result.get("success"):
                print_success(f"Generated SVG for {test_smiles}")
                print_info(f"SVG URL: {result.get('svg_url')}")
                print_info(f"Cached: {result.get('cached')}")

                # Wait a moment for file to be written
                time.sleep(0.5)

                # Check cache count
                final_count = count_cache_files(cache_path)
                print_info(f"Final cache files: {final_count}")

                if final_count > initial_count:
                    print_success("Cache file was created successfully!")
                else:
                    print_info("File was already cached (no new file created)")

                # Get cache stats
                cache_stats = requests.get("http://localhost:5001/api/cache/stats").json()
                print_info(f"Cached entries: {cache_stats.get('cached_entries', 0)}")
                print_info(f"Total files: {cache_stats.get('total_files', 0)}")
                print_info(f"Storage size: {cache_stats.get('storage_size_mb', 0)} MB")

                return True
            else:
                print_error(f"Generation failed: {result.get('error')}")
                return False
        else:
            print_error(f"Failed to generate SVG: Status {response.status_code}")
            return False

    except Exception as e:
        print_error(f"Error testing mol2chemfig: {str(e)}")
        return False

def verify_separation():
    """Verify that caches are properly separated"""
    print_header("CACHE SEPARATION VERIFICATION")

    mv_cache = Path("MoleculeViewer/cache/moleculeviewer")
    m2c_cache = Path("cache/mol2chemfig")

    issues = []

    # Check if both caches exist
    if not mv_cache.exists():
        print_warning("MoleculeViewer cache doesn't exist yet")
    else:
        mv_files = list(mv_cache.glob("*"))
        print_info(f"MoleculeViewer cache: {len(mv_files)} files")

        # Check for mol2chemfig-style files in MoleculeViewer cache
        for f in mv_files:
            if "_m2cf_" in f.name or "default_m2cf" in f.name:
                issues.append(f"mol2chemfig file found in MoleculeViewer cache: {f.name}")

    if not m2c_cache.exists():
        print_warning("mol2chemfig cache doesn't exist yet")
    else:
        m2c_files = list(m2c_cache.glob("*"))
        print_info(f"mol2chemfig cache: {len(m2c_files)} files")

        # Check for MoleculeViewer-style files in mol2chemfig cache
        for f in m2c_files:
            # MoleculeViewer uses MD5 hashes (32 chars) without descriptive names
            if len(f.stem) == 32 and not any(x in f.name for x in ["_", "-"]):
                issues.append(f"Potential MoleculeViewer file in mol2chemfig cache: {f.name}")

    if issues:
        print_error("Cache separation issues found:")
        for issue in issues:
            print_error(f"  - {issue}")
        return False
    else:
        print_success("Caches are properly separated!")
        return True

def check_old_caches():
    """Check for old cache directories that should be migrated/removed"""
    print_header("OLD CACHE DIRECTORY CHECK")

    old_caches = [
        ("MoleculeViewer/svg-cache", "MoleculeViewer old cache"),
        ("mol2chemfig_storage", "mol2chemfig old storage")
    ]

    found_old = False
    for path, name in old_caches:
        if os.path.exists(path):
            file_count = count_cache_files(path)
            if file_count > 0:
                found_old = True
                print_warning(f"{name} still exists with {file_count} files: {path}")
                print_info("Consider backing up and removing after verifying new cache works")
            else:
                print_info(f"{name} exists but is empty: {path}")
        else:
            print_success(f"{name} not found (good!): {path}")

    if not found_old:
        print_success("No old cache directories with files found!")

    return not found_old

def main():
    """Main test runner"""
    print_header("CACHE SEPARATION TEST SUITE")
    print_info("Testing separated cache folders for MoleculeViewer and mol2chemfig")
    print_info("Project: Chemparser")

    results = {
        "MoleculeViewer": False,
        "mol2chemfig": False,
        "Separation": False,
        "OldCaches": False
    }

    # Test MoleculeViewer
    results["MoleculeViewer"] = test_moleculeviewer()

    # Test mol2chemfig
    results["mol2chemfig"] = test_mol2chemfig()

    # Verify separation
    results["Separation"] = verify_separation()

    # Check old caches
    results["OldCaches"] = check_old_caches()

    # Summary
    print_header("TEST SUMMARY")

    for name, passed in results.items():
        if passed:
            print_success(f"{name}: PASSED")
        else:
            print_error(f"{name}: FAILED")

    all_passed = all(results.values())

    print()
    if all_passed:
        print_success("All tests passed! Cache separation is working correctly.")
        return 0
    else:
        print_error("Some tests failed. Please review the output above.")
        return 1

if __name__ == "__main__":
    sys.exit(main())
