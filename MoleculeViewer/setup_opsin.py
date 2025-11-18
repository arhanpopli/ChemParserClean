#!/usr/bin/env python
"""
Setup OPSIN parser for MoleculeViewer (Linux/Mac/Windows compatible)

OPSIN is a Java-based IUPAC nomenclature to SMILES converter.
This script downloads and configures it for MoleculeViewer.
"""

import subprocess
import os
import sys
import urllib.request
from pathlib import Path

def check_java():
    """Check if Java is installed"""
    try:
        result = subprocess.run(['java', '-version'], capture_output=True, text=True)
        if result.returncode == 0:
            print("✓ Java found")
            # Print version
            version_line = (result.stderr or result.stdout).split('\n')[0]
            print(f"  {version_line}")
            return True
    except FileNotFoundError:
        pass
    return False

def opsin_exists():
    """Check if OPSIN JAR already exists"""
    return Path('opsin-cli.jar').exists()

def download_opsin():
    """Download OPSIN JAR file"""
    print("\nDownloading OPSIN parser...")
    
    urls = [
        'https://github.com/dan2097/opsin/releases/download/v2.8.1/opsin-cli.jar',
        'https://github.com/opsin/opsin/releases/download/v2.8.0/opsin-cli-2.8.0.jar',
    ]
    
    for url in urls:
        try:
            print(f"  Trying: {url}")
            urllib.request.urlretrieve(url, 'opsin-cli.jar', reporthook=_download_hook)
            print("\n✓ Downloaded successfully!")
            return True
        except Exception as e:
            print(f"\n  ✗ Failed: {e}")
            continue
    
    return False

def _download_hook(count, block_size, total_size):
    """Progress indicator for download"""
    if total_size > 0:
        percent = min((count * block_size * 100) / total_size, 100)
        sys.stdout.write(f"\r  Progress: {percent:.1f}%")
        sys.stdout.flush()

def test_opsin():
    """Test if OPSIN works"""
    try:
        result = subprocess.run(
            ['java', '-jar', 'opsin-cli.jar', 'ethanol'],
            capture_output=True,
            text=True,
            timeout=5
        )
        return result.returncode == 0
    except Exception:
        return False

def main():
    print("=" * 60)
    print("OPSIN Parser Setup for MoleculeViewer")
    print("=" * 60)
    
    # Check Java
    print("\nChecking Java installation...")
    if not check_java():
        print("\n✗ ERROR: Java not found!")
        print("\nOPSIN requires Java Runtime Environment (JRE)")
        print("Install from: https://www.java.com/en/download/")
        print("\nAfter installing Java, run this script again.")
        return 1
    
    # Check if already installed
    print("\nChecking for existing OPSIN installation...")
    if opsin_exists():
        print("✓ OPSIN already installed (opsin-cli.jar exists)")
        if test_opsin():
            print("✓ OPSIN is working correctly")
        print("\nSetup complete!")
        return 0
    
    # Download OPSIN
    print("\nDownloading OPSIN...")
    if not download_opsin():
        print("\n✗ Download failed!")
        print("\nManual setup:")
        print("1. Visit: https://github.com/dan2097/opsin/releases")
        print("2. Download: opsin-cli-*.jar")
        print("3. Copy to this directory")
        print("4. Rename to: opsin-cli.jar")
        return 1
    
    # Test OPSIN
    print("\nTesting OPSIN...")
    if test_opsin():
        print("✓ OPSIN test successful")
    else:
        print("⚠ Warning: OPSIN test failed (may still work)")
    
    # Success
    print("\n" + "=" * 60)
    print("OPSIN Setup Complete!")
    print("=" * 60)
    print("\nYou can now use IUPAC nomenclature parsing.")
    print("\nConfigure in app/config.py:")
    print("  NOMENCLATURE_PARSER = 'auto'  (or 'opsin' for IUPAC only)")
    print("  ENABLE_OPSIN = True")
    
    return 0

if __name__ == '__main__':
    sys.exit(main())
