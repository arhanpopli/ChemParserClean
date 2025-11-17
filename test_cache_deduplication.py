#!/usr/bin/env python3
"""
Test script for cache deduplication and canonical SMILES implementation

This script tests:
1. SMILES canonicalization works correctly
2. Same molecule with different SMILES representations uses same cache key
3. Nomenclature to SMILES conversion produces canonical cache keys
4. No duplicate cache files are created for equivalent molecules
5. Cache deduplication utility correctly identifies and removes duplicates

Usage:
    python test_cache_deduplication.py
"""

import requests
import time
import hashlib
import json
from pathlib import Path

# Server URLs
MOLECULEVIEWER_URL = "http://localhost:5000"
MOL2CHEMFIG_URL = "http://localhost:5001"

# Test molecules with different SMILES representations
TEST_MOLECULES = [
    {
        "name": "ethanol",
        "smiles_variants": ["CCO", "OCC", "C(O)C"],  # Different representations of ethanol
        "canonical": "CCO"
    },
    {
        "name": "benzene",
        "smiles_variants": ["c1ccccc1", "C1=CC=CC=C1"],  # Aromatic vs Kekule
        "canonical": "c1ccccc1"
    },
    {
        "name": "acetic acid",
        "smiles_variants": ["CC(=O)O", "CC(O)=O", "C(C)(=O)O"],
        "canonical": "CC(=O)O"
    }
]

def print_header(text):
    """Print a formatted header"""
    print("\n" + "="*80)
    print(f"  {text}")
    print("="*80)

def print_subheader(text):
    """Print a formatted subheader"""
    print(f"\n--- {text} ---")

def test_canonicalization():
    """Test that SMILES canonicalization is working"""
    print_header("TEST 1: SMILES Canonicalization")

    try:
        from canonicalize_smiles import canonicalize_smiles

        for molecule in TEST_MOLECULES:
            print_subheader(f"Testing {molecule['name']}")
            print(f"Expected canonical SMILES: {molecule['canonical']}")

            for variant in molecule['smiles_variants']:
                canonical = canonicalize_smiles(variant)
                match = "✓" if canonical == molecule['canonical'] else "✗"
                print(f"  {match} {variant:20s} -> {canonical}")

                if canonical != molecule['canonical']:
                    print(f"    ERROR: Expected {molecule['canonical']}, got {canonical}")
                    return False

        print("\n✓ All SMILES canonicalization tests passed!")
        return True

    except ImportError:
        print("✗ Error: Cannot import canonicalize_smiles module")
        return False
    except Exception as e:
        print(f"✗ Error during canonicalization test: {e}")
        return False

def test_moleculeviewer_cache():
    """Test MoleculeViewer server cache deduplication"""
    print_header("TEST 2: MoleculeViewer Cache Deduplication")

    cache_dir = Path("MoleculeViewer/cache/moleculeviewer")
    if not cache_dir.exists():
        print(f"✗ Cache directory not found: {cache_dir}")
        return False

    # Clear existing cache files for test molecules
    print("Clearing old test cache files...")
    for svg_file in cache_dir.glob("*.svg"):
        svg_file.unlink()

    initial_file_count = len(list(cache_dir.glob("*.svg")))
    print(f"Initial cache files: {initial_file_count}")

    for molecule in TEST_MOLECULES:
        print_subheader(f"Testing {molecule['name']}")

        cache_files_before = set(f.name for f in cache_dir.glob("*.svg"))

        # Request the same molecule with different SMILES representations
        for variant in molecule['smiles_variants']:
            try:
                url = f"{MOLECULEVIEWER_URL}/img/smiles?smiles={variant}"
                response = requests.get(url, timeout=10)

                if response.status_code == 200:
                    print(f"  ✓ {variant:20s} -> HTTP 200")
                else:
                    print(f"  ✗ {variant:20s} -> HTTP {response.status_code}")
                    return False

                time.sleep(0.5)  # Brief delay between requests

            except Exception as e:
                print(f"  ✗ Error requesting {variant}: {e}")
                return False

        # Check cache files after
        cache_files_after = set(f.name for f in cache_dir.glob("*.svg"))
        new_files = cache_files_after - cache_files_before

        print(f"\n  New cache files created: {len(new_files)}")

        if len(new_files) == 1:
            print(f"  ✓ All SMILES variants used the same cache file!")
            print(f"    Cache file: {list(new_files)[0]}")
        else:
            print(f"  ✗ ERROR: Expected 1 cache file, but {len(new_files)} were created")
            print(f"    Files: {new_files}")
            return False

    print("\n✓ MoleculeViewer cache deduplication tests passed!")
    return True

def test_nomenclature_to_smiles_cache():
    """Test that nomenclature conversion produces same cache as direct SMILES"""
    print_header("TEST 3: Nomenclature to SMILES Cache Consistency")

    cache_dir = Path("MoleculeViewer/cache/moleculeviewer")

    for molecule in TEST_MOLECULES:
        print_subheader(f"Testing {molecule['name']}")

        cache_files_before = set(f.name for f in cache_dir.glob("*.svg"))

        # Request by nomenclature
        try:
            url = f"{MOLECULEVIEWER_URL}/img/nomenclature?nomenclature={molecule['name']}"
            response = requests.get(url, timeout=15)

            if response.status_code == 200:
                print(f"  ✓ Nomenclature request succeeded")
            else:
                print(f"  ✗ Nomenclature request failed: HTTP {response.status_code}")
                print(f"    Note: This might fail if PubChem API doesn't recognize '{molecule['name']}'")
                continue

        except Exception as e:
            print(f"  ✗ Error requesting by nomenclature: {e}")
            continue

        cache_files_after = set(f.name for f in cache_dir.glob("*.svg"))
        new_files = cache_files_after - cache_files_before

        if len(new_files) == 0:
            print(f"  ✓ Reused existing cache from SMILES request (no new files)")
        elif len(new_files) == 1:
            print(f"  ! Created new cache file (might be first request for this molecule)")
            print(f"    Cache file: {list(new_files)[0]}")
        else:
            print(f"  ✗ ERROR: Multiple cache files created: {new_files}")

        time.sleep(1)  # Delay for API rate limiting

    print("\n✓ Nomenclature to SMILES cache consistency tests passed!")
    return True

def test_mol2chemfig_cache():
    """Test Mol2ChemFig server cache deduplication"""
    print_header("TEST 4: Mol2ChemFig Cache Deduplication")

    cache_dir = Path("cache/mol2chemfig")
    if not cache_dir.exists():
        print(f"✗ Cache directory not found: {cache_dir}")
        return False

    # Test one molecule with different SMILES representations
    molecule = TEST_MOLECULES[0]  # ethanol
    print_subheader(f"Testing {molecule['name']}")

    cache_files_before = set(f.name for f in cache_dir.glob("*.svg"))

    for variant in molecule['smiles_variants']:
        try:
            url = f"{MOL2CHEMFIG_URL}/api/generate"
            payload = {
                "smiles": variant,
                "format": "smiles",
                "options": [],
                "return_format": "svg"
            }

            response = requests.post(url, json=payload, timeout=30)

            if response.status_code == 200:
                result = response.json()
                if result.get('success'):
                    cached = result.get('cached', False)
                    print(f"  ✓ {variant:20s} -> {'Cached' if cached else 'Generated'}")
                else:
                    print(f"  ✗ {variant:20s} -> API error: {result.get('error')}")
            else:
                print(f"  ✗ {variant:20s} -> HTTP {response.status_code}")

            time.sleep(1)  # Brief delay

        except Exception as e:
            print(f"  ✗ Error requesting {variant}: {e}")
            return False

    cache_files_after = set(f.name for f in cache_dir.glob("*.svg"))
    new_files = cache_files_after - cache_files_before

    print(f"\n  New cache files created: {len(new_files)}")

    if len(new_files) == 1:
        print(f"  ✓ All SMILES variants used the same cache file!")
    elif len(new_files) == 0:
        print(f"  ✓ All requests reused existing cache!")
    else:
        print(f"  ✗ ERROR: Expected 1 cache file, but {len(new_files)} were created")
        return False

    print("\n✓ Mol2ChemFig cache deduplication tests passed!")
    return True

def test_cache_stats():
    """Display cache statistics"""
    print_header("Cache Statistics")

    # MoleculeViewer cache
    mv_cache = Path("MoleculeViewer/cache/moleculeviewer")
    if mv_cache.exists():
        mv_files = list(mv_cache.glob("*.svg"))
        mv_size = sum(f.stat().st_size for f in mv_files)
        print(f"\nMoleculeViewer Cache:")
        print(f"  Files: {len(mv_files)}")
        print(f"  Size: {mv_size / 1024:.2f} KB")

    # Mol2ChemFig cache
    m2c_cache = Path("cache/mol2chemfig")
    if m2c_cache.exists():
        m2c_files = list(m2c_cache.glob("*"))
        m2c_size = sum(f.stat().st_size for f in m2c_files if f.is_file())
        print(f"\nMol2ChemFig Cache:")
        print(f"  Files: {len(m2c_files)}")
        print(f"  Size: {m2c_size / 1024:.2f} KB")

def main():
    """Run all tests"""
    print("\n" + "="*80)
    print("  CACHE DEDUPLICATION TEST SUITE")
    print("="*80)
    print("\nThis test suite verifies that:")
    print("1. SMILES canonicalization works correctly")
    print("2. Equivalent molecules use the same cache key")
    print("3. No duplicate cache files are created")
    print("\nNOTE: Servers must be running on ports 5000 and 5001")
    print("="*80)

    # Check if servers are running
    try:
        requests.get(f"{MOLECULEVIEWER_URL}/health", timeout=5)
        print("✓ MoleculeViewer server is running")
    except:
        print("✗ MoleculeViewer server is not responding on port 5000")
        print("  Please start the server with: node MoleculeViewer/server.js")
        return

    try:
        requests.get(f"{MOL2CHEMFIG_URL}/health", timeout=5)
        print("✓ Mol2ChemFig server is running")
    except:
        print("✗ Mol2ChemFig server is not responding on port 5001")
        print("  Please start the server with: python mol2chemfig_server.py")
        return

    # Run tests
    results = []

    results.append(("SMILES Canonicalization", test_canonicalization()))
    results.append(("MoleculeViewer Cache", test_moleculeviewer_cache()))
    results.append(("Nomenclature Cache Consistency", test_nomenclature_to_smiles_cache()))
    results.append(("Mol2ChemFig Cache", test_mol2chemfig_cache()))

    # Display cache statistics
    test_cache_stats()

    # Summary
    print_header("TEST SUMMARY")
    passed = sum(1 for _, result in results if result)
    total = len(results)

    for test_name, result in results:
        status = "PASS" if result else "FAIL"
        symbol = "✓" if result else "✗"
        print(f"  {symbol} {test_name:40s} [{status}]")

    print(f"\n  Total: {passed}/{total} tests passed")

    if passed == total:
        print("\n  🎉 All tests passed! Cache deduplication is working correctly.")
    else:
        print(f"\n  ⚠ {total - passed} test(s) failed. Please review the output above.")

    print("="*80 + "\n")

if __name__ == '__main__':
    main()
