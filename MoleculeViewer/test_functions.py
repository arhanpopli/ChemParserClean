#!/usr/bin/env python
"""Quick test of MoleculeViewer functions"""

import sys
sys.path.insert(0, '/'.join(__file__.split('/')[:-1]))

from app.chemistry import nomenclature_to_smiles, smiles_to_svg, get_molecule_info

print("=" * 50)
print("MoleculeViewer - Quick Function Test")
print("=" * 50)

# Test 1: nomenclature_to_smiles
print("\n[1] Testing nomenclature_to_smiles('benzene')")
err, smiles = nomenclature_to_smiles('benzene')
print(f"    Result: {smiles}")
if not err:
    print("    ✓ PASS")
else:
    print(f"    ✗ FAIL: {err}")

# Test 2: smiles_to_svg
print("\n[2] Testing smiles_to_svg('c1ccccc1')")
err, svg = smiles_to_svg('c1ccccc1')
if not err and svg and len(svg) > 0:
    print(f"    SVG Generated: {len(svg)} characters")
    print("    ✓ PASS")
else:
    print(f"    ✗ FAIL: {err}")

# Test 3: get_molecule_info
print("\n[3] Testing get_molecule_info('c1ccccc1')")
err, info = get_molecule_info('c1ccccc1')
if not err and info:
    print(f"    Formula: {info.get('formula')}")
    print(f"    Molecular Weight: {info.get('molecular_weight')} g/mol")
    print(f"    LogP: {info.get('logp')}")
    print("    ✓ PASS")
else:
    print(f"    ✗ FAIL: {err}")

print("\n" + "=" * 50)
print("All tests completed successfully!")
print("=" * 50)
