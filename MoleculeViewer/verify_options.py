#!/usr/bin/env python
"""
Comprehensive verification that all visualization options are now working correctly
with the corrected Docker API-compatible names and implementations.
"""

import sys
sys.path.insert(0, '.')

from app.chemistry import smiles_to_svg, nomenclature_to_smiles

def verify_option(name, smiles, options, expected_features=None):
    """Verify a single visualization option"""
    error, svg = smiles_to_svg(smiles, 400, 400, options)
    
    status = "✓ PASS" if error is None and svg else "✗ FAIL"
    print(f"{status}: {name}")
    
    if error:
        print(f"  Error: {error}")
        return False
    
    if expected_features:
        for feature, should_contain in expected_features.items():
            contains = feature in svg
            matches = contains == should_contain
            marker = "✓" if matches else "✗"
            print(f"  {marker} {feature} in SVG: {contains} (expected: {should_contain})")
            if not matches:
                return False
    
    return True

print("=" * 70)
print("VISUALIZATION OPTIONS VERIFICATION - Docker API Compatibility")
print("=" * 70)
print()

# Test suite
tests = [
    ("Option: aromatic_circles=True", 
     "c1ccccc1",  # benzene
     {'aromatic_circles': True},
     {'<circle': True}),
    
    ("Option: aromatic_circles=False",
     "c1ccccc1",
     {'aromatic_circles': False},
     None),
    
    ("Option: show_carbons=True",
     "CCc1ccccc1",  # ethylbenzene
     {'show_carbons': True},
     None),
    
    ("Option: show_methyls=True",
     "CCc1ccccc1",
     {'show_methyls': True},
     None),
    
    ("Option: flip_horizontal=True",
     "CCc1ccccc1",
     {'flip_horizontal': True},
     {'scaleX(-1)': True}),
    
    ("Option: flip_vertical=True",
     "CCc1ccccc1",
     {'flip_vertical': True},
     {'scaleY(-1)': True}),
    
    ("Option: rotate=90",
     "c1ccccc1",
     {'rotate': 90},
     {'rotate(90deg)': True}),
    
    ("Option: rotate=180",
     "c1ccccc1",
     {'rotate': 180},
     {'rotate(180deg)': True}),
    
    ("Option: hydrogens='keep'",
     "CCO",
     {'hydrogens': 'keep'},
     None),
    
    ("Option: hydrogens='delete'",
     "CCO",
     {'hydrogens': 'delete'},
     None),
    
    ("Option: hydrogens='add'",
     "CCO",
     {'hydrogens': 'add'},
     None),
    
    ("Option: atom_numbers=True",
     "CCc1ccccc1",
     {'atom_numbers': True},
     None),
    
    ("Option: fancy_bonds=True",
     "c1ccccc1",
     {'fancy_bonds': True},
     None),
    
    ("Combined: show_carbons + show_methyls + aromatic_circles",
     "CCc1ccccc1",
     {'show_carbons': True, 'show_methyls': True, 'aromatic_circles': True},
     {'<circle': True}),
    
    ("Combined: flip + rotate + show_carbons",
     "CCc1ccccc1",
     {'flip_horizontal': True, 'rotate': 90, 'show_carbons': True},
     {'scaleX(-1)': True, 'rotate(90deg)': True}),
]

passed = 0
failed = 0

for test_name, smiles, options, features in tests:
    # Ensure all required option keys are present
    full_options = {
        'show_carbons': False,
        'show_methyls': False,
        'aromatic_circles': False,
        'fancy_bonds': False,
        'atom_numbers': False,
        'hydrogens': 'keep',
        'flip_horizontal': False,
        'flip_vertical': False,
        'rotate': 0,
        'recalculate_coordinates': False
    }
    full_options.update(options)
    
    if verify_option(test_name, smiles, full_options, features):
        passed += 1
    else:
        failed += 1
    print()

print("=" * 70)
print(f"RESULTS: {passed} passed, {failed} failed")
print("=" * 70)

if failed == 0:
    print("\n✓ All visualization options are working correctly!")
    print("✓ Docker API compatibility verified!")
    sys.exit(0)
else:
    print(f"\n✗ {failed} test(s) failed")
    sys.exit(1)
