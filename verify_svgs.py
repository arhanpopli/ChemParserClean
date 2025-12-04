"""
Inspect SVG files to verify they contain specific text elements.
"""
import os
import re

def check_svg_features(filename, expected_patterns):
    print(f"\n🔍 Inspecting {os.path.basename(filename)}...")
    
    if not os.path.exists(filename):
        print("❌ File not found")
        return

    with open(filename, 'r', encoding='utf-8') as f:
        content = f.read()
    
    # Simple check for text elements
    # RDKit usually puts text in <text> tags or paths, but let's look for the strings first
    
    all_found = True
    for pattern in expected_patterns:
        if pattern in content:
            print(f"  ✅ Found '{pattern}'")
        else:
            print(f"  ❌ Missing '{pattern}'")
            all_found = False
            
    # Also look for <circle> tags for aromaticity
    if "aromatic" in filename or "circles" in filename:
        if "<circle" in content:
             print(f"  ✅ Found aromatic circle (<circle> tag)")
        else:
             print(f"  ⚠️ No <circle> tag found (might use dashed paths instead)")

    return all_found

# Define checks
checks = [
    ("rdkit_test_output/toluene_show_methyls.svg", ["CH", "3"]), # RDKit might split CH3
    ("rdkit_test_output/ethanol_show_carbons_hydrogens.svg", ["H", "C", "O"]),
    ("rdkit_test_output/benzene_aromatic_circles.svg", []), # Just check for circle
]

print("Verifying RDKit SVG Features:")
for filename, patterns in checks:
    check_svg_features(filename, patterns)
