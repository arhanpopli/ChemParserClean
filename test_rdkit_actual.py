"""
Test RDKit rendering capabilities with various options.
Generate actual SVG files to see what RDKit can produce.
"""
import sys
import os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), 'MoleculeViewer'))

from app.chemistry import smiles_to_svg

# Create output directory
OUTPUT_DIR = "rdkit_test_output"
os.makedirs(OUTPUT_DIR, exist_ok=True)

test_cases = [
    ("benzene_default", "c1ccccc1", {}),
    ("benzene_aromatic_circles", "c1ccccc1", {"aromatic_circles": True}),
    ("benzene_show_carbons", "c1ccccc1", {"show_carbons": True}),
    ("benzene_both", "c1ccccc1", {"aromatic_circles": True, "show_carbons": True}),
    ("toluene_show_methyls", "Cc1ccccc1", {"show_methyls": True, "aromatic_circles": True}),
    ("ethanol_show_carbons_hydrogens", "CCO", {"show_carbons": True, "hydrogens": "add"}),
]

print("🧪 Testing RDKit rendering capabilities...\n")

for name, smiles, options in test_cases:
    print(f"Testing: {name}")
    print(f"  SMILES: {smiles}")
    print(f"  Options: {options}")
    
    error, svg = smiles_to_svg(smiles, 400, 300, options)
    
    if error:
        print(f"  ❌ Error: {error}\n")
    else:
        filename = os.path.join(OUTPUT_DIR, f"{name}.svg")
        with open(filename, "w", encoding="utf-8") as f:
            f.write(svg)
        print(f"  ✅ Saved: {filename}\n")

print(f"\n✨ Done! Check {OUTPUT_DIR}/ for SVG files")
