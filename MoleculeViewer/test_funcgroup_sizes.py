"""
Test to see what font size RDKit uses for functional groups
"""
import sys
sys.path.insert(0, '.')

from app.chemistry import smiles_to_svg

# Test molecules with functional groups
test_cases = [
    ('CC(=O)O', 'Acetic acid (COOH)'),
    ('CCO', 'Ethanol (OH)'),
    ('c1ccc(cc1)N(=O)=O', 'Nitrobenzene (NO2)'),
    ('CC(C)O', 'Isopropanol'),
]

print("Testing functional group font sizes:")
print("=" * 70)

for smiles, name in test_cases:
    print(f"\n{name} ({smiles}):")
    
    error, svg = smiles_to_svg(smiles, 400, 300, {})
    
    if error:
        print(f"  ❌ Error: {error}")
    else:
        # Save SVG
        filename = f"test_funcgroup_{name.split()[0].lower()}.svg"
        with open(filename, 'w', encoding='utf-8') as f:
            f.write(svg)
        print(f"  ✓ Saved: {filename}")

print("\n" + "=" * 70)
print("✓ Check the SVG files to see font sizes used by RDKit")
