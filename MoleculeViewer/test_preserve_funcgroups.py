"""
Test to ensure real functional groups (OH, NO2, etc.) are preserved
"""
import sys
sys.path.insert(0, '.')

from app.chemistry import smiles_to_svg

test_cases = [
    ('CCO', 'Ethanol (OH group)', {'show_methyls': True}),
    ('CC(=O)O', 'Acetic acid (COOH)', {'show_methyls': True}),
    ('c1ccc(cc1)[N+](=O)[O-]', 'Nitrobenzene (NO2)', {'show_methyls': True}),
]

print("Testing to ensure real functional groups are NOT replaced:")
print("=" * 70)

for smiles, name, opts in test_cases:
    print(f"\n{name} ({smiles}):")
    
    error, svg = smiles_to_svg(smiles, 400, 300, opts)
    
    if error:
        print(f"  ❌ Error: {error}")
    else:
        # Check what's in the SVG
        has_ch3 = 'CH₃' in svg
        has_oh = '>O<' in svg or 'fill=\'#FF0000\'' in svg  # O is red
        has_n = '>N<' in svg or '#0000FF' in svg  # N is blue
        
        print(f"  ✓ SVG generated")
        print(f"    Contains CH₃: {has_ch3}")
        print(f"    Contains O (red): {has_oh}")
        print(f"    Contains N (blue): {has_n}")
        
        # Save for inspection
        filename = f"test_funcgroup_{name.split()[0].lower()}.svg"
        with open(filename, 'w', encoding='utf-8') as f:
            f.write(svg)
        print(f"    Saved: {filename}")

print("\n" + "=" * 70)
print("✓ Tests complete! OH, NO2, etc. should be preserved.")
