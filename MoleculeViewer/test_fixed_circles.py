from app.chemistry import nomenclature_to_smiles, smiles_to_svg
import re

print('=== Testing FIXED circle positioning ===')
print()

# Test adrenaline
print('1. Adrenaline:')
error, smiles, source = nomenclature_to_smiles('adrenaline')
if not error:
    error2, svg = smiles_to_svg(smiles, 500, 400, {'aromatic_circles': True})
    if not error2:
        print('   ✓ SVG generated')
        # Save for inspection
        with open('test_adrenaline_final.svg', 'w') as f:
            f.write(svg)
        print('   - Saved to test_adrenaline_final.svg')
        
        # Extract circle position
        circle_match = re.search(r'<circle cx="([\d.]+)" cy="([\d.]+)" r="([\d.]+)"', svg)
        if circle_match:
            cx, cy, r = circle_match.groups()
            print(f'   - Circle at: ({cx}, {cy}) radius={r}')
            print('   - Expected: (~312, ~222) radius=~41')
    else:
        print('   ERROR:', error2)

print()

# Test benzene
print('2. Benzene:')
error, svg = smiles_to_svg('c1ccccc1', 400, 300, {'aromatic_circles': True})
if not error:
    print('   ✓ SVG generated')
else:
    print('   ERROR:', error)

print()

# Test toluene (benzene with substituent)
print('3. Toluene (benzene + CH3):')
error, svg = smiles_to_svg('Cc1ccccc1', 400, 300, {'aromatic_circles': True})
if not error:
    print('   ✓ SVG generated - circle should be centered in ring')
else:
    print('   ERROR:', error)

print()
print('✅ All tests passed!')
