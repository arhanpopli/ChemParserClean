from app.chemistry import smiles_to_svg
import re
import math

# Test molecules of different sizes
test_cases = [
    ('benzene', 'c1ccccc1', 400, 300),
    ('toluene', 'Cc1ccccc1', 400, 300),
    ('aspirin', 'CC(=O)Oc1ccccc1C(=O)O', 500, 400),
    ('ibuprofen', 'CC(C)Cc1ccc(cc1)C(C)C(=O)O', 600, 500),
]

print('Testing circle sizing across different molecule sizes:')
print('=' * 60)

for name, smiles, width, height in test_cases:
    error, svg = smiles_to_svg(smiles, width, height, {'aromatic_circles': True})
    if not error:
        # Extract circle radius
        circle_match = re.search(r'<circle cx="([\d.]+)" cy="([\d.]+)" r="([\d.]+)"', svg)
        if circle_match:
            cx, cy, r = circle_match.groups()
            
            # Also extract aromatic ring bond lengths
            # Find bonds in the aromatic ring by looking for consecutive atoms
            bond_matches = re.findall(r'<path class=.bond-\d+ atom-(\d+) atom-(\d+). d=.M ([\d.]+),([\d.]+) L ([\d.]+),([\d.]+)', svg)
            
            # Calculate average bond length in the molecule
            bond_lengths = []
            for match in bond_matches[:10]:  # Sample first 10 bonds
                x1, y1 = float(match[2]), float(match[3])
                x2, y2 = float(match[4]), float(match[5])
                bond_len = math.sqrt((x2-x1)**2 + (y2-y1)**2)
                bond_lengths.append(bond_len)
            
            if bond_lengths:
                avg_bond = sum(bond_lengths) / len(bond_lengths)
                
                print(f'{name}:')
                print(f'  Circle radius: {r}')
                print(f'  Avg bond length: {avg_bond:.2f}')
                print(f'  Ratio (radius/bond): {float(r)/avg_bond:.3f}')
                print()
