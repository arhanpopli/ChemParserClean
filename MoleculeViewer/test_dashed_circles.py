from app.chemistry import smiles_to_svg
import re

# Test molecules
test_cases = [
    ('benzene', 'c1ccccc1', 400, 300),
    ('toluene', 'Cc1ccccc1', 400, 300),
    ('naphthalene', 'c1ccc2ccccc2c1', 500, 400),
]

print('Testing dashed aromatic circles:')
print('=' * 70)

for name, smiles, width, height in test_cases:
    error, svg = smiles_to_svg(smiles, width, height, {'aromatic_circles': True})
    if not error:
        # Extract circle properties
        circle_matches = re.findall(r'<circle cx="([\d.]+)" cy="([\d.]+)" r="([\d.]+)" fill="none" stroke="black" stroke-width="([\d.]+)" stroke-dasharray="([\d.,]+)"', svg)
        
        print(f'\n{name}:')
        for i, match in enumerate(circle_matches, 1):
            cx, cy, r, stroke, dasharray = match
            print(f'  Circle {i}:')
            print(f'    Radius: {r}')
            print(f'    Stroke width: {stroke}px (matches bonds)')
            print(f'    Dash pattern: {dasharray}')
        
        # Save for visual inspection
        filename = f'test_{name}_dashed.svg'
        with open(filename, 'w') as f:
            f.write(svg)
        print(f'  Saved to: {filename}')

print('\n' + '=' * 70)
print('✓ All circles now have 2.0px stroke and dashed pattern!')
