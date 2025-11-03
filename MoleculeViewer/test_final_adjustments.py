import sys
sys.path.insert(0, '.')
from app.chemistry import smiles_to_svg

# Test benzene with aromatic circles and methyls
test_cases = [
    ('c1ccccc1', 'Benzene with aromatic circle', {'aromatic_circles': True}),
    ('Cc1ccccc1', 'Toluene with CH3 and circle', {'aromatic_circles': True, 'show_methyls': True}),
    ('CC(C)C', 'Isobutane with methyls', {'show_methyls': True}),
]

for smiles, desc, options in test_cases:
    print(f'\nTesting: {desc} ({smiles})')
    error, svg = smiles_to_svg(smiles, 400, 300, options)
    
    if error:
        print(f'  ✗ Error: {error}')
    else:
        print(f'  ✓ SVG generated')
        
        if options.get('aromatic_circles'):
            has_circle = '<circle' in svg
            print(f'    Has aromatic circle: {has_circle}')
        
        if options.get('show_methyls'):
            ch3_count = svg.count('CH₃')
            print(f'    CH₃ labels: {ch3_count}')
            # Check font size
            if 'font-size:32px' in svg:
                print(f'    Font size: 32px ✓')
            elif 'font-size:' in svg:
                import re
                sizes = re.findall(r'font-size:(\d+)px', svg)
                print(f'    Font sizes found: {set(sizes)}')
        
        # Save SVG
        filename = f'test_final_{smiles.replace("(", "").replace(")", "").replace("1", "").replace("2", "").replace("3", "").replace("4", "").replace("5", "").replace("6", "")[:20]}.svg'
        with open(filename, 'w', encoding='utf-8') as f:
            f.write(svg)
        print(f'    Saved: {filename}')

print('\n' + '='*60)
print('✓ All tests complete!')
print('Check the SVG files to verify:')
print('  1. CH₃ labels are 32px (slightly smaller)')
print('  2. Aromatic circles are better centered')
