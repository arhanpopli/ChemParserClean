from app.chemistry import smiles_to_svg

print('Testing show_carbons and show_methyls options:')
print('=' * 70)

# Test 1: Hexane (simple chain) - should show C at each vertex
print('\n1. Hexane (CCCCCC) - Testing show_carbons:')
error, svg = smiles_to_svg('CCCCCC', 500, 300, {'show_carbons': True})
if not error:
    print('   ✓ SVG generated')
    c_count = svg.count('CH₃') + svg.count('>C<')
    print(f'   - Carbon labels found: {c_count}')
    with open('test_hexane_carbons.svg', 'w', encoding='utf-8') as f:
        f.write(svg)
    print('   - Saved to: test_hexane_carbons.svg')
else:
    print(f'   ERROR: {error}')

# Test 2: Propane - should show CH3 at ends
print('\n2. Propane (CCC) - Testing show_methyls:')
error, svg = smiles_to_svg('CCC', 400, 300, {'show_methyls': True})
if not error:
    print('   ✓ SVG generated')
    ch3_count = svg.count('CH₃') + svg.count('CH3')
    print(f'   - CH₃ labels found: {ch3_count}')
    with open('test_propane_methyls.svg', 'w', encoding='utf-8') as f:
        f.write(svg)
    print('   - Saved to: test_propane_methyls.svg')
else:
    print(f'   ERROR: {error}')

# Test 3: Toluene (benzene + CH3) - test methyls on aromatic
print('\n3. Toluene (Cc1ccccc1) - Testing show_methyls:')
error, svg = smiles_to_svg('Cc1ccccc1', 400, 300, {'show_methyls': True})
if not error:
    print('   ✓ SVG generated')
    ch3_count = svg.count('CH₃') + svg.count('CH3')
    print(f'   - CH₃ labels found: {ch3_count}')
    with open('test_toluene_methyls.svg', 'w', encoding='utf-8') as f:
        f.write(svg)
    print('   - Saved to: test_toluene_methyls.svg')
else:
    print(f'   ERROR: {error}')

# Test 4: Butane - both options
print('\n4. Butane (CCCC) - Testing BOTH options:')
error, svg = smiles_to_svg('CCCC', 500, 300, {'show_carbons': True, 'show_methyls': True})
if not error:
    print('   ✓ SVG generated')
    with open('test_butane_both.svg', 'w', encoding='utf-8') as f:
        f.write(svg)
    print('   - Saved to: test_butane_both.svg')
else:
    print(f'   ERROR: {error}')

print('\n' + '=' * 70)
print('✓ Tests complete! Check the SVG files.')
