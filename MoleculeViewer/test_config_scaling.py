"""
Test script demonstrating the configurable CH₃ label sizing system.

This script shows how the auto-scaling works with different molecule sizes
and how to adjust the config settings for custom sizing.
"""

import sys
sys.path.insert(0, '.')

# Import and temporarily modify config for testing
from app import config
from app.chemistry import smiles_to_svg

print("="*70)
print("CARBON LABEL CONFIGURATION TEST")
print("="*70)

# Test molecules of different sizes
test_cases = [
    ('C', 'Methane (very small)', 200, 150),
    ('CCC', 'Propane (small)', 300, 200),
    ('CCCCCCCC', 'Octane (medium)', 500, 300),
    ('CC(C)C', 'Isobutane (branched)', 400, 300),
    ('Cc1ccccc1', 'Toluene (aromatic)', 400, 300),
]

print(f"\nCurrent Configuration:")
print(f"  Scaling Mode: {config.CARBON_LABEL_SCALING}")
print(f"  Base Font Size: {config.CARBON_LABEL_FONT_SIZE}px")
print(f"  Scale Factor: {config.CARBON_LABEL_SCALE_FACTOR}")
print(f"  Min Size: {config.CARBON_LABEL_MIN_SIZE}px")
print(f"  Max Size: {config.CARBON_LABEL_MAX_SIZE}px")
print(f"  Offsets: X={config.CARBON_LABEL_OFFSET_X}, Y={config.CARBON_LABEL_OFFSET_Y}")
print(f"  Scale Offsets: {config.CARBON_LABEL_SCALE_OFFSETS}")
print()

for smiles, desc, width, height in test_cases:
    print(f"\nTesting: {desc}")
    print(f"  SMILES: {smiles}")
    print(f"  Canvas: {width}x{height}px")
    
    error, svg = smiles_to_svg(smiles, width, height, {'show_methyls': True})
    
    if error:
        print(f"  ✗ Error: {error}")
    else:
        print(f"  ✓ SVG generated")
        
        # Extract font sizes from SVG
        import re
        font_sizes = re.findall(r'font-size:([\d.]+)px', svg)
        if font_sizes:
            sizes = [float(s) for s in font_sizes]
            print(f"  CH₃ Font Size(s): {', '.join(f'{s:.1f}px' for s in set(sizes))}")
        
        ch3_count = svg.count('CH₃')
        print(f"  CH₃ labels: {ch3_count}")
        
        # Save for inspection
        filename = f'test_config_{smiles.replace("(", "").replace(")", "")[:10]}.svg'
        with open(filename, 'w', encoding='utf-8') as f:
            f.write(svg)
        print(f"  Saved: {filename}")

print("\n" + "="*70)
print("CONFIGURATION EXAMPLES")
print("="*70)

print("""
To customize CH₃ label sizes, edit app/config.py:

Example 1: Fixed size (always 28px)
-----------------------------------
CARBON_LABEL_SCALING = 'fixed'
CARBON_LABEL_FONT_SIZE = 28

Example 2: Auto-scaling (current default)
------------------------------------------
CARBON_LABEL_SCALING = 'auto'
CARBON_LABEL_SCALE_FACTOR = 0.55  # Adjust this for size
CARBON_LABEL_MIN_SIZE = 18
CARBON_LABEL_MAX_SIZE = 48

Example 3: Larger labels for all molecules
-------------------------------------------
CARBON_LABEL_FONT_SIZE = 40  # Or adjust SCALE_FACTOR to 0.65

Example 4: Smaller labels for all molecules
--------------------------------------------
CARBON_LABEL_FONT_SIZE = 24  # Or adjust SCALE_FACTOR to 0.45

Example 5: Fine-tune positioning
---------------------------------
CARBON_LABEL_OFFSET_X = 10  # Move right
CARBON_LABEL_OFFSET_Y = 12  # Move down
CARBON_LABEL_SCALE_OFFSETS = True  # Scale with font size

""")

print("="*70)
print("✓ Test complete!")
print("="*70)
