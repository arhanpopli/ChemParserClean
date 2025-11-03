"""
Visual comparison of different CH₃ label sizes.
Shows same molecule with different config settings.
"""

import sys
sys.path.insert(0, '.')
from app import config
from app.chemistry import smiles_to_svg

# Test molecule
SMILES = 'CC(C)C'  # Isobutane
DESC = 'Isobutane (3 methyls)'

print("="*70)
print(f"VISUAL COMPARISON TEST: {DESC}")
print("="*70)

# Save original settings
orig_scaling = config.CARBON_LABEL_SCALING
orig_size = config.CARBON_LABEL_FONT_SIZE
orig_factor = config.CARBON_LABEL_SCALE_FACTOR

test_configs = [
    ('Current (auto, 0.55)', 'auto', 32, 0.55),
    ('Larger (auto, 0.65)', 'auto', 36, 0.65),
    ('Smaller (auto, 0.45)', 'auto', 24, 0.45),
    ('Fixed 28px', 'fixed', 28, 0.55),
    ('Fixed 36px', 'fixed', 36, 0.55),
]

for name, scaling, size, factor in test_configs:
    # Temporarily modify config
    config.CARBON_LABEL_SCALING = scaling
    config.CARBON_LABEL_FONT_SIZE = size
    config.CARBON_LABEL_SCALE_FACTOR = factor
    
    print(f"\n{name}")
    print(f"  Settings: scaling={scaling}, size={size}, factor={factor}")
    
    error, svg = smiles_to_svg(SMILES, 400, 300, {'show_methyls': True})
    
    if error:
        print(f"  ✗ Error: {error}")
    else:
        import re
        font_sizes = re.findall(r'font-size:([\d.]+)px', svg)
        if font_sizes:
            actual_size = float(font_sizes[0])
            print(f"  Actual font size: {actual_size:.1f}px")
        
        filename = f'compare_{name.replace(" ", "_").replace("(", "").replace(")", "").replace(",", "")}.svg'
        with open(filename, 'w', encoding='utf-8') as f:
            f.write(svg)
        print(f"  Saved: {filename}")

# Restore original settings
config.CARBON_LABEL_SCALING = orig_scaling
config.CARBON_LABEL_FONT_SIZE = orig_size
config.CARBON_LABEL_SCALE_FACTOR = orig_factor

print("\n" + "="*70)
print("✓ Comparison complete!")
print("\nOpen the generated SVG files side-by-side to compare sizes.")
print("Choose the setting that looks best, then update app/config.py")
print("="*70)
