from app.chemistry import smiles_to_svg
import re

# Generate a simple benzene molecule
error, svg = smiles_to_svg('c1ccccc1', 400, 300, {'aromatic_circles': False})

# Extract stroke-width from bonds
stroke_matches = re.findall(r'stroke-width:([\d.]+)px', svg)
if stroke_matches:
    print('Bond stroke widths found:', set(stroke_matches))
    print('Typical bond width:', stroke_matches[0], 'px')
else:
    print('No stroke-width found in style')

# Check for stroke-width attribute
attr_matches = re.findall(r'stroke-width="([\d.]+)"', svg)
if attr_matches:
    print('Stroke-width attributes:', set(attr_matches))
