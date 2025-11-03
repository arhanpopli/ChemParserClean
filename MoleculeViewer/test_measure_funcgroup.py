import sys
sys.path.insert(0, '.')
from rdkit import Chem
from rdkit.Chem import Draw
import re

# Generate SVG with NO2 and OH groups
mol = Chem.MolFromSmiles('CCO')  # Ethanol - simple OH
drawer = Draw.MolDraw2DSVG(400, 300)
drawer.DrawMolecule(mol)
drawer.FinishDrawing()
svg = drawer.GetDrawingText()

# RDKit renders atoms as vector paths
# Measure the dimensions
path_pattern = r"<path class='atom-(\d+)' d='M ([\d.]+) ([\d.]+)([^']+)' fill='[^']+'/>"
matches = list(re.finditer(path_pattern, svg))

print(f'Found {len(matches)} atom labels in ethanol SVG')
print()

for i, match in enumerate(matches):
    atom_idx = match.group(1)
    start_x = float(match.group(2))
    start_y = float(match.group(3))
    full_path = match.group(4)
    
    # Extract all coordinate pairs from the path
    coord_pairs = re.findall(r'L?\s*([\d.]+)\s+([\d.]+)', full_path)
    if coord_pairs:
        x_coords = [float(x) for x, y in coord_pairs]
        y_coords = [float(y) for x, y in coord_pairs]
        
        # Add starting point
        x_coords.append(start_x)
        y_coords.append(start_y)
        
        width = max(x_coords) - min(x_coords)
        height = max(y_coords) - min(y_coords)
        
        print(f'Atom {atom_idx}: ~{width:.1f}px × {height:.1f}px')

# Save for visual inspection
with open('test_rdkit_funcgroup_size.svg', 'w') as f:
    f.write(svg)
print('\nSaved: test_rdkit_funcgroup_size.svg')
print('RDKit functional groups are typically 28-35px tall')
