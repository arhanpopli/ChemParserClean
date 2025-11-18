import re
import math

# Read the fixed SVG
with open('test_adrenaline_fixed.svg', 'r') as f:
    svg = f.read()

# Extract all bonds with their atom indices and coordinates
bond_pattern = r'<path class=.bond-\d+ atom-(\d+) atom-(\d+). d=.M ([\d.]+),([\d.]+) L ([\d.]+),([\d.]+)'
matches = re.findall(bond_pattern, svg)

# Aromatic ring atoms
ring_atoms = {5, 6, 7, 9, 11, 12}

print('Bonds in the aromatic ring:')
ring_coords = {}
for match in matches:
    atom1, atom2, x1, y1, x2, y2 = match
    atom1, atom2 = int(atom1), int(atom2)
    
    # Check if both atoms are in the ring
    if atom1 in ring_atoms and atom2 in ring_atoms:
        print(f'  Bond: atom {atom1} <-> atom {atom2}')
        print(f'    Atom {atom1} at: ({x1}, {y1})')
        print(f'    Atom {atom2} at: ({x2}, {y2})')
        
        # Store positions (use first occurrence)
        if atom1 not in ring_coords:
            ring_coords[atom1] = (float(x1), float(y1))
        if atom2 not in ring_coords:
            ring_coords[atom2] = (float(x2), float(y2))

print()
print('Ring atom positions extracted:')
for atom_idx in sorted(ring_coords.keys()):
    x, y = ring_coords[atom_idx]
    print(f'  Atom {atom_idx}: ({x:.1f}, {y:.1f})')

# Calculate center
if ring_coords:
    center_x = sum(x for x, y in ring_coords.values()) / len(ring_coords)
    center_y = sum(y for x, y in ring_coords.values()) / len(ring_coords)
    print()
    print(f'Calculated ring center: ({center_x:.2f}, {center_y:.2f})')
    
    # Calculate radius
    min_dist = min(math.sqrt((x - center_x)**2 + (y - center_y)**2) 
                   for x, y in ring_coords.values())
    radius = min_dist * 0.6  # 60% for 6-membered ring
    print(f'Calculated radius (60% of inradius): {radius:.2f}')
    
    # Extract what's currently in the SVG
    circle_match = re.search(r'<circle cx="([\d.]+)" cy="([\d.]+)" r="([\d.]+)"', svg)
    if circle_match:
        current_x, current_y, current_r = circle_match.groups()
        print()
        print(f'Current circle in SVG: ({current_x}, {current_y}) r={current_r}')
        print()
        print(f'PROBLEM: Circle should be at ({center_x:.2f}, {center_y:.2f}) r={radius:.2f}')
        print(f'         But it is at ({current_x}, {current_y}) r={current_r}')
