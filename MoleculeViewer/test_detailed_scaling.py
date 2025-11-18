from app.chemistry import smiles_to_svg
from rdkit import Chem
import re
import math

# Test molecules of different sizes
test_cases = [
    ('benzene', 'c1ccccc1', 400, 300),
    ('toluene', 'Cc1ccccc1', 400, 300),
    ('aspirin', 'CC(=O)Oc1ccccc1C(=O)O', 500, 400),
    ('adrenaline', 'CNCC(O)C1=CC(O)=C(O)C=C1', 500, 400),
    ('ibuprofen', 'CC(C)Cc1ccc(cc1)C(C)C(=O)O', 600, 500),
]

print('Detailed circle sizing analysis:')
print('=' * 70)

for name, smiles, width, height in test_cases:
    print(f'\n{name.upper()}:')
    print('-' * 70)
    
    # Get aromatic ring atoms
    mol = Chem.MolFromSmiles(smiles)
    aromatic_rings = []
    for ring in mol.GetRingInfo().AtomRings():
        ring_atoms = [mol.GetAtomWithIdx(i) for i in ring]
        if all(atom.GetIsAromatic() for atom in ring_atoms):
            aromatic_rings.append(ring)
            print(f'Aromatic ring atoms: {ring}')
    
    # Generate SVG
    error, svg = smiles_to_svg(smiles, width, height, {'aromatic_circles': True})
    if not error:
        # Extract circle radius
        circle_match = re.search(r'<circle cx="([\d.]+)" cy="([\d.]+)" r="([\d.]+)"', svg)
        if circle_match:
            cx, cy, r = circle_match.groups()
            print(f'Circle: center=({cx}, {cy}) radius={r}')
            
            # Extract all bond coordinates
            bond_pattern = r'<path class=.bond-\d+ atom-(\d+) atom-(\d+). d=.M ([\d.]+),([\d.]+) L ([\d.]+),([\d.]+)'
            bond_matches = re.findall(bond_pattern, svg)
            
            # Find bonds within the aromatic ring
            if aromatic_rings:
                ring = set(aromatic_rings[0])
                ring_bonds = []
                
                for match in bond_matches:
                    atom1, atom2 = int(match[0]), int(match[1])
                    # Check if this bond is in the aromatic ring
                    if atom1 in ring and atom2 in ring:
                        x1, y1 = float(match[2]), float(match[3])
                        x2, y2 = float(match[4]), float(match[5])
                        bond_len = math.sqrt((x2-x1)**2 + (y2-y1)**2)
                        ring_bonds.append(bond_len)
                
                if ring_bonds:
                    avg_ring_bond = sum(ring_bonds) / len(ring_bonds)
                    min_ring_bond = min(ring_bonds)
                    max_ring_bond = max(ring_bonds)
                    
                    print(f'Aromatic ring bonds:')
                    print(f'  Count: {len(ring_bonds)}')
                    print(f'  Average: {avg_ring_bond:.2f}')
                    print(f'  Min: {min_ring_bond:.2f}')
                    print(f'  Max: {max_ring_bond:.2f}')
                    print(f'  Circle radius / Avg ring bond: {float(r) / avg_ring_bond:.3f}')
                    print(f'  Circle should fit nicely in ring!')

print('\n' + '=' * 70)
