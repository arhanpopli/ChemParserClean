def _add_carbon_labels(svg, mol, show_carbons=False, show_methyls=False):
    """Post-process SVG to add carbon and methyl labels."""
    try:
        import re
        bond_pattern = r'<path class=[\'"]bond-\d+ atom-(\d+) atom-(\d+)[\'"] d=[\'"]M ([\d.]+),([\d.]+) L ([\d.]+),([\d.]+)'
        bond_matches = re.findall(bond_pattern, svg)
        atom_coords = {}
        for match in bond_matches:
            atom1_idx = int(match[0])
            atom2_idx = int(match[1])
            x1, y1 = float(match[2]), float(match[3])
            x2, y2 = float(match[4]), float(match[5])
            if atom1_idx not in atom_coords:
                atom_coords[atom1_idx] = []
            atom_coords[atom1_idx].append((x1, y1))
            if atom2_idx not in atom_coords:
                atom_coords[atom2_idx] = []
            atom_coords[atom2_idx].append((x2, y2))
        atom_positions = {}
        for atom_idx, coords in atom_coords.items():
            avg_x = sum(x for x, y in coords) / len(coords)
            avg_y = sum(y for x, y in coords) / len(coords)
            atom_positions[atom_idx] = (avg_x, avg_y)
        labels_svg = []
        for atom_idx, (x, y) in atom_positions.items():
            if atom_idx >= mol.GetNumAtoms():
                continue
            atom = mol.GetAtomWithIdx(atom_idx)
            if atom.GetSymbol() != 'C':
                continue
            label = None
            num_h = atom.GetTotalNumHs()
            if show_methyls and num_h == 3:
                label = 'CH₃'
            elif show_methyls and num_h == 2:
                label = 'CH₂'
            elif show_methyls and num_h == 1:
                label = 'CH'
            elif show_carbons:
                label = 'C'
            if label:
                label_svg = f'<text x="{x:.2f}" y="{y:.2f}" class="atom-{atom_idx}" style="font-size:14px;font-family:sans-serif;text-anchor:middle;dominant-baseline:middle;fill:black">{label}</text>'
                labels_svg.append(label_svg)
        if labels_svg:
            labels_group = '\n'.join(labels_svg)
            svg = svg.replace('</svg>', labels_group + '\n</svg>', 1)
        return svg
    except Exception as e:
        return svg
