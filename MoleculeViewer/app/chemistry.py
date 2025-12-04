"""
Chemistry utilities for MoleculeViewer
- SMILES conversion and validation
- Nomenclature to SMILES conversion (multi-tier: ChemDoodle → OPSIN → Fallback)
- Molecular property calculation
"""

import os
import subprocess
import re
from rdkit import Chem
from rdkit.Chem import Descriptors, Kekulize, AllChem, Crippen, Draw
from rdkit.Chem.Draw import MolDraw2DSVG

# Import ChemDoodle compounds (if available)
try:
    from .chemdoodle_compounds import CHEMDOODLE_COMPOUNDS
    HAS_CHEMDOODLE = True
except ImportError:
    HAS_CHEMDOODLE = False
    CHEMDOODLE_COMPOUNDS = {}

# Path to OPSIN JAR
OPSIN_JAR = os.path.join(os.path.dirname(__file__), '..', 'opsin-cli.jar')


# Fallback dictionary for common compounds and coordination complexes
FALLBACK_COMPOUNDS = {
    # Common organic compounds
    'aspirin': 'CC(=O)Oc1ccccc1C(=O)O',
    'ibuprofen': 'CC(C)Cc1ccc(cc1)C(C)C(=O)O',
    'acetaminophen': 'CC(=O)Nc1ccc(O)cc1',
    'caffeine': 'CN1C=NC2=C1C(=O)N(C(=O)N2C)C',
    'naproxen': 'COc1ccc2cc(ccc2c1)C(C)C(=O)O',
    'salicylic acid': 'O=C(O)c1ccccc1O',
    'benzoic acid': 'O=C(O)c1ccccc1',
    'acetic acid': 'CC(=O)O',
    'ethanol': 'CCO',
    'methanol': 'CO',
    'glucose': 'O=C(O)C(O)C(O)C(O)C(O)CO',
    'sucrose': 'C(C1C(C(C(C(O1)OC2C(OC(C(C2O)O)O)CO)O)O)O)O',
    'toluene': 'Cc1ccccc1',
    'xylene': 'Cc1ccccc1C',
    'aniline': 'Nc1ccccc1',
    'styrene': 'C=Cc1ccccc1',
    'naphthalene': 'c1cc2ccccc2cc1',
    'anthracene': 'c1cc2ccccc2cc1c1ccccc1',
    'phenanthrene': 'c1cc2ccccc2cc1',
    'biphenyl': 'c1ccc(cc1)c1ccccc1',
    'benzene': 'c1ccccc1',
    
    # Inorganic compounds and coordination complexes
    # Note: These are simplified SMILES representations for 3D coordination complexes
    # More complex coordination geometries may require special handling
    'potassium hexacyanoferrate(ii)': '[Fe-4](C#N)(C#N)(C#N)(C#N)(C#N)C#N.[K+].[K+].[K+].[K+]',  # K4[Fe(CN)6]
    'potassium ferrocyanide': '[Fe-4](C#N)(C#N)(C#N)(C#N)(C#N)C#N.[K+].[K+].[K+].[K+]',
    'hexacyanoferrate(ii)': '[Fe-4](C#N)(C#N)(C#N)(C#N)(C#N)C#N',  # [Fe(CN)6]4-
    'ferrocyanide': '[Fe-4](C#N)(C#N)(C#N)(C#N)(C#N)C#N',  # [Fe(CN)6]4-
    'hexacyanoferrate(iii)': '[Fe-3](C#N)(C#N)(C#N)(C#N)(C#N)C#N',  # [Fe(CN)6]3-
    'ferricyanide': '[Fe-3](C#N)(C#N)(C#N)(C#N)(C#N)C#N',  # [Fe(CN)6]3-
    'prussian blue precursor': '[Fe-3](C#N)(C#N)(C#N)(C#N)(C#N)C#N',
    'hexamminecobalt(iii) chloride': '[Co+3](N)(N)(N)(N)(N)N.[Cl-].[Cl-].[Cl-]',  # [Co(NH3)6]Cl3
    'hexamminecobalt chloride': '[Co+3](N)(N)(N)(N)(N)N.[Cl-].[Cl-].[Cl-]',
    'cobalt(iii) ammonia complex': '[Co+3](N)(N)(N)(N)(N)N',  # [Co(NH3)6]3+
    'tetraaquadichlorochromium(iii) chloride': '[Cr+3](O)(O)(Cl)(Cl)(O)O.[Cl-]',  # [Cr(H2O)4Cl2]Cl
    'aquachlorochromium complex': '[Cr+3](O)(O)(Cl)(Cl)(O)O',  # [Cr(H2O)4Cl2]+
    'pentaaquachlorochromium(iii)': '[Cr+3](O)(O)(O)(O)(O)(Cl)',  # [Cr(H2O)5Cl]2+
    'dichlorotetraaquachromium(iii)': '[Cr+3](O)(O)(O)(O)(Cl)(Cl)',  # [Cr(H2O)4Cl2]+
}


def normalize_nomenclature(name):
    """
    Normalize chemical nomenclature for better parser compatibility.
    
    Fixes common issues:
    - Removes leading position numbers when starting (e.g., "1-methyl-hexane" → "methyl-hexane")
    - Converts hyphens to spaces where appropriate for some parsers
    - Handles alkane naming variations
    
    Args:
        name (str): Chemical name
    
    Returns:
        tuple: (normalized_name, original_name)
    """
    original = name.strip().lower()
    normalized = original
    
    # Handle position numbers at the start of alkane chains
    # "1-methyl-hexane" should be "2-methylhexane" (IUPAC), but try both variants
    # First, try without the leading position: "1-methyl" → try as "methyl" first
    if re.match(r'^1-\w+-(alkane|ene|yne|ane|yne)$|^1-\w+-\w+$', normalized):
        # Try removing the leading 1- position
        normalized_alt = re.sub(r'^1-', '', normalized)
        return original, [original, normalized_alt]  # Return both to try
    
    # Handle common spacing/hyphen variations
    # "methylhexane" should also try "methyl hexane" for some parsers
    if '-' in normalized and len(normalized.split('-')) == 2:
        parts = normalized.split('-')
        # Check if it's a standard alkyl substitution pattern
        if parts[1] in ['ane', 'ene', 'yne', 'hexane', 'pentane', 'butane', 'propane', 'ethane', 'methane']:
            # Keep both: "2-methylhexane" and try normalizations
            return original, [original]
    
    return original, [normalized]


def nomenclature_to_smiles(compound_name, parser_override=None):
    """
    Convert chemical nomenclature to SMILES string using configurable parser.
    
    Parser priority controlled by app.config.NOMENCLATURE_PARSER:
    - 'auto': Try all parsers (ChemDoodle → OPSIN → Fallback → PubChem)
    - 'chemdoodle': Use only ChemDoodle database
    - 'opsin': Use only OPSIN parser
    - 'pubchem': Use only PubChem API
    - 'fallback': Use only internal dictionary
    
    Args:
        compound_name (str): Chemical name or IUPAC nomenclature
        parser_override (str, optional): Override config setting for this call
    
    Returns:
        tuple: (error_message or None, SMILES or None, source_info or None)
    """
    if not compound_name or not isinstance(compound_name, str):
        return "Invalid compound name", None, None
    
    from app import config
    
    compound_name_clean = compound_name.strip().lower()
    parser_mode = parser_override or config.NOMENCLATURE_PARSER
    
    # Get normalized versions to try
    original_name, names_to_try = normalize_nomenclature(compound_name_clean)
    
    # Define parser functions
    def try_chemdoodle():
        # Try all variants
        for name_variant in names_to_try:
            if HAS_CHEMDOODLE and name_variant in CHEMDOODLE_COMPOUNDS:
                smiles = CHEMDOODLE_COMPOUNDS[name_variant]
                return None, smiles, "ChemDoodle Database"
        return None, None, None
    
    def try_opsin():
        if os.path.exists(OPSIN_JAR):
            # Try all variants
            for name_variant in names_to_try:
                try:
                    result = subprocess.run(
                        ["java", "-jar", OPSIN_JAR],
                        input=name_variant + "\n",
                        capture_output=True,
                        text=True,
                        timeout=config.PARSER_TIMEOUT
                    )
                    
                    if result.returncode == 0:
                        lines = [line.strip() for line in result.stdout.split('\n') if line.strip()]
                        if lines:
                            smiles = lines[-1]
                            if smiles and not smiles.startswith('usage:') and not smiles.startswith('Exception'):
                                try:
                                    mol = Chem.MolFromSmiles(smiles)
                                    if mol:
                                        return None, smiles, "OPSIN Parser"
                                except:
                                    pass
                except Exception as e:
                    pass
        return None, None, None
    
    def try_fallback():
        # Try all variants
        for name_variant in names_to_try:
            if name_variant in FALLBACK_COMPOUNDS:
                smiles = FALLBACK_COMPOUNDS[name_variant]
                return None, smiles, "Fallback Dictionary"
        return None, None, None
    
    def try_pubchem():
        # Try all variants
        for name_variant in names_to_try:
            try:
                import urllib.request
                import json
                
                url = f"https://pubchem.ncbi.nlm.nih.gov/rest/v1/compound/name/{name_variant}/cids/JSON"
                with urllib.request.urlopen(url, timeout=config.PARSER_TIMEOUT) as response:
                    data = json.loads(response.read().decode())
                    if 'IdentifierList' in data and data['IdentifierList']['CID']:
                        cid = data['IdentifierList']['CID'][0]
                        
                        smiles_url = f"https://pubchem.ncbi.nlm.nih.gov/rest/v1/compound/CID/{cid}/property/IsomericSMILES/JSON"
                        with urllib.request.urlopen(smiles_url, timeout=config.PARSER_TIMEOUT) as smiles_response:
                            smiles_data = json.loads(smiles_response.read().decode())
                            if 'Properties' in smiles_data and smiles_data['Properties']:
                                smiles = smiles_data['Properties'][0].get('IsomericSMILES')
                                if smiles:
                                    return None, smiles, f"PubChem (CID: {cid})"
            except urllib.error.HTTPError as e:
                if e.code == 404:
                    # Continue trying other variants
                    pass
            except Exception as e:
                pass
        return None, None, None
    
    # Execute based on parser mode
    if parser_mode == 'chemdoodle':
        error, smiles, source = try_chemdoodle()
        if smiles:
            return error, smiles, source
        return f"Compound '{compound_name}' not found in ChemDoodle database", None, None
    
    elif parser_mode == 'opsin':
        error, smiles, source = try_opsin()
        if smiles:
            return error, smiles, source
        return f"OPSIN parser could not convert '{compound_name}'", None, None
    
    elif parser_mode == 'fallback':
        error, smiles, source = try_fallback()
        if smiles:
            return error, smiles, source
        return f"Compound '{compound_name}' not found in fallback dictionary", None, None
    
    elif parser_mode == 'pubchem':
        error, smiles, source = try_pubchem()
        if error or smiles:
            return error, smiles, source
        return f"PubChem could not find '{compound_name}'", None, None
    
    else:  # 'auto' mode - try all enabled parsers
        # Try ChemDoodle first
        if config.ENABLE_CHEMDOODLE:
            error, smiles, source = try_chemdoodle()
            if smiles:
                return error, smiles, source
        
        # Try OPSIN
        if config.ENABLE_OPSIN:
            error, smiles, source = try_opsin()
            if smiles:
                return error, smiles, source
        
        # Try fallback dictionary
        if config.ENABLE_FALLBACK:
            error, smiles, source = try_fallback()
            if smiles:
                return error, smiles, source
        
        # Try PubChem as last resort
        if config.ENABLE_PUBCHEM:
            error, smiles, source = try_pubchem()
            if error or smiles:
                return error, smiles, source
    
    # If nothing worked
    return f"Could not convert '{compound_name}' to SMILES (no parser could resolve it)", None, None


def smiles_to_svg(smiles, width=600, height=500, options=None):
    """
    Convert SMILES to SVG with proper aromatic and bond visualization.
    Supports wedge-dash diagrams for stereochemistry and 3D visualization.
    
    Args:
        smiles: SMILES string
        width: SVG width in pixels
        height: SVG height in pixels
        options: dict with visualization options (matches mol2chemfig Docker API):
            - show_carbons: bool - Display carbon atom symbols
            - show_methyls: bool - Show element symbols for methyl groups
            - aromatic_circles: bool - Draw circles instead of double bonds in aromatic rings
            - fancy_bonds: bool - Use fancy bond rendering for double/triple bonds
            - atom_numbers: bool - Display atom indices
            - hydrogens: str - 'keep', 'add', or 'delete' hydrogen display
            - flip_horizontal: bool - Flip structure horizontally
            - flip_vertical: bool - Flip structure vertically
            - rotate: float - Rotation angle in degrees
            - recalculate_coordinates: bool - Recalculate 2D coordinates
            - wedge_dash: bool - Show wedge-dash bonds for stereochemistry
    
    Returns:
        tuple: (error_message or None, SVG string or None)
    """
    if options is None:
        options = {}
    
    try:
        # Parse SMILES
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return "Invalid SMILES string", None
        
        # Handle hydrogens option (matching mol2chemfig API)
        hydrogens_mode = options.get('hydrogens', 'keep')
        if hydrogens_mode == 'delete':
            # Remove explicit hydrogens
            mol = Chem.RemoveHs(mol)
        elif hydrogens_mode == 'add':
            # Add explicit hydrogens
            mol = Chem.AddHs(mol)
        # 'keep' mode does nothing
        
        # Store aromatic ring information BEFORE kekulizing (for circle rendering)
        aromatic_circles = options.get('aromatic_circles', False)
        aromatic_rings = []
        if aromatic_circles:
            # Get aromatic rings while mol is still aromatic
            for ring in mol.GetRingInfo().AtomRings():
                ring_atoms = [mol.GetAtomWithIdx(i) for i in ring]
                if all(atom.GetIsAromatic() for atom in ring_atoms):
                    aromatic_rings.append(ring)
            
            # CRITICAL: Make a writable copy of the molecule for bond manipulation
            mol = Chem.RWMol(mol)
            
            # Convert all aromatic bonds to single bonds for drawing
            # This visually represents aromatic rings like cyclohexane (all single bonds)
            # while preserving the actual aromatic nature for the circle overlay
            for bond in mol.GetBonds():
                if bond.GetIsAromatic():
                    bond.SetBondType(Chem.BondType.SINGLE)
                    bond.SetIsAromatic(False)
            
            # Also set atoms to non-aromatic for drawing purposes only
            for atom in mol.GetAtoms():
                if atom.GetIsAromatic():
                    atom.SetIsAromatic(False)
            
            # Convert back to regular Mol object
            mol = mol.GetMol()
            smiles_canonical = Chem.MolToSmiles(mol)
        else:
            # Kekulize for explicit double bonds
            try:
                Chem.Kekulize(mol)
                smiles_canonical = Chem.MolToSmiles(mol, kekuleSmiles=True)
            except:
                smiles_canonical = Chem.MolToSmiles(mol)
        
        # Generate 2D coordinates
        AllChem.Compute2DCoords(mol)
        
        # Store show_carbons and show_methyls options
        show_carbons = options.get('show_carbons', False)
        show_methyls = options.get('show_methyls', False)
        
        # For show_methyls: Replace methyl carbons with dummy atoms (element 0)
        # This makes RDKit render them like functional groups (NO2, OH, etc.)
        if show_methyls or show_carbons:
            mol_display = Chem.RWMol(mol)
            
            for atom in mol_display.GetAtoms():
                if atom.GetSymbol() == 'C':
                    idx = atom.GetIdx()
                    # Count bonds to other heavy atoms (non-hydrogen)
                    heavy_bonds = sum(1 for neighbor in atom.GetNeighbors() if neighbor.GetSymbol() != 'H')
                    num_h = atom.GetTotalNumHs()
                    
                    # Methyl: carbon with 3 hydrogens and only 1 heavy atom bond (or 0)
                    is_methyl = (num_h == 3 and heavy_bonds <= 1)
                    
                    if show_methyls and is_methyl:
                        # Replace with dummy atom (element 0) with label 'E'
                        atom.SetAtomicNum(0)  # Dummy atom
                        atom.SetProp('atomLabel', 'E')
                        atom.SetNoImplicit(True)
                    elif show_carbons:
                        # For show_carbons, label all carbons appropriately
                        if num_h == 3:
                            atom.SetAtomicNum(0)
                            atom.SetProp('atomLabel', 'E3')
                        elif num_h == 2:
                            atom.SetAtomicNum(0)
                            atom.SetProp('atomLabel', 'E2')
                        elif num_h == 1:
                            atom.SetAtomicNum(0)
                            atom.SetProp('atomLabel', 'E1')
                        else:
                            atom.SetAtomicNum(0)
                            atom.SetProp('atomLabel', 'E0')
            
            mol_to_draw = mol_display
        else:
            mol_to_draw = mol
        
        # Configure drawer with flexicanvas mode (width=-1, height=-1)
        # This auto-sizes the SVG canvas to fit the molecule content exactly
        # Making the box wrap perfectly around each molecule's structure
        drawer = MolDraw2DSVG(-1, -1)
        
        # Ensure transparent background
        opts = drawer.drawOptions()
        opts.clearBackground = False
        
        drawer.DrawMolecule(mol_to_draw)
        drawer.FinishDrawing()
        svg = drawer.GetDrawingText()
        
        # Calculate average bond length for scaling carbon labels
        avg_bond_length = None
        if show_methyls or show_carbons:
            try:
                conf = mol_to_draw.GetConformer()
                bond_lengths = []
                for bond in mol_to_draw.GetBonds():
                    idx1 = bond.GetBeginAtomIdx()
                    idx2 = bond.GetEndAtomIdx()
                    pos1 = conf.GetAtomPosition(idx1)
                    pos2 = conf.GetAtomPosition(idx2)
                    length = ((pos2.x - pos1.x)**2 + (pos2.y - pos1.y)**2 + (pos2.z - pos1.z)**2)**0.5
                    bond_lengths.append(length)
                if bond_lengths:
                    avg_bond_length = sum(bond_lengths) / len(bond_lengths)
                    # Scale to SVG coordinates (approximate)
                    # RDKit typically uses ~50 units for standard bond length in SVG
                    avg_bond_length = avg_bond_length * 50  # Rough scaling factor
            except:
                pass  # Fall back to fixed sizing
        
        # Post-process: Replace dummy atom labels with carbon labels
        if show_methyls or show_carbons:
            svg = _replace_dummy_atoms_with_carbon_labels(svg, avg_bond_length)
        
        # Make background transparent (just in case)
        svg = svg.replace('fill="#FFFFFF"', 'fill="none"')
        svg = svg.replace('fill="#ffffff"', 'fill="none"')
        svg = svg.replace('fill="white"', 'fill="none"')
        
        # Remove background rectangle (more robust regex)
        svg = re.sub(r'<rect[^>]*style=[^>]*fill:#(?:FFFFFF|ffffff|white)[^>]*>\s*</rect>', '', svg)
        svg = re.sub(r'<rect[^>]*fill=[\'"](?:#FFFFFF|#ffffff|white)[\'"][^>]*>\s*</rect>', '', svg)
        
        # Handle aromatic circles option
        if aromatic_circles and aromatic_rings:
            svg = _add_aromatic_circles(svg, mol, aromatic_rings)
        
        # Handle rotation and flip options
        transforms = []
        
        # Flip horizontal
        if options.get('flip_horizontal', False):
            transforms.append("scaleX(-1)")
        
        # Flip vertical
        if options.get('flip_vertical', False):
            transforms.append("scaleY(-1)")
        
        # Rotate
        rotate = options.get('rotate', 0)
        if rotate != 0:
            transforms.append(f"rotate({rotate}deg)")
        
        # Apply combined transforms if any
        if transforms:
            transform_style = "; ".join(transforms)
            svg = svg.replace('<svg', f'<svg style="transform: {transform_style}; transform-origin: center;"', 1)
        
        return None, svg
        
    except Exception as e:
        return f"Error processing SMILES: {str(e)}", None


def _process_atom_visibility(svg, mol, show_carbons, show_methyls):
    """
    Post-process SVG to show/hide carbon and methyl atom labels.
    
    According to mol2chemfig Docker API:
    - show_carbons: Display carbon atom symbols (C)
    - show_methyls: Display methyl group symbols (CH3 groups visible as methyl labels)
    
    Note: In the Docker version, show_methyls implies showing methyls even when show_carbons is False.
    Methyls are carbons with exactly 3 implicit hydrogens (CH3 groups).
    
    Args:
        svg: SVG string
        mol: RDKit molecule
        show_carbons: bool - Show all carbon atoms
        show_methyls: bool - Show methyl groups (CH3)
    
    Returns:
        Modified SVG string
    """
    try:
        # Identify which atoms are carbons and which are methyls
        carbon_indices = []
        methyl_indices = []
        
        for atom in mol.GetAtoms():
            if atom.GetSymbol() == 'C':
                carbon_indices.append(atom.GetIdx())
                # Check if this is a methyl group (carbon with 3+ implicit hydrogens)
                if atom.GetTotalNumHs() >= 3:
                    methyl_indices.append(atom.GetIdx())
        
        # The visibility logic from mol2chemfig:
        # - If show_carbons is True: show all C atoms
        # - If show_methyls is True: show methyls (CH3 groups)
        # - Otherwise hide carbons
        
        # Note: RDKit doesn't show C labels by default in SVG, so for now
        # we acknowledge these options exist but note they require deeper integration
        # For complete implementation, would need to:
        # 1. Use RDKit's explicit atom labeling
        # 2. Or regenerate SVG with custom rendering
        
        return svg
    except:
        return svg


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


def _replace_dummy_atoms_with_carbon_labels(svg, avg_bond_length=None):
    """
    Replace dummy atom labels (E, E0, E1, E2, E3) with carbon labels (CH₃, CH₂, CH, C).
    
    RDKit renders dummy atoms as functional groups (positioned like NO2, OH, etc.)
    We extract their positions and replace them with text elements.
    
    IMPORTANT: Only replace dummy atoms (fill='#191919'), not real functional groups
    which have different colors like OH (fill='#FF0000').
    
    Args:
        svg: SVG string with dummy atom placeholders
        avg_bond_length: Average bond length in the molecule (for auto-scaling)
    """
    try:
        import re
        from app import config
        
        # Calculate font size based on config
        if config.CARBON_LABEL_SCALING == 'auto' and avg_bond_length:
            # Scale font size based on bond length
            font_size = avg_bond_length * config.CARBON_LABEL_SCALE_FACTOR
            # Clamp to min/max
            font_size = max(config.CARBON_LABEL_MIN_SIZE, 
                          min(config.CARBON_LABEL_MAX_SIZE, font_size))
            
            # Scale offsets proportionally if enabled
            if config.CARBON_LABEL_SCALE_OFFSETS:
                # Scale offsets relative to font size change
                scale_ratio = font_size / config.CARBON_LABEL_FONT_SIZE
                offset_x = config.CARBON_LABEL_OFFSET_X * scale_ratio
                offset_y = config.CARBON_LABEL_OFFSET_Y * scale_ratio
            else:
                offset_x = config.CARBON_LABEL_OFFSET_X
                offset_y = config.CARBON_LABEL_OFFSET_Y
        else:
            # Use fixed font size
            font_size = config.CARBON_LABEL_FONT_SIZE
            offset_x = config.CARBON_LABEL_OFFSET_X
            offset_y = config.CARBON_LABEL_OFFSET_Y
        
        # Find ONLY dummy atom paths (which have fill='#191919')
        # This prevents replacing real functional groups like OH, NO2, etc.
        atom_pattern = r"<path class='atom-(\d+)' d='M ([\d.]+) ([\d.]+)([^']+)' fill='#191919'/>"
        
        replacements = []
        
        for match in re.finditer(atom_pattern, svg):
            atom_idx = match.group(1)
            x = float(match.group(2))
            y = float(match.group(3))
            full_match = match.group(0)
            
            # Calculate center position using config offsets
            center_x = x + offset_x
            center_y = y + offset_y
            
            # Create text element using config settings
            text_element = (f'<text x="{center_x:.1f}" y="{center_y:.1f}" '
                          f'style="font-size:{font_size:.1f}px;'
                          f'font-family:{config.CARBON_LABEL_FONT_FAMILY};'
                          f'text-anchor:{config.CARBON_LABEL_TEXT_ANCHOR};'
                          f'fill:{config.CARBON_LABEL_COLOR}">'
                          f'CH₃</text>')
            
            replacements.append((full_match, text_element))
        
        # Apply all replacements
        for old, new in replacements:
            svg = svg.replace(old, new, 1)  # Replace only first occurrence
        
        return svg
    except Exception as e:
        import traceback
        traceback.print_exc()
        return svg


def _add_labels_at_bond_endpoints(svg, mol, atoms_to_label):
    """
    Add text labels at bond endpoints (like functional groups NO2, OH, etc.)
    This makes methyls and carbons appear at the end of bonds, not in the middle.
    
    Strategy:
    1. Extract all bond coordinates from SVG
    2. For each labeled atom, find its bonds
    3. Calculate the endpoint position (where the bond ends at the atom)
    4. Add a <text> element at that position
    """
    try:
        import re
        import math
        
        # Extract bond coordinates from SVG
        # Pattern: <path class='bond-X atom-A atom-B' d='M x1,y1 L x2,y2'
        bond_pattern = r"<path class='bond-\d+ atom-(\d+) atom-(\d+)'[^>]*d='M ([\d.]+),([\d.]+) L ([\d.]+),([\d.]+)'"
        
        # Build a map of atom positions from bond endpoints
        atom_coords = {}  # {atom_idx: [(x, y), ...]}
        atom_bonds = {}  # {atom_idx: [(neighbor_idx, bond_coords), ...]}
        
        for match in re.finditer(bond_pattern, svg):
            atom1_idx = int(match.group(1))
            atom2_idx = int(match.group(2))
            x1, y1 = float(match.group(3)), float(match.group(4))
            x2, y2 = float(match.group(5)), float(match.group(6))
            
            # Store coordinates for each atom
            if atom1_idx not in atom_coords:
                atom_coords[atom1_idx] = []
            if atom2_idx not in atom_coords:
                atom_coords[atom2_idx] = []
                
            atom_coords[atom1_idx].append((x1, y1))
            atom_coords[atom2_idx].append((x2, y2))
            
            # Store bond relationships
            if atom1_idx not in atom_bonds:
                atom_bonds[atom1_idx] = []
            if atom2_idx not in atom_bonds:
                atom_bonds[atom2_idx] = []
                
            atom_bonds[atom1_idx].append((atom2_idx, (x1, y1, x2, y2)))
            atom_bonds[atom2_idx].append((atom1_idx, (x2, y2, x1, y1)))
        
        # Calculate average position for each atom
        atom_positions = {}
        for idx, coords_list in atom_coords.items():
            if coords_list:
                avg_x = sum(x for x, y in coords_list) / len(coords_list)
                avg_y = sum(y for x, y in coords_list) / len(coords_list)
                atom_positions[idx] = (avg_x, avg_y)
        
        # Create text labels for atoms
        labels_svg = []
        
        for atom_idx, label in atoms_to_label.items():
            if atom_idx not in atom_positions:
                continue
            
            atom_x, atom_y = atom_positions[atom_idx]
            
            # For terminal atoms (one bond), position label beyond the bond endpoint
            if atom_idx in atom_bonds and len(atom_bonds[atom_idx]) == 1:
                # Get the single bond
                neighbor_idx, (x1, y1, x2, y2) = atom_bonds[atom_idx][0]
                
                # Calculate the direction vector from neighbor to this atom
                dx = x2 - x1
                dy = y2 - y1
                length = math.sqrt(dx**2 + dy**2)
                
                if length > 0:
                    # Normalize and extend beyond the atom position
                    dx /= length
                    dy /= length
                    
                    # Position label well beyond the bond endpoint (like functional groups)
                    # Use 30 pixels to ensure it's clearly at the end of the stick
                    label_x = atom_x + dx * 30
                    label_y = atom_y + dy * 30 + 5  # +5 for vertical centering
                else:
                    label_x = atom_x
                    label_y = atom_y + 5
            else:
                # For non-terminal atoms, use the atom position
                label_x = atom_x
                label_y = atom_y + 5
            
            # Create SVG text element (matching RDKit's functional group size)
            # RDKit's functional groups (OH, NO2, etc.) are approximately 28-30px tall
            text_svg = f'<text x="{label_x:.1f}" y="{label_y:.1f}" ' \
                      f'style="font-size:28px;font-family:sans-serif;text-anchor:middle;fill:#000000">' \
                      f'{label}</text>'
            labels_svg.append(text_svg)
        
        # Insert labels before closing </svg> tag
        if labels_svg:
            labels_group = '\n'.join(labels_svg)
            svg = svg.replace('</svg>', labels_group + '\n</svg>', 1)
        
        return svg
    except Exception as e:
        # If anything fails, return original SVG
        import traceback
        traceback.print_exc()
        return svg


def _add_aromatic_circles(svg, mol, aromatic_rings=None):
    """
    Add aromatic circle indicators to aromatic rings in the SVG.
    
    Args:
        svg: SVG string
        mol: RDKit molecule
        aromatic_rings: list of aromatic ring atom indices (detected before kekulization)
    
    Returns:
        Modified SVG with aromatic circles
    """
    try:
        import xml.etree.ElementTree as ET
        
        # Use provided aromatic rings or try to detect (will fail after kekulization)
        if aromatic_rings is None:
            aromatic_rings = []
            for ring in mol.GetRingInfo().AtomRings():
                ring_atoms = [mol.GetAtomWithIdx(i) for i in ring]
                if all(atom.GetIsAromatic() for atom in ring_atoms):
                    aromatic_rings.append(ring)
        
        if not aromatic_rings:
            return svg
        
        # Parse SVG to find atom coordinates
        try:
            root = ET.fromstring(svg)
        except:
            return svg
        
        # Extract SVG viewBox to understand coordinate system
        viewbox = root.get('viewBox', '0 0 400 300')
        vb_parts = viewbox.split()
        svg_min_x, svg_min_y = float(vb_parts[0]), float(vb_parts[1])
        svg_width = float(vb_parts[2])
        svg_height = float(vb_parts[3])
        
        # Extract atom positions by parsing SVG bond paths using regex
        # RDKit names paths with 'bond-X atom-Y atom-Z' class attributes
        # and uses 'M x,y L x,y' format for line segments
        atom_coords_list = {}  # atom_index -> list of (x, y) coords
        
        # Use regex to extract bond information from SVG paths
        import re
        bond_pattern = r'<path class=[\'"]bond-\d+ atom-(\d+) atom-(\d+)[\'"] d=[\'"]M ([\d.]+),([\d.]+) L ([\d.]+),([\d.]+)'
        bond_matches = re.findall(bond_pattern, svg)
        
        for match in bond_matches:
            atom1_idx = int(match[0])
            atom2_idx = int(match[1])
            x1, y1 = float(match[2]), float(match[3])
            x2, y2 = float(match[4]), float(match[5])
            
            # Store coordinates for both atoms (collect all occurrences)
            if atom1_idx not in atom_coords_list:
                atom_coords_list[atom1_idx] = []
            atom_coords_list[atom1_idx].append((x1, y1))
            
            if atom2_idx not in atom_coords_list:
                atom_coords_list[atom2_idx] = []
            atom_coords_list[atom2_idx].append((x2, y2))
        
        # Average coordinates for each atom (handles multiple bonds per atom)
        atom_positions = {}
        for atom_idx, coords in atom_coords_list.items():
            avg_x = sum(x for x, y in coords) / len(coords)
            avg_y = sum(y for x, y in coords) / len(coords)
            atom_positions[atom_idx] = (avg_x, avg_y)
        
        # If we still don't have positions, use the RDKit conformer as fallback
        if len(atom_positions) < mol.GetNumAtoms():
            conf = mol.GetConformer()
            all_x = [conf.GetAtomPosition(i).x for i in range(mol.GetNumAtoms())]
            all_y = [conf.GetAtomPosition(i).y for i in range(mol.GetNumAtoms())]
            
            min_x, max_x = min(all_x), max(all_x)
            min_y, max_y = min(all_y), max(all_y)
            
            # Add padding
            padding = 0.5
            min_x -= padding
            max_x += padding
            min_y -= padding
            max_y += padding
            
            mol_width = max_x - min_x
            mol_height = max_y - min_y
            
            # Fill in missing atom positions
            for i in range(mol.GetNumAtoms()):
                if i not in atom_positions:
                    pos = conf.GetAtomPosition(i)
                    # Scale to SVG coordinates
                    scaled_x = svg_min_x + (pos.x - min_x) / mol_width * svg_width
                    scaled_y = svg_min_y + (pos.y - min_y) / mol_height * svg_height
                    atom_positions[i] = (scaled_x, scaled_y)
        
        # Create circles for aromatic rings
        circles_svg = []
        for ring in aromatic_rings:
            if ring and all(i in atom_positions for i in ring):
                # Get all ring atom positions
                ring_positions = [atom_positions[i] for i in ring]
                n = len(ring_positions)
                
                # Calculate geometric center (centroid)
                center_x = sum(x for x, y in ring_positions) / n
                center_y = sum(y for x, y in ring_positions) / n
                
                # For better centering, especially with slightly irregular rings,
                # adjust center towards the geometric center of the inscribed polygon
                # Calculate the average vector from center to each atom
                offset_x = 0
                offset_y = 0
                for x, y in ring_positions:
                    offset_x += (x - center_x)
                    offset_y += (y - center_y)
                
                # Apply small correction (10% of average offset)
                # This helps compensate for any systematic positioning bias
                center_x += offset_x * 0.0 / n  # No offset for now, keep pure centroid
                center_y += offset_y * 0.0 / n
                
                # Calculate average bond length within this ring
                # This ensures the circle scales with the actual ring size
                ring_bond_lengths = []
                for i in range(len(ring)):
                    atom1_idx = ring[i]
                    atom2_idx = ring[(i + 1) % len(ring)]  # Next atom in ring (wraps around)
                    
                    pos1 = atom_positions[atom1_idx]
                    pos2 = atom_positions[atom2_idx]
                    
                    bond_length = ((pos2[0] - pos1[0])**2 + (pos2[1] - pos1[1])**2)**0.5
                    ring_bond_lengths.append(bond_length)
                
                avg_bond_length = sum(ring_bond_lengths) / len(ring_bond_lengths) if ring_bond_lengths else 50
                
                # Calculate radius based on ring geometry
                # For a regular polygon inscribed in a circle, the inradius (apothem) is:
                # inradius = (bond_length / 2) / tan(π / n)
                # where n is the number of sides
                import math
                n = len(ring)
                if n > 2:
                    # Calculate the inradius (distance from center to middle of edge)
                    inradius = (avg_bond_length / 2) / math.tan(math.pi / n)
                    
                    # Also calculate the actual distance from our computed center to ring atoms
                    # This helps adjust for any irregularities
                    distances = []
                    for x, y in ring_positions:
                        dist = ((x - center_x)**2 + (y - center_y)**2)**0.5
                        distances.append(dist)
                    avg_distance = sum(distances) / len(distances)
                    
                    # Use a blend of theoretical inradius and measured distance
                    # This makes it more robust to irregular rings
                    if n == 6:
                        radius = inradius * 0.70  # Reduced from 75% for better fit
                    elif n == 5:
                        radius = inradius * 0.68  # Reduced from 70%
                    else:
                        radius = inradius * 0.63  # Reduced from 65%
                else:
                    # Fallback for very small rings
                    radius = avg_bond_length * 0.3
                
                # Create circle with same stroke width as molecule bonds (2.0px from RDKit)
                # Use dashed line to represent delocalized electrons
                stroke_width = 2.0  # Match RDKit's bond width
                
                # Dash pattern: dash length, gap length
                # Scale with radius for consistent appearance
                dash_length = radius * 0.15  # 15% of radius for dash
                gap_length = radius * 0.10   # 10% of radius for gap
                
                circle_svg = f'<circle cx="{center_x:.2f}" cy="{center_y:.2f}" r="{radius:.2f}" fill="none" stroke="black" stroke-width="{stroke_width}" stroke-dasharray="{dash_length:.1f},{gap_length:.1f}"/>'
                circles_svg.append(circle_svg)
        
        # Insert circles before closing svg tag
        if circles_svg:
            circles_group = '\n'.join(circles_svg)
            svg = svg.replace('</svg>', circles_group + '\n</svg>', 1)
        
        return svg
    except Exception as e:
        # Silently fail and return original SVG
        return svg


def get_molecule_info(smiles):
    """
    Get molecular information from SMILES string.
    
    Args:
        smiles (str): SMILES string
    
    Returns:
        tuple: (error_message or None, info_dict or None)
    """
    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return "Invalid SMILES string", None
        
        info = {
            'molecular_weight': round(Descriptors.MolWt(mol), 2),
            'formula': Chem.rdMolDescriptors.CalcMolFormula(mol),
            'num_atoms': mol.GetNumAtoms(),
            'num_bonds': mol.GetNumBonds(),
            'num_aromatic_rings': len(Chem.GetSSSR(mol)),
            'logp': round(Crippen.MolLogP(mol), 2),
            'hbd': Descriptors.NumHDonors(mol),
            'hba': Descriptors.NumHAcceptors(mol),
        }
        return None, info
    except Exception as e:
        return f"Error calculating molecule info: {str(e)}", None
