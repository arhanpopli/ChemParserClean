import json
import rdkit
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem.Draw import rdMolDraw2D

# List of molecules (Name, SMILES)
molecules_list = [
    # Simple
    ("ethane", "CC"),
    ("propane", "CCC"),
    ("methanol", "CO"),
    ("ethanol", "CCO"),
    ("dimethylEther", "COC"),
    ("acetone", "CC(=O)C"),
    ("acetaldehyde", "CC=O"),
    ("isobutane", "CC(C)C"),
    
    # Rings
    ("benzene", "c1ccccc1"),
    ("cyclohexane", "C1CCCCC1"),
    ("pyridine", "c1ccncc1"),
    ("THF", "C1CCOC1"),
    ("phenol", "Oc1ccccc1"),
    ("toluene", "Cc1ccccc1"),
    ("catechol", "Oc1c(O)cccc1"),
    ("naphthalene", "c1ccc2ccccc2c1"),
    ("anisole", "COc1ccccc1"),
    ("dioxane", "C1COCCO1"),
    
    # Functional / Medium
    ("aceticAcid", "CC(=O)O"),
    ("histamine", "NCCc1c[nH]cn1"),
    ("amphetamine", "CC(N)Cc1ccccc1"),
    ("diethylEther", "CCOCC"),
    ("diisopropylEther", "CC(C)OC(C)C"),
    ("MTBE", "COC(C)(C)C"),
    ("aspirin", "CC(=O)Oc1ccccc1C(=O)O"),
    ("paracetamol", "CC(=O)Nc1ccc(O)cc1"),
    ("ibuprofen", "CC(C)Cc1ccc(C(C)C(=O)O)cc1"),
    ("dopamine", "NCCc1ccc(O)c(O)c1"),
    ("serotonin", "NCCc1c[nH]c2ccc(O)cc12"),
    
    # Larger / Requested
    ("caffeine", "Cn1cnc2c1c(=O)n(C)c(=O)n2C"),
    ("methylphenidate", "COC(=O)C(C1CCCCN1)c2ccccc2"),
    ("atomoxetine", "CNCCC(Oc1ccccc1)c2ccccc2C"),
    ("adrenaline", "CNc1cc(O)c(O)cc1C(O)CN"),
    ("morphine", "CN1CCC23C4C1Cc5ccc(O)c(O5)c2C3(O)C=C4"),
]

def generate_svg(mol, name, style="rdkit"):
    try:
        # Compute 2D coordinates
        AllChem.Compute2DCoords(mol)
        
        # Get conformer to calculate bounds
        conf = mol.GetConformer()
        min_x = min_y = float('inf')
        max_x = max_y = float('-inf')
        
        for i in range(mol.GetNumAtoms()):
            pos = conf.GetAtomPosition(i)
            min_x = min(min_x, pos.x)
            min_y = min(min_y, pos.y)
            max_x = max(max_x, pos.x)
            max_y = max(max_y, pos.y)
            
        # Calculate dimensions based on a fixed scale (pixels per Angstrom)
        scale = 30.0 # pixels per unit
        
        width = (max_x - min_x) * scale + 60 # Add padding
        height = (max_y - min_y) * scale + 60
        
        # Ensure minimum size
        width = max(width, 100)
        height = max(height, 100)
        
        d = rdMolDraw2D.MolDraw2DSVG(int(width), int(height))
        opts = d.drawOptions()
        
        # Style configuration
        if style == "mol2chemfig":
            # Black and white, clean
            if hasattr(opts, 'useBWAtomPalette'):
                opts.useBWAtomPalette()
            if hasattr(opts, 'bondLineWidth'):
                opts.bondLineWidth = 2
            opts.clearBackground = False 
        else:
            # RDKit style (colored)
            if hasattr(opts, 'bondLineWidth'):
                opts.bondLineWidth = 2
            opts.clearBackground = False
            
        d.DrawMolecule(mol)
        d.FinishDrawing()
        svg = d.GetDrawingText()
        
        # Clean up SVG string if needed (remove newlines)
        svg = svg.replace('\n', '')
        
        return svg
    except Exception as e:
        print(f"Error in generate_svg for {name}: {e}")
        return ""

output_data = {
    "rdkit": {},
    "mol2chemfig": {}
}

for name, smiles in molecules_list:
    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol:
            # Generate RDKit style
            svg_rdkit = generate_svg(mol, name, "rdkit")
            if svg_rdkit:
                output_data["rdkit"][name] = svg_rdkit
            
            # Generate Mol2ChemFig style
            svg_m2cf = generate_svg(mol, name, "mol2chemfig")
            if svg_m2cf:
                output_data["mol2chemfig"][name] = svg_m2cf
            
            print(f"Generated {name}")
        else:
            print(f"Failed to parse {name}")
    except Exception as e:
        print(f"Error generating {name}: {e}")

with open("generated_molecules.json", "w") as f:
    json.dump(output_data, f, indent=2)
