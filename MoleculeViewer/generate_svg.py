#!/usr/bin/env python3
"""
Generate SVG from SMILES using RDKit
Called by Node.js server via subprocess
"""

import sys
import json
from io import StringIO

try:
    from rdkit import Chem
    from rdkit.Chem import Draw
    import re

    def smiles_to_svg(smiles, width=300, height=200, options=None):
        """Convert SMILES to SVG using RDKit"""
        if options is None:
            options = {}

        try:
            # Sanitize SMILES
            smiles = smiles.strip()
            if not smiles:
                return {"error": "Empty SMILES string"}

            # Create molecule
            mol = Chem.MolFromSmiles(smiles)
            if mol is None:
                return {"error": f"Invalid SMILES: {smiles}"}

            # Generate 2D coordinates
            from rdkit.Chem import AllChem
            AllChem.Compute2DCoords(mol)

            # Draw molecule to SVG
            drawer = Draw.MolDraw2DSVG(width, height)
            drawer.DrawMolecule(mol)
            drawer.FinishDrawing()
            svg = drawer.GetDrawingText()

            return {"svg": svg}

        except Exception as e:
            return {"error": str(e)}

    # Main execution
    if __name__ == "__main__":
        input_data = json.loads(sys.argv[1])
        result = smiles_to_svg(
            input_data['smiles'],
            input_data.get('width', 300),
            input_data.get('height', 200),
            input_data.get('options', {})
        )
        print(json.dumps(result))

except ImportError:
    print(json.dumps({"error": "RDKit not installed. Run: pip install rdkit"}))
    sys.exit(1)
