#!/usr/bin/env python3
"""
Generate SVG from SMILES using RDKit
Called by Node.js server via subprocess
"""

import sys
import json
from io import StringIO
from datetime import datetime

try:
    from rdkit import Chem
    from rdkit.Chem import Draw
    import re

    def smiles_to_svg(smiles, options=None):
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

            # Draw molecule to SVG with options
            drawer = Draw.MolDraw2DSVG(-1, -1)

            # Apply drawing options if available
            try:
                draw_options = drawer.drawOptions()
                draw_options.clearBackground = False
                # Note: RDKit doesn't have a direct "aromatic circles" option
                # It draws aromatic bonds as Kekulé structure by default
                # The aromaticCircles option is more relevant for mol2chemfig/chemfig
            except Exception:
                pass  # Older RDKit versions may not have drawOptions

            drawer.DrawMolecule(mol)
            drawer.FinishDrawing()
            svg = drawer.GetDrawingText()

            # Make background transparent by removing/modifying the background rect
            # RDKit generates a white rect as background - we want it transparent
            # RDKit format: <rect style='opacity:1.0;fill:#FFFFFF;stroke:none' width='300.0' height='200.0' x='0.0' y='0.0'> </rect>

            # Remove rect with style containing fill:#FFFFFF - handles both /> and > </rect> endings
            svg = re.sub(
                r"<rect[^>]*style='[^']*fill:#FFFFFF[^']*'[^>]*>[\s]*</rect>",
                "",
                svg
            )
            svg = re.sub(
                r"<rect[^>]*style='[^']*fill:#FFFFFF[^']*'[^>]*/>",
                "",
                svg
            )

            # Remove rect with fill='#FFFFFF' attribute
            svg = re.sub(
                r"<rect[^>]*fill='#FFFFFF'[^>]*>[\s]*</rect>",
                "",
                svg
            )
            svg = re.sub(
                r"<rect[^>]*fill='#FFFFFF'[^>]*/>",
                "",
                svg
            )

            # Remove any rect that has white fill (case insensitive variations)
            svg = re.sub(
                r"<rect[^>]*fill='white'[^>]*>[\s]*</rect>",
                "",
                svg,
                flags=re.IGNORECASE
            )

            # Final fallback: Remove the first rect element after the header comment
            # This catches any format RDKit might use
            svg = re.sub(
                r"(<!-- END OF HEADER -->)\s*<rect[^>]*>[\s]*</rect>",
                r"\1",
                svg
            )

            # Embed SMILES metadata in SVG for deduplication
            # Get canonical SMILES
            canonical_smiles = Chem.MolToSmiles(mol, canonical=True)

            # Insert metadata comment after SVG declaration
            svg_lines = svg.split('\n')
            if len(svg_lines) > 1:
                svg_lines.insert(1, f'<!-- SMILES: {canonical_smiles} -->')
                svg_lines.insert(2, f'<!-- Generated: {datetime.now().isoformat()} -->')
                svg = '\n'.join(svg_lines)

            return {"svg": svg}

        except Exception as e:
            return {"error": str(e)}

    # Main execution
    if __name__ == "__main__":
        input_data = json.loads(sys.argv[1])
        result = smiles_to_svg(
            input_data['smiles'],
            input_data.get('options', {})
        )
        print(json.dumps(result))

except ImportError:
    print(json.dumps({"error": "RDKit not installed. Run: pip install rdkit"}))
    sys.exit(1)
