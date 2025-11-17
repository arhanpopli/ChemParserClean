#!/usr/bin/env python3
"""
Canonicalize SMILES strings using RDKit
This ensures the same molecule always produces the same cache key
Called by Node.js server via subprocess
"""

import sys
import json

try:
    from rdkit import Chem

    def canonicalize_smiles(smiles):
        """
        Convert any SMILES representation to canonical SMILES

        Args:
            smiles: SMILES string (can be non-canonical)

        Returns:
            dict with 'canonical_smiles' or 'error'
        """
        try:
            smiles = smiles.strip()
            if not smiles:
                return {"error": "Empty SMILES string"}

            # Parse SMILES
            mol = Chem.MolFromSmiles(smiles)
            if mol is None:
                return {"error": f"Invalid SMILES: {smiles}"}

            # Convert to canonical SMILES
            canonical = Chem.MolToSmiles(mol, canonical=True)

            return {"canonical_smiles": canonical}

        except Exception as e:
            return {"error": str(e)}

    # Main execution
    if __name__ == "__main__":
        smiles_input = sys.argv[1] if len(sys.argv) > 1 else ""
        result = canonicalize_smiles(smiles_input)
        print(json.dumps(result))

except ImportError:
    print(json.dumps({"error": "RDKit not installed. Run: pip install rdkit"}))
    sys.exit(1)
