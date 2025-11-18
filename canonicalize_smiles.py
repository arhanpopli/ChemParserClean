#!/usr/bin/env python3
"""
Canonicalize SMILES strings using RDKit
This ensures the same molecule always produces the same cache key
Shared utility for mol2chemfig_server.py
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
            Canonical SMILES string or None if invalid
        """
        try:
            smiles = smiles.strip()
            if not smiles:
                return None

            # Parse SMILES
            mol = Chem.MolFromSmiles(smiles)
            if mol is None:
                return None

            # Convert to canonical SMILES
            canonical = Chem.MolToSmiles(mol, canonical=True)

            return canonical

        except Exception as e:
            print(f"Error canonicalizing SMILES: {e}", file=sys.stderr)
            return None

    # Main execution (for command-line usage)
    if __name__ == "__main__":
        smiles_input = sys.argv[1] if len(sys.argv) > 1 else ""
        result = canonicalize_smiles(smiles_input)
        if result:
            print(json.dumps({"canonical_smiles": result}))
        else:
            print(json.dumps({"error": "Could not canonicalize SMILES"}))

except ImportError:
    print(json.dumps({"error": "RDKit not installed. Run: pip install rdkit"}))
    sys.exit(1)
