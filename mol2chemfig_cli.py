#!/usr/bin/env python
"""
Command-line wrapper for mol2chemfig processor.
Reads SMILES from stdin and outputs LaTeX chemfig code.
"""
import sys
import os

# Add mol2chemfig directory to path so Python 2-style imports work
mol2chemfig_dir = os.path.join(os.path.dirname(__file__), 'mol2chemfig')
sys.path.insert(0, mol2chemfig_dir)

# Now import processor (this will find common, options, molecule in the same dir)
import processor

if __name__ == "__main__":
    # Get SMILES from stdin
    smiles = sys.stdin.read().strip()
    
    # Process with mol2chemfig
    success, result = processor.process(rawargs=sys.argv[1:], data=smiles, rpc=True)
    
    if success:
        # Output LaTeX code
        print(result.render_server())
        sys.exit(0)
    else:
        # Output error
        print(result, file=sys.stderr)
        sys.exit(1)
