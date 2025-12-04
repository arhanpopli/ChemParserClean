"""
Simple test to verify mol2chemfig works with RPC mode.
This bypasses the file reading logic.
"""
import sys
import os

# Add mol2chemfig to path
sys.path.insert(0, os.path.join(os.path.dirname(__file__), 'mol2chemfig'))

# Monkey-patch the imports before importing processor
import mol2chemfig.common as common
import mol2chemfig.options as options
import mol2chemfig.molecule as molecule

# Now import processor
from mol2chemfig import processor

print("Testing native mol2chemfig with rpc=True...")
print("SMILES: CCO")

try:
    success, result = processor.process(data='CCO', rpc=True)
    
    if success:
        print("✅ Success!")
        print(f"Result type: {type(result)}")
        chemfig_code = result.render_server()
        print(f"LaTeX code: {chemfig_code}")
    else:
        print(f"❌ Failed: {result}")
except Exception as e:
    print(f"❌ Exception: {e}")
    import traceback
    traceback.print_exc()
