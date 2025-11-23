from rdkit import Chem
from rdkit.Chem import AllChem

def add_hydrogens_to_structure(text_data):
    try:
        mol = Chem.MolFromSmiles(text_data)
        if mol is None:
            return "Error: Invalid SMILES"
        
        mol_h = Chem.AddHs(mol)
        return Chem.MolToSmiles(mol_h, isomericSmiles=True, allHsExplicit=True)
    except Exception as e:
        return f"Error: {e}"

smiles = "OC1=CC=CC=C1"
explicit = add_hydrogens_to_structure(smiles)
print(f"Original: {smiles}")
print(f"Explicit: {explicit}")
