# MolView Search API Enhancement - Chemical Data Fields

## Overview
Updated the MolView Search API (`search-server.js` on port 8001) to provide comprehensive chemical data for all molecule types: compounds, minerals, and biomolecules.

## New Fields Added

### 1. `chemical_name` (String)
- **Purpose**: Proper chemical name (IUPAC or systematic name)
- **Compounds**: IUPAC name from PubChem (e.g., "Acetylsalicylic acid" for aspirin)
- **Minerals**: Chemical name from CIF file (e.g., "Silicon oxide" for quartz)
- **Biomolecules**: Currently null (can be enhanced later with protein names)

### 2. `formula` (String)
- **Purpose**: Chemical formula in standard notation
- **Compounds**: Molecular formula from PubChem (e.g., "C9H8O4" for aspirin)
- **Minerals**: Formula from CIF file (e.g., "O2 Si" for quartz, "Pt Sn" for niggliite)
- **Biomolecules**: Currently null (can be enhanced later)

### 3. Enhanced SMILES for Minerals
- **`canonical_smiles`**: Extracted from CIF `_chemical_smiles_canonical` tag
- **`isomeric_smiles`**: Extracted from CIF `_chemical_smiles_isomeric` tag
- **Fallback**: When SMILES is not available in CIF, the formula field provides an alternative identifier

## Implementation Details

### CIF Parser Function
Added `parseCIFData(cifContent)` helper function that extracts:
- Chemical name from `_chemical_name_common` or `_chemical_name_mineral`
- Formula from `_chemical_formula_sum`
- Canonical SMILES from `_chemical_smiles_canonical`
- Isomeric SMILES from `_chemical_smiles_isomeric`

### PubChem API Enhancement
Updated all PubChem property requests to include:
- `MolecularFormula` - For the formula field
- `IUPACName` - For the chemical_name field
- Existing SMILES properties (CanonicalSMILES, IsomericSMILES, etc.)

## Example API Responses

### Compound (Aspirin)
```json
{
  "name": "Aspirin",
  "chemical_name": "2-acetyloxybenzoic acid",
  "formula": "C9H8O4",
  "canonical_smiles": "CC(=O)OC1=CC=CC=C1C(=O)O",
  "isomeric_smiles": "CC(=O)OC1=CC=CC=C1C(=O)O"
}
```

### Mineral with SMILES (Quartz - COD 5000035)
```json
{
  "name": "Quartz",
  "chemical_name": "Silicon oxide",
  "formula": "O2 Si",
  "canonical_smiles": "[Si]1([O])([O])O[Si]([O])([O])O[Si]([O])([O])O[Si]([O])([O])O[Si]([O])([O])O1",
  "isomeric_smiles": null,
  "codid": "5000035"
}
```

### Mineral without SMILES (Niggliite - COD 9008913)
```json
{
  "name": "Niggliite",
  "chemical_name": null,
  "formula": "Pt Sn",
  "canonical_smiles": null,
  "isomeric_smiles": null,
  "codid": "9008913"
}
```

## Usage in Extension

The extension's `content.js` can now access these fields from the search API response:
- Display proper chemical names in tooltips or labels
- Use formulas as fallback identifiers when SMILES aren't available
- Show both common names and IUPAC names for better clarity

## Testing

To test the changes:
1. Restart the search server: `node search-server.js`
2. Query examples:
   - Compound: `http://localhost:8001/search?q=aspirin`
   - Mineral with SMILES: `http://localhost:8001/search?q=quartz`
   - Mineral without SMILES: `http://localhost:8001/search?q=niggliite`

## Future Enhancements

- Add protein/biomolecule chemical names and formulas
- Support additional CIF tags for more comprehensive mineral data
- Cache parsed CIF data to improve performance
- Add validation for extracted SMILES strings
