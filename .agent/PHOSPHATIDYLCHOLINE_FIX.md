# 3D Viewer Fix for Phosphatidylcholine

## The Issue
The user reported that "phosphatidylcholine" was not showing a 3D model, even though "Phosphatidylcholine_1" works on the server and in PubChem.

## The Root Cause
1. "phosphatidylcholine" (generic name) resolves to a CID that often lacks 3D data or fails lookup.
2. "Phosphatidylcholine_1" (CID 10425706) is the representative compound that HAS 3D data.
3. The extension was falling back to a 2D structure (flat 3D) which the user explicitly disliked.

## The Fix

### 1. Smart CID Lookup (`content.js`)
Updated `getPubChemCID` to automatically try appending `_1` if the direct name lookup fails.
- Input: "phosphatidylcholine"
- Try 1: "phosphatidylcholine" -> Fails or returns CID without 3D
- Try 2: "phosphatidylcholine_1" -> **Success** (Returns CID 10425706)

### 2. Strict 3D Fetching (`3dmol-viewer.js`)
Updated `loadMolecule` to:
- **Remove** the 2D structure fallback (no more flat molecules).
- **Keep** the 3D Conformer fetch (Priority 1).
- **Keep** the Computed 3D fetch (Priority 2).
- Throw a clear error if no 3D data is found.

## Verification
- **Server behavior**: Uses MolView with "Phosphatidylcholine_1".
- **Extension behavior**: Now resolves "phosphatidylcholine" to "Phosphatidylcholine_1" CID and fetches the real 3D SDF.

## Files Modified
- `chem-extension/content.js`
- `chem-extension/3dmol-viewer.js`
