# Quantum Chemistry Implementation Plan v2

## Current Status
- Pattern detection works ✅ (chem:quan=H+1s: is detected)
- Placeholder creation works ✅ (img element created)
- Lazy loading NOT working ❌ (error: "Cannot read properties of undefined (reading 'substring')")
- 3D iframe NOT rendering ❌

## Architecture (Mirroring Biomolecule Flow)

### How Biomolecules Work (For Reference)
```
1. chem:biomol=hemoglobin: detected
2. Server lookup → returns pdbid (e.g., "1HHO")
3. 2D Image: https://cdn.rcsb.org/images/structures/1hho_assembly-1.jpeg
4. 3D Viewer (on click): MolView or Mol* iframe with pdbid
```

### Proposed Quantum Orbital Flow
```
1. chem:quan=C+2p: detected
2. Parse: element=C, n=2, orbital=p
3. 2D Image: Winter Group or generated orbital image URL
4. 3D Viewer (on click): Falstad iframe with parameters
```

## Static Orbital Images Sources

### 1. Winter Group (Sheffield University)
- URL: https://winter.group.shef.ac.uk/orbitron/
- Examples:
  - 1s: https://winter.group.shef.ac.uk/orbitron/AOs/1s/wave_function/1s.html (has images)
  - 2p: https://winter.group.shef.ac.uk/orbitron/AOs/2p/wave_function/2p.html
- **Problem**: No direct image API, would need to scrape or host locally

### 2. Pre-rendered Images (Best Option)
- Create a small set of orbital images (1s, 2s, 2p, 3s, 3p, 3d, 4s, 4p, 4d, 4f)
- Host on your server or CDN
- ~20-30 images total covers 99% of chemistry use cases

### 3. Wikipedia Commons
- Has some orbital images but inconsistent quality

## 3D Viewers

### 1. Falstad Quantum Applet ✅ (Chosen)
- URL: https://falstad.com/qmatom/qmatom.html
- Deep linking: ?vc=0&n=0&l=0&m=0
- Supports: Atomic orbitals, hybrids, superpositions
- Embeddable: Yes (iframe)

### 2. Orbital Viewer (Alternative)
- URL: https://www.falstad.com/qm1d/
- 1D wave function visualizer

## Molecular Orbitals (Future - Complex)

### The Challenge
- Atomic orbitals are element-independent (hydrogen-like)
- Molecular orbitals depend on:
  - Molecular geometry
  - Electron configuration
  - Computation (Hartree-Fock, DFT)

### Solutions
1. **Pre-computed common molecules** (bonus, ~100 molecules)
2. **Server-side computation** (Psi4/PySCF) - expensive to host
3. **External API** (none exist publicly)
4. **Cube file hosting** (very large files, not practical)

### Recommendation: Focus on Atomic Orbitals First
- 99% of educational use cases
- ChatGPT can explain "carbon uses sp3 hybridization" with chem:hybrid=sp3:
- Molecular orbitals = Phase 3 or later

## Implementation Phases

### Phase 1: Fix Current Bug + Basic Display (TODAY)
1. Fix the "substring" error in lazy loading
2. Get Falstad iframe displaying for quantum orbitals
3. Test with: chem:quan=H+1s:, chem:quan=C+2p:, chem:hybrid=sp3:

### Phase 2: 2D Images + 3D Toggle (NEXT)
1. Create/source static orbital images (2D)
2. Add to server or CDN
3. Show 2D by default, 3D on click (like biomolecules)
4. Add popup settings for quantum display

### Phase 3: Enhanced Features
1. Orbital superposition visualization
2. Phase coloring options
3. Slice views (xy, xz, yz planes)
4. Animation support

### Phase 4: Molecular Orbitals (Future)
1. Pre-computed common molecules
2. Server-side computation (if hosting Psi4)
3. Cube file support

## Syntax Reference

### Atomic Orbitals
```
chem:quan=H+1s:       → Hydrogen 1s (default display)
chem:quan=C+2p:       → Carbon 2p orbital
chem:quan=C+2px:      → Carbon 2px specifically
chem:quan=Fe+3dxy:    → Iron 3dxy orbital
chem:quan=U+5f:       → Uranium 5f orbital
```

### Hybrid Orbitals
```
chem:hybrid=sp:       → Linear (180°)
chem:hybrid=sp2:      → Trigonal planar (120°)
chem:hybrid=sp3:      → Tetrahedral (109.5°)
chem:hybrid=sp3d:     → Trigonal bipyramidal
chem:hybrid=sp3d2:    → Octahedral
```

### Future: Molecular Orbitals
```
chem:mo=water+homo:   → Water's HOMO
chem:mo=benzene+pi:   → Benzene's π system
chem:mo=ethene+pi*:   → Ethene's antibonding π*
```

## Popup Settings (Phase 2)

### Quantum Chemistry Section
```
☑ Enable quantum orbital visualization
  ○ Real orbitals (chemistry style)
  ○ Complex orbitals (physics style)
  
☑ Show phase colors
  Color scheme: [Red/Blue ▼]
  
☑ Auto-show 3D viewer for orbitals
  Viewer size: [Normal ▼]
```

## ChatGPT Prompt Addition (Phase 2)

```
For atomic orbitals: chem:quan=ELEMENT+ORBITAL:
  Examples: chem:quan=H+1s: chem:quan=C+2p: chem:quan=Fe+3dxy:

For hybrid orbitals: chem:hybrid=TYPE:
  Examples: chem:hybrid=sp3: chem:hybrid=sp2: chem:hybrid=sp:

When explaining electron configuration or bonding, use these tags to visualize orbitals.
```

## Immediate Next Steps

1. **Debug the substring error** - check what's undefined in lazy loader
2. **Test iframe directly** - try the Falstad URL in browser first
3. **Verify data flow** - add console.log in DEV build (not production)
4. **Get basic 3D working before adding 2D images**

## Notes

### About Cube Files
- Cube files contain volumetric electron density data
- Very large (10-100MB per molecule)
- Not practical for web delivery
- Better to use pre-rendered images or real-time calculation

### About Python Computation
- Libraries like Psi4, PySCF can compute orbitals
- Would require server-side infrastructure
- Good for Phase 4 if you want full molecular orbital support
- For now, atomic orbitals from Falstad are sufficient

### About ChatGPT Hallucination
- You're right that ChatGPT shouldn't guess SMILES for complex molecules
- For orbitals, the syntax is simpler (just orbital names)
- ChatGPT can reliably say "carbon has 2p orbitals" even if it can't draw complex SMILES
