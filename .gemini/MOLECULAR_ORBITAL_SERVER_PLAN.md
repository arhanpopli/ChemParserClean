# Molecular Orbital Viewer - Server Implementation Plan

## Architecture Overview

```
┌─────────────────┐     ┌─────────────────────────────────────┐
│   Extension     │     │         ChemTex Server              │
│                 │     │                                     │
│  chem:mo=water  │────▶│  1. LCAO Computation (Node.js)     │
│     +homo:      │     │  2. Volume Data Generation          │
│                 │     │  3. Orbital Viewer (hosted page)    │
│  ┌───────────┐  │     │                                     │
│  │  iframe   │◀─┼─────│  Returns: Interactive 3D viewer     │
│  │ (viewer)  │  │     │                                     │
│  └───────────┘  │     └─────────────────────────────────────┘
└─────────────────┘
```

---

## Open Source 3D Viewers for Orbitals

### Option 1: 3Dmol.js (Self-Hosted) ⭐ RECOMMENDED
- **License:** BSD-3-Clause (free for commercial use)
- **Size:** ~500KB
- **Features:** Volumetric rendering, isosurfaces, WebGL
- **Orbital Support:** Yes (can render cube data)
- **URL:** https://github.com/3dmol/3Dmol.js

**Why recommended:** Already supports orbital visualization, lightweight, easy to embed.

### Option 2: NGL Viewer
- **License:** MIT
- **Size:** ~1MB
- **Features:** Molecular structures + volumetric rendering
- **Orbital Support:** Yes (via DX format)
- **URL:** https://github.com/nglviewer/ngl

### Option 3: Speck (Minimalist)
- **License:** MIT
- **Size:** ~100KB
- **Features:** Beautiful molecular rendering
- **Orbital Support:** Limited (would need extension)
- **URL:** https://github.com/wwwtyro/speck

### Option 4: Jmol (Java-based)
- **License:** LGPL
- **Features:** Most complete orbital support
- **Downside:** Java applet, heavy
- **URL:** http://jmol.sourceforge.net/

---

## Recommended Stack

```
Server Stack:
├── Node.js (already have for ChemTex server)
├── Express.js (API endpoints)
├── LCAO Engine (new - computes orbital values)
├── Volume Generator (new - creates grid data)
└── Static Hosting (3Dmol.js viewer page)
```

---

## Implementation Phases

### Phase 1: LCAO Computation Engine

**File:** `server/orbital/lcao-engine.js`

```javascript
// Atomic orbital formulas
const atomicOrbitals = {
  '1s': (x, y, z, zeff) => {
    const r = Math.sqrt(x*x + y*y + z*z);
    return Math.exp(-zeff * r);
  },
  '2s': (x, y, z, zeff) => {
    const r = Math.sqrt(x*x + y*y + z*z);
    return (1 - zeff*r/2) * Math.exp(-zeff * r / 2);
  },
  '2pz': (x, y, z, zeff) => {
    const r = Math.sqrt(x*x + y*y + z*z);
    return z * Math.exp(-zeff * r / 2);
  },
  '2px': (x, y, z, zeff) => {
    const r = Math.sqrt(x*x + y*y + z*z);
    return x * Math.exp(-zeff * r / 2);
  },
  '2py': (x, y, z, zeff) => {
    const r = Math.sqrt(x*x + y*y + z*z);
    return y * Math.exp(-zeff * r / 2);
  },
  // ... more orbitals
};

// Evaluate molecular orbital at a point
function evaluateMO(x, y, z, molecule, orbitalName) {
  const orbital = molecule.orbitals[orbitalName];
  let value = 0;
  
  for (const [atomOrbital, coeff] of Object.entries(orbital.coeffs)) {
    const [atomIdx, orbType] = parseAtomOrbital(atomOrbital);
    const atom = molecule.atoms[atomIdx];
    const dx = x - atom.pos[0];
    const dy = y - atom.pos[1];
    const dz = z - atom.pos[2];
    const zeff = atom.zeff[orbType];
    
    value += coeff * atomicOrbitals[orbType](dx, dy, dz, zeff);
  }
  
  return value;
}
```

### Phase 2: Volume Data Generator

**File:** `server/orbital/volume-generator.js`

```javascript
function generateVolumeData(molecule, orbitalName, gridSize = 50) {
  // Determine bounding box
  const bounds = calculateBounds(molecule);
  const padding = 3.0; // Angstroms
  
  const minX = bounds.minX - padding;
  const maxX = bounds.maxX + padding;
  // ... similar for Y, Z
  
  const dx = (maxX - minX) / gridSize;
  const dy = (maxY - minY) / gridSize;
  const dz = (maxZ - minZ) / gridSize;
  
  // Generate 3D grid of values
  const data = new Float32Array(gridSize * gridSize * gridSize);
  let idx = 0;
  
  for (let iz = 0; iz < gridSize; iz++) {
    for (let iy = 0; iy < gridSize; iy++) {
      for (let ix = 0; ix < gridSize; ix++) {
        const x = minX + ix * dx;
        const y = minY + iy * dy;
        const z = minZ + iz * dz;
        
        data[idx++] = evaluateMO(x, y, z, molecule, orbitalName);
      }
    }
  }
  
  return {
    data: data,
    origin: [minX, minY, minZ],
    dimensions: [gridSize, gridSize, gridSize],
    spacing: [dx, dy, dz]
  };
}
```

### Phase 3: Orbital Viewer Page

**File:** `server/public/orbital-viewer.html`

```html
<!DOCTYPE html>
<html>
<head>
  <title>Orbital Viewer</title>
  <script src="https://3dmol.org/build/3Dmol-min.js"></script>
  <style>
    body { margin: 0; overflow: hidden; background: #000; }
    #viewer { width: 100vw; height: 100vh; }
    .label {
      position: absolute;
      bottom: 10px;
      left: 10px;
      color: white;
      font-family: system-ui;
      font-size: 14px;
      background: rgba(0,0,0,0.5);
      padding: 5px 10px;
      border-radius: 4px;
    }
  </style>
</head>
<body>
  <div id="viewer"></div>
  <div class="label" id="label"></div>
  
  <script>
    // Get parameters from URL
    const params = new URLSearchParams(window.location.search);
    const molecule = params.get('mol') || 'water';
    const orbital = params.get('orbital') || 'homo';
    
    // Initialize viewer
    const viewer = $3Dmol.createViewer('viewer', {
      backgroundColor: 'black'
    });
    
    // Fetch volume data from server
    fetch(`/api/orbital/${molecule}/${orbital}`)
      .then(res => res.json())
      .then(data => {
        // Add positive lobe (blue)
        viewer.addVolumetricData(data.volumeData, "cube", {
          isoval: 0.02,
          color: '#4488ff',
          opacity: 0.8
        });
        
        // Add negative lobe (red)
        viewer.addVolumetricData(data.volumeData, "cube", {
          isoval: -0.02,
          color: '#ff4444',
          opacity: 0.8
        });
        
        // Add molecule structure
        viewer.addModel(data.xyz, "xyz");
        viewer.setStyle({}, {stick: {radius: 0.1}, sphere: {scale: 0.25}});
        
        viewer.zoomTo();
        viewer.render();
        
        document.getElementById('label').textContent = 
          `${data.moleculeName} - ${orbital.toUpperCase()}`;
      });
  </script>
</body>
</html>
```

### Phase 4: API Endpoint

**File:** `server/api/orbital.js`

```javascript
const express = require('express');
const router = express.Router();
const { generateVolumeData } = require('../orbital/volume-generator');
const molecules = require('../orbital/molecules.json');

// GET /api/orbital/water/homo
router.get('/:molecule/:orbital', (req, res) => {
  const { molecule, orbital } = req.params;
  
  const molData = molecules[molecule.toLowerCase()];
  if (!molData) {
    return res.status(404).json({ error: 'Molecule not found' });
  }
  
  const orbitalData = molData.orbitals[orbital.toLowerCase()];
  if (!orbitalData) {
    return res.status(404).json({ error: 'Orbital not found' });
  }
  
  // Generate volume data (cached after first request)
  const volumeData = generateVolumeData(molData, orbital);
  
  // Generate XYZ format for molecule display
  const xyz = generateXYZ(molData);
  
  res.json({
    moleculeName: molData.name,
    orbital: orbital,
    volumeData: volumeData,
    xyz: xyz
  });
});

module.exports = router;
```

### Phase 5: Molecule Database

**File:** `server/orbital/molecules.json`

```json
{
  "water": {
    "name": "Water (H₂O)",
    "atoms": [
      {"element": "O", "pos": [0, 0, 0.117], "zeff": {"2s": 2.275, "2p": 2.275}},
      {"element": "H", "pos": [0.756, 0, -0.469], "zeff": {"1s": 1.24}},
      {"element": "H", "pos": [-0.756, 0, -0.469], "zeff": {"1s": 1.24}}
    ],
    "orbitals": {
      "homo": {
        "name": "HOMO (1b₁)",
        "energy": -0.496,
        "coeffs": {
          "0_2py": 1.0
        }
      },
      "homo-1": {
        "name": "HOMO-1 (3a₁)",
        "energy": -0.558,
        "coeffs": {
          "0_2pz": 0.78,
          "1_1s": 0.21,
          "2_1s": 0.21
        }
      },
      "lumo": {
        "name": "LUMO (4a₁)",
        "energy": 0.21,
        "coeffs": {
          "0_3s": 0.6,
          "1_1s": -0.3,
          "2_1s": -0.3
        }
      }
    }
  },
  "methane": {
    "name": "Methane (CH₄)",
    "atoms": [
      {"element": "C", "pos": [0, 0, 0], "zeff": {"2s": 1.625, "2p": 1.625}},
      {"element": "H", "pos": [0.629, 0.629, 0.629], "zeff": {"1s": 1.24}},
      {"element": "H", "pos": [-0.629, -0.629, 0.629], "zeff": {"1s": 1.24}},
      {"element": "H", "pos": [-0.629, 0.629, -0.629], "zeff": {"1s": 1.24}},
      {"element": "H", "pos": [0.629, -0.629, -0.629], "zeff": {"1s": 1.24}}
    ],
    "orbitals": {
      "homo": {
        "name": "HOMO (1t₂)",
        "coeffs": {
          "0_2px": 0.5,
          "1_1s": 0.25,
          "2_1s": -0.25,
          "3_1s": -0.25,
          "4_1s": 0.25
        }
      }
    }
  }
}
```

---

## Extension Integration

**In `content.js`:**

```javascript
// Handle chem:mo=water+homo:
if (flags.compoundType === 'molecular-orbital' || flags.isMolecularOrbital) {
  const [molecule, orbital] = lookupValue.split('+');
  const viewerUrl = `${CHEMTEX_SERVER}/orbital-viewer.html?mol=${molecule}&orbital=${orbital}`;
  
  // Create iframe
  const iframe = document.createElement('iframe');
  iframe.src = viewerUrl;
  iframe.style.cssText = 'width: 400px; height: 350px; border: none; border-radius: 8px;';
  
  container.appendChild(iframe);
}
```

---

## File Structure

```
chem-extension-server/
├── server/
│   ├── api/
│   │   ├── render.js (existing)
│   │   └── orbital.js (NEW)
│   ├── orbital/
│   │   ├── lcao-engine.js (NEW)
│   │   ├── volume-generator.js (NEW)
│   │   └── molecules.json (NEW)
│   └── public/
│       └── orbital-viewer.html (NEW)
└── package.json
```

---

## Performance Considerations

### Caching Strategy
1. **Volume data caching:** Pre-compute on first request, cache in memory
2. **Static file caching:** Viewer HTML + 3Dmol.js cached by browser
3. **CDN optional:** Volume data for common molecules could be edge-cached

### Grid Size Optimization
| Grid Size | Points | Compute Time | Quality |
|-----------|--------|--------------|---------|
| 30³ | 27,000 | ~5ms | Low |
| 50³ | 125,000 | ~20ms | Good |
| 80³ | 512,000 | ~100ms | High |

50³ recommended - fast enough for real-time, good visual quality.

---

## Syntax for Extension

```
chem:mo=water+homo:        → Water's HOMO orbital
chem:mo=water+lumo:        → Water's LUMO orbital
chem:mo=methane+homo:      → Methane's HOMO
chem:mo=benzene+pi:        → Benzene's π orbital
chem:mo=ethene+pi*:        → Ethene's π* antibonding
```

---

## Initial Molecule Set (50 molecules)

### Priority 1 (Most Common)
- Water, Methane, Ammonia, Carbon dioxide, Hydrogen
- Ethane, Ethene, Ethyne
- Benzene, Toluene
- Formaldehyde, Acetaldehyde

### Priority 2 (Organic)
- Methanol, Ethanol
- Acetic acid
- Acetone
- Chloromethane, Dichloromethane

### Priority 3 (Inorganic)
- Sulfur hexafluoride
- Phosphine
- Hydrogen sulfide
- Borane, Diborane

---

## Timeline Estimate

| Phase | Task | Time |
|-------|------|------|
| 1 | LCAO engine | 2 hours |
| 2 | Volume generator | 1 hour |
| 3 | Viewer page | 1 hour |
| 4 | API endpoint | 1 hour |
| 5 | Extension integration | 1 hour |
| 6 | 10 molecules data | 2 hours |
| **Total** | **MVP** | **~8 hours** |

---

## Summary

**What we're building:**
- Server computes LCAO → volume data (~20ms)
- Self-hosted 3Dmol.js page renders it
- Extension embeds as iframe
- ~50 pre-defined molecules
- Syntax: `chem:mo=water+homo:`

**Why this approach:**
- No client-side JS bloat
- Fast (cached after first compute)
- Interactive 3D (not just images)
- Same pattern as MolView integration

Ready to implement?
