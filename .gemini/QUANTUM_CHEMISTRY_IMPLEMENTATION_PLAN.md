# Quantum Chemistry Visualization - Implementation Plan

## 🎯 Overview

This document outlines the implementation plan for adding **quantum chemistry visualization** to ChemTex, enabling the display of:
- Atomic orbitals (s, p, d, f)
- Molecular orbitals (HOMO, LUMO, bonding/antibonding)
- Electron density surfaces
- Wave function visualizations
- Orbital overlap demonstrations

---

## 📊 Current State Analysis

### What ChemTex Already Supports:
| Category | Tool Used | Data Source |
|----------|-----------|-------------|
| 2D Molecular Structures | SmilesDrawer (custom) | PubChem SMILES |
| 3D Compounds | MolView (embedded) | PubChem 3D conformers |
| Biomolecules/Proteins | Mol* + MolView | RCSB PDB |
| Minerals/Crystals | COD + custom renderer | COD database |
| Static Images | Server-generated SVGs | RCSB images |

### Gap Analysis:
**Quantum chemistry is fundamentally different** - it requires volumetric data (cube files) or mathematical computation, not just coordinate files.

---

## 🔬 Research Findings: Available Tools & Approaches

### Option 1: 3Dmol.js (⭐ RECOMMENDED for Molecular Orbitals)
**Best for:** Molecular orbital visualization of real compounds

| Feature | Support |
|---------|---------|
| Cube file rendering | ✅ Native support |
| Isosurface generation | ✅ Built-in |
| WebGL acceleration | ✅ Yes |
| Browser embedding | ✅ Easy iframe/div |
| License | BSD-3 (open source) |

**Pros:**
- Already focused on chemistry
- Can load `.cube` files directly
- Supports electron density, molecular orbitals
- Works well in iframes (CSP friendly)
- API for programmatic control

**Cons:**
- Requires pre-computed cube files (server needs to provide them)
- No way to compute orbitals on-the-fly in browser

**Resources:**
- https://3dmol.csb.pitt.edu/
- GitHub: 3dmol/3Dmol.js

---

### Option 2: Mol* (Already in ChemTex!)
**Best for:** Electron density from experimental data (cryo-EM, X-ray)

Mol* already has a **Density panel** and supports volumetric data visualization. The example on molstar.org shows "Alpha orbitals and density of Atorvastatin".

| Feature | Support |
|---------|---------|
| Electron density | ✅ Yes |
| Volumetric data | ✅ Via extensions |
| Already integrated | ✅ In ChemTex |
| Quantum orbitals | ⚠️ Limited (needs cube files) |

**Action:** Check if your existing Mol* wrapper can load `.map` or `.cube` files.

---

### Option 3: Mathematical Orbital Generation (⭐ RECOMMENDED for Atomic Orbitals)
**Best for:** Pure atomic orbital visualization (hydrogen-like orbitals)

This approach generates orbitals mathematically using the hydrogen wave function formulas. No external data needed!

| Feature | Possibility |
|---------|-------------|
| s, p, d, f orbitals | ✅ Can compute |
| Any element | ⚠️ Hydrogen-like only (approximation) |
| Real-time rendering | ✅ WebGL volume rendering |
| No server needed | ✅ All client-side |

**Existing implementations to reference:**
1. **Falstad's Hydrogen Atom Viewer**: https://falstad.com/qmatom/
   - Open source, JavaScript
   - Can display n=1 to n=7 orbitals
   - Supports combinations of orbitals
   
2. **Lisyarus WebGL Orbitals**: https://lisyarus.github.io/webgl/orbitals.html
   - Pure WebGL, very visual
   
3. **The Orbitron** (static images): https://winter.group.shef.ac.uk/orbitron/
   - Pre-rendered images of all orbitals
   - Could be used as a fallback

**Mathematical basis:**
```
Ψ(n,l,m) = R(n,l) × Y(l,m)
```
Where R is the radial function and Y is the spherical harmonic.

---

### Option 4: Pre-computed Orbital Database (Hybrid Approach)
**Best for:** Fast loading, consistent quality

Create/curate a database of:
- Pre-rendered orbital images (PNG/SVG)
- Pre-computed cube files for common orbitals
- Hosted on your server

| Orbital Type | Elements | Files Needed |
|--------------|----------|--------------|
| 1s | H, He, Li, etc. | ~10 cube files |
| 2s, 2p | Li→Ne | ~40 cube files |
| 3s, 3p, 3d | Na→Ar | ~90 cube files |
| 4s, 4p, 4d, 4f | K→Kr | ~160 cube files |

---

## 🏗️ Proposed Architecture

### Tier 1: Atomic Orbitals (Individual Atoms)
```
chem:orbital=carbon+2p:     → Shows carbon's 2p orbital
chem:orbital=hydrogen+1s:   → Shows hydrogen's 1s orbital
chem:orbital=oxygen+2px:    → Shows specific 2px orbital
chem:orbital=iron+3d:       → Shows iron's 3d orbitals
```

**Implementation:** Client-side mathematical generation OR pre-computed images

### Tier 2: Interactive Orbital Explorer
```
chem:quantexplore=hydrogen:          → Interactive orbital selector
chem:quantexplore=carbon+config:     → Show electron configuration
```

**Implementation:** Embedded viewer with controls (like Falstad)

### Tier 3: Molecular Orbitals (Compounds)
```
chem:mo=H2+bond:           → Shows σ bonding orbital
chem:mo=O2+homo:           → Shows HOMO of oxygen
chem:mo=benzene+pi:        → Shows π system
chem:morbital=caffeine:    → All MOs for a compound (needs server)
```

**Implementation:** 3Dmol.js with pre-computed cube files from server

### Tier 4: Advanced Demonstrations
```
chem:overlap=H+H:           → Animated orbital overlap → H₂ formation
chem:hybrid=sp3:            → Shows sp³ hybridization
chem:resonance=benzene:     → Animated π electron delocalization
```

**Implementation:** Custom animations using Three.js or similar

---

## 🎨 Proposed Syntax

### Basic Atomic Orbital Syntax
```
chem:quan=<element>+<orbital>:
chem:quan=carbon+2p:
chem:quan=nitrogen+2s:
chem:quan=iron+3dxy:
chem:quan=hydrogen+3dz2:
```

### Quantum Numbers Syntax
```
chem:quann=1,0,0:           → n=1, l=0, m=0 (1s)
chem:quann=2,1,0:           → n=2, l=1, m=0 (2pz)
chem:quann=3,2,1:           → n=3, l=2, m=1 (3d)
```

### Named Orbital Syntax
```
chem:Hydrogen 1s Orbitalquan=H+1s:
chem:Carbon p Orbitalsquan=C+2p:
```

### Flags for Quantum Visualizations
```
+phase     → Show positive/negative phases with colors
+node      → Highlight nodal planes
+radial    → Show radial distribution function
+angular   → Show angular distribution
+wave      → Animate wave function
+prob      → Show probability density (|Ψ|²)
+3d        → Force 3D view (vs 2D cross-section)
+slice=xy  → Show XY plane slice
```

---

## 📦 Data Sources

### For Atomic Orbitals:
1. **Generate mathematically** (recommended)
   - Hydrogen wave function formulas are well-known
   - Works for all elements (hydrogen-like approximation)
   - No external dependencies

2. **Pre-compute and host**
   - Use PySCF, Gaussian, or ORCA to generate cube files
   - Host on your server
   - Cache aggressively

### For Molecular Orbitals:
1. **PubChem** - Does NOT provide orbital data ❌
2. **Quantum chemistry APIs:**
   - Molecular Data API (RapidAPI) - Limited
   - ChemSpider - No MO data
   
3. **Option: Server-side computation**
   - Install PySCF or Psi4 on server
   - Compute MOs on-demand
   - Cache results
   - ⚠️ CPU intensive!

4. **Option: Pre-computed for common molecules**
   - Benzene, H₂, O₂, CO₂, H₂O, etc.
   - Store cube files
   - Cover educational use cases

---

## 🛠️ Implementation Phases

### Phase 1: Foundational Research & Prototype (2-3 weeks)
- [ ] Fork/embed Falstad's orbital viewer as proof-of-concept
- [ ] Test 3Dmol.js cube file rendering
- [ ] Design final syntax for `chem:quan=` tags
- [ ] Create wrapper page for quantum viewer (like your Mol* wrapper)

### Phase 2: Atomic Orbital Viewer (3-4 weeks)
- [ ] Build standalone orbital viewer page
  - Mathematical generation of s, p, d, f orbitals
  - WebGL rendering with Three.js or Babylon.js
  - Interactive rotation, zoom
  - Orbital selector (n, l, m controls)
- [ ] Integrate into ChemTex content.js
- [ ] Add `chem:quan=` tag parsing
- [ ] Server endpoint for orbital metadata

### Phase 3: Pre-computed Orbital Gallery (2 weeks)
- [ ] Generate high-quality orbital images for documentation
- [ ] Create cube files for common orbitals
- [ ] Host on server with CDN caching
- [ ] Fallback to static images for unsupported browsers

### Phase 4: Molecular Orbital Support (4-6 weeks)
- [ ] 3Dmol.js integration for cube file rendering
- [ ] Server-side MO computation setup (PySCF)
- [ ] Pre-compute MOs for 20-50 common molecules
- [ ] API endpoint: `/api/orbital/{compound}`
- [ ] Caching layer for computed orbitals

### Phase 5: Educational Demonstrations (3-4 weeks)
- [ ] Orbital overlap animations
- [ ] Hybridization visualizations
- [ ] Electron configuration diagrams
- [ ] Interactive tutorials

### Phase 6: ChatGPT Integration (1-2 weeks)
- [ ] Update ChatGPT prompt with quantum syntax
- [ ] Document best practices for AI-generated quantum tags
- [ ] Test with various chemistry queries

---

## 🖥️ Technical Details

### Client-Side (Extension)
```javascript
// New orbital detection in content.js
const quantumPattern = /chem:(?:(\w+))?quan(?:n)?=([^:]+):/g;

// Example match: "chem:Hydrogen 1s Orbitalquan=H+1s:"
// Group 1: "Hydrogen 1s Orbital" (optional name)
// Group 2: "H+1s" (element + orbital)
```

### Server-Side (New Endpoints)
```
GET /api/orbital/atomic/{element}/{orbital}
    → Returns orbital visualization config/data

GET /api/orbital/molecular/{compound}
    → Returns cube file or 3Dmol.js viewer URL

GET /api/orbital/cube/{id}
    → Returns pre-computed cube file
```

### New Files Needed
```
ChemistryLaTeX-source/
├── orbital-viewer.html      # Embedded quantum viewer
├── orbital.js               # Orbital computation logic
└── assets/
    └── orbitals/            # Pre-computed cube files

chem-extension-server/
└── server/api/
    └── orbital.js           # Orbital API endpoints
```

---

## 🎯 MVP Definition

For a **Minimum Viable Product**, focus on:

1. **Atomic orbitals only** (no molecular orbitals initially)
2. **Embed Falstad viewer** as first pass (it's open source)
3. **Simple syntax**: `chem:quan=H+1s:`, `chem:quan=C+2p:`
4. **Hydrogen-like approximation** for all elements
5. **3D interactive viewer** with rotation

This can be done in **2-3 weeks** and provides immediate value for educational content.

---

## 📚 Resources & References

### JavaScript Libraries
- **3Dmol.js**: https://3dmol.csb.pitt.edu/
- **Three.js**: https://threejs.org/
- **Babylon.js**: https://babylonjs.com/

### Existing Orbital Viewers
- **Falstad Hydrogen Orbital Viewer**: https://falstad.com/qmatom/
- **The Orbitron (Sheffield)**: https://winter.group.shef.ac.uk/orbitron/
- **Lisyarus WebGL Orbitals**: https://lisyarus.github.io/webgl/orbitals.html

### Quantum Chemistry Computation
- **PySCF** (Python): https://pyscf.org/
- **Psi4** (Python): https://psicode.org/
- **ORCA** (Free for academics): https://orcaforum.kofo.mpg.de/

### Cube File Format
- Gaussian cube format specification
- Used by all major quantum chemistry packages

---

## ❓ Open Questions

1. **Server compute budget** - Can your server handle on-demand MO calculations?
2. **Scope prioritization** - Do you want atomic orbitals first, or molecular orbitals?
3. **ChatGPT integration** - Should AI be able to request specific orbitals?
4. **Educational focus** - Do you want animated demonstrations (like orbital mixing)?
5. **Fallback strategy** - Static images vs. error message for unsupported cases?

---

## 📊 Effort Estimate Summary

| Phase | Effort | Priority |
|-------|--------|----------|
| Atomic Orbital Viewer (MVP) | 2-3 weeks | 🔴 High |
| Pre-computed Gallery | 1-2 weeks | 🟡 Medium |
| Molecular Orbitals | 4-6 weeks | 🟡 Medium |
| Animations/Demos | 3-4 weeks | 🟢 Low |
| ChatGPT Integration | 1 week | 🟢 Low |

**Total for full implementation:** ~12-16 weeks

**MVP (atomic orbitals only):** ~3-4 weeks

---

## 🚀 Recommended Next Steps

1. **Try embedding Falstad viewer** in an iframe as a quick prototype
2. **Test 3Dmol.js** with a sample cube file
3. **Decide on scope** - atomic vs molecular orbitals first
4. **Design the final syntax** you want for `chem:quan=` tags
5. **Create a demo page** outside the extension first

---

*Document created: January 14, 2026*
*Last updated: January 14, 2026*
