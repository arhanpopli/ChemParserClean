# QCArchive Integration - Detailed Implementation Plan

## Project Overview

**Goal**: Integrate QCArchive as a real-time quantum chemistry data source for molecular orbital visualization, enabling visualization of orbitals for any molecule in their database.

**Data Sources Stack (After Integration)**:
| Data Type | Source | API |
|-----------|--------|-----|
| 2D Molecules | PubChem | REST |
| 3D Structures | RCSB PDB | REST |
| Minerals | COD | REST |
| **Molecular Orbitals** | **QCArchive** | **REST** |

---

# Phase 1: QCArchive API Research & Connection
**Estimated Time: 3-4 hours**
**Status: 🟡 IN PROGRESS - Critical Limitation Found**

## 1.1 API Exploration
- [x] Test QCArchive demo server connectivity - ✅ Working
- [x] Understand authentication requirements - ✅ No auth for read-only
- [ ] Document rate limits and usage policies
- [x] Identify available datasets - ✅ Many datasets available

## 1.2 Molecule Query Endpoint
- [x] Test `/api/v1/molecules/query` endpoint - ✅ Working
- [x] Query by molecular formula (e.g., "H2O") - ✅ Returns molecule IDs
- [ ] Query by SMILES string - TODO
- [ ] Query by identifiers - TODO
- [x] Document response format - ✅ Returns array of molecule IDs

## 1.3 Singlepoint Records Endpoint
- [x] Test `/api/v1/records/singlepoint/query` endpoint - ✅ Working
- [x] Understand relationship between molecules and records - ✅ Documented
- [x] Identify which records contain wavefunction data - ⚠️ **CRITICAL ISSUE**
- [x] Document computation methods available - ✅ b3lyp, psi4, etc.

## 1.4 Wavefunction Data Structure
- [x] Retrieve a sample wavefunction - ⚠️ **RETURNED NULL**
- [ ] Document QCSchema format for wavefunctions
- [ ] Identify location of orbital coefficients
- [ ] Identify location of orbital energies

---

## ⚠️ CRITICAL FINDING

### Problem: Wavefunction Data Not Stored

**Most QCArchive records do NOT store wavefunction (orbital) data.**

When I queried for a water molecule record:
```
GET /api/v1/records/singlepoint/144541839/wavefunction
Response: null
```

**Reason**: Wavefunction data is large (MBs per molecule) and is only stored when:
1. Explicitly requested when submitting the calculation
2. The dataset requires it for ML training

### What IS Available in QCArchive:
- ✅ Molecular geometries
- ✅ Energies (total, SCF, XC)
- ✅ Dipole moments
- ✅ Gradients
- ✅ Hessians
- ❌ Orbital coefficients (not by default)
- ❌ Orbital energies (not directly)

### Implications:

**QCArchive cannot be used as a real-time orbital coefficient source** without:
1. Finding specific datasets that DO store wavefunctions, OR
2. Running our own calculations and submitting them, OR
3. Using a different approach entirely

---

## Alternative Discovery: Orbital Energies Available

While full orbital coefficients are rare, I noticed `calcinfo_nmo` (number of MOs) and orbital energy references exist. Some records MAY have partial orbital info.

## Recommended Pivot Options:

### Option A: Search for Wavefunction-Enabled Datasets
Some ML-focused datasets might store wavefunctions. I can search for these.

### Option B: Use NIST CCCBDB Instead
NIST provides `.wfn` files that DO contain orbital coefficients for ~2000 molecules.

### Option C: Run Our Own Calculations
Set up a Psi4 server that computes ON DEMAND and caches results.

### Option D: Hybrid Approach
1. Use QCArchive for geometry + energy
2. Use our LCAO approximation for orbital visualization
3. The orbitals won't be "exact" but will be educationally correct

---

## Next Steps - Decision Required

**User input needed**: Given this limitation, which path should we take?

1. **Continue searching QCArchive** for wavefunction-enabled datasets
2. **Pivot to NIST CCCBDB** for orbital data
3. **Build Psi4 calculation server** (most powerful but complex)
4. **Use hybrid approach** (QCArchive geometry + LCAO orbitals)

---

# Phase 2: QCArchive Client Implementation
**Estimated Time: 4-5 hours**

## 2.1 Create QCArchive Client Module
**File**: `server/orbital/qcarchive-client.js`

- [ ] Create base HTTP client with proper headers
- [ ] Implement retry logic for failed requests
- [ ] Add request timeout handling (QCArchive can be slow)
- [ ] Add error handling for API errors

## 2.2 Molecule Lookup Functions
- [ ] `findMoleculeBySmiles(smiles)` - Primary lookup method
- [ ] `findMoleculeByFormula(formula)` - Backup lookup
- [ ] `findMoleculeByName(name)` - Via PubChem → SMILES → QCArchive
- [ ] Handle "molecule not found" gracefully

## 2.3 Record Retrieval Functions
- [ ] `getSinglepointRecords(moleculeId)` - Get calculation records
- [ ] `getBestRecord(records)` - Select highest quality calculation
- [ ] `getWavefunction(recordId)` - Fetch wavefunction data
- [ ] Implement quality ranking:
  1. ccsd(t) > mp2 > dft > hf
  2. Larger basis sets preferred

## 2.4 Caching Layer
- [ ] In-memory cache for current session
- [ ] Consider file-based cache for persistence
- [ ] Cache structure: `{ smiles: { moleculeId, recordId, wavefunction } }`
- [ ] Cache expiry (24 hours for records, indefinite for orbitals)

### Phase 2 Deliverables:
- [ ] Complete `qcarchive-client.js` module
- [ ] All lookup functions working
- [ ] Caching implemented
- [ ] Unit tests passing

---

# Phase 3: Basis Set Implementation
**Estimated Time: 6-8 hours**

## 3.1 Understand Gaussian Basis Sets
QCArchive uses Gaussian-type orbitals (GTOs), not Slater-type (STOs).

**Gaussian orbital formula**:
```
χ = N * x^l * y^m * z^n * exp(-α * r²)
```

vs our current Slater:
```
ψ = N * r^(n-1) * exp(-ζ * r)
```

- [ ] Document basis set format in QCArchive
- [ ] Understand contracted vs primitive Gaussians
- [ ] Map basis function indices to atom positions

## 3.2 Implement Gaussian Basis Functions
**File**: `server/orbital/gaussian-basis.js`

- [ ] Implement primitive Gaussian evaluation
- [ ] Implement contracted Gaussian evaluation
- [ ] Support s-type functions (l=m=n=0)
- [ ] Support p-type functions (l+m+n=1)
- [ ] Support d-type functions (l+m+n=2)
- [ ] Support f-type functions (l+m+n=3) [if needed]

## 3.3 Basis Set Library
- [ ] Parse basis set definitions from QCArchive
- [ ] Implement common basis sets:
  - [ ] STO-3G (minimal, fast)
  - [ ] 6-31G* (common for organics)
  - [ ] cc-pVDZ (correlation consistent)
  - [ ] def2-SVP (Ahlrichs basis)
- [ ] Store basis set parameters locally for speed

## 3.4 Shell and Contraction Handling
- [ ] Map QCArchive shell format to our functions
- [ ] Handle normalization constants
- [ ] Handle spherical vs Cartesian d/f functions

### Phase 3 Deliverables:
- [ ] Complete `gaussian-basis.js` module
- [ ] All supported basis sets working
- [ ] Verification against known values
- [ ] Performance benchmarks

---

# Phase 4: Molecular Orbital Reconstruction
**Estimated Time: 5-6 hours**

## 4.1 Orbital Coefficient Processing
**File**: `server/orbital/mo-builder.js`

- [ ] Parse orbital coefficient matrix from QCArchive
- [ ] Identify HOMO, LUMO, and other orbitals by energy
- [ ] Handle alpha/beta spin orbitals (if different)
- [ ] Map coefficients to basis functions

## 4.2 MO Evaluation Function
```javascript
function evaluateMO(x, y, z, orbitalIndex, coefficients, basisSet, atoms) {
  let value = 0;
  for (let i = 0; i < basisSet.length; i++) {
    value += coefficients[orbitalIndex][i] * evaluateBasisFunction(x, y, z, basisSet[i], atoms);
  }
  return value;
}
```

- [ ] Implement MO evaluation at a point
- [ ] Optimize for performance (this runs 64,000+ times per orbital)
- [ ] Add SIMD/vectorization if possible

## 4.3 Volume Data Generation
- [ ] Adapt existing `generateVolumeData` for Gaussian basis
- [ ] Determine appropriate grid spacing for GTO resolution
- [ ] Implement adaptive grid (finer near atoms)
- [ ] Profile and optimize grid generation

## 4.4 Cube File Generation
- [ ] Update `generateCubeFile` for GTO-based data
- [ ] Ensure proper unit conversion (Bohr/Angstrom)
- [ ] Verify cube file validity with 3Dmol.js

### Phase 4 Deliverables:
- [ ] Complete `mo-builder.js` module
- [ ] MO evaluation working
- [ ] Cube file generation from QCArchive data
- [ ] Performance within acceptable limits (<5s per orbital)

---

# Phase 5: API Endpoint Integration
**Estimated Time: 3-4 hours**

## 5.1 Update Orbital API
**File**: `server/api/orbital.js`

Current flow:
```
Request → Local molecules.json → LCAO render → Response
```

New flow:
```
Request → Local cache → QCArchive lookup → GTO render → Cache → Response
```

- [ ] Add QCArchive lookup path
- [ ] Implement fallback to local LCAO for unsupported molecules
- [ ] Add query parameter for data source selection
- [ ] Handle loading states for slow QCArchive queries

## 5.2 Request Handling
- [ ] Accept molecule by: name, SMILES, formula, CID
- [ ] Accept orbital by: "homo", "lumo", "homo-1", or index
- [ ] Return appropriate error for missing data
- [ ] Add metadata about data source in response

## 5.3 Response Format
```json
{
  "success": true,
  "source": "qcarchive",  // or "local"
  "qualityLevel": "b3lyp/6-31g*",
  "molecule": { "name": "Water", "formula": "H2O" },
  "orbital": { "name": "HOMO", "energy": "-0.496 Ha" },
  "cubeFile": "...",
  "xyz": "..."
}
```

- [ ] Add source metadata
- [ ] Add computation method info
- [ ] Add timing/cache info for debugging

## 5.4 Rate Limiting & Queueing
- [ ] Implement request queue for QCArchive calls
- [ ] Add rate limiting to avoid API abuse
- [ ] Add timeout for long-running queries (30s max)
- [ ] Return cached partial results while fetching new data

### Phase 5 Deliverables:
- [ ] Updated `/api/orbital` endpoint
- [ ] Full QCArchive integration working
- [ ] Fallback to local data working
- [ ] Error handling complete

---

# Phase 6: Viewer Updates
**Estimated Time: 2-3 hours**

## 6.1 Loading States
**File**: `server/orbital-viewer.html`

- [ ] Add better loading indicator (QCArchive is slower)
- [ ] Add progress messages ("Fetching from QCArchive...", "Generating orbital...")
- [ ] Add timeout handling with user-friendly message

## 6.2 Quality Indicators
- [ ] Show computation method used
- [ ] Show data source (QCArchive vs local)
- [ ] Add "quality badge" (approximate vs high-quality)

## 6.3 Error Handling
- [ ] Handle "molecule not found in QCArchive"
- [ ] Handle "no orbital data available"
- [ ] Suggest alternatives or show atomic orbitals instead

## 6.4 Performance Optimization
- [ ] Lazy load 3Dmol.js
- [ ] Add WebWorker for cube parsing (if needed)
- [ ] Optimize initial render time

### Phase 6 Deliverables:
- [ ] Updated viewer with loading states
- [ ] Quality indicators visible
- [ ] All error cases handled gracefully

---

# Phase 7: Extension Integration
**Estimated Time: 2-3 hours**

## 7.1 New Syntax Support
**File**: `ChemistryLaTeX-source/content.js`

Current:
- `chem:quan=H+1s:` → Atomic orbital (Falstad)
- `chem:hybrid=sp3:` → Hybrid orbital (Falstad)

New:
- `chem:mo=water+homo:` → Molecular orbital (QCArchive/local)
- `chem:mo=CCO+homo:` → MO by SMILES
- `chem:mo=benzene+pi:` → Named orbital

- [ ] Parse `chem:mo=` syntax
- [ ] Extract molecule identifier and orbital name
- [ ] Generate iframe URL for orbital viewer
- [ ] Handle loading state in extension

## 7.2 Settings
- [ ] Add "Molecular Orbital Quality" setting (fast/accurate)
- [ ] Add "Use QCArchive" toggle (for fallback control)
- [ ] Show warning for slow connections

## 7.3 ChatGPT Prompt Updates
- [ ] Document `chem:mo=` syntax for AI
- [ ] Add examples of common molecular orbitals
- [ ] Suggest when to use MO vs atomic orbital syntax

### Phase 7 Deliverables:
- [ ] `chem:mo=` syntax working in extension
- [ ] Settings UI updated
- [ ] ChatGPT prompt updated

---

# Phase 8: Testing & Documentation
**Estimated Time: 2-3 hours**

## 8.1 Test Cases
**File**: `server/tests/orbital-qcarchive.test.js`

- [ ] Test QCArchive connection
- [ ] Test molecule lookup (water, benzene, ethanol)
- [ ] Test orbital retrieval
- [ ] Test cube file generation
- [ ] Test fallback to local data
- [ ] Test error handling

## 8.2 Performance Testing
- [ ] Benchmark cube generation time
- [ ] Measure QCArchive response times
- [ ] Identify bottlenecks
- [ ] Optimize if needed

## 8.3 Documentation
- [ ] Update README with MO feature
- [ ] Document API endpoints
- [ ] Document supported molecules/orbitals
- [ ] Add troubleshooting guide

### Phase 8 Deliverables:
- [ ] All tests passing
- [ ] Performance benchmarks documented
- [ ] User documentation complete

---

# Phase 9: Deployment & Monitoring
**Estimated Time: 1-2 hours**

## 9.1 Vercel Configuration
- [ ] Update `vercel.json` for new files
- [ ] Ensure orbital folder is included
- [ ] Set appropriate memory limits (may need more for GTO)
- [ ] Configure increased timeout if needed

## 9.2 Monitoring
- [ ] Add logging for QCArchive calls
- [ ] Track cache hit/miss rates
- [ ] Monitor QCArchive API health
- [ ] Alert on high error rates

## 9.3 Production Deployment
- [ ] Deploy to staging
- [ ] Test all features
- [ ] Deploy to production
- [ ] Verify with real users

### Phase 9 Deliverables:
- [ ] Production deployment complete
- [ ] Monitoring in place
- [ ] Feature live for users

---

# Summary Timeline

| Phase | Description | Hours | Status |
|-------|-------------|-------|--------|
| 1 | API Research & Connection | 3-4 | ⏳ Not Started |
| 2 | QCArchive Client | 4-5 | ⏳ Not Started |
| 3 | Basis Set Implementation | 6-8 | ⏳ Not Started |
| 4 | MO Reconstruction | 5-6 | ⏳ Not Started |
| 5 | API Endpoint Integration | 3-4 | ⏳ Not Started |
| 6 | Viewer Updates | 2-3 | ⏳ Not Started |
| 7 | Extension Integration | 2-3 | ⏳ Not Started |
| 8 | Testing & Documentation | 2-3 | ⏳ Not Started |
| 9 | Deployment & Monitoring | 1-2 | ⏳ Not Started |
| **Total** | | **28-38 hours** | |

---

# File Structure (After Implementation)

```
server/
├── api/
│   ├── orbital.js              # Main API (updated)
│   └── ...
├── orbital/
│   ├── lcao-engine.js          # Existing Slater-type orbitals
│   ├── molecules.json          # Local fallback data
│   ├── qcarchive-client.js     # NEW: QCArchive API client
│   ├── gaussian-basis.js       # NEW: Gaussian basis functions
│   ├── mo-builder.js           # NEW: MO from coefficients
│   ├── basis-sets/             # NEW: Basis set parameters
│   │   ├── 6-31g.json
│   │   ├── cc-pvdz.json
│   │   └── sto-3g.json
│   └── cache/                  # NEW: Cached QCArchive data
│       └── .gitkeep
├── public/
│   └── orbital-viewer.html     # Viewer (updated)
└── tests/
    └── orbital-qcarchive.test.js  # NEW: Tests
```

---

# Starting Point

**Ready to begin Phase 1?**

I'll start by testing the QCArchive API connection and understanding their data format.
