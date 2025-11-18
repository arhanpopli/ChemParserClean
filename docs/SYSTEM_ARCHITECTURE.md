# ChemParser System Architecture

## System Overview

```
┌─────────────────────────────────────────────────────────────────┐
│                    ChemParser System                              │
│                                                                   │
│  ┌─────────────┐  ┌─────────────┐  ┌─────────────┐  ┌─────────┐ │
│  │ Chrome      │  │   Master    │  │   Test      │  │  Cache  │ │
│  │ Extension   │  │  Launcher   │  │   Pages     │  │ Storage │ │
│  └──────┬──────┘  └──────┬──────┘  └──────┬──────┘  └────┬────┘ │
│         │                │                │               │       │
│         └────────────────┴────────────────┴───────────────┘       │
│                              ▼                                    │
│         ┌────────────────────────────────────────────────┐        │
│         │         HTTP API Layer (localhost)             │        │
│         └────────────────────────────────────────────────┘        │
│                              ▼                                    │
│  ┌─────────────┬─────────────┬──────────────┬──────────────┐     │
│  │             │             │              │              │     │
│  ▼             ▼             ▼              ▼              │     │
│ Port 5000    Port 8000    Port 5001      Port 5002        │     │
│ ┌─────────┐  ┌─────────┐  ┌──────────┐  ┌──────────┐     │     │
│ │Molecule │  │  Docker │  │Mol2Chemfig│ │ PubChem  │     │     │
│ │ Viewer  │  │ Backend │  │  Server   │ │  Server  │     │     │
│ │(Node.js)│  │(Docker) │  │ (Flask)   │ │ (Flask)  │     │     │
│ └─────────┘  └─────────┘  └──────────┘  └──────────┘     │     │
│                                                            │     │
└────────────────────────────────────────────────────────────┘     │
```

## Server Architecture

### 1. MoleculeViewer Server (Port 5000)

```
┌──────────────────────────────────────────────────────┐
│           MoleculeViewer (Node.js)                    │
├──────────────────────────────────────────────────────┤
│                                                       │
│  HTTP Endpoints:                                     │
│  ├─ GET /img/smiles?smiles=CCO                       │
│  ├─ GET /img/nomenclature?nomenclature=benzene       │
│  ├─ GET /health                                      │
│  └─ GET /cache-info                                  │
│                                                       │
│  Processing:                                         │
│  ├─ SMILES → Python (RDKit) → SVG                    │
│  ├─ Name → OPSIN → SMILES → SVG                      │
│  └─ Cache (MD5 hash)                                 │
│                                                       │
│  Cache: MoleculeViewer/cache/moleculeviewer/         │
└──────────────────────────────────────────────────────┘
```

### 2. Mol2ChemFig Docker Backend (Port 8000)

```
┌──────────────────────────────────────────────────────┐
│        Mol2ChemFig Docker Container                   │
├──────────────────────────────────────────────────────┤
│                                                       │
│  HTTP Endpoints:                                     │
│  ├─ POST /m2cf/submit     (Simple conversion)        │
│  ├─ POST /m2cf/apply      (With options)             │
│  ├─ POST /m2cf/layers     (Layered SVG)              │
│  └─ GET  /m2cf/search     (Name lookup)              │
│                                                       │
│  Processing:                                         │
│  ├─ SMILES → mol2chemfig → ChemFig code              │
│  ├─ ChemFig → LaTeX → PDF/SVG                        │
│  └─ Multiple rendering options                       │
│                                                       │
│  Technology:                                         │
│  └─ Python + LaTeX + mol2chemfig library             │
└──────────────────────────────────────────────────────┘
```

### 3. Mol2ChemFig Server (Port 5001)

```
┌──────────────────────────────────────────────────────┐
│       Mol2ChemFig Flask Wrapper                       │
├──────────────────────────────────────────────────────┤
│                                                       │
│  HTTP Endpoints:                                     │
│  ├─ POST /api/generate       (Generate image)        │
│  ├─ POST /api/generate-3d    (3D SMILES via OPSIN)   │
│  ├─ GET  /api/search         (Name search)           │
│  ├─ GET  /api/opsin          (Name → 3D SMILES)      │
│  ├─ POST /api/options/save   (Save defaults)         │
│  └─ GET  /images/{hash}.svg  (Serve cached)          │
│                                                       │
│  Features:                                           │
│  ├─ Persistent caching (SHA256 hash)                 │
│  ├─ SMILES canonicalization (dedupe)                 │
│  ├─ Options persistence                              │
│  ├─ OPSIN 3D integration                             │
│  └─ Docker backend wrapper                           │
│                                                       │
│  Dependencies:                                       │
│  ├─ Calls Docker backend (port 8000)                 │
│  ├─ Calls OPSIN API (external)                       │
│  └─ Uses RDKit for canonicalization                  │
│                                                       │
│  Cache: cache/mol2chemfig/                            │
└──────────────────────────────────────────────────────┘
```

### 4. PubChem Server (Port 5002)

```
┌──────────────────────────────────────────────────────┐
│            PubChem Integration Server                 │
├──────────────────────────────────────────────────────┤
│                                                       │
│  HTTP Endpoints:                                     │
│  ├─ GET /pubchem/img/{name}       (Direct PNG)       │
│  ├─ GET /pubchem/image?name=X     (With metadata)    │
│  ├─ GET /pubchem/3d-model?name=X  (SDF format)       │
│  ├─ GET /pubchem/3d-viewer?name=X (Interactive)      │
│  └─ GET /pubchem/info?name=X      (Compound info)    │
│                                                       │
│  Processing:                                         │
│  ├─ Name → PubChem CID lookup                        │
│  ├─ CID → Image/3D model fetch                       │
│  ├─ Cache locally                                    │
│  └─ Serve from cache or PubChem                      │
│                                                       │
│  External API:                                       │
│  └─ PubChem REST API (pubchem.ncbi.nlm.nih.gov)      │
│                                                       │
│  Cache: pubchem-cache/                                │
│     ├─ images/  (PNG files)                          │
│     └─ sdf/     (3D models)                          │
└──────────────────────────────────────────────────────┘
```

## Data Flow

### Simple SMILES Rendering

```
User Request
    │
    ▼
┌────────────────────┐
│  Chrome Extension  │  (detects "chem:benzene:")
│  or Direct API     │
└────────┬───────────┘
         │
         ▼
┌────────────────────────────────────────┐
│  Try MoleculeViewer (Port 5000)        │
│  GET /img/nomenclature?name=benzene    │
└────────┬───────────────────────────────┘
         │
         ├─ Success → Return SVG
         │
         └─ Fail ──→ Fallback to Mol2ChemFig
                     │
                     ▼
              ┌─────────────────────────┐
              │ Mol2ChemFig (Port 5001) │
              │ POST /api/generate      │
              └────────┬────────────────┘
                       │
                       ├──→ Check cache → Found → Return
                       │
                       └──→ Not cached
                            │
                            ▼
                     ┌──────────────────┐
                     │ Docker (Port 8000)│
                     │ POST /m2cf/submit │
                     └────────┬──────────┘
                              │
                              ▼
                        Generate SVG
                              │
                              ▼
                        Cache & Return
```

### 3D SMILES with OPSIN

```
User Request: "glucose" (with 3D)
    │
    ▼
┌────────────────────────┐
│ Mol2ChemFig Server     │
│ POST /api/generate-3d  │
└────────┬───────────────┘
         │
         ▼
┌────────────────────────┐
│  OPSIN External API    │
│  Name → 3D SMILES      │
│  "glucose" →           │
│  "O=C[C@H](O)..."      │
└────────┬───────────────┘
         │
         ▼
┌────────────────────────┐
│  Canonicalize SMILES   │
│  (RDKit)               │
└────────┬───────────────┘
         │
         ├─→ Check cache (canonical SMILES)
         │   │
         │   ├─ Found → Return
         │   │
         │   └─ Not found
         │       │
         ▼       ▼
┌────────────────────────┐
│ Docker Backend         │
│ POST /m2cf/apply       │
│ (with options)         │
└────────┬───────────────┘
         │
         ▼
┌────────────────────────┐
│  Generate SVG          │
│  Cache & Return        │
└────────────────────────┘
```

### PubChem 3D Viewer

```
User Request: 3D view of "histamine"
    │
    ▼
┌────────────────────────┐
│  PubChem Server        │
│  GET /3d-viewer?       │
│  name=histamine        │
└────────┬───────────────┘
         │
         ▼
┌────────────────────────┐
│  Name → CID Lookup     │
│  "histamine" → 774     │
│  (PubChem API)         │
└────────┬───────────────┘
         │
         ▼
┌────────────────────────┐
│  Fetch 3D Model (SDF)  │
│  CID 774 → 3D coords   │
│  (PubChem API)         │
└────────┬───────────────┘
         │
         ▼
┌────────────────────────┐
│  Cache SDF locally     │
│  pubchem-cache/sdf/    │
└────────┬───────────────┘
         │
         ▼
┌────────────────────────┐
│  Return HTML viewer    │
│  with embedded 3D      │
│  visualization         │
└────────────────────────┘
```

## Cache Strategy

### Multi-Level Caching

```
┌─────────────────────────────────────────────────────┐
│                  Cache Levels                        │
├─────────────────────────────────────────────────────┤
│                                                      │
│  Level 1: Server In-Memory Cache                    │
│  ├─ Mol2ChemFig: image_cache dict                   │
│  └─ Fast lookup for recent requests                 │
│                                                      │
│  Level 2: Disk Cache (Separated by Server)          │
│  ├─ MoleculeViewer:                                 │
│  │   MoleculeViewer/cache/moleculeviewer/           │
│  │   Key: MD5(type:value:widthxheight)              │
│  │                                                   │
│  ├─ Mol2ChemFig:                                    │
│  │   cache/mol2chemfig/                             │
│  │   Key: SHA256(canonical_smiles:sorted_options)   │
│  │                                                   │
│  └─ PubChem:                                        │
│      pubchem-cache/                                 │
│      Key: MD5(img_cid_size_type)                    │
│                                                      │
│  Level 3: Browser Cache                             │
│  └─ HTTP Cache-Control headers                      │
│                                                      │
└─────────────────────────────────────────────────────┘
```

### Deduplication Strategy

```
Problem: Same molecule, different cache entries
  - "CCO" and "C(C)O" are the same (ethanol)
  - Should have ONE cache entry

Solution: Canonicalization
  ┌──────────────┐
  │ Input SMILES │
  └──────┬───────┘
         │
         ▼
  ┌──────────────────┐
  │ Canonicalize     │
  │ (RDKit)          │
  │ "CCO" → "CCO"    │
  │ "C(C)O" → "CCO"  │
  └──────┬───────────┘
         │
         ▼
  ┌──────────────────┐
  │ Generate hash    │
  │ from CANONICAL   │
  │ SMILES only      │
  └──────┬───────────┘
         │
         ▼
  Single cache entry
```

## Master Launcher

### Control Flow

```
┌─────────────────────────────────────────────────────┐
│              launcher.html (Browser)                 │
├─────────────────────────────────────────────────────┤
│                                                      │
│  Every 5 seconds:                                   │
│  ┌────────────────────────────────────┐             │
│  │ Health Check Loop                  │             │
│  │ ├─ fetch(localhost:5000/health)    │             │
│  │ ├─ fetch(localhost:8000)           │             │
│  │ ├─ fetch(localhost:5001/health)    │             │
│  │ └─ fetch(localhost:5002/health)    │             │
│  └────────┬───────────────────────────┘             │
│           │                                          │
│           ├─ Success → Show green badge             │
│           └─ Fail → Show red badge                  │
│                                                      │
│  User Actions:                                      │
│  ├─ Start All → Alert to run start_all.bat         │
│  ├─ Stop All → Alert to run stop_all.bat           │
│  └─ Individual → Alert to run specific .bat        │
│                                                      │
└─────────────────────────────────────────────────────┘
```

### Batch Script Flow

```
start_all.bat
│
├─ Check prerequisites
│  ├─ Python installed?
│  ├─ Node.js installed?
│  └─ Docker available?
│
├─ Create .env if missing
│
├─ Stop existing processes
│
├─ Start Docker containers
│  └─ docker-compose up -d
│
├─ Start MoleculeViewer
│  └─ new window: node server.js
│
├─ Start Mol2ChemFig Server
│  └─ new window: python mol2chemfig_server.py
│
├─ Start PubChem Server
│  └─ new window: python pubchem_server.py
│
└─ Open launcher.html in browser
```

## File Structure

```
ChemParser/
│
├── Launcher System
│   ├── launcher.html              # Master control panel
│   ├── start_all.bat             # Start all servers
│   ├── stop_all.bat              # Stop all servers
│   ├── status.bat                # Check status
│   ├── start_moleculeviewer.bat  # Individual server start
│   ├── start_docker.bat
│   ├── start_mol2chemfig.bat
│   ├── start_pubchem.bat
│   ├── LAUNCHER_README.md        # Full documentation
│   ├── QUICK_START_LAUNCHER.md   # Quick start guide
│   └── SYSTEM_ARCHITECTURE.md    # This file
│
├── MoleculeViewer/ (Port 5000)
│   ├── server.js                 # Node.js server
│   ├── generate_svg.py           # Python SVG generator
│   ├── canonicalize_smiles.py    # SMILES canonicalization
│   └── cache/
│       └── moleculeviewer/       # SVG cache
│
├── Mol2ChemFig Server (Port 5001)
│   ├── mol2chemfig_server.py     # Flask server
│   ├── canonicalize_smiles.py    # SMILES canonicalization
│   └── cache/
│       └── mol2chemfig/          # SVG/PDF cache
│
├── PubChem Server (Port 5002)
│   ├── pubchem_server.py         # Flask server
│   └── pubchem-cache/
│       ├── images/               # PNG cache
│       └── sdf/                  # 3D model cache
│
├── Docker Backend (Port 8000)
│   ├── docker-compose.yml        # Docker configuration
│   ├── .env                      # Environment variables
│   └── m2cf_fixed.py             # Custom mol2chemfig
│
├── Chrome Extension
│   └── chem-extension/
│       ├── manifest.json
│       ├── content.js            # Text replacement
│       ├── popup.html
│       └── popup.js
│
└── Test Files
    ├── test_m2cf_full.html
    ├── test_pubchem_3d.html
    ├── test_3d_smiles.html
    ├── test_opsin_3d.html
    └── TEST_EXTENSION.html
```

## Port Summary

| Port | Server | Type | Protocol |
|------|--------|------|----------|
| 5000 | MoleculeViewer | Node.js | HTTP |
| 8000 | Mol2ChemFig Backend | Docker | HTTP |
| 5001 | Mol2ChemFig Server | Flask | HTTP |
| 5002 | PubChem Server | Flask | HTTP |

## Technology Stack

```
┌─────────────────────────────────────────┐
│          Technology Breakdown            │
├─────────────────────────────────────────┤
│                                          │
│  Backend Languages:                     │
│  ├─ Python 3.8+ (Flask servers)         │
│  ├─ Node.js 14+ (MoleculeViewer)        │
│  └─ Bash/Batch (Launcher scripts)       │
│                                          │
│  Chemistry Libraries:                   │
│  ├─ RDKit (SMILES processing)           │
│  ├─ mol2chemfig (LaTeX rendering)       │
│  └─ OPSIN (Name → SMILES)               │
│                                          │
│  Web Technologies:                      │
│  ├─ HTML5/CSS3/JavaScript (Launcher)    │
│  ├─ Flask + Flask-CORS (API servers)    │
│  ├─ Express + CORS (Node server)        │
│  └─ Chrome Extension API                │
│                                          │
│  Infrastructure:                        │
│  ├─ Docker + Docker Compose             │
│  ├─ LaTeX + tikz (SVG generation)       │
│  └─ Windows Batch Scripts               │
│                                          │
│  External APIs:                         │
│  ├─ PubChem REST API                    │
│  ├─ OPSIN Web Service                   │
│  └─ Docker Hub (images)                 │
│                                          │
└─────────────────────────────────────────┘
```

## Security Considerations

```
┌─────────────────────────────────────────┐
│             Security Notes               │
├─────────────────────────────────────────┤
│                                          │
│  ✓ All servers run on localhost only    │
│  ✓ CORS enabled for local development   │
│  ✓ No authentication (local only)       │
│  ✓ Input sanitization for SMILES        │
│  ✓ Cache isolated by server type        │
│                                          │
│  ⚠ Not for production use as-is         │
│  ⚠ No HTTPS (local development)         │
│  ⚠ No rate limiting                     │
│  ⚠ No user authentication               │
│                                          │
└─────────────────────────────────────────┘
```

## Performance Characteristics

```
┌─────────────────────────────────────────────────────┐
│              Performance Metrics                     │
├─────────────────────────────────────────────────────┤
│                                                      │
│  Cold Start (First Request):                        │
│  ├─ MoleculeViewer: ~500ms - 2s                     │
│  ├─ Mol2ChemFig: ~2s - 5s (LaTeX compilation)       │
│  └─ PubChem: ~1s - 3s (external API)                │
│                                                      │
│  Cached Request:                                    │
│  ├─ MoleculeViewer: <50ms                           │
│  ├─ Mol2ChemFig: <50ms                              │
│  └─ PubChem: <50ms                                  │
│                                                      │
│  Cache Size (Typical):                              │
│  ├─ MoleculeViewer: ~100KB per SVG                  │
│  ├─ Mol2ChemFig: ~50-200KB per SVG                  │
│  └─ PubChem: ~20-50KB per PNG, ~10KB per SDF        │
│                                                      │
└─────────────────────────────────────────────────────┘
```

## Future Enhancements

- [ ] Real process control from launcher (via Node.js backend)
- [ ] Real-time log streaming in launcher UI
- [ ] Health check with detailed metrics
- [ ] Auto-restart on failure
- [ ] Configuration UI for server settings
- [ ] Cache management UI with preview
- [ ] Batch processing interface
- [ ] API usage statistics
- [ ] Export/import cache functionality
- [ ] Linux/Mac support for batch scripts
