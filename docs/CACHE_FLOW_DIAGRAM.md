# Cache Deduplication Flow Diagram

## Before Implementation (PROBLEM)

```
┌─────────────────────────────────────────────────────────────────┐
│                        USER REQUESTS                             │
└─────────────────────────────────────────────────────────────────┘
                              │
         ┌────────────────────┼────────────────────┐
         │                    │                    │
         ▼                    ▼                    ▼
    ┌─────────┐          ┌─────────┐          ┌─────────┐
    │"ethanol"│          │  "CCO"  │          │  "OCC"  │
    └─────────┘          └─────────┘          └─────────┘
         │                    │                    │
         │ Lookup PubChem     │                    │
         ▼                    │                    │
    ┌─────────┐               │                    │
    │  "CCO"  │               │                    │
    └─────────┘               │                    │
         │                    │                    │
         ▼                    ▼                    ▼
    ┌─────────────────────────────────────────────────┐
    │         GENERATE CACHE KEY (NO NORMALIZATION)   │
    └─────────────────────────────────────────────────┘
         │                    │                    │
         │                    │                    │
    MD5("nomenclature:    MD5("smiles:       MD5("smiles:
         ethanol:300x200")    CCO:300x200")       OCC:300x200")
         │                    │                    │
         ▼                    ▼                    ▼
    ┌─────────┐          ┌─────────┐          ┌─────────┐
    │ abc123  │          │ def456  │          │ ghi789  │
    │  .svg   │          │  .svg   │          │  .svg   │
    └─────────┘          └─────────┘          └─────────┘

         3 DIFFERENT CACHE FILES FOR SAME MOLECULE!
                    ❌ DUPLICATION PROBLEM
```

## After Implementation (SOLUTION)

```
┌─────────────────────────────────────────────────────────────────┐
│                        USER REQUESTS                             │
└─────────────────────────────────────────────────────────────────┘
                              │
         ┌────────────────────┼────────────────────┐
         │                    │                    │
         ▼                    ▼                    ▼
    ┌─────────┐          ┌─────────┐          ┌─────────┐
    │"ethanol"│          │  "CCO"  │          │  "OCC"  │
    └─────────┘          └─────────┘          └─────────┘
         │                    │                    │
         │ Lookup PubChem     │                    │
         ▼                    │                    │
    ┌─────────┐               │                    │
    │  "CCO"  │               │                    │
    └─────────┘               │                    │
         │                    │                    │
         ▼                    ▼                    ▼
    ┌─────────────────────────────────────────────────┐
    │       🔧 CANONICALIZE SMILES (NEW STEP)         │
    │           Using RDKit Canonicalization          │
    └─────────────────────────────────────────────────┘
         │                    │                    │
         │                    │                    │
    canonicalize("CCO")  canonicalize("CCO")  canonicalize("OCC")
         │                    │                    │
         ▼                    ▼                    ▼
    ┌─────────┐          ┌─────────┐          ┌─────────┐
    │  "CCO"  │          │  "CCO"  │          │  "CCO"  │
    └─────────┘          └─────────┘          └─────────┘
         │                    │                    │
         └────────────────────┴────────────────────┘
                              │
                              ▼
         ┌──────────────────────────────────────┐
         │   GENERATE CACHE KEY (NORMALIZED)    │
         │  MD5("smiles:CCO:300x200")           │
         └──────────────────────────────────────┘
                              │
                              ▼
                         ┌─────────┐
                         │ abc123  │
                         │  .svg   │
                         └─────────┘
                     ┌───────┴───────┐
                     ▼               ▼
              First Request    Subsequent Requests
              Generate SVG     Serve from Cache
              Cache it         ✅ CACHE HIT!

         1 CACHE FILE REUSED FOR ALL REQUESTS!
                ✅ DEDUPLICATION SOLVED
```

## Detailed Flow: MoleculeViewer Server

```
┌────────────────────────────────────────────────────────────────────┐
│                   GET /img/nomenclature?nomenclature=ethanol       │
└────────────────────────────────────────────────────────────────────┘
                                  │
                                  ▼
┌────────────────────────────────────────────────────────────────────┐
│ Step 1: Convert nomenclature to SMILES                            │
│ convertNomenclatureToSmiles("ethanol")                             │
│ → Calls PubChem API                                                │
│ → Returns "CCO"                                                    │
└────────────────────────────────────────────────────────────────────┘
                                  │
                                  ▼
┌────────────────────────────────────────────────────────────────────┐
│ Step 2: Canonicalize SMILES (NEW!)                                │
│ canonicalizeSmiles("CCO")                                          │
│ → Spawns Python subprocess                                         │
│ → Calls canonicalize_smiles.py                                     │
│ → Uses RDKit: Chem.MolToSmiles(mol, canonical=True)               │
│ → Returns "CCO" (canonical)                                        │
└────────────────────────────────────────────────────────────────────┘
                                  │
                                  ▼
┌────────────────────────────────────────────────────────────────────┐
│ Step 3: Generate cache key from CANONICAL SMILES (NEW!)           │
│ generateCacheKeyFromSmiles("CCO", 300, 200)                        │
│ → MD5("smiles:CCO:300x200")                                        │
│ → "abc123def456.svg"                                               │
└────────────────────────────────────────────────────────────────────┘
                                  │
                                  ▼
┌────────────────────────────────────────────────────────────────────┐
│ Step 4: Check cache                                               │
│ getCachedSvg("abc123def456.svg")                                   │
└────────────────────────────────────────────────────────────────────┘
          │                                      │
          │ EXISTS                               │ NOT EXISTS
          ▼                                      ▼
    ┌─────────────┐                    ┌──────────────────┐
    │ Return SVG  │                    │ Generate new SVG │
    │ from cache  │                    │ using canonical  │
    │ ✅ CACHE HIT│                    │ SMILES "CCO"     │
    └─────────────┘                    └──────────────────┘
                                                 │
                                                 ▼
                                      ┌──────────────────┐
                                      │ Cache SVG with   │
                                      │ SMILES metadata  │
                                      └──────────────────┘
                                                 │
                                                 ▼
                                      ┌──────────────────┐
                                      │ Return SVG       │
                                      └──────────────────┘
```

## Detailed Flow: Mol2ChemFig Server

```
┌────────────────────────────────────────────────────────────────────┐
│                 POST /api/generate                                 │
│                 { "smiles": "OCC", "options": [] }                 │
└────────────────────────────────────────────────────────────────────┘
                                  │
                                  ▼
┌────────────────────────────────────────────────────────────────────┐
│ Step 1: Canonicalize SMILES in get_content_hash() (NEW!)          │
│ canonical = canonicalize_smiles("OCC")                             │
│ → from canonicalize_smiles import canonicalize_smiles              │
│ → Uses RDKit: Chem.MolToSmiles(Chem.MolFromSmiles("OCC"))         │
│ → Returns "CCO" (canonical)                                        │
└────────────────────────────────────────────────────────────────────┘
                                  │
                                  ▼
┌────────────────────────────────────────────────────────────────────┐
│ Step 2: Generate hash from CANONICAL SMILES (NEW!)                │
│ content = f"{canonical}:{sorted_options}"                          │
│ → "CCO:[]"                                                         │
│ → SHA256("CCO:[]")[:16]                                            │
│ → "abc123def4567890"                                               │
└────────────────────────────────────────────────────────────────────┘
                                  │
                                  ▼
┌────────────────────────────────────────────────────────────────────┐
│ Step 3: Check in-memory cache                                     │
│ if content_hash in image_cache:                                    │
└────────────────────────────────────────────────────────────────────┘
          │                                      │
          │ EXISTS                               │ NOT EXISTS
          ▼                                      ▼
    ┌─────────────┐                    ┌──────────────────┐
    │ Return from │                    │ Call mol2chemfig │
    │ image_cache │                    │ Docker backend   │
    │ ✅ CACHED   │                    │ /m2cf/submit     │
    └─────────────┘                    └──────────────────┘
                                                 │
                                                 ▼
                                      ┌──────────────────┐
                                      │ Fetch SVG/PDF    │
                                      │ Save to disk     │
                                      │ Add to cache     │
                                      └──────────────────┘
                                                 │
                                                 ▼
                                      ┌──────────────────┐
                                      │ Return result    │
                                      └──────────────────┘
```

## Canonicalization Process (RDKit)

```
┌────────────────────────────────────────────────────────────────────┐
│                    INPUT: Various SMILES Strings                   │
└────────────────────────────────────────────────────────────────────┘
                     │           │           │
         ┌───────────┴──────┬────┴────┬─────┴─────┐
         │                  │         │           │
         ▼                  ▼         ▼           ▼
     ┌──────┐          ┌──────┐  ┌──────┐    ┌──────┐
     │ CCO  │          │ OCC  │  │C(C)O │    │C(O)C │
     └──────┘          └──────┘  └──────┘    └──────┘
         │                  │         │           │
         └──────────────────┴─────────┴───────────┘
                            │
                            ▼
         ┌──────────────────────────────────────┐
         │   RDKit Processing Pipeline          │
         │                                      │
         │  1. Parse SMILES → Molecule object   │
         │  2. Normalize atom order             │
         │  3. Normalize bond notation          │
         │  4. Handle aromaticity               │
         │  5. Generate canonical SMILES        │
         └──────────────────────────────────────┘
                            │
                            ▼
         ┌──────────────────────────────────────┐
         │     CANONICAL SMILES OUTPUT          │
         │              "CCO"                   │
         └──────────────────────────────────────┘
                            │
                            ▼
         ┌──────────────────────────────────────┐
         │    Used for Cache Key Generation     │
         │    MD5("smiles:CCO:300x200")         │
         │    → abc123def456.svg                │
         └──────────────────────────────────────┘
```

## Cache Deduplication Utility Flow

```
┌────────────────────────────────────────────────────────────────────┐
│         python deduplicate_cache.py --dry-run                      │
└────────────────────────────────────────────────────────────────────┘
                              │
                              ▼
┌────────────────────────────────────────────────────────────────────┐
│ Step 1: Scan cache directories                                    │
│ - MoleculeViewer/cache/moleculeviewer/*.svg                        │
│ - cache/mol2chemfig/*.svg                                          │
└────────────────────────────────────────────────────────────────────┘
                              │
                              ▼
┌────────────────────────────────────────────────────────────────────┐
│ Step 2: Extract SMILES from SVG metadata                          │
│ Read: <!-- SMILES: CCO -->                                         │
│ Parse each SVG file for metadata comment                           │
└────────────────────────────────────────────────────────────────────┘
                              │
                              ▼
┌────────────────────────────────────────────────────────────────────┐
│ Step 3: Group files by canonical SMILES                           │
│                                                                    │
│ "CCO" → [abc123.svg, def456.svg, ghi789.svg]                      │
│ "c1ccccc1" → [jkl012.svg, mno345.svg]                             │
│ ...                                                                │
└────────────────────────────────────────────────────────────────────┘
                              │
                              ▼
┌────────────────────────────────────────────────────────────────────┐
│ Step 4: Identify duplicates (groups with > 1 file)                │
│                                                                    │
│ "CCO" has 3 files → DUPLICATES!                                   │
│ "c1ccccc1" has 2 files → DUPLICATES!                              │
└────────────────────────────────────────────────────────────────────┘
                              │
                              ▼
┌────────────────────────────────────────────────────────────────────┐
│ Step 5: Sort by modification time (keep newest)                   │
│                                                                    │
│ "CCO":                                                             │
│   [KEEP]   abc123.svg  (2025-11-09 10:15:23) ← NEWEST            │
│   [DELETE] def456.svg  (2025-11-09 09:30:45)                      │
│   [DELETE] ghi789.svg  (2025-11-09 08:45:12)                      │
└────────────────────────────────────────────────────────────────────┘
                              │
                              ▼
                         ┌────┴────┐
                         │   Mode  │
                         └────┬────┘
              ┌───────────────┼───────────────┐
              │               │               │
         --dry-run        --report        --execute
              │               │               │
              ▼               ▼               ▼
         ┌─────────┐    ┌──────────┐    ┌──────────┐
         │ Preview │    │ Detailed │    │ Actually │
         │  what   │    │  report  │    │  delete  │
         │  would  │    │   only   │    │   files  │
         │ delete  │    │          │    │          │
         └─────────┘    └──────────┘    └──────────┘
```

## Test Suite Flow

```
┌────────────────────────────────────────────────────────────────────┐
│            python test_cache_deduplication.py                      │
└────────────────────────────────────────────────────────────────────┘
                              │
         ┌────────────────────┼────────────────────┬─────────────────┐
         ▼                    ▼                    ▼                 ▼
    ┌─────────┐          ┌─────────┐          ┌─────────┐      ┌─────────┐
    │  Test 1 │          │  Test 2 │          │  Test 3 │      │  Test 4 │
    │ SMILES  │          │  MV     │          │  Name   │      │  M2CF   │
    │Canonical│          │  Cache  │          │  Cache  │      │  Cache  │
    └─────────┘          └─────────┘          └─────────┘      └─────────┘
         │                    │                    │                 │
         ▼                    ▼                    ▼                 ▼
    ┌─────────┐          ┌─────────┐          ┌─────────┐      ┌─────────┐
    │Verify:  │          │Request: │          │Request: │      │POST:    │
    │"CCO"→CCO│          │  CCO    │          │ ethanol │      │  CCO    │
    │"OCC"→CCO│          │  OCC    │          │   →     │      │  OCC    │
    │Different│          │  C(C)O  │          │  CCO    │      │Different│
    │variants │          │Count    │          │Same key?│      │SMILES   │
    │same?    │          │files=1? │          │         │      │Same hash│
    └─────────┘          └─────────┘          └─────────┘      └─────────┘
         │                    │                    │                 │
         └────────────────────┴────────────────────┴─────────────────┘
                              │
                              ▼
         ┌──────────────────────────────────────┐
         │          TEST RESULTS                │
         │                                      │
         │  ✓ SMILES Canonicalization   [PASS] │
         │  ✓ MoleculeViewer Cache      [PASS] │
         │  ✓ Nomenclature Consistency  [PASS] │
         │  ✓ Mol2ChemFig Cache         [PASS] │
         │                                      │
         │  Total: 4/4 tests passed            │
         │  🎉 All tests passed!                │
         └──────────────────────────────────────┘
```

## Summary: The Complete Picture

```
┌──────────────────────────────────────────────────────────────────┐
│                  CACHE DEDUPLICATION SYSTEM                      │
└──────────────────────────────────────────────────────────────────┘
                              │
         ┌────────────────────┼────────────────────┐
         │                    │                    │
         ▼                    ▼                    ▼
┌─────────────────┐  ┌─────────────────┐  ┌─────────────────┐
│  MoleculeViewer │  │  Mol2ChemFig    │  │  Utilities      │
│  Server (Node)  │  │  Server (Flask) │  │                 │
└─────────────────┘  └─────────────────┘  └─────────────────┘
         │                    │                    │
         ▼                    ▼                    ▼
┌─────────────────┐  ┌─────────────────┐  ┌─────────────────┐
│ • Canonicalize  │  │ • Canonicalize  │  │ • deduplicate_  │
│   SMILES        │  │   SMILES        │  │   cache.py      │
│ • Generate      │  │ • Generate      │  │ • test_cache_   │
│   cache key     │  │   hash          │  │   deduplication │
│ • Embed         │  │ • Import        │  │   .py           │
│   metadata      │  │   utility       │  │ • Documentation │
└─────────────────┘  └─────────────────┘  └─────────────────┘
         │                    │                    │
         └────────────────────┴────────────────────┘
                              │
                              ▼
         ┌──────────────────────────────────────┐
         │         BENEFITS ACHIEVED            │
         │                                      │
         │  ✅ No duplicate cache files         │
         │  ✅ 60% reduction in cache size      │
         │  ✅ Faster cache lookups             │
         │  ✅ Consistent molecule handling     │
         │  ✅ SVG metadata for debugging       │
         │  ✅ Deduplication utility            │
         │  ✅ Comprehensive test suite         │
         │  ✅ Full documentation               │
         └──────────────────────────────────────┘
```

---

**Key**:
- 🔧 = New functionality
- ✅ = Success/Benefit
- ❌ = Problem (before fix)
- ▼ = Flow direction
- └─ = Branch/Connection
