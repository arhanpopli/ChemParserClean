# 🎯 MoleculeViewer Architecture & Flow Diagrams

## 🏗️ System Architecture

```
┌─────────────────────────────────────────────────────────────────┐
│                         ChatGPT Interface                         │
│                     (Browser - Any Website)                       │
└────────────────────┬──────────────────────────────────────────────┘
                     │
                     ▼
┌─────────────────────────────────────────────────────────────────┐
│           Chrome Extension (content.js)                          │
│  ┌──────────────────────────────────────────────────────────┐  │
│  │ 1. Detects: chem:acetone or chem:CCO                   │  │
│  │ 2. Identifies: nomenclature or SMILES                  │  │
│  │ 3. Generates: Direct image URL                         │  │
│  │ 4. Creates: <img src="http://localhost:5000/...">    │  │
│  └──────────────────────────────────────────────────────────┘  │
└────────────────────┬──────────────────────────────────────────────┘
                     │
                     ▼
          [Browser Requests Image URL]
                     │
        ┌────────────┴────────────┐
        │                         │
        ▼                         ▼
  /img/smiles?           /img/nomenclature?
  smiles=CCO             nomenclature=acetone
        │                         │
        └────────────┬────────────┘
                     │
                     ▼
┌─────────────────────────────────────────────────────────────────┐
│      Node.js Server (server.js) - Port 5000                    │
│  ┌──────────────────────────────────────────────────────────┐  │
│  │ Step 1: Parse query parameters                          │  │
│  │ Step 2: Generate cache key (MD5 hash)                   │  │
│  │ Step 3: Check if cached: svg-cache/{hash}.svg           │  │
│  │         - IF YES: Return cached SVG (50ms) ✅           │  │
│  │         - IF NO: Continue to step 4                     │  │
│  │ Step 4: Call Python helpers via subprocess             │  │
│  │         - For nomenclature: name → SMILES              │  │
│  │         - For SMILES: generate SVG                     │  │
│  │ Step 5: Cache the SVG                                  │  │
│  │ Step 6: Return SVG with Content-Type: image/svg+xml   │  │
│  └──────────────────────────────────────────────────────────┘  │
└────────────────────┬──────────────────────────────────────────────┘
                     │
            ┌────────┴────────┐
            │                 │
            ▼                 ▼
    [RDKit Helper]    [PubChem API Helper]
    generate_svg.py   nomenclature_to_smiles.py
            │                 │
            ▼                 ▼
      SMILES→SVG          Name→SMILES
                     │
                     ▼
┌─────────────────────────────────────────────────────────────────┐
│           SVG Cache Directory                                   │
│         svg-cache/ {hash}.svg                                   │
└────────────────────┬──────────────────────────────────────────────┘
                     │
                     ▼
          [SVG Response Sent to Browser]
                     │
                     ▼
┌─────────────────────────────────────────────────────────────────┐
│                    Browser Renders                              │
│              Inline Molecule Image ✅                           │
│                    In ChatGPT                                   │
└─────────────────────────────────────────────────────────────────┘
```

---

## 📊 Request Flow Comparison

### Old Way (Flask + JSON + Blob)
```
Extension ──POST──> Flask /api/render-smiles
                           │
                           ▼
                       Return JSON:
                       {svg: "<svg>...", smiles: "CCO"}
                           │
                           ▼
                    Parse JSON in JavaScript
                           │
                           ▼
                    Convert SVG to Blob
                           │
                           ▼
                    Create blob:// URL
                           │
                           ▼
                    Set img.src = "blob://..."
                           │
                           ▼
                    Browser loads blob URL
                           │
                           ▼
                    Display image
```
❌ **Complicated**: 7 steps  
❌ **Fragile**: Multiple conversions  
❌ **No sharing**: Blob URLs are temporary  

---

### New Way (Node.js + Direct URL)
```
Extension ──GET──> Node.js /img/smiles?smiles=CCO
                           │
                           ▼
                       Check cache
                           │
                ┌──────────┴──────────┐
                │                    │
         Cache HIT ❌        Cache MISS ✅
                │                    │
                ▼                    ▼
           Return cached        Process:
           SVG (50ms) ✅        1. Generate SVG
                │               2. Save cache
                │               3. Return SVG
                │               │
                └───────┬───────┘
                        │
                        ▼
              Set img.src = URL
                        │
                        ▼
          Browser fetches SVG directly
                        │
                        ▼
              Display image ✅
```
✅ **Simple**: 3 steps  
✅ **Robust**: Direct URL handling  
✅ **Shareable**: URLs work forever  

---

## 🔄 Pattern Detection Flowchart

```
User types: chem:acetone (or chem:CCO, etc.)
             │
             ▼
   Extension detects "chem:"
             │
             ▼
  Extract content after "chem:"
     (e.g., "acetone" or "CCO")
             │
             ▼
   Check for special chemistry characters:
   = [ ] ( ) @ + # \
             │
    ┌────────┴────────┐
    │                 │
 YES ▼                 ▼ NO
 Contains            No special
 special chars       characters
    │                 │
    ▼                 ▼
 SMILES          NOMENCLATURE
    │                 │
    ▼                 ▼
/img/smiles?   /img/nomenclature?
smiles=CCO     nomenclature=acetone
    │                 │
    └────────┬────────┘
             │
             ▼
    [Send to MoleculeViewer]
```

---

## 💾 Cache Key Generation

```
Input: nomenclature=acetone (width=300, height=200)
       │
       ▼
Create cache key string:
"nomenclature:acetone:300x200"
       │
       ▼
MD5 Hash:
e4f2a1b7c3d9f1a5...
       │
       ▼
Filename:
e4f2a1b7.svg
       │
       ▼
Full path:
svg-cache/e4f2a1b7.svg
```

### Cache Key Examples:
```
SMILES:CCO:300x200          → a7f3b2c1.svg
SMILES:c1ccccc1:300x200     → d2e8f1a9.svg
nomenclature:acetone:300x200 → e4f2a1b7.svg
nomenclature:benzene:300x200 → f5c2b8d3.svg
```

---

## ⏱️ Timeline: First vs Cached Request

### First Request (No Cache)
```
T=0ms   Browser: GET /img/smiles?smiles=CCO
        │
T=10ms  Server: Check cache → NOT FOUND
        │
T=50ms  Server: Call Python generate_svg.py
        │
T=150ms Python: Create molecule from SMILES
        │
T=200ms Python: Render to SVG
        │
T=250ms Python: Return SVG to Node.js
        │
T=260ms Node.js: Save to svg-cache/a7f3b2c1.svg
        │
T=270ms Node.js: Send SVG to browser
        │
T=280ms Browser: Receive SVG
        │
T=300ms Browser: Render SVG inline ✅
        
Total: ~300ms
```

### Second Request (Cached)
```
T=0ms   Browser: GET /img/smiles?smiles=CCO
        │
T=5ms   Server: Check cache → FOUND ✅
        │
T=50ms  Server: Read from svg-cache/a7f3b2c1.svg
        │
T=100ms Node.js: Send cached SVG to browser
        │
T=150ms Browser: Receive SVG
        │
T=170ms Browser: Render SVG inline ✅
        
Total: ~170ms (5x faster!)
```

---

## 🎯 Routing Logic

```
Request arrives:
GET /img/smiles?smiles=CCO

            ▼

Check parameters:
- smiles: ✅ Present
- nomenclature: ❌ Absent

            ▼

Route to SMILES handler:
1. smiles = "CCO"
2. width = 300 (default)
3. height = 200 (default)

            ▼

Generate cache key:
"smiles:CCO:300x200"

            ▼

Check cache:
- svg-cache/a7f3b2c1.svg exists?
  - YES: Return cached ✅
  - NO: Generate new

            ▼

Call: generate_svg.py CCO
Returns: <svg>...</svg>

            ▼

Save to cache:
svg-cache/a7f3b2c1.svg

            ▼

Return to browser:
Content-Type: image/svg+xml
<svg>...</svg>
```

---

## 🐍 Python Helper Execution

### generate_svg.py Execution
```
Node.js spawns process:
spawn('python', [
  'generate_svg.py',
  JSON.stringify({
    smiles: 'CCO',
    width: 300,
    height: 200,
    options: {...}
  })
])

            ▼

Python receives JSON input

            ▼

Parse JSON (sys.argv[1])

            ▼

Create RDKit molecule:
Chem.MolFromSmiles('CCO')

            ▼

Generate 2D coordinates:
AllChem.Compute2DCoords(mol)

            ▼

Draw to SVG:
Draw.MolDraw2DSVG(300, 200)

            ▼

Output JSON:
{
  "svg": "<svg>...</svg>"
}

            ▼

Node.js receives stdout

            ▼

Parse JSON, extract SVG
```

### nomenclature_to_smiles.py Execution
```
Node.js spawns process:
spawn('python', ['nomenclature_to_smiles.py', 'acetone'])

            ▼

Python receives nomenclature

            ▼

Query PubChem API:
GET /rest/pug/compound/name/acetone/property/CanonicalSMILES/JSON

            ▼

PubChem responds:
{
  "properties": [{
    "CanonicalSMILES": "CC(=O)C"
  }]
}

            ▼

Extract SMILES:
"CC(=O)C"

            ▼

Output JSON:
{
  "smiles": "CC(=O)C"
}

            ▼

Node.js receives stdout

            ▼

Parse JSON, get SMILES

            ▼

Call generate_svg.py with SMILES
```

---

## 📈 Performance Profile

```
Request/Response Waterfall:

1. Extension detects chem:acetone
   └─ Time: ~10ms

2. Extension creates URL
   └─ Time: ~5ms

3. Browser fetches URL
   └─ Time: ~50ms (network)

4. Server receives request
   └─ Time: ~1ms

5. Server checks cache
   ├─ Cache HIT: ~5ms ✅
   └─ Cache MISS: Processes...

6. For cache MISS:
   ├─ nomenclature_to_smiles.py
   │  └─ Time: ~1000-2000ms (PubChem API)
   │
   └─ generate_svg.py
      └─ Time: ~100-200ms (RDKit)

7. Server saves cache
   └─ Time: ~10ms

8. Server returns SVG
   └─ Time: ~50ms (network)

9. Browser renders
   └─ Time: ~20ms

Total First Load:  ~1200-2300ms
Total Cached Load: ~140-200ms
Cache Speedup:     8-16x faster
```

---

## 🎯 Decision Tree: Which Endpoint?

```
User input detected: "chem:something"

            ▼

Extract: "something"

            ▼

Does "something" contain
any of: = [ ] ( ) @ + # \ ?

    ┌─────────┴─────────┐
    │                   │
   YES                  NO
    │                   │
    ▼                   ▼
 SMILES           NOMENCLATURE
    │                   │
    ▼                   ▼
/img/smiles      /img/nomenclature
?smiles=...      ?nomenclature=...
    │                   │
    ▼                   ▼
RDKit renders      PubChem converts
SMILES to SVG      name to SMILES
    │               then RDKit
    │               renders SVG
    └───┬───────────┘
        │
        ▼
    Return SVG
    to browser
```

---

## 🔗 URL Structure

### SMILES Endpoint URL
```
http://localhost:5000/img/smiles?smiles=CCO&width=300&height=200
│      │    │          │        │   │  │   │
│      │    │          │        │   │  │   └─ Height (optional)
│      │    │          │        │   │  └───── Width (optional)
│      │    │          │        │   └──────── Value encoded
│      │    │          │        └──────────── Parameter name
│      │    │          └───────────────────── Endpoint
│      │    └──────────────────────────────── Port
│      └─────────────────────────────────── Localhost
└──────────────────────────────────────────── Protocol
```

### Nomenclature Endpoint URL
```
http://localhost:5000/img/nomenclature?nomenclature=acetone&width=300&height=200
```

---

## 📊 File Organization

```
MoleculeViewer/
│
├── 🚀 Node.js Server
│   ├── server.js (450 lines)
│   │   ├── Express app setup
│   │   ├── GET /img/smiles handler
│   │   ├── GET /img/nomenclature handler
│   │   ├── Cache management
│   │   └── Error handling
│   │
│   └── package.json (dependencies)
│       ├── express
│       ├── cors
│       └── axios
│
├── 🐍 Python Helpers
│   ├── generate_svg.py (70 lines)
│   │   ├── Import RDKit
│   │   ├── Parse SMILES
│   │   └── Render to SVG
│   │
│   └── nomenclature_to_smiles.py (60 lines)
│       ├── Query PubChem API
│       ├── Extract SMILES
│       └── Return JSON
│
├── 💾 Generated Files
│   └── svg-cache/
│       ├── a7f3b2c1.svg
│       ├── d2e8f1a9.svg
│       ├── e4f2a1b7.svg
│       └── ...
│
└── 📖 Documentation
    ├── README.md
    ├── SETUP_GUIDE.md
    ├── QUICK_REF.md
    ├── CONVERSION_SUMMARY.md
    └── ARCHITECTURE.md (this file)
```

---

## 🎊 Summary

✅ **Simple**: Direct URL → SVG image  
✅ **Fast**: Caching makes repeats instant  
✅ **Smart**: Automatic SMILES/nomenclature detection  
✅ **Reliable**: Error handling for edge cases  
✅ **Scalable**: Can handle many concurrent requests  
✅ **Shareable**: URLs persist for cached molecules  

**That's how your molecule viewer works!** 🧪⚗️
