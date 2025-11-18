# Quick Start: 3D SMILES with OPSIN

## 30-Second Setup

### 1. Start Servers
```bash
# Terminal 1: Docker backend
cd mol2chemfig-docker && docker-compose up

# Terminal 2: mol2chemfig server
python mol2chemfig_server.py
```

### 2. Test OPSIN
```bash
# Get 3D SMILES
curl "http://localhost:5001/api/opsin?name=glucose"

# Generate 3D structure
curl -X POST http://localhost:5001/api/generate-3d \
  -H "Content-Type: application/json" \
  -d '{"name":"D-glucose"}'
```

### 3. Open Test Page
```
Open: test_3d_smiles.html
URL: http://localhost:5001 (or file:// path)
```

## Quick API Reference

### OPSIN Conversion
```javascript
// GET
fetch('http://localhost:5001/api/opsin?name=glucose')
  .then(r => r.json())
  .then(d => console.log(d.smiles));
// Output: O=C[C@H](O)[C@@H](O)[C@H](O)[C@H](O)CO

// POST
fetch('http://localhost:5001/api/opsin', {
  method: 'POST',
  headers: { 'Content-Type': 'application/json' },
  body: JSON.stringify({ name: 'glucose' })
});
```

### 3D Structure Generation
```javascript
fetch('http://localhost:5001/api/generate-3d', {
  method: 'POST',
  headers: { 'Content-Type': 'application/json' },
  body: JSON.stringify({
    name: 'D-glucose',
    options: ['-o', '-m'],  // Optional: chemfig options
    return_format: 'svg'    // 'svg', 'pdf', or 'both'
  })
})
.then(r => r.json())
.then(d => {
  console.log('3D SMILES:', d.smiles);
  console.log('SVG URL:', d.svg_url);
});
```

## Extension Usage

1. Open extension popup
2. Select "mol2chemfig" engine
3. Enable "3D Stereochemistry (OPSIN)"
4. On webpage: `chem:D-glucose:`
5. Extension uses OPSIN for 3D SMILES

## Test Cases

### Simple Molecules
- `glucose` - Basic sugar
- `D-glucose` - With stereochemistry
- `L-alanine` - Amino acid

### Complex Molecules
- `cholesterol` - Steroid with rings
- `sucrose` - Disaccharide
- `caffeine` - Alkaloid

### Stereochemistry Examples
- `D-glucose` vs `L-glucose` (enantiomers)
- `D-fructose` (ketose)
- `aspirin` (simple drug)

## 2D vs 3D Comparison

| Molecule | 2D SMILES | 3D SMILES (OPSIN) |
|----------|-----------|-------------------|
| Glucose | `C(C1C(C(C(C(O1)O)O)O)O)O` | `O=C[C@H](O)[C@@H](O)[C@H](O)[C@H](O)CO` |
| Alanine | `CC(C(=O)O)N` | `N[C@@H](C)C(=O)O` |

Notice the `@` and `@@` markers in 3D SMILES indicating stereochemistry.

## Troubleshooting

### OPSIN Returns Error
```bash
# Test OPSIN directly
curl "https://opsin.ch.cam.ac.uk/opsin/glucose.json"
# If this fails, check internet connection
```

### No SVG Generated
```bash
# Check mol2chemfig server
curl http://localhost:5001/health

# Check Docker backend
curl http://localhost:8000/
```

### Extension Not Using 3D
1. Verify toggle is enabled in popup
2. Check Chrome DevTools console for errors
3. Ensure "mol2chemfig" engine is selected

## Files Reference

| File | Purpose |
|------|---------|
| `test_3d_smiles.html` | Test interface |
| `mol2chemfig_server.py` | Server with OPSIN endpoints |
| `chem-extension/popup.html` | UI with 3D toggle |
| `OPSIN_3D_IMPLEMENTATION.md` | Full documentation |

## Quick Examples

### HTML Test Page
```html
<!DOCTYPE html>
<html>
<body>
  <div id="result"></div>
  <script>
    fetch('http://localhost:5001/api/generate-3d', {
      method: 'POST',
      headers: { 'Content-Type': 'application/json' },
      body: JSON.stringify({ name: 'glucose' })
    })
    .then(r => r.json())
    .then(d => {
      document.getElementById('result').innerHTML =
        `<img src="http://localhost:5001${d.svg_url}">`;
    });
  </script>
</body>
</html>
```

### Python Script
```python
import requests

# Get 3D SMILES
response = requests.get('http://localhost:5001/api/opsin',
                       params={'name': 'glucose'})
data = response.json()
print(f"3D SMILES: {data['smiles']}")

# Generate structure
response = requests.post('http://localhost:5001/api/generate-3d',
                        json={'name': 'D-glucose'})
data = response.json()
print(f"SVG URL: {data['svg_url']}")
```

### Node.js Script
```javascript
const fetch = require('node-fetch');

async function test3D() {
  // Get 3D SMILES
  const opsin = await fetch('http://localhost:5001/api/opsin?name=glucose');
  const opsinData = await opsin.json();
  console.log('3D SMILES:', opsinData.smiles);

  // Generate structure
  const gen = await fetch('http://localhost:5001/api/generate-3d', {
    method: 'POST',
    headers: { 'Content-Type': 'application/json' },
    body: JSON.stringify({ name: 'D-glucose' })
  });
  const genData = await gen.json();
  console.log('SVG URL:', genData.svg_url);
}

test3D();
```

## Server Endpoints Summary

```
http://localhost:5001
├─ /api/opsin           GET/POST  Name → 3D SMILES
├─ /api/generate-3d     POST      Name → 3D structure
├─ /api/generate        POST      SMILES → structure
├─ /api/search          GET/POST  Name → structure
└─ /images/<hash>.svg   GET       Serve cached images
```

## Next Steps

1. ✅ Test with `test_3d_smiles.html`
2. ✅ Try different molecules
3. ✅ Compare 2D vs 3D rendering
4. ✅ Enable in Chrome extension
5. ✅ Read full docs: `OPSIN_3D_IMPLEMENTATION.md`

---

**Ready to use!** Start with `test_3d_smiles.html` for interactive testing.
