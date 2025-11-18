# MOL2CHEMFIG Docker - API Documentation

## Status
✅ **RUNNING** on localhost:8000 (Docker)

## Starting the System

```bash
cd C:\Users\Kapil\Personal\PROJECTS\Mol2chemfig
docker-compose up -d
```

## API Endpoint

### POST /m2cf/submit

Converts SMILES, InChI, or chemical names to chemfig LaTeX code.

**Request:**
```json
{
  "textAreaData": "CCO"
}
```

**Response:**
```json
{
  "chem_data": "CCO",
  "chem_format": "smiles",
  "chemfig": "\\chemfig{...}",
  "pdflink": "data:application/pdf;...",
  "svglink": "data:image/svg+xml;...",
  "error": null
}
```

## Response Fields You Can Extract

| Field | Type | Purpose |
|-------|------|---------|
| `chemfig` | string | LaTeX chemfig code for the molecule structure |
| `chem_format` | string | Input format detected (`smiles`, `name`, `inchi`) |
| `chem_data` | string | The normalized input data |
| `pdflink` | string | Base64-encoded PDF of the molecule (data URI) |
| `svglink` | string | Base64-encoded SVG of the molecule (data URI) |
| `error` | null or string | Error message if conversion failed |

## Examples

### Example 1: Ethanol (SMILES)
```bash
curl -X POST http://localhost:8000/m2cf/submit \
  -H "Content-Type: application/json" \
  -d '{"textAreaData": "CCO"}'
```

**Output:**
```json
{
  "chemfig": "\\chemfig{\n              % 1\n       -[:330]% 2\n    -[:30,,,1]OH% 3\n}",
  "chem_format": "smiles",
  "pdflink": "data:application/pdf;base64,JVBERi0xLjQK...",
  "svglink": "data:image/svg+xml;base64,PD94bWwgdmVyc2...",
  "error": null
}
```

### Example 2: Benzene (SMILES)
```bash
curl -X POST http://localhost:8000/m2cf/submit \
  -H "Content-Type: application/json" \
  -d '{"textAreaData": "c1ccccc1"}'
```

### Example 3: Acetone (SMILES)
```bash
curl -X POST http://localhost:8000/m2cf/submit \
  -H "Content-Type: application/json" \
  -d '{"textAreaData": "CC(=O)C"}'
```

### Example 4: Aspirin (Name lookup)
```bash
curl -X POST http://localhost:8000/m2cf/submit \
  -H "Content-Type: application/json" \
  -d '{"textAreaData": "aspirin"}'
```

## Python Example

```python
import requests
import json
import base64

url = "http://localhost:8000/m2cf/submit"

# Convert molecule
payload = {"textAreaData": "CCO"}
response = requests.post(url, json=payload)
data = response.json()

# Extract data
print("Chemfig LaTeX:", data['chemfig'])
print("Format detected:", data['chem_format'])

# Extract SVG from base64
if 'svglink' in data and data['svglink']:
    svg_base64 = data['svglink'].replace('data:image/svg+xml;base64,', '')
    svg_content = base64.b64decode(svg_base64).decode('utf-8')
    print("SVG Content:", svg_content)

# Extract PDF from base64
if 'pdflink' in data and data['pdflink']:
    pdf_base64 = data['pdflink'].replace('data:application/pdf;base64,', '')
    pdf_bytes = base64.b64decode(pdf_base64)
    # Save to file
    with open('molecule.pdf', 'wb') as f:
        f.write(pdf_bytes)
```

## What You Can Do With This Data

### 1. **Extract Chemfig LaTeX**
Use the `chemfig` field to create LaTeX documents:
```latex
\documentclass{article}
\usepackage{chemfig}
\begin{document}
\chemfig{...}
\end{document}
```

### 2. **Extract SVG Images**
Decode the `svglink` Base64 to get molecule structure images

### 3. **Extract PDF Images**
Decode the `pdflink` Base64 to get high-quality PDFs

### 4. **Store in Database**
Save the chemfig, SVG, and PDF for later use

### 5. **Convert to Other Formats**
Parse the chemfig data to generate other formats

## Containers Status

```
CONTAINER ID   IMAGE                                COMMAND              PORTS
ed8531110315   pychemist/m2cf_web_frontend:latest   "docker-entrypoint"  0.0.0.0:8080->9000/tcp
9b25072e2cf5   pychemist/m2cf_web_backend:latest    "python src/api.py"  0.0.0.0:8000->5000/tcp
```

## Common SMILES Examples

| Molecule | SMILES |
|----------|--------|
| Ethanol | `CCO` |
| Ethane | `CC` |
| Methane | `C` |
| Benzene | `c1ccccc1` |
| Toluene | `Cc1ccccc1` |
| Acetone | `CC(=O)C` |
| Acetic Acid | `CC(=O)O` |
| Ethyl Acetate | `CC(=O)OCC` |

## Troubleshooting

**Docker not responding:**
```bash
docker-compose restart
```

**Want to see logs:**
```bash
docker logs m2cf_backend --tail 20
```

**Stop containers:**
```bash
docker-compose down
```

---

**All set! The mol2chemfig API is ready to use for extracting molecule structures.**
