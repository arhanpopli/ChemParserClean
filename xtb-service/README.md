# xTB Orbital Service

Microservice that calculates molecular orbitals using xTB and returns cube files.

## API Endpoints

### `GET /orbital?smiles=CCO&orbital=homo`

Generates a molecular orbital cube file.

**Parameters:**
- `smiles` (required) - SMILES string of the molecule
- `orbital` (optional) - `homo`, `lumo`, `homo-1`, etc. Default: `homo`
- `gridSize` (optional) - Grid resolution. Default: `50`

**Response:**
```json
{
  "success": true,
  "smiles": "CCO",
  "orbital": "homo",
  "cubeFile": "...",
  "atoms": [...]
}
```

## Deployment

This service is deployed on Railway and called by the main Vercel API.

Railway URL: `https://xtb-service-production.up.railway.app`
