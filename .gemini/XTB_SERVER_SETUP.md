# xTB Server Setup Plan

## Overview

Set up xTB on a VPS to run molecular orbital calculations on demand.

---

## Server Requirements

### Minimum Specs
- **OS**: Ubuntu 22.04 LTS
- **RAM**: 1 GB (2 GB recommended)
- **Storage**: 20 GB
- **Cost**: ~$5/month (DigitalOcean, Linode, Vultr)

### Software Stack
- xTB (tight-binding calculations)
- Node.js (API server)
- nginx (reverse proxy)

---

## Architecture

```
┌──────────────┐     ┌────────────────────────────────────┐
│ ChemTex      │     │         xTB VPS Server             │
│ Extension    │     │                                    │
│              │     │  ┌─────────────┐   ┌────────────┐  │
│ chem:mo=     │────▶│  │  Node.js    │──▶│   xTB      │  │
│ caffeine+    │     │  │  API        │   │   Binary   │  │
│ homo:        │     │  │             │◀──│            │  │
│              │◀────│  │  Returns    │   │  Molden    │  │
│   ┌──────┐   │     │  │  Cube File  │   │  Output    │  │
│   │Viewer│   │     │  └─────────────┘   └────────────┘  │
│   │      │   │     │                                    │
│   └──────┘   │     │  + Cache Layer (saves results)    │
└──────────────┘     └────────────────────────────────────┘
```

---

## API Design

### Endpoint: POST /calculate

**Request:**
```json
{
  "smiles": "Cn1cnc2c1c(=O)n(c(=O)n2C)C",
  "orbital": "homo",
  "gridSize": 40
}
```

**Response:**
```json
{
  "success": true,
  "molecule": "caffeine",
  "orbital": "HOMO",
  "energy": "-0.234 Ha",
  "cubeFile": "... cube data ...",
  "atoms": [
    {"element": "C", "x": 0.0, "y": 0.0, "z": 0.0},
    ...
  ],
  "computeTime": "1.2s"
}
```

---

## Installation Script

```bash
#!/bin/bash
# install-xtb.sh

# Update system
apt update && apt upgrade -y

# Install dependencies
apt install -y build-essential gfortran cmake wget unzip

# Download xTB binary
wget https://github.com/grimme-lab/xtb/releases/download/v6.6.1/xtb-6.6.1-linux-x86_64.tar.xz
tar -xf xtb-6.6.1-linux-x86_64.tar.xz
mv xtb-6.6.1 /opt/xtb

# Add to PATH
echo 'export PATH=/opt/xtb/bin:$PATH' >> /etc/profile.d/xtb.sh
source /etc/profile.d/xtb.sh

# Verify installation
xtb --version

# Install Node.js
curl -fsSL https://deb.nodesource.com/setup_20.x | bash -
apt install -y nodejs

# Create app directory
mkdir -p /opt/xtb-server
cd /opt/xtb-server
npm init -y
npm install express cors

echo "xTB server installed successfully!"
```

---

## API Server Code

```javascript
// /opt/xtb-server/server.js

const express = require('express');
const cors = require('cors');
const { exec } = require('child_process');
const fs = require('fs');
const path = require('path');
const crypto = require('crypto');

const app = express();
app.use(cors());
app.use(express.json());

const CACHE_DIR = '/opt/xtb-server/cache';
const WORK_DIR = '/opt/xtb-server/work';

// Ensure directories exist
fs.mkdirSync(CACHE_DIR, { recursive: true });
fs.mkdirSync(WORK_DIR, { recursive: true });

// Convert SMILES to XYZ using Open Babel (or fetch from PubChem)
async function smilesToXyz(smiles) {
  // Use PubChem to get 3D structure
  const url = `https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/smiles/${encodeURIComponent(smiles)}/SDF?record_type=3d`;
  const response = await fetch(url);
  if (!response.ok) throw new Error('Could not get 3D structure');
  
  const sdf = await response.text();
  // Convert SDF to XYZ (simple parser)
  return sdfToXyz(sdf);
}

// Run xTB calculation
async function runXtb(xyzFile, workDir) {
  return new Promise((resolve, reject) => {
    const cmd = `cd ${workDir} && xtb ${xyzFile} --gfn 2 --molden --sp`;
    
    exec(cmd, { timeout: 60000 }, (error, stdout, stderr) => {
      if (error) {
        reject(new Error(`xTB failed: ${stderr}`));
        return;
      }
      resolve({
        stdout,
        moldenFile: path.join(workDir, 'molden.input')
      });
    });
  });
}

// Parse Molden file to extract orbital coefficients
function parseMolden(moldenPath) {
  const content = fs.readFileSync(moldenPath, 'utf8');
  // Parse [MO] section for orbital coefficients
  // ... parsing logic ...
  return { orbitals: [], energies: [] };
}

// Generate cube file from orbital data
function generateCube(atoms, orbital, gridSize) {
  // ... cube generation logic ...
}

// Main calculation endpoint
app.post('/calculate', async (req, res) => {
  try {
    const { smiles, orbital = 'homo', gridSize = 40 } = req.body;
    
    if (!smiles) {
      return res.status(400).json({ error: 'SMILES required' });
    }
    
    // Create unique work directory
    const hash = crypto.createHash('md5').update(smiles).digest('hex').slice(0, 8);
    const workDir = path.join(WORK_DIR, hash);
    fs.mkdirSync(workDir, { recursive: true });
    
    // Check cache first
    const cacheKey = `${hash}_${orbital}_${gridSize}`;
    const cachePath = path.join(CACHE_DIR, `${cacheKey}.json`);
    
    if (fs.existsSync(cachePath)) {
      console.log(`Cache hit: ${cacheKey}`);
      const cached = JSON.parse(fs.readFileSync(cachePath, 'utf8'));
      return res.json(cached);
    }
    
    console.log(`Computing: ${smiles} (${orbital})`);
    const startTime = Date.now();
    
    // Get 3D structure
    const xyz = await smilesToXyz(smiles);
    const xyzPath = path.join(workDir, 'molecule.xyz');
    fs.writeFileSync(xyzPath, xyz);
    
    // Run xTB
    const result = await runXtb(xyzPath, workDir);
    
    // Parse Molden output
    const moldenData = parseMolden(result.moldenFile);
    
    // Generate cube file
    const cubeFile = generateCube(moldenData.atoms, moldenData.orbitals[orbital], gridSize);
    
    const computeTime = ((Date.now() - startTime) / 1000).toFixed(2);
    
    const response = {
      success: true,
      smiles,
      orbital,
      energy: moldenData.energies[orbital],
      cubeFile,
      atoms: moldenData.atoms,
      computeTime: `${computeTime}s`
    };
    
    // Cache result
    fs.writeFileSync(cachePath, JSON.stringify(response));
    
    res.json(response);
    
  } catch (error) {
    console.error('Calculation error:', error);
    res.status(500).json({ error: error.message });
  }
});

// Health check
app.get('/health', (req, res) => {
  res.json({ status: 'ok', version: '1.0.0' });
});

const PORT = process.env.PORT || 3001;
app.listen(PORT, () => {
  console.log(`xTB server running on port ${PORT}`);
});
```

---

## Deployment Steps

### Step 1: Provision VPS
- Create Ubuntu 22.04 droplet on DigitalOcean ($5/month)
- Note the IP address

### Step 2: SSH and Install
```bash
ssh root@YOUR_VPS_IP
wget https://raw.githubusercontent.com/.../install-xtb.sh
chmod +x install-xtb.sh
./install-xtb.sh
```

### Step 3: Start Server
```bash
cd /opt/xtb-server
node server.js
```

### Step 4: Set Up Process Manager
```bash
npm install -g pm2
pm2 start server.js --name xtb-server
pm2 startup
pm2 save
```

### Step 5: Configure nginx (optional)
```nginx
server {
    listen 80;
    server_name xtb.yourdomain.com;
    
    location / {
        proxy_pass http://localhost:3001;
        proxy_http_version 1.1;
    }
}
```

---

## Testing

```bash
# Test the API
curl -X POST http://YOUR_VPS_IP:3001/calculate \
  -H "Content-Type: application/json" \
  -d '{"smiles": "O", "orbital": "homo"}'
```

---

## Cost Analysis

| Item | Monthly Cost |
|------|--------------|
| DigitalOcean Basic Droplet | $4-6 |
| Domain (optional) | $1 |
| **Total** | **~$5-7/month** |

---

## Next Steps

1. [ ] Set up VPS
2. [ ] Install xTB
3. [ ] Deploy API server
4. [ ] Test with sample molecules
5. [ ] Integrate with ChemTex extension
6. [ ] Modify Falstad for multi-atom display

---

## Alternative: Local Testing First

Before deploying to VPS, we can test locally on Windows using WSL:

```bash
# In WSL/Ubuntu
sudo apt install xTB
# OR download binary
```

Want me to start with local testing in WSL first?
