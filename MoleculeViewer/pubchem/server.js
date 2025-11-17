/**
 * PubChem Integration Server
 * Lightweight Node.js server that provides direct links to PubChem images and 3D viewers
 * Port: 5002
 */

const express = require('express');
const cors = require('cors');
const axios = require('axios');
const fs = require('fs');
const path = require('path');
const crypto = require('crypto');

const app = express();
const PORT = 5002;

// Middleware
app.use(cors());
app.use(express.json());

// Serve static files (for 3D library)
app.use('/static', express.static(path.join(__dirname, 'static')));

// Cache directory for CID lookups (NOT images - those are direct from PubChem)
const CACHE_DIR = path.join(__dirname, '..', 'cache', 'pubchem');
if (!fs.existsSync(CACHE_DIR)) {
  fs.mkdirSync(CACHE_DIR, { recursive: true });
}

const CID_CACHE_FILE = path.join(CACHE_DIR, 'cid_cache.json');

// In-memory CID cache (name -> CID mapping)
let cidCache = {};

// Load CID cache from file
function loadCidCache() {
  try {
    if (fs.existsSync(CID_CACHE_FILE)) {
      const data = fs.readFileSync(CID_CACHE_FILE, 'utf-8');
      cidCache = JSON.parse(data);
      console.log(`   Loaded ${Object.keys(cidCache).length} cached CID lookups`);
    }
  } catch (e) {
    console.error('   Error loading CID cache:', e.message);
    cidCache = {};
  }
}

// Save CID cache to file
function saveCidCache() {
  try {
    fs.writeFileSync(CID_CACHE_FILE, JSON.stringify(cidCache, null, 2));
  } catch (e) {
    console.error('   Error saving CID cache:', e.message);
  }
}

console.log('='.repeat(70));
console.log('PubChem Integration Server Starting...');
console.log(`Cache directory: ${CACHE_DIR}`);
console.log('='.repeat(70));

// Load CID cache on startup
loadCidCache();

// ============================================================
// HELPER FUNCTIONS
// ============================================================

/**
 * Get PubChem CID from compound name or SMILES
 */
async function getCompoundCid(name) {
  const normalizedName = name.toLowerCase().trim();

  // Check cache first
  if (cidCache[normalizedName]) {
    return cidCache[normalizedName];
  }

  try {
    // Try by name first
    let url = `https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/${encodeURIComponent(name)}/cids/JSON`;
    let response = await axios.get(url, { timeout: 10000 });

    if (response.status === 200 && response.data.IdentifierList && response.data.IdentifierList.CID) {
      const cid = response.data.IdentifierList.CID[0];
      cidCache[normalizedName] = cid;
      saveCidCache();
      return cid;
    }
  } catch (e) {
    // Name lookup failed, try SMILES
    try {
      let url = `https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/smiles/${encodeURIComponent(name)}/cids/JSON`;
      let response = await axios.get(url, { timeout: 10000 });

      if (response.status === 200 && response.data.IdentifierList && response.data.IdentifierList.CID) {
        const cid = response.data.IdentifierList.CID[0];
        cidCache[normalizedName] = cid;
        saveCidCache();
        return cid;
      }
    } catch (e2) {
      console.error(`   Error getting CID for ${name}:`, e2.message);
    }
  }

  return null;
}

/**
 * Generate PubChem image URL
 */
function getPubChemImageUrl(cid, imageSize = 'large', recordType = '2d') {
  const baseUrl = `https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/${cid}/PNG`;
  const params = new URLSearchParams();

  if (imageSize) {
    params.append('image_size', imageSize);
  }
  if (recordType && recordType !== '2d') {
    params.append('record_type', recordType);
  }

  const queryString = params.toString();
  return queryString ? `${baseUrl}?${queryString}` : baseUrl;
}

/**
 * Generate PubChem 3D viewer URL
 */
function getPubChem3DViewerUrl(cid) {
  return `https://pubchem.ncbi.nlm.nih.gov/compound/${cid}#section=3D-Conformer`;
}

/**
 * Generate PubChem SDF URL
 */
function getPubChemSdfUrl(cid, recordType = '3d') {
  return `https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/${cid}/record/SDF?record_type=${recordType}`;
}

// ============================================================
// ROUTES - Server Info
// ============================================================

/**
 * GET / - Server info
 */
app.get('/', (req, res) => {
  res.json({
    name: 'PubChem Integration Server',
    version: '1.0.0',
    port: PORT,
    description: 'Lightweight server providing direct links to PubChem images and 3D viewers',
    endpoints: {
      health: 'GET /health',
      image: 'GET /img/:name',
      image_3d: 'GET /img/:name?type=3d',
      viewer_3d: 'GET /3d/:name',
      info: 'GET /info/:name',
      cache_info: 'GET /cache-info',
      clear_cache: 'DELETE /clear-cache'
    },
    examples: {
      image_2d: `http://localhost:${PORT}/img/histamine`,
      image_3d: `http://localhost:${PORT}/img/histamine?type=3d`,
      viewer_3d: `http://localhost:${PORT}/3d/histamine`,
      info: `http://localhost:${PORT}/info/histamine`
    },
    note: 'This server provides DIRECT links to PubChem. Images are not cached locally - only CID lookups are cached.'
  });
});

/**
 * GET /health - Health check
 */
app.get('/health', (req, res) => {
  res.json({
    status: 'ok',
    uptime: process.uptime(),
    timestamp: new Date().toISOString(),
    cached_cids: Object.keys(cidCache).length
  });
});

// ============================================================
// ROUTES - PubChem Image API
// ============================================================

/**
 * GET /img/:name - Direct image redirect/info
 * Returns the PubChem image URL or redirects to it
 * Query params:
 *   - type: '2d' or '3d' (default: '2d')
 *   - size: 'small', 'large', or custom like '500x500' (default: 'large')
 *   - redirect: 'true' to redirect to PubChem, 'false' to return JSON (default: 'true')
 */
app.get('/img/:name', async (req, res) => {
  try {
    const name = req.params.name.trim();
    const recordType = req.query.type || '2d';
    const imageSize = req.query.size || 'large';
    const shouldRedirect = req.query.redirect !== 'false';

    console.log('\n' + '='.repeat(70));
    console.log(`[PubChem] GET /img/${name}`);
    console.log(`   Type: ${recordType}, Size: ${imageSize}`);

    // Get CID
    const cid = await getCompoundCid(name);

    if (!cid) {
      console.log(`   Compound not found: ${name}`);
      console.log('='.repeat(70) + '\n');
      return res.status(404).json({
        error: `Compound not found: ${name}`,
        name: name
      });
    }

    console.log(`   Found CID: ${cid}`);

    // Generate PubChem URL
    const imageUrl = getPubChemImageUrl(cid, imageSize, recordType);
    console.log(`   URL: ${imageUrl}`);
    console.log('='.repeat(70) + '\n');

    if (shouldRedirect) {
      // Redirect directly to PubChem image
      res.redirect(imageUrl);
    } else {
      // Return JSON with URL
      res.json({
        cid: cid,
        name: name,
        image_url: imageUrl,
        direct_access: true
      });
    }
  } catch (error) {
    console.error(`   Error: ${error.message}`);
    console.log('='.repeat(70) + '\n');
    res.status(500).json({ error: error.message });
  }
});

// ============================================================
// ROUTES - PubChem 3D Viewer
// ============================================================

/**
 * GET /3d/:name - 3D viewer redirect or embeddable HTML
 * Query params:
 *   - embed: 'true' to return embeddable HTML, 'false' to redirect (default: 'false')
 *   - format: 'json' to return URLs in JSON format
 */
app.get('/3d/:name', async (req, res) => {
  try {
    const name = req.params.name.trim();
    const embed = req.query.embed === 'true';
    const format = req.query.format;

    console.log('\n' + '='.repeat(70));
    console.log(`[PubChem] GET /3d/${name}`);
    console.log(`   Embed: ${embed}, Format: ${format || 'html'}`);

    // Get CID
    const cid = await getCompoundCid(name);

    if (!cid) {
      console.log(`   Compound not found: ${name}`);
      console.log('='.repeat(70) + '\n');
      return res.status(404).json({
        error: `Compound not found: ${name}`,
        name: name
      });
    }

    console.log(`   Found CID: ${cid}`);

    const viewer3dUrl = getPubChem3DViewerUrl(cid);
    const sdfUrl = getPubChemSdfUrl(cid);

    console.log(`   3D Viewer URL: ${viewer3dUrl}`);
    console.log('='.repeat(70) + '\n');

    // JSON format
    if (format === 'json') {
      return res.json({
        cid: cid,
        name: name,
        viewer_3d_url: viewer3dUrl,
        sdf_url: sdfUrl,
        pubchem_url: `https://pubchem.ncbi.nlm.nih.gov/compound/${cid}`
      });
    }

    // Embed format (HTML widget)
    if (embed) {
      const html = `
<!DOCTYPE html>
<html>
<head>
  <title>3D Viewer - ${name}</title>
  <style>
    body {
      font-family: Arial, sans-serif;
      margin: 0;
      padding: 20px;
      background: #f5f5f5;
    }
    .container {
      max-width: 900px;
      margin: 0 auto;
      background: white;
      padding: 20px;
      border-radius: 8px;
      box-shadow: 0 2px 8px rgba(0,0,0,0.1);
    }
    h1 {
      color: #333;
      margin-top: 0;
    }
    .info {
      background: #e7f3ff;
      padding: 15px;
      border-radius: 4px;
      margin: 15px 0;
    }
    .info p {
      margin: 5px 0;
    }
    .button {
      display: inline-block;
      padding: 12px 24px;
      background: #007bff;
      color: white;
      text-decoration: none;
      border-radius: 4px;
      margin: 5px;
      font-weight: 500;
    }
    .button:hover {
      background: #0056b3;
    }
    .button.secondary {
      background: #6c757d;
    }
    .button.secondary:hover {
      background: #545b62;
    }
    .viewer-note {
      background: #fff3cd;
      padding: 15px;
      border-radius: 4px;
      margin: 20px 0;
      border-left: 4px solid #ffc107;
    }
    .features {
      background: #f8f9fa;
      padding: 15px;
      border-radius: 4px;
      margin: 15px 0;
    }
    .features ul {
      margin: 10px 0;
      padding-left: 20px;
    }
    .features li {
      margin: 5px 0;
    }
  </style>
</head>
<body>
  <div class="container">
    <h1>3D Molecular Viewer</h1>

    <div class="info">
      <p><strong>Compound:</strong> ${name}</p>
      <p><strong>PubChem CID:</strong> ${cid}</p>
    </div>

    <div class="features">
      <strong>3D Viewer Features:</strong>
      <ul>
        <li>Interactive 3D rotation (drag to rotate)</li>
        <li>Zoom controls (scroll wheel)</li>
        <li>Multiple display modes (ball-and-stick, space-filling, stick)</li>
        <li>Show/hide hydrogen atoms</li>
        <li>Animated rotation</li>
        <li>Download 3D structure (SDF format)</li>
      </ul>
    </div>

    <div>
      <a href="${viewer3dUrl}" target="_blank" class="button">Open 3D Viewer on PubChem</a>
      <a href="${sdfUrl}" class="button secondary">Download SDF File</a>
    </div>

    <div class="viewer-note">
      <strong>Note:</strong> The interactive 3D viewer will open on PubChem's website.
      It uses WebGL technology for real-time 3D rendering. Make sure your browser supports WebGL.
    </div>
  </div>
</body>
</html>`;
      return res.send(html);
    }

    // Default: redirect to PubChem 3D viewer
    res.redirect(viewer3dUrl);
  } catch (error) {
    console.error(`   Error: ${error.message}`);
    console.log('='.repeat(70) + '\n');
    res.status(500).json({ error: error.message });
  }
});

// ============================================================
// ROUTES - Compound Information
// ============================================================

/**
 * GET /info/:name - Get all URLs and info for a compound
 */
app.get('/info/:name', async (req, res) => {
  try {
    const name = req.params.name.trim();

    console.log('\n' + '='.repeat(70));
    console.log(`[PubChem] GET /info/${name}`);

    // Get CID
    const cid = await getCompoundCid(name);

    if (!cid) {
      console.log(`   Compound not found: ${name}`);
      console.log('='.repeat(70) + '\n');
      return res.status(404).json({
        error: `Compound not found: ${name}`,
        name: name
      });
    }

    console.log(`   Found CID: ${cid}`);
    console.log('='.repeat(70) + '\n');

    res.json({
      cid: cid,
      name: name,
      pubchem_url: `https://pubchem.ncbi.nlm.nih.gov/compound/${cid}`,
      images: {
        '2d_small': getPubChemImageUrl(cid, 'small', '2d'),
        '2d_large': getPubChemImageUrl(cid, 'large', '2d'),
        '3d_small': getPubChemImageUrl(cid, 'small', '3d'),
        '3d_large': getPubChemImageUrl(cid, 'large', '3d')
      },
      viewer_3d: {
        url: getPubChem3DViewerUrl(cid),
        embed_html: `http://localhost:${PORT}/3d/${name}?embed=true`
      },
      data_files: {
        sdf_2d: getPubChemSdfUrl(cid, '2d'),
        sdf_3d: getPubChemSdfUrl(cid, '3d')
      },
      local_endpoints: {
        image_2d: `http://localhost:${PORT}/img/${name}`,
        image_3d: `http://localhost:${PORT}/img/${name}?type=3d`,
        viewer_3d: `http://localhost:${PORT}/3d/${name}`
      }
    });
  } catch (error) {
    console.error(`   Error: ${error.message}`);
    console.log('='.repeat(70) + '\n');
    res.status(500).json({ error: error.message });
  }
});

// ============================================================
// ROUTES - Cache Management
// ============================================================

/**
 * GET /cache-info - Get cache statistics
 */
app.get('/cache-info', (req, res) => {
  res.json({
    cache_directory: CACHE_DIR,
    cached_cid_lookups: Object.keys(cidCache).length,
    cache_file: CID_CACHE_FILE,
    note: 'Only CID lookups are cached. Images are served directly from PubChem.'
  });
});

/**
 * DELETE /clear-cache - Clear CID cache
 */
app.delete('/clear-cache', (req, res) => {
  try {
    const count = Object.keys(cidCache).length;
    cidCache = {};
    saveCidCache();

    console.log(`Cache cleared - ${count} CID lookups removed`);

    res.json({
      message: `Cleared ${count} cached CID lookups`,
      status: 'success'
    });
  } catch (error) {
    res.status(500).json({ error: error.message });
  }
});

// ============================================================
// ROUTES - 3D Viewer with Canvas
// ============================================================

/**
 * GET /viewer-3d/:name - Serve the local 3D viewer HTML
 */
app.get('/viewer-3d/:name', async (req, res) => {
  try {
    const name = req.params.name.trim();
    
    console.log('\n' + '='.repeat(70));
    console.log(`[PubChem] GET /viewer-3d/${name}`);
    
    // Get CID
    const cid = await getCompoundCid(name);
    
    if (!cid) {
      console.log(`   Compound not found: ${name}`);
      console.log('='.repeat(70) + '\n');
      return res.status(404).json({
        error: `Compound not found: ${name}`,
        name: name
      });
    }
    
    console.log(`   Found CID: ${cid}, serving viewer HTML`);
    console.log('='.repeat(70) + '\n');
    
    // Redirect to static viewer HTML with query params
    res.redirect(`/static/viewer-3d.html?name=${encodeURIComponent(name)}&cid=${cid}`);
  } catch (error) {
    res.status(500).json({ error: error.message });
  }
});

/**
 * GET /3d-model/:cidOrName - Get SDF file data
 */
app.get('/3d-model/:cidOrName', async (req, res) => {
  try {
    const cidOrName = req.params.cidOrName.trim();
    const format = req.query.format || 'sdf';
    
    console.log('\n' + '='.repeat(70));
    console.log(`[PubChem] GET /3d-model/${cidOrName}`);
    console.log(`   Format: ${format}`);
    
    // Get CID
    let cid = cidOrName;
    if (isNaN(cidOrName)) {
      cid = await getCompoundCid(cidOrName);
      if (!cid) {
        console.log(`   Compound not found: ${cidOrName}`);
        console.log('='.repeat(70) + '\n');
        return res.status(404).json({
          error: `Compound not found: ${cidOrName}`,
          name: cidOrName
        });
      }
    }
    
    console.log(`   Found CID: ${cid}, fetching SDF data`);
    
    // Fetch SDF from PubChem
    const sdfUrl = getPubChemSdfUrl(cid);
    const response = await axios.get(sdfUrl);
    
    if (response.status !== 200) {
      throw new Error('Failed to fetch SDF from PubChem');
    }
    
    const sdfData = response.data;
    console.log(`   SDF data fetched, length: ${sdfData.length}`);
    console.log('='.repeat(70) + '\n');
    
    // Return based on format
    if (format === 'json') {
      return res.json({
        cid: parseInt(cid),
        name: cidOrName,
        format: 'sdf',
        data: sdfData,
        viewer_url: `http://localhost:${PORT}/viewer-3d/${cidOrName}`
      });
    } else {
      // Return as SDF file
      res.setHeader('Content-Type', 'chemical/x-mdl-sdfile');
      res.setHeader('Content-Disposition', `attachment; filename="${cidOrName}_3d.sdf"`);
      res.send(sdfData);
    }
  } catch (error) {
    console.error('Error fetching SDF:', error.message);
    res.status(500).json({ error: error.message });
  }
});

// ============================================================
// ERROR HANDLING
// ============================================================

app.use((err, req, res, next) => {
  console.error('Server error:', err);
  res.status(500).json({
    error: err.message,
    timestamp: new Date().toISOString()
  });
});

// ============================================================
// START SERVER
// ============================================================

app.listen(PORT, () => {
  console.log('\n' + '='.repeat(70));
  console.log(`PubChem Server running on http://localhost:${PORT}`);
  console.log('='.repeat(70));
  console.log('\nAPI Endpoints:');
  console.log(`   Image (2D):       http://localhost:${PORT}/img/histamine`);
  console.log(`   Image (3D):       http://localhost:${PORT}/img/histamine?type=3d`);
  console.log(`   3D Viewer:        http://localhost:${PORT}/3d/histamine`);
  console.log(`   Compound Info:    http://localhost:${PORT}/info/histamine`);
  console.log(`   Health Check:     http://localhost:${PORT}/health`);
  console.log(`   Cache Info:       http://localhost:${PORT}/cache-info`);
  console.log(`\nCache Directory: ${CACHE_DIR}`);
  console.log(`Cached CID Lookups: ${Object.keys(cidCache).length}`);
  console.log('='.repeat(70) + '\n');
});

// Handle shutdown gracefully
process.on('SIGINT', () => {
  console.log('\nServer shutting down...');
  saveCidCache();
  process.exit(0);
});

process.on('SIGTERM', () => {
  console.log('\nServer shutting down...');
  saveCidCache();
  process.exit(0);
});
