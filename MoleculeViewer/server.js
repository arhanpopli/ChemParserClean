/**
 * MoleculeViewer Node.js Server
 * Hosts molecule visualization as a URL-based image service (like CodeCogs)
 * Generates SVG images from SMILES and chemical nomenclature
 */

const express = require('express');
const cors = require('cors');
const axios = require('axios');
const fs = require('fs');
const path = require('path');
const { spawn } = require('child_process');

const app = express();
const PORT = 5000;

// Middleware
app.use(cors());
app.use(express.json({ limit: '10mb' }));
app.use(express.urlencoded({ limit: '10mb' }));

// Serve static files from root directory
app.use(express.static(path.join(__dirname, '..')));
app.use(express.static(path.join(__dirname)));

// Cache directory for generated SVGs - MoleculeViewer specific
const CACHE_DIR = path.join(__dirname, 'cache', 'moleculeviewer');
if (!fs.existsSync(CACHE_DIR)) {
  fs.mkdirSync(CACHE_DIR, { recursive: true });
}

console.log('🚀 MoleculeViewer Node.js Server Starting...');
console.log(`📁 Cache directory: ${CACHE_DIR}`);

// ============================================================
// HELPER FUNCTIONS
// ============================================================

/**
 * Call Python chemistry module to generate SVG
 */
async function generateSvgViaPython(smiles, options = {}) {
  return new Promise((resolve, reject) => {
    const pythonProcess = spawn('python', [
      path.join(__dirname, 'generate_svg.py'),
      JSON.stringify({
        smiles: smiles,
        options: options
      })
    ]);

    let stdout = '';
    let stderr = '';

    pythonProcess.stdout.on('data', (data) => {
      stdout += data.toString();
    });

    pythonProcess.stderr.on('data', (data) => {
      stderr += data.toString();
    });

    pythonProcess.on('close', (code) => {
      if (code !== 0) {
        console.error(`❌ Python error: ${stderr}`);
        reject(new Error(stderr || 'SVG generation failed'));
        return;
      }

      try {
        const result = JSON.parse(stdout);
        if (result.error) {
          reject(new Error(result.error));
        } else {
          resolve(result.svg);
        }
      } catch (e) {
        reject(new Error(`Failed to parse SVG output: ${e.message}`));
      }
    });
  });
}

/**
 * Convert nomenclature to SMILES using Python
 */
async function convertNomenclatureToSmiles(nomenclature) {
  return new Promise((resolve, reject) => {
    const pythonProcess = spawn('python', [
      path.join(__dirname, 'nomenclature_to_smiles.py'),
      nomenclature
    ]);

    let stdout = '';
    let stderr = '';

    pythonProcess.stdout.on('data', (data) => {
      stdout += data.toString();
    });

    pythonProcess.stderr.on('data', (data) => {
      stderr += data.toString();
    });

    pythonProcess.on('close', (code) => {
      if (code !== 0) {
        console.error(`❌ Nomenclature conversion error: ${stderr}`);
        reject(new Error(stderr || 'Nomenclature conversion failed'));
        return;
      }

      try {
        const result = JSON.parse(stdout);
        if (result.error) {
          reject(new Error(result.error));
        } else {
          resolve(result.smiles);
        }
      } catch (e) {
        reject(new Error(`Failed to parse conversion output: ${e.message}`));
      }
    });
  });
}

/**
 * Canonicalize SMILES using Python/RDKit
 */
async function canonicalizeSmiles(smiles) {
  return new Promise((resolve, reject) => {
    const pythonProcess = spawn('python', [
      path.join(__dirname, 'canonicalize_smiles.py'),
      smiles
    ]);

    let stdout = '';
    let stderr = '';

    pythonProcess.stdout.on('data', (data) => {
      stdout += data.toString();
    });

    pythonProcess.stderr.on('data', (data) => {
      stderr += data.toString();
    });

    pythonProcess.on('close', (code) => {
      if (code !== 0) {
        console.error(`❌ SMILES canonicalization error: ${stderr}`);
        reject(new Error(stderr || 'SMILES canonicalization failed'));
        return;
      }

      try {
        const result = JSON.parse(stdout);
        if (result.error) {
          reject(new Error(result.error));
        } else {
          resolve(result.canonical_smiles);
        }
      } catch (e) {
        reject(new Error(`Failed to parse canonicalization output: ${e.message}`));
      }
    });
  });
}

/**
 * Generate cache key for SVG
 * NOTE: For SMILES-based keys, use generateCacheKeyFromSmiles instead
 */
function generateCacheKey(type, value) {
  const crypto = require('crypto');
  const key = `${type}:${value}`;
  return crypto.createHash('md5').update(key).digest('hex') + '.svg';
}

/**
 * Generate cache key from canonical SMILES
 * This ensures the same molecule always gets the same cache key
 */
function generateCacheKeyFromSmiles(canonicalSmiles) {
  const crypto = require('crypto');
  const key = `smiles:${canonicalSmiles}`;
  return crypto.createHash('md5').update(key).digest('hex') + '.svg';
}

/**
 * Check if SVG exists in cache
 */
function getCachedSvg(cacheKey) {
  const filePath = path.join(CACHE_DIR, cacheKey);
  if (fs.existsSync(filePath)) {
    return fs.readFileSync(filePath, 'utf-8');
  }
  return null;
}

/**
 * Cache SVG file
 */
function cacheSvg(cacheKey, svgContent) {
  const filePath = path.join(CACHE_DIR, cacheKey);
  fs.writeFileSync(filePath, svgContent);
}

// ============================================================
// ROUTES - Direct Image URLs (Like CodeCogs)
// ============================================================

/**
 * GET /img/smiles - Render SMILES as SVG image
 * Usage: http://localhost:5000/img/smiles?smiles=CCO&width=300&height=200
 */
app.get('/img/smiles', async (req, res) => {
  try {
    const smiles = req.query.smiles?.trim();
    const wantJson = req.query.json === 'true';

    if (!smiles) {
      if (wantJson) {
        return res.json({ success: false, error: 'SMILES required' });
      }
      return res.status(400).send(createErrorSvg('SMILES required'));
    }

    console.log(`\n${'='.repeat(70)}`);
    console.log(`📥 [MoleculeViewer] GET /img/smiles`);
    console.log(`   SMILES: ${smiles}`);
    console.log(`   JSON mode: ${wantJson}`);

    // Canonicalize SMILES to prevent duplicate cache entries
    let canonicalSmiles;
    try {
      canonicalSmiles = await canonicalizeSmiles(smiles);
      console.log(`   Canonical SMILES: ${canonicalSmiles}`);
    } catch (e) {
      console.log(`   Warning: Could not canonicalize SMILES, using original: ${e.message}`);
      canonicalSmiles = smiles;
    }

    // Check cache using canonical SMILES
    const cacheKey = generateCacheKeyFromSmiles(canonicalSmiles);
    let svg = getCachedSvg(cacheKey);

    if (svg) {
      console.log(`✅ Served from cache: ${cacheKey}`);
      if (wantJson) {
        return res.json({
          success: true,
          svg: svg,
          cache_url: `http://localhost:${PORT}/cache/moleculeviewer/${cacheKey}`,
          cached: true
        });
      }
      res.set('Content-Type', 'image/svg+xml');
      res.set('Cache-Control', 'public, max-age=86400');
      res.set('Access-Control-Allow-Origin', '*');
      return res.send(svg);
    }

    // Generate SVG using canonical SMILES
    svg = await generateSvgViaPython(canonicalSmiles, {
      aromaticCircles: true,
      fancyBonds: true
    });

    // Cache it
    cacheSvg(cacheKey, svg);

    console.log(`✅ Generated and cached: ${cacheKey}`);
    console.log(`${'='.repeat(70)}\n`);

    if (wantJson) {
      return res.json({
        success: true,
        svg: svg,
        cache_url: `http://localhost:${PORT}/cache/moleculeviewer/${cacheKey}`,
        cached: false
      });
    }

    res.set('Content-Type', 'image/svg+xml');
    res.set('Cache-Control', 'public, max-age=86400');
    res.set('Access-Control-Allow-Origin', '*');
    res.send(svg);
  } catch (error) {
    console.error(`❌ Error: ${error.message}`);
    if (req.query.json === 'true') {
      return res.json({ success: false, error: error.message });
    }
    res.set('Content-Type', 'image/svg+xml');
    res.send(createErrorSvg(error.message));
  }
});

/**
 * GET /img/nomenclature - Render chemical name as SVG image
 * Usage: http://localhost:5000/img/nomenclature?nomenclature=acetone&width=300&height=200
 * Add json=true for JSON response with SVG data
 */
app.get('/img/nomenclature', async (req, res) => {
  try {
    const nomenclature = req.query.nomenclature?.trim();
    const wantJson = req.query.json === 'true';

    if (!nomenclature) {
      if (wantJson) {
        return res.json({ success: false, error: 'Nomenclature required' });
      }
      return res.status(400).send(createErrorSvg('Nomenclature required'));
    }

    console.log(`\n${'='.repeat(70)}`);
    console.log(`📥 [MoleculeViewer] GET /img/nomenclature`);
    console.log(`   Nomenclature: ${nomenclature}`);
    console.log(`   JSON mode: ${wantJson}`);

    // Step 1: Convert nomenclature to SMILES
    console.log(`   Converting to SMILES...`);
    const smiles = await convertNomenclatureToSmiles(nomenclature);
    console.log(`   ✓ SMILES: ${smiles}`);

    // Step 2: Canonicalize SMILES to prevent duplicate cache entries
    let canonicalSmiles;
    try {
      canonicalSmiles = await canonicalizeSmiles(smiles);
      console.log(`   ✓ Canonical SMILES: ${canonicalSmiles}`);
    } catch (e) {
      console.log(`   Warning: Could not canonicalize SMILES, using original: ${e.message}`);
      canonicalSmiles = smiles;
    }

    // Step 3: Check cache using canonical SMILES
    const cacheKey = generateCacheKeyFromSmiles(canonicalSmiles);
    let svg = getCachedSvg(cacheKey);

    if (svg) {
      console.log(`✅ Served from cache: ${cacheKey}`);
      if (wantJson) {
        return res.json({
          success: true,
          svg: svg,
          cache_url: `http://localhost:${PORT}/cache/moleculeviewer/${cacheKey}`,
          cached: true
        });
      }
      res.set('Content-Type', 'image/svg+xml');
      res.set('Cache-Control', 'public, max-age=86400');
      res.set('Access-Control-Allow-Origin', '*');
      return res.send(svg);
    }

    // Step 4: Generate SVG from canonical SMILES
    svg = await generateSvgViaPython(canonicalSmiles, {
      aromaticCircles: true,
      fancyBonds: true
    });

    // Cache it
    cacheSvg(cacheKey, svg);

    console.log(`✅ Generated and cached: ${cacheKey}`);
    console.log(`${'='.repeat(70)}\n`);

    if (wantJson) {
      return res.json({
        success: true,
        svg: svg,
        cache_url: `http://localhost:${PORT}/cache/moleculeviewer/${cacheKey}`,
        cached: false
      });
    }

    res.set('Content-Type', 'image/svg+xml');
    res.set('Cache-Control', 'public, max-age=86400');
    res.set('Access-Control-Allow-Origin', '*');
    res.send(svg);
  } catch (error) {
    console.error(`❌ Error: ${error.message}`);
    if (req.query.json === 'true') {
      return res.json({ success: false, error: error.message });
    }
    res.set('Content-Type', 'image/svg+xml');
    res.send(createErrorSvg(`Not found: ${error.message}`));
  }
});

// ============================================================
// ROUTES - Helper Routes
// ============================================================

/**
 * GET / - Info page
 */
app.get('/', (req, res) => {
  res.json({
    name: 'MoleculeViewer Node.js Server',
    version: '1.0.0',
    endpoints: {
      smiles: 'GET /img/smiles?smiles=CCO&width=300&height=200',
      nomenclature: 'GET /img/nomenclature?nomenclature=acetone&width=300&height=200',
      health: 'GET /health'
    },
    examples: {
      smiles: 'http://localhost:5000/img/smiles?smiles=c1ccccc1',
      nomenclature: 'http://localhost:5000/img/nomenclature?nomenclature=benzene'
    }
  });
});

/**
 * POST /api/nomenclature-to-smiles - Convert chemical name to SMILES
 * Used by mol2chemfig extension to get SMILES for rendering
 */
app.post('/api/nomenclature-to-smiles', async (req, res) => {
  try {
    const { nomenclature } = req.body;

    if (!nomenclature) {
      return res.status(400).json({ error: 'Nomenclature required' });
    }

    console.log(`\n${'='.repeat(70)}`);
    console.log(`📥 [MoleculeViewer] POST /api/nomenclature-to-smiles`);
    console.log(`   Nomenclature: ${nomenclature}`);

    const smiles = await convertNomenclatureToSmiles(nomenclature);

    console.log(`   ✓ SMILES: ${smiles}`);
    console.log(`${'='.repeat(70)}\n`);

    res.json({
      success: true,
      nomenclature: nomenclature,
      smiles: smiles,
      source: 'MoleculeViewer/PubChem'
    });
  } catch (error) {
    console.error(`❌ Nomenclature conversion error: ${error.message}`);
    res.status(404).json({
      success: false,
      error: error.message,
      nomenclature: req.body?.nomenclature
    });
  }
});

/**
 * GET /health - Health check
 */
app.get('/health', (req, res) => {
  res.json({
    status: 'ok',
    uptime: process.uptime(),
    timestamp: new Date().toISOString()
  });
});

/**
 * GET /cache-info - Show cache statistics
 */
app.get('/cache-info', (req, res) => {
  const files = fs.readdirSync(CACHE_DIR);
  const totalSize = files.reduce((sum, file) => {
    const filePath = path.join(CACHE_DIR, file);
    return sum + fs.statSync(filePath).size;
  }, 0);

  res.json({
    cacheDirectory: CACHE_DIR,
    cachedSvgs: files.length,
    totalCacheSize: `${(totalSize / 1024).toFixed(2)} KB`,
    files: files.slice(0, 10) // Show first 10
  });
});

/**
 * DELETE /clear-cache - Clear all cached SVGs
 */
app.delete('/clear-cache', (req, res) => {
  try {
    const files = fs.readdirSync(CACHE_DIR);
    files.forEach(file => {
      fs.unlinkSync(path.join(CACHE_DIR, file));
    });
    res.json({
      message: `Cleared ${files.length} cached SVGs`,
      cacheSize: '0 KB'
    });
    console.log(`🗑️  Cache cleared - ${files.length} files deleted`);
  } catch (error) {
    res.status(500).json({ error: error.message });
  }
});

// ============================================================
// UTILITY FUNCTIONS
// ============================================================

/**
 * Create error SVG for display
 */
function createErrorSvg(message) {
  const width = 300;
  const height = 200;
  return `<svg width="${width}" height="${height}" xmlns="http://www.w3.org/2000/svg">
    <rect width="${width}" height="${height}" fill="#f0f0f0"/>
    <text x="${width / 2}" y="50" text-anchor="middle" font-family="monospace" font-size="14" fill="#d00" font-weight="bold">Error</text>
    <text x="${width / 2}" y="80" text-anchor="middle" font-family="monospace" font-size="11" fill="#666">${escapeXml(message.substring(0, 60))}</text>
  </svg>`;
}

/**
 * Escape XML special characters
 */
function escapeXml(str) {
  return String(str)
    .replace(/&/g, '&amp;')
    .replace(/</g, '&lt;')
    .replace(/>/g, '&gt;')
    .replace(/"/g, '&quot;')
    .replace(/'/g, '&#39;');
}

// ============================================================
// ERROR HANDLING
// ============================================================

app.use((err, req, res, next) => {
  console.error('❌ Server error:', err);
  res.status(500).json({
    error: err.message,
    timestamp: new Date().toISOString()
  });
});

// ============================================================
// START SERVER
// ============================================================

app.listen(PORT, () => {
  console.log(`\n${'='.repeat(70)}`);
  console.log(`✅ MoleculeViewer Server running on http://localhost:${PORT}`);
  console.log(`${'='.repeat(70)}`);
  console.log('\n📍 API Endpoints:');
  console.log(`   SMILES:       http://localhost:${PORT}/img/smiles?smiles=CCO`);
  console.log(`   Nomenclature: http://localhost:${PORT}/img/nomenclature?nomenclature=acetone`);
  console.log(`   Health:       http://localhost:${PORT}/health`);
  console.log(`   Cache Info:   http://localhost:${PORT}/cache-info`);
  console.log(`\n💾 Cache Directory: ${CACHE_DIR}`);
  console.log(`${'='.repeat(70)}\n`);
});

// Handle shutdown gracefully
process.on('SIGINT', () => {
  console.log('\n🛑 Server shutting down...');
  process.exit(0);
});
