/**
 * Vercel serverless function for SMILES to SVG conversion
 * Endpoint: /api/img/smiles?smiles=CCO&width=300&height=200
 */

const { spawn } = require('child_process');
const path = require('path');

/**
 * Generate error SVG
 */
function createErrorSvg(message, width = 300, height = 200) {
  const escapeXml = (str) => String(str)
    .replace(/&/g, '&amp;')
    .replace(/</g, '&lt;')
    .replace(/>/g, '&gt;')
    .replace(/"/g, '&quot;')
    .replace(/'/g, '&#39;');

  return `<svg width="${width}" height="${height}" xmlns="http://www.w3.org/2000/svg">
    <rect width="${width}" height="${height}" fill="#f0f0f0"/>
    <text x="${width / 2}" y="50" text-anchor="middle" font-family="monospace" font-size="14" fill="#d00" font-weight="bold">Error</text>
    <text x="${width / 2}" y="80" text-anchor="middle" font-family="monospace" font-size="11" fill="#666">${escapeXml(message.substring(0, 60))}</text>
  </svg>`;
}

/**
 * Generate simple SVG using basic drawing (fallback when RDKit unavailable)
 */
function generateFallbackSvg(smiles, width = 300, height = 200) {
  return `<svg width="${width}" height="${height}" xmlns="http://www.w3.org/2000/svg">
    <rect width="${width}" height="${height}" fill="#ffffff"/>
    <text x="${width / 2}" y="${height / 2 - 20}" text-anchor="middle" font-family="Arial" font-size="14" fill="#333">ChemParser</text>
    <text x="${width / 2}" y="${height / 2}" text-anchor="middle" font-family="monospace" font-size="12" fill="#666">${smiles}</text>
    <text x="${width / 2}" y="${height / 2 + 20}" text-anchor="middle" font-family="Arial" font-size="10" fill="#999">RDKit not available - showing SMILES</text>
  </svg>`;
}

/**
 * Try to generate SVG via Python RDKit (if available)
 */
async function tryGenerateSvgViaPython(smiles, width, height) {
  return new Promise((resolve, reject) => {
    const pythonScript = path.join(process.cwd(), 'MoleculeViewer', 'generate_svg.py');

    const pythonProcess = spawn('python3', [
      pythonScript,
      JSON.stringify({
        smiles: smiles,
        width: width,
        height: height,
        options: {
          aromaticCircles: true,
          fancyBonds: true
        }
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

    // Timeout after 8 seconds (Vercel has 10s limit)
    setTimeout(() => {
      pythonProcess.kill();
      reject(new Error('SVG generation timeout'));
    }, 8000);
  });
}

module.exports = async (req, res) => {
  res.setHeader('Access-Control-Allow-Origin', '*');
  res.setHeader('Access-Control-Allow-Methods', 'GET');
  res.setHeader('Content-Type', 'image/svg+xml');
  res.setHeader('Cache-Control', 'public, max-age=86400');

  const smiles = req.query.smiles?.trim();
  const width = parseInt(req.query.width) || 300;
  const height = parseInt(req.query.height) || 200;

  if (!smiles) {
    return res.status(400).send(createErrorSvg('SMILES required', width, height));
  }

  try {
    // Try to generate via Python/RDKit
    const svg = await tryGenerateSvgViaPython(smiles, width, height);
    res.status(200).send(svg);
  } catch (error) {
    console.error('Python generation failed:', error.message);
    // Fallback to simple SVG
    const fallbackSvg = generateFallbackSvg(smiles, width, height);
    res.status(200).send(fallbackSvg);
  }
};
