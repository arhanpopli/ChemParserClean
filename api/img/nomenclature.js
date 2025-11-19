/**
 * Vercel serverless function for nomenclature to SVG conversion
 * Endpoint: /api/img/nomenclature?nomenclature=ethanol&width=300&height=200
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
 * Generate simple SVG showing nomenclature
 */
function generateFallbackSvg(nomenclature, width = 300, height = 200) {
  return `<svg width="${width}" height="${height}" xmlns="http://www.w3.org/2000/svg">
    <rect width="${width}" height="${height}" fill="#ffffff"/>
    <text x="${width / 2}" y="${height / 2 - 20}" text-anchor="middle" font-family="Arial" font-size="14" fill="#333">ChemParser</text>
    <text x="${width / 2}" y="${height / 2}" text-anchor="middle" font-family="monospace" font-size="12" fill="#666">${nomenclature}</text>
    <text x="${width / 2}" y="${height / 2 + 20}" text-anchor="middle" font-family="Arial" font-size="10" fill="#999">Conversion service not available</text>
  </svg>`;
}

module.exports = async (req, res) => {
  res.setHeader('Access-Control-Allow-Origin', '*');
  res.setHeader('Access-Control-Allow-Methods', 'GET');
  res.setHeader('Content-Type', 'image/svg+xml');
  res.setHeader('Cache-Control', 'public, max-age=86400');

  const nomenclature = req.query.nomenclature?.trim();
  const width = parseInt(req.query.width) || 300;
  const height = parseInt(req.query.height) || 200;

  if (!nomenclature) {
    return res.status(400).send(createErrorSvg('Nomenclature required', width, height));
  }

  try {
    // For now, return fallback SVG
    // In production, you would integrate with OPSIN or similar service
    const fallbackSvg = generateFallbackSvg(nomenclature, width, height);
    res.status(200).send(fallbackSvg);
  } catch (error) {
    console.error('Error:', error.message);
    res.status(500).send(createErrorSvg(error.message, width, height));
  }
};
