/**
 * Generate static icon files for ChemistryLaTeX extension
 * Run once to create icons that will be included with all builds
 */

const sharp = require('sharp');
const fs = require('fs');
const path = require('path');

const assetsDir = path.join(__dirname, 'assets');

// Create assets directory if it doesn't exist
if (!fs.existsSync(assetsDir)) {
    fs.mkdirSync(assetsDir, { recursive: true });
}

// Blue beaker SVG - visible on both light and dark backgrounds
const beakerSvg = (size) => `<svg xmlns="http://www.w3.org/2000/svg" width="${size}" height="${size}" viewBox="0 0 24 24">
  <path fill="#2196F3" d="M3 3h2v11c0 2.21 1.79 4 4 4s4-1.79 4-4V3h2c.55 0 1-.45 1-1s-.45-1-1-1H3c-.55 0-1 .45-1 1s.45 1 1 1m4 11V3h4v11c0 1.1-.9 2-2 2s-2-.9-2-2m11.5 0c-.28 0-.5.22-.5.5s.22.5.5.5.5-.22.5-.5-.22-.5-.5-.5m-1 3c-.28 0-.5.22-.5.5s.22.5.5.5.5-.22.5-.5-.22-.5-.5-.5m3-3c-.28 0-.5.22-.5.5s.22.5.5.5.5-.22.5-.5-.22-.5-.5-.5m-1.38 5.38c-1.34 1.34-3.5 1.34-4.85 0l4.85-4.85c1.35 1.34 1.35 3.51 0 4.85z"/>
</svg>`;

async function generateIcons() {
    console.log('🎨 Generating static icon files...\n');

    const sizes = [16, 48, 128];

    for (const size of sizes) {
        const pngPath = path.join(assetsDir, `icon${size}.png`);
        await sharp(Buffer.from(beakerSvg(size)))
            .resize(size, size)
            .png()
            .toFile(pngPath);

        const stats = fs.statSync(pngPath);
        console.log(`   ✅ icon${size}.png (${stats.size} bytes)`);
    }

    console.log('\n✅ Static icons generated in assets/ folder');
    console.log('   These icons will be copied to all builds (Chrome, Dev, Firefox)');
}

generateIcons().catch(console.error);
