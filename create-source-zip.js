/**
 * Create source ZIP with proper forward slashes for Mozilla
 */
const fs = require('fs');
const path = require('path');
const archiver = require('archiver');

const sourceDir = path.join(__dirname, 'ChemistryLaTeX-source');
const outputPath = path.join(__dirname, 'ChemistryLaTeX-source.zip');

// Delete existing ZIP
if (fs.existsSync(outputPath)) {
    fs.unlinkSync(outputPath);
}

const output = fs.createWriteStream(outputPath);
const archive = archiver('zip', { zlib: { level: 9 } });

output.on('close', () => {
    console.log(`✅ Created ${outputPath}`);
    console.log(`   Size: ${(archive.pointer() / 1024).toFixed(1)} KB`);
});

archive.on('error', (err) => { throw err; });
archive.pipe(output);

// Files to add (without node_modules or dist)
const filesToAdd = [
    'README.md',
    'USAGE.md',
    'package.json',
    'package-lock.json',
    'build-firefox.js',
    'content.js',
    'popup.js',
    'background.js',
    'size-controls.js',
    'popup.html',
    'manifest.json',
];

for (const file of filesToAdd) {
    const filePath = path.join(sourceDir, file);
    if (fs.existsSync(filePath)) {
        archive.file(filePath, { name: file });
        console.log(`   + ${file}`);
    }
}

// Add directories
const dirs = ['assets', 'physics'];
for (const dir of dirs) {
    const dirPath = path.join(sourceDir, dir);
    if (fs.existsSync(dirPath)) {
        archive.directory(dirPath, dir);
        console.log(`   + ${dir}/`);
    }
}

archive.finalize();
