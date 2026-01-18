/**
 * xTB Orbital Service
 * 
 * Calculates molecular orbitals using xTB and returns cube files.
 * Deployed on Railway, called by Vercel API.
 */

const express = require('express');
const cors = require('cors');
const { execSync, exec } = require('child_process');
const fs = require('fs');
const path = require('path');
const crypto = require('crypto');

const app = express();
app.use(cors());
app.use(express.json());

const PORT = process.env.PORT || 3000;
const WORK_DIR = '/tmp/xtb-work';

// Ensure work directory exists
if (!fs.existsSync(WORK_DIR)) {
    fs.mkdirSync(WORK_DIR, { recursive: true });
}

// Simple in-memory cache
const cache = new Map();
const MAX_CACHE_SIZE = 100;

/**
 * Generate a unique job ID from SMILES
 */
function getJobId(smiles) {
    return crypto.createHash('md5').update(smiles.toLowerCase()).digest('hex').slice(0, 12);
}

/**
 * Convert SMILES to XYZ using Open Babel
 */
function smilesToXYZ(smiles, outputPath) {
    try {
        // Use Open Babel to convert SMILES to 3D XYZ
        execSync(`echo "${smiles}" | obabel -ismi -oxyz --gen3d -O "${outputPath}"`, {
            timeout: 30000,
            stdio: 'pipe'
        });
        return true;
    } catch (error) {
        console.error('SMILES to XYZ conversion failed:', error.message);
        return false;
    }
}

/**
 * Run xTB calculation
 */
function runXTB(xyzPath, workDir) {
    try {
        execSync(`cd "${workDir}" && xtb "${xyzPath}" --gfn 2 --sp --molden`, {
            timeout: 120000, // 2 minute timeout
            stdio: 'pipe'
        });
        return true;
    } catch (error) {
        console.error('xTB calculation failed:', error.message);
        return false;
    }
}

/**
 * Parse Molden file and generate cube data
 */
function parseMoldenAndGenerateCube(moldenPath, orbitalType, gridSize) {
    const content = fs.readFileSync(moldenPath, 'utf8');
    const lines = content.split('\n');

    // Parse atoms
    const atoms = [];
    let section = null;

    for (let i = 0; i < lines.length; i++) {
        const line = lines[i].trim();
        if (line.startsWith('[')) {
            section = line.toLowerCase();
            continue;
        }
        if (!line) continue;

        if (section && section.includes('[atoms]')) {
            const parts = line.split(/\s+/);
            if (parts.length >= 6 && !isNaN(parseInt(parts[1]))) {
                const isBohr = section.includes('au');
                const scale = isBohr ? 0.529177 : 1.0;
                atoms.push({
                    element: parts[0],
                    atomicNumber: parseInt(parts[2]),
                    x: parseFloat(parts[3]) * scale,
                    y: parseFloat(parts[4]) * scale,
                    z: parseFloat(parts[5]) * scale
                });
            }
        }
    }

    // Parse orbital energies and find HOMO/LUMO
    const orbitals = [];
    let currentOrbital = null;

    for (let i = 0; i < lines.length; i++) {
        const line = lines[i].trim();
        if (line.startsWith('Sym=')) {
            if (currentOrbital) orbitals.push(currentOrbital);
            currentOrbital = { symmetry: line.split('=')[1], coefficients: [] };
        }
        if (line.startsWith('Ene=') && currentOrbital) {
            currentOrbital.energy = parseFloat(line.split('=')[1]);
        }
        if (line.startsWith('Occup=') && currentOrbital) {
            currentOrbital.occupation = parseFloat(line.split('=')[1]);
        }
        const coeffMatch = line.match(/^\s*(\d+)\s+([\d.E+-]+)/i);
        if (coeffMatch && currentOrbital) {
            currentOrbital.coefficients[parseInt(coeffMatch[1]) - 1] = parseFloat(coeffMatch[2]);
        }
    }
    if (currentOrbital) orbitals.push(currentOrbital);

    // Find the requested orbital
    const occupied = orbitals.filter(o => o.occupation > 0.5).sort((a, b) => b.energy - a.energy);
    const virtual = orbitals.filter(o => o.occupation <= 0.5).sort((a, b) => a.energy - b.energy);

    let targetOrbital;
    switch (orbitalType.toLowerCase()) {
        case 'homo': targetOrbital = occupied[0]; break;
        case 'homo-1': targetOrbital = occupied[1]; break;
        case 'homo-2': targetOrbital = occupied[2]; break;
        case 'lumo': targetOrbital = virtual[0]; break;
        case 'lumo+1': targetOrbital = virtual[1]; break;
        case 'lumo+2': targetOrbital = virtual[2]; break;
        default: targetOrbital = occupied[0];
    }

    if (!targetOrbital) {
        throw new Error(`Orbital ${orbitalType} not found`);
    }

    // Generate a simple cube file based on atom positions
    // This is a simplified version - for full accuracy we'd need basis functions
    const cubeFile = generateSimpleCube(atoms, targetOrbital, gridSize);

    return {
        atoms,
        orbitalEnergy: targetOrbital.energy,
        orbitalOccupation: targetOrbital.occupation,
        cubeFile
    };
}

/**
 * Generate a simple cube file (simplified orbital representation)
 */
function generateSimpleCube(atoms, orbital, gridSize = 50) {
    // Calculate bounding box
    let minX = Infinity, maxX = -Infinity;
    let minY = Infinity, maxY = -Infinity;
    let minZ = Infinity, maxZ = -Infinity;

    for (const atom of atoms) {
        minX = Math.min(minX, atom.x); maxX = Math.max(maxX, atom.x);
        minY = Math.min(minY, atom.y); maxY = Math.max(maxY, atom.y);
        minZ = Math.min(minZ, atom.z); maxZ = Math.max(maxZ, atom.z);
    }

    // Add padding
    const padding = 3.0;
    minX -= padding; maxX += padding;
    minY -= padding; maxY += padding;
    minZ -= padding; maxZ += padding;

    // Step sizes (convert to Bohr)
    const ANG_TO_BOHR = 1.889726;
    const dx = (maxX - minX) / (gridSize - 1);
    const dy = (maxY - minY) / (gridSize - 1);
    const dz = (maxZ - minZ) / (gridSize - 1);

    // Generate simplified orbital density
    // This uses a simple Gaussian approximation for visualization
    const data = [];
    let maxVal = 0;

    for (let iz = 0; iz < gridSize; iz++) {
        for (let iy = 0; iy < gridSize; iy++) {
            for (let ix = 0; ix < gridSize; ix++) {
                const x = minX + ix * dx;
                const y = minY + iy * dy;
                const z = minZ + iz * dz;

                let val = 0;
                for (let a = 0; a < atoms.length && a < orbital.coefficients.length; a++) {
                    const atom = atoms[a];
                    const coeff = orbital.coefficients[a] || 0;
                    const r2 = (x - atom.x) ** 2 + (y - atom.y) ** 2 + (z - atom.z) ** 2;
                    val += coeff * Math.exp(-0.5 * r2);
                }
                data.push(val);
                maxVal = Math.max(maxVal, Math.abs(val));
            }
        }
    }

    // Build cube file
    const lines = [];
    lines.push('Generated by xTB Orbital Service');
    lines.push(`Orbital: ${orbital.symmetry} E=${orbital.energy?.toFixed(6)} Ha`);

    // Origin (in Bohr)
    const originBohr = [minX * ANG_TO_BOHR, minY * ANG_TO_BOHR, minZ * ANG_TO_BOHR];
    lines.push(`${atoms.length} ${originBohr[0].toFixed(6)} ${originBohr[1].toFixed(6)} ${originBohr[2].toFixed(6)}`);

    // Grid dimensions
    lines.push(`${gridSize} ${(dx * ANG_TO_BOHR).toFixed(6)} 0.000000 0.000000`);
    lines.push(`${gridSize} 0.000000 ${(dy * ANG_TO_BOHR).toFixed(6)} 0.000000`);
    lines.push(`${gridSize} 0.000000 0.000000 ${(dz * ANG_TO_BOHR).toFixed(6)}`);

    // Atoms
    for (const atom of atoms) {
        const posBohr = [atom.x * ANG_TO_BOHR, atom.y * ANG_TO_BOHR, atom.z * ANG_TO_BOHR];
        lines.push(`${atom.atomicNumber} ${atom.atomicNumber}.000000 ${posBohr[0].toFixed(6)} ${posBohr[1].toFixed(6)} ${posBohr[2].toFixed(6)}`);
    }

    // Volume data
    let dataLine = '';
    for (let i = 0; i < data.length; i++) {
        dataLine += ` ${data[i].toExponential(5)}`;
        if ((i + 1) % 6 === 0) {
            lines.push(dataLine.trim());
            dataLine = '';
        }
    }
    if (dataLine.trim()) lines.push(dataLine.trim());

    return lines.join('\n');
}

/**
 * Health check endpoint
 */
app.get('/', (req, res) => {
    res.json({
        service: 'xTB Orbital Service',
        status: 'running',
        endpoints: {
            orbital: 'GET /orbital?smiles=CCO&orbital=homo'
        }
    });
});

/**
 * Main orbital calculation endpoint
 */
app.get('/orbital', async (req, res) => {
    const { smiles, orbital = 'homo', gridSize = 50 } = req.query;

    if (!smiles) {
        return res.status(400).json({ error: 'Missing smiles parameter' });
    }

    const jobId = getJobId(smiles);
    const cacheKey = `${jobId}_${orbital}_${gridSize}`;

    // Check cache
    if (cache.has(cacheKey)) {
        console.log(`Cache hit: ${cacheKey}`);
        return res.json({ ...cache.get(cacheKey), cached: true });
    }

    const jobDir = path.join(WORK_DIR, jobId);

    try {
        // Create job directory
        if (!fs.existsSync(jobDir)) {
            fs.mkdirSync(jobDir, { recursive: true });
        }

        const xyzPath = path.join(jobDir, 'molecule.xyz');
        const moldenPath = path.join(jobDir, 'molden.input');

        console.log(`Processing: ${smiles} (${orbital})`);

        // Step 1: Convert SMILES to XYZ
        if (!smilesToXYZ(smiles, xyzPath)) {
            throw new Error('Failed to convert SMILES to XYZ');
        }

        // Step 2: Run xTB
        if (!runXTB(xyzPath, jobDir)) {
            throw new Error('xTB calculation failed');
        }

        // Step 3: Parse Molden and generate cube
        if (!fs.existsSync(moldenPath)) {
            throw new Error('Molden file not generated');
        }

        const result = parseMoldenAndGenerateCube(moldenPath, orbital, parseInt(gridSize));

        const response = {
            success: true,
            smiles,
            orbital,
            gridSize: parseInt(gridSize),
            atoms: result.atoms,
            orbitalEnergy: result.orbitalEnergy,
            cubeFile: result.cubeFile
        };

        // Cache result
        if (cache.size >= MAX_CACHE_SIZE) {
            const firstKey = cache.keys().next().value;
            cache.delete(firstKey);
        }
        cache.set(cacheKey, response);

        console.log(`Completed: ${smiles}`);
        res.json(response);

    } catch (error) {
        console.error(`Error processing ${smiles}:`, error);
        res.status(500).json({
            error: error.message,
            smiles
        });
    } finally {
        // Cleanup job directory
        try {
            if (fs.existsSync(jobDir)) {
                fs.rmSync(jobDir, { recursive: true, force: true });
            }
        } catch (e) {
            // Ignore cleanup errors
        }
    }
});

app.listen(PORT, () => {
    console.log(`xTB Orbital Service running on port ${PORT}`);

    // Verify xTB is installed
    try {
        const version = execSync('xtb --version 2>&1').toString();
        console.log('xTB version:', version.split('\n')[0]);
    } catch (e) {
        console.warn('Warning: xTB not found in PATH');
    }

    // Verify Open Babel is installed
    try {
        const version = execSync('obabel -V 2>&1').toString();
        console.log('Open Babel:', version.trim());
    } catch (e) {
        console.warn('Warning: Open Babel not found');
    }
});
