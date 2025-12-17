// Get CID from URL parameters
const urlParams = new URLSearchParams(window.location.search);
const cid = urlParams.get('cid');
const codid = urlParams.get('codid');
const pdbid = urlParams.get('pdbid');  // PDB ID for biomolecules
const bioassembly = urlParams.get('bioassembly') === 'true';  // Show biological assembly
const name = urlParams.get('name');
const style = urlParams.get('style') || 'stick:sphere';
const stickRadius = parseFloat(urlParams.get('stickRadius') || '0.15');
const sphereRadius = parseFloat(urlParams.get('sphereRadius') || '0.3');
const autoRotate = urlParams.get('autoRotate') !== 'false';
const bgColor = decodeURIComponent(urlParams.get('bgColor') || '#1a1a2e');
// Protein-specific styling parameters
const chainType = urlParams.get('chainType') || 'cartoon';  // cartoon, ribbon, trace, tube, stick
const chainColor = urlParams.get('chainColor') || 'spectrum';  // spectrum, chain, ss, residue

// Log all URL parameters for debugging
console.log('%c🔍 3DMOL VIEWER INITIALIZATION', 'background: #2196F3; color: white; font-weight: bold; padding: 8px; font-size: 16px;');
console.log('URL Parameters:', {
    cid,
    codid,
    pdbid,
    bioassembly,
    name,
    style,
    chainType,
    chainColor,
    bgColor,
    autoRotate
});

// Set body background color from URL parameter
document.documentElement.style.setProperty('--bg-color', bgColor);

// External MolView embed URL (used as fallback for minerals/crystals)
const MOLVIEW_EMBED = 'https://embed.molview.org/v1';


// Van der Waals radii (in Angstroms) - scientifically accurate atomic radii
// These are the actual physical sizes of atoms used in chemistry
// Source: Bondi (1964), Rowland & Taylor (1996), Alvarez (2013)
const vanDerWaalsRadii = {
    // ===== ORGANIC ATOMS (Most Common) =====
    'H': 1.20,   // Hydrogen
    'C': 1.70,   // Carbon
    'N': 1.55,   // Nitrogen
    'O': 1.52,   // Oxygen
    'F': 1.47,   // Fluorine
    'P': 1.80,   // Phosphorus
    'S': 1.80,   // Sulfur
    'Cl': 1.75,  // Chlorine
    'Br': 1.85,  // Bromine
    'I': 1.98,   // Iodine

    // ===== HALOGENS (Complete) =====
    'At': 2.02,  // Astatine

    // ===== NOBLE GASES =====
    'He': 1.40,  // Helium
    'Ne': 1.54,  // Neon
    'Ar': 1.88,  // Argon
    'Kr': 2.02,  // Krypton
    'Xe': 2.16,  // Xenon
    'Rn': 2.20,  // Radon

    // ===== ALKALI METALS (Group 1) =====
    'Li': 1.82,  // Lithium
    'Na': 2.27,  // Sodium
    'K': 2.75,   // Potassium
    'Rb': 3.03,  // Rubidium
    'Cs': 3.43,  // Cesium

    // ===== ALKALINE EARTH METALS (Group 2) =====
    'Be': 1.53,  // Beryllium
    'Mg': 1.73,  // Magnesium
    'Ca': 2.31,  // Calcium
    'Sr': 2.49,  // Strontium
    'Ba': 2.68,  // Barium

    // ===== TRANSITION METALS (Common) =====
    'Sc': 2.15,  // Scandium
    'Ti': 2.11,  // Titanium
    'V': 2.07,   // Vanadium
    'Cr': 2.06,  // Chromium
    'Mn': 2.05,  // Manganese
    'Fe': 2.04,  // Iron
    'Co': 2.00,  // Cobalt
    'Ni': 1.63,  // Nickel
    'Cu': 1.40,  // Copper
    'Zn': 1.39,  // Zinc
    'Y': 2.32,   // Yttrium
    'Zr': 2.23,  // Zirconium
    'Nb': 2.18,  // Niobium
    'Mo': 2.17,  // Molybdenum
    'Tc': 2.16,  // Technetium
    'Ru': 2.13,  // Ruthenium
    'Rh': 2.10,  // Rhodium
    'Pd': 1.63,  // Palladium
    'Ag': 1.72,  // Silver
    'Cd': 1.58,  // Cadmium
    'Hf': 2.23,  // Hafnium
    'Ta': 2.22,  // Tantalum
    'W': 2.18,   // Tungsten
    'Re': 2.16,  // Rhenium
    'Os': 2.16,  // Osmium
    'Ir': 2.13,  // Iridium
    'Pt': 1.75,  // Platinum
    'Au': 1.66,  // Gold
    'Hg': 1.55,  // Mercury

    // ===== POST-TRANSITION METALS =====
    'Al': 1.84,  // Aluminum
    'Ga': 1.87,  // Gallium
    'In': 1.93,  // Indium
    'Tl': 1.96,  // Thallium
    'Sn': 2.17,  // Tin
    'Pb': 2.02,  // Lead
    'Bi': 2.07,  // Bismuth

    // ===== METALLOIDS =====
    'B': 1.92,   // Boron
    'Si': 2.10,  // Silicon
    'Ge': 2.11,  // Germanium
    'As': 1.85,  // Arsenic
    'Sb': 2.06,  // Antimony
    'Te': 2.06,  // Tellurium
    'Po': 1.97,  // Polonium

    // ===== OTHER NON-METALS =====
    'Se': 1.90,  // Selenium

    // ===== LANTHANIDES (Rare Earths) =====
    'La': 2.43,  // Lanthanum
    'Ce': 2.42,  // Cerium
    'Pr': 2.40,  // Praseodymium
    'Nd': 2.39,  // Neodymium
    'Pm': 2.38,  // Promethium
    'Sm': 2.36,  // Samarium
    'Eu': 2.35,  // Europium
    'Gd': 2.34,  // Gadolinium
    'Tb': 2.33,  // Terbium
    'Dy': 2.31,  // Dysprosium
    'Ho': 2.30,  // Holmium
    'Er': 2.29,  // Erbium
    'Tm': 2.27,  // Thulium
    'Yb': 2.26,  // Ytterbium
    'Lu': 2.24,  // Lutetium

    // ===== ACTINIDES =====
    'Ac': 2.47,  // Actinium
    'Th': 2.45,  // Thorium
    'Pa': 2.43,  // Protactinium
    'U': 1.86,   // Uranium
    'Np': 2.39,  // Neptunium
    'Pu': 2.43,  // Plutonium
    'Am': 2.44,  // Americium
    'Cm': 2.45,  // Curium
};

if (pdbid) {
    // Load PDB structure (biomolecule/protein)
    loadPDBStructure(pdbid, bioassembly);
} else if (cid) {
    loadMolecule(cid);
} else if (codid) {
    fallbackToMolView(codid, 'codid');
} else if (name) {
    resolveAndLoad(name);
} else {
    showError('No molecule ID or Name provided');
}

/**
 * Load biological assembly using RCSB ModelServer API
 * This fetches ONLY C-alpha backbone atoms for the full biological assembly
 * Result: ~15,000 atoms instead of ~500,000 - perfect for blob visualization!
 * @param {string} pdbId - The PDB ID
 * @param {Object} viewer - 3Dmol viewer instance
 */
async function loadBiologicalAssembly(pdbId, viewer) {
    console.log(`%c🔬 LOADING BIOLOGICAL ASSEMBLY (LOW-POLY METHOD)`, 'background: #9C27B0; color: white; font-weight: bold; padding: 4px;');
    console.log(`   Using RCSB ModelServer to get C-alpha backbone only`);
    console.log(`   This reduces ~500k atoms to ~15k atoms for smooth blob visualization!`);

    // Update loading message
    const loadingEl = document.getElementById('loading');
    if (loadingEl) {
        loadingEl.innerHTML = `<div style="text-align: center;">
            <div>🧬 Loading Biological Assembly...</div>
            <div style="font-size: 12px; margin-top: 8px;">Fetching simplified backbone structure</div>
            <div style="font-size: 11px; color: #888; margin-top: 4px;">Creating low-poly visualization</div>
        </div>`;
    }

    // Use RCSB ModelServer API to get C-alpha atoms for biological assembly
    // This endpoint computes the assembly AND filters to only CA atoms!
    // Format: https://models.rcsb.org/v1/{id}/assembly?name=1&label_atom_id=CA&encoding=cif
    const modelServerUrl = `https://models.rcsb.org/v1/${pdbId.toUpperCase()}/assembly?name=1&label_atom_id=CA&encoding=cif`;
    console.log(`   Fetching from ModelServer: ${modelServerUrl}`);

    const response = await fetch(modelServerUrl);
    if (!response.ok) {
        console.warn(`   ModelServer failed (${response.status}), falling back to standard method...`);
        return await loadBiologicalAssemblyFallback(pdbId);
    }

    let cifData = await response.text();
    const sizeKB = (cifData.length / 1024).toFixed(1);
    console.log(`   Downloaded ${sizeKB} KB of C-alpha backbone data`);

    // Log first 500 chars of data for debugging
    console.log(`   Data preview: ${cifData.substring(0, 300)}...`);

    // For now, don't subsample - just return the data to test if it works
    // Subsampling CIF format is complex and may break the file structure
    console.log(`   Skipping subsampling - loading full C-alpha dataset`);

    return { data: cifData, format: 'cif', isSimplified: true };
}

/**
 * Fallback: Download .pdb1 file and extract C-alpha atoms manually
 * Used if ModelServer is unavailable
 */
async function loadBiologicalAssemblyFallback(pdbId) {
    console.log(`   ⚠️ Using fallback: downloading .pdb1 and extracting C-alpha atoms`);

    const bioAssemblyUrl = `https://files.rcsb.org/download/${pdbId.toUpperCase()}.pdb1`;
    const response = await fetch(bioAssemblyUrl);

    if (!response.ok) {
        throw new Error(`Failed to fetch biological assembly: ${response.status}`);
    }

    const pdbData = await response.text();
    const sizeMB = (pdbData.length / 1024 / 1024).toFixed(2);
    console.log(`   Downloaded ${sizeMB} MB of full assembly data`);

    // Extract only C-alpha atoms
    const lines = pdbData.split('\n');
    const caLines = [];

    for (const line of lines) {
        // Keep ATOM lines for CA (C-alpha) atoms only
        if ((line.startsWith('ATOM') || line.startsWith('HETATM'))) {
            const atomName = line.substring(12, 16).trim();
            if (atomName === 'CA') {
                caLines.push(line);
            }
        }
        // Keep MODEL/ENDMDL for multi-model support
        else if (line.startsWith('MODEL') || line.startsWith('ENDMDL') ||
            line.startsWith('TER') || line.startsWith('END')) {
            caLines.push(line);
        }
    }

    const caData = caLines.join('\n');
    const reducedSizeKB = (caData.length / 1024).toFixed(1);
    console.log(`   Extracted ${caLines.length} lines, reduced to ${reducedSizeKB} KB`);

    return { data: caData, format: 'pdb', isSimplified: true };
}


/**
 * Load PDB structure from RCSB with optional biological assembly
 * @param {string} pdbId - The PDB ID (e.g., "1HHO" for hemoglobin)
 * @param {boolean} loadBioAssembly - Whether to load biological assembly (true) or asymmetric unit (false)
 */
async function loadPDBStructure(pdbId, loadBioAssembly = false) {
    console.log(`%c🧬 LOADING PDB STRUCTURE`, 'background: #E91E63; color: white; font-weight: bold; padding: 8px; font-size: 16px;');
    console.log(`   PDB ID: ${pdbId}`);
    console.log(`   Biological Assembly: ${loadBioAssembly}`);
    console.log(`   Chain Type: ${chainType}`);
    console.log(`   Chain Color: ${chainColor}`);

    try {
        const viewer = $3Dmol.createViewer('viewer', {
            backgroundColor: bgColor
        });

        // For biological assemblies, use low-poly C-alpha backbone visualization
        // This downloads only C-alpha atoms and renders them as overlapping spheres
        if (loadBioAssembly) {
            console.log(`%c📦 BIOLOGICAL ASSEMBLY MODE (LOW-POLY)`, 'background: #4CAF50; color: white; font-weight: bold; padding: 4px;');
            console.log(`   Using C-alpha backbone for blob-like visualization`);

            try {
                const assemblyData = await loadBiologicalAssembly(pdbId, viewer);
                renderLowPolyAssembly(viewer, assemblyData.data, pdbId, assemblyData.format);
                return;
            } catch (error) {
                console.error('Failed to load biological assembly:', error);
                showError(`Failed to load biological assembly for ${pdbId}`);
                return;
            }
        }

        // Standard asymmetric unit - load full detail
        const pdbUrl = `https://files.rcsb.org/download/${pdbId.toUpperCase()}.pdb`;
        console.log(`%c📦 Fetching ASYMMETRIC UNIT`, 'background: #FF9800; color: white; font-weight: bold; padding: 4px;');
        console.log(`   URL: ${pdbUrl}`);

        // Fetch the PDB file
        console.log(`⏳ Fetching PDB file...`);
        const response = await fetch(pdbUrl);
        console.log(`   Response status: ${response.status} ${response.statusText}`);

        if (!response.ok) {
            throw new Error(`Failed to fetch PDB: ${pdbId} (Status: ${response.status})`);
        }

        const pdbData = await response.text();
        console.log(`%c✅ Successfully loaded ${pdbData.length} bytes of PDB data`, 'color: #4CAF50; font-weight: bold;');
        renderPDBStructure(viewer, pdbData, pdbId, loadBioAssembly);

    } catch (error) {
        console.error('%c❌ ERROR LOADING PDB:', 'color: #f44336; font-weight: bold; font-size: 14px;', error);
        showError(`Failed to load PDB ${pdbId}: ${error.message}`);
    }
}

/**
 * Render PDB structure with protein-appropriate styling
 */
function renderPDBStructure(viewer, pdbData, pdbId, isBioAssembly) {
    console.log(`🎨 Rendering PDB structure: ${pdbId} (Bio Assembly: ${isBioAssembly})`);

    // Add the model
    viewer.addModel(pdbData, 'pdb');

    // Get atom count for logging
    const atoms = viewer.getModel().selectedAtoms({});
    console.log(`   Total atoms in structure: ${atoms.length}`);

    // Apply chain coloring based on chainColor parameter
    let colorScheme;
    switch (chainColor) {
        case 'spectrum':
            colorScheme = { color: 'spectrum' };
            break;
        case 'chain':
            colorScheme = { colorscheme: 'chainHetatm' };
            break;
        case 'ss':
            // Secondary structure coloring
            colorScheme = { colorscheme: 'ssJmol' };
            break;
        case 'residue':
            colorScheme = { colorscheme: 'amino' };
            break;
        case 'bfactor':
            colorScheme = { colorscheme: { prop: 'b', gradient: 'roygb', min: 0, max: 100 } };
            break;
        default:
            colorScheme = { color: 'spectrum' };
    }

    // FIRST: Hide everything by default
    viewer.setStyle({}, { sphere: { hidden: true } });

    // THEN: Show only protein/nucleic acid chains with proper styling
    // Apply chain representation based on chainType parameter
    switch (chainType) {
        case 'cartoon':
            viewer.setStyle({ hetflag: false }, { cartoon: { ...colorScheme } });
            break;
        case 'ribbon':
            viewer.setStyle({ hetflag: false }, { cartoon: { style: 'ribbon', ...colorScheme } });
            break;
        case 'trace':
            viewer.setStyle({ hetflag: false }, { cartoon: { style: 'trace', ...colorScheme } });
            break;
        case 'tube':
            viewer.setStyle({ hetflag: false }, { cartoon: { tubes: true, ...colorScheme } });
            break;
        case 'stick':
            viewer.setStyle({ hetflag: false }, { stick: { radius: 0.2, ...colorScheme } });
            break;
        case 'line':
            viewer.setStyle({ hetflag: false }, { line: { ...colorScheme } });
            break;
        case 'sphere':
            viewer.setStyle({ hetflag: false }, { sphere: { scale: 0.3, ...colorScheme } });
            break;
        default:
            viewer.setStyle({ hetflag: false }, { cartoon: { ...colorScheme } });
    }

    // OPTIONALLY: Show heteroatoms (ligands) as small sticks (commented out for cleaner view of large assemblies)
    // viewer.setStyle({ hetflag: true }, { stick: { radius: 0.15, colorscheme: 'Jmol' } });

    // Hide water molecules (they clutter the view)
    viewer.setStyle({ resn: 'HOH' }, { sphere: { hidden: true } });

    console.log(`   Applying style: ${chainType} with ${chainColor} coloring`);

    // Center and zoom
    viewer.zoomTo();
    viewer.zoom(0.8);  // Zoom out more for large assemblies
    viewer.render();

    console.log(`   Rendered and centered on structure`);

    // Add slow rotation for biomolecules if autoRotate is enabled
    if (autoRotate) {
        viewer.spin('y', 0.5);  // Slower spin for large structures
    }

    // Custom zoom control
    const viewerElement = document.getElementById('viewer');
    viewerElement.addEventListener('wheel', (event) => {
        event.preventDefault();
        event.stopPropagation();
        const zoomFactor = event.deltaY > 0 ? 0.97 : 1.03;
        viewer.zoom(zoomFactor);
        viewer.render();
    }, { passive: false });

    // Hide loading
    document.getElementById('loading').style.display = 'none';

    const assemblyType = isBioAssembly ? 'Biological Assembly' : 'Asymmetric Unit';
    console.log(`✅ PDB ${pdbId} (${assemblyType}) loaded successfully`);
    document.title = `${pdbId.toUpperCase()} - ${assemblyType}`;
}

/**
 * Render biological assembly as colored point cloud
 * Simple, reliable visualization using small spheres
 */
function renderLowPolyAssembly(viewer, structureData, pdbId, format) {
    console.log(`%c🎨 RENDERING BIOLOGICAL ASSEMBLY`, 'background: #9C27B0; color: white; font-weight: bold; padding: 4px;');
    console.log(`   PDB ID: ${pdbId}`);
    console.log(`   Data format: ${format}`);

    // Add the model
    viewer.addModel(structureData, format);

    // Get all atoms
    const atoms = viewer.getModel().selectedAtoms({});
    console.log(`   Atoms loaded: ${atoms.length}`);

    if (atoms.length === 0) {
        console.error('   No atoms found in structure!');
        showError('Failed to parse biological assembly structure');
        return;
    }

    // Simple, reliable visualization: small colored spheres
    // Not trying to be fancy - just show the structure clearly
    viewer.setStyle({}, {
        sphere: {
            radius: 3.0,              // Small spheres
            colorscheme: 'chainHetatm' // Color by chain - shows symmetry
        }
    });

    console.log(`   Applied sphere visualization (${atoms.length} atoms, radius 3.0Å)`);

    // Center and zoom
    viewer.zoomTo();
    viewer.zoom(0.5);
    viewer.render();

    console.log(`   Rendered and centered`);

    // Slow rotation
    if (autoRotate) {
        viewer.spin('y', 0.15);
    }

    // Zoom control
    const viewerElement = document.getElementById('viewer');
    viewerElement.addEventListener('wheel', (event) => {
        event.preventDefault();
        event.stopPropagation();
        const zoomFactor = event.deltaY > 0 ? 0.95 : 1.05;
        viewer.zoom(zoomFactor);
        viewer.render();
    }, { passive: false });

    // Hide loading
    document.getElementById('loading').style.display = 'none';

    console.log(`%c✅ ASSEMBLY LOADED`, 'background: #4CAF50; color: white; font-weight: bold; padding: 4px;');
    console.log(`   ${pdbId.toUpperCase()} - Biological Assembly (${atoms.length} atoms)`);
    document.title = `${pdbId.toUpperCase()} - Biological Assembly`;
}

async function resolveAndLoad(name) {
    console.log(`🔍 Resolving molecule: ${name}`);

    // STEP 1: Use Integrated Search (no server needed!)
    if (window.IntegratedSearch) {
        try {
            console.log(`🔍 [IntegratedSearch] Querying for: "${name}"`);
            const data = await window.IntegratedSearch.search(name, { format: 'full' });

            if (!data.error) {
                console.log('📦 IntegratedSearch response:', data);
                const category = data.primary_type?.toLowerCase();

                // CASE 1: Compound → Use CID or SMILES
                if (category === 'compound') {
                    console.log(`✅ Found compound: ${data.name}`);

                    if (data.cid) {
                        loadMolecule(data.cid);
                        return;
                    } else if (data.canonical_smiles) {
                        const cid = await getCIDFromSMILES(data.canonical_smiles);
                        if (cid) {
                            loadMolecule(cid);
                            return;
                        }
                    }
                }

                // CASE 2: Mineral → Use COD ID with MolView
                else if (category === 'mineral' && data.codid) {
                    console.log(`💎 Found mineral: ${data.name} (COD ID: ${data.codid})`);
                    fallbackToMolView(data.codid, 'codid');
                    return;
                }

                // CASE 3: Biomolecule/Protein → Use PDB ID with MolView embed
                else if ((category === 'biomolecule' || category === 'protein') && data.pdbid) {
                    console.log(`🧬 Found biomolecule: ${data.name} (PDB ID: ${data.pdbid})`);
                    const molviewUrl = `${MOLVIEW_EMBED}/?pdbid=${data.pdbid}`;
                    embedMolView(molviewUrl, data.name);
                    return;
                }
            }
        } catch (e) {
            console.warn(`⚠️ IntegratedSearch failed: ${e.message}`);
        }
    }

    // STEP 2: Fallback to external Search API (port 8001)
    const SEARCH_API = 'http://localhost:8001';
    try {
        console.log(`🌐 Querying Search API: ${SEARCH_API}/search?q=${encodeURIComponent(name)}`);
        const response = await fetch(`${SEARCH_API}/search?q=${encodeURIComponent(name)}`);

        if (response.ok) {
            const data = await response.json();
            console.log('📦 Search API response:', data);

            if (data.found) {
                const category = data.category?.toLowerCase();

                // CASE 1: Compound or Mineral → Use SMILES/SDF with local 3Dmol.js
                if (category === 'compound' || category === 'mineral') {
                    console.log(`✅ Found ${category}: ${name}`);

                    // Check for mineral ID (codid)
                    if (category === 'mineral' && (data.codid || data.id)) {
                        const mineralId = data.codid || data.id;
                        console.log('💎 Mineral detected with ID:', mineralId);

                        // If no direct 3D data, use MolView with codid
                        if (!data.sdf && !data.smiles) {
                            fallbackToMolView(mineralId, 'codid');
                            return;
                        }
                    }

                    if (data.sdf) {
                        // Render SDF directly with 3Dmol.js
                        console.log('📄 Using SDF data from MolView');
                        loadMoleculeFromSDF(data.sdf, name);
                        return;
                    } else if (data.smiles) {
                        // Convert SMILES to SDF via PubChem, then render
                        console.log('🧪 Using SMILES from MolView:', data.smiles);
                        const cid = await getCIDFromSMILES(data.smiles);
                        if (cid) {
                            loadMolecule(cid);
                            return;
                        }
                    }
                }

                // CASE 2: Protein/Biomolecule → Use MolView embed
                else if (category === 'protein' || category === 'biomolecule' || category === 'pdb') {
                    console.log(`🧬 Found ${category}: ${name}`);

                    // Use MolView embed with the query
                    const molviewUrl = `${MOLVIEW_EMBED}/?q=${encodeURIComponent(name)}`;
                    embedMolView(molviewUrl, name);
                    return;
                }
            }
        }
    } catch (e) {
        console.warn(`⚠️ Search API failed: ${e.message}`);
    }

    // STEP 3: Fallback to PubChem for small molecules
    console.log('🔄 Falling back to PubChem lookup...');
    try {
        const response = await fetch(`https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/${encodeURIComponent(name)}/cids/JSON`);
        if (response.ok) {
            const data = await response.json();
            if (data.IdentifierList && data.IdentifierList.CID && data.IdentifierList.CID.length > 0) {
                const cid = data.IdentifierList.CID[0];
                console.log(`✅ Resolved name "${name}" to CID ${cid} via PubChem`);
                loadMolecule(cid);
                return;
            }
        }
    } catch (e) {
        console.warn(`⚠️ PubChem lookup failed: ${e.message}`);
    }

    // STEP 4: Final fallback to online MolView embed
    console.warn(`⚠️ Could not resolve "${name}" locally. Using online MolView.`);
    fallbackToMolView(name, false);
}

// Helper: Get CID from SMILES
async function getCIDFromSMILES(smiles) {
    try {
        const url = `https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/smiles/${encodeURIComponent(smiles)}/cids/JSON`;
        const response = await fetch(url);
        if (response.ok) {
            const data = await response.json();
            if (data.IdentifierList && data.IdentifierList.CID && data.IdentifierList.CID.length > 0) {
                return data.IdentifierList.CID[0];
            }
        }
    } catch (e) {
        console.warn('Failed to get CID from SMILES:', e);
    }
    return null;
}

// Helper: Load molecule from SDF data
function loadMoleculeFromSDF(sdfData, name) {
    console.log('🎨 Rendering molecule from SDF data');

    // Hide loading
    document.getElementById('loading').style.display = 'none';

    // Create 3Dmol viewer
    const bgColor = decodeURIComponent(urlParams.get('bgColor') || '#1a1a2e');
    let viewer = $3Dmol.createViewer('viewer', {
        backgroundColor: bgColor,
        antialias: true
    });

    // Add model from SDF
    const model = viewer.addModel(sdfData, 'sdf');

    // Apply default styling
    viewer.setStyle({}, {
        stick: { radius: 0.15, colorscheme: 'Jmol' },
        sphere: { colorscheme: 'Jmol' }
    });

    viewer.zoomTo();
    viewer.zoom(1.2);
    viewer.render();

    // Enable rotation if configured
    const autoRotate = urlParams.get('autoRotate') !== 'false';
    if (autoRotate) {
        viewer.spin(true);
    }

    console.log('✅ Molecule loaded from SDF:', name);
}

// Helper: Embed MolView iframe
function embedMolView(url, name) {
    console.log(`🖼️ Embedding MolView: ${url}`);

    // Clear the viewer container
    document.body.innerHTML = '';
    document.body.style.margin = '0';
    document.body.style.padding = '0';
    document.body.style.overflow = 'hidden';

    // Create MolView iframe
    const iframe = document.createElement('iframe');
    iframe.src = url;
    iframe.style.width = '100vw';
    iframe.style.height = '100vh';
    iframe.style.border = 'none';
    iframe.style.display = 'block';
    iframe.title = `MolView: ${name}`;

    document.body.appendChild(iframe);

    console.log(`✅ MolView embedded for: ${name}`);
}

function fallbackToMolView(query, type = 'q') {
    console.warn(`⚠️ Falling back to MolView for ${type}: ${query}`);

    // Clear the viewer container
    document.body.innerHTML = '';
    document.body.style.margin = '0';
    document.body.style.padding = '0';
    document.body.style.overflow = 'hidden';

    // Create MolView iframe pointing to external MolView embed
    const iframe = document.createElement('iframe');

    // Use external MolView embed endpoint with query parameter
    let param;
    if (type === 'cid' || type === true) {
        param = `cid=${query}`;
    } else if (type === 'codid') {
        param = `codid=${query}`;
    } else if (type === 'pdbid') {
        param = `pdbid=${query}`;
    } else {
        param = `smiles=${encodeURIComponent(query)}`;
    }
    iframe.src = `${MOLVIEW_EMBED}/?${param}`;

    console.log('%c🔗 3D VIEWER IFRAME SRC:', 'background: #FF5722; color: white; font-weight: bold; padding: 4px; font-size: 14px;');
    console.log(`   ${iframe.src}`);

    iframe.style.width = '100vw';
    iframe.style.height = '100vh';
    iframe.style.border = 'none';
    iframe.style.display = 'block';

    document.body.appendChild(iframe);
}


function showError(message) {
    document.getElementById('loading').style.display = 'none';
    document.getElementById('error').style.display = 'block';
    document.getElementById('errorMessage').textContent = message;
}

async function loadMolecule(cid) {
    try {
        // Create 3Dmol viewer with user-selected background color
        // For transparent, we set it to match the body background dynamically
        const effectiveBgColor = bgColor === 'transparent' ?
            getComputedStyle(document.body).backgroundColor : bgColor;

        let viewer = $3Dmol.createViewer('viewer', {
            backgroundColor: effectiveBgColor,
            antialias: true
        });

        // Ensure the canvas element has no borders (wait a moment for canvas to be created)
        setTimeout(() => {
            const canvas = document.querySelector('#viewer canvas');
            if (canvas) {
                canvas.style.border = 'none';
                canvas.style.outline = 'none';
                canvas.style.margin = '0';
                canvas.style.padding = '0';
            }
        }, 100);

        // Try multiple PubChem endpoints for 3D structure data
        // Priority 1: 3D conformer (pre-computed 3D structure)
        // Priority 2: Computed 3D structure
        // Note: 2D structure fallback removed as per user request (strict 3D only)

        let sdfData = null;
        let fetchMethod = '';

        // Try 3D conformer first
        console.log(`🔍 Trying to fetch 3D conformer for CID ${cid}...`);
        try {
            const url3d = `https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/${cid}/record/SDF?record_type=3d`;
            const response3d = await fetch(url3d);

            if (response3d.ok) {
                sdfData = await response3d.text();
                fetchMethod = '3D conformer';
                console.log(`✅ Successfully fetched 3D conformer for CID ${cid}`);
            } else {
                console.warn(`⚠️ 3D conformer not available (status ${response3d.status})`);
            }
        } catch (e) {
            console.warn(`⚠️ 3D conformer fetch failed: ${e.message}`);
        }

        // Fallback: Try computed 3D
        if (!sdfData) {
            console.log(`🔍 Trying to fetch computed 3D structure for CID ${cid}...`);
            try {
                const urlComputed = `https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/${cid}/SDF?record_type=3d`;
                const responseComputed = await fetch(urlComputed);

                if (responseComputed.ok) {
                    sdfData = await responseComputed.text();
                    fetchMethod = 'computed 3D structure';
                    console.log(`✅ Successfully fetched computed 3D structure for CID ${cid}`);
                }
            } catch (e) {
                console.warn(`⚠️ Computed 3D fetch failed: ${e.message}`);
            }
        }

        // If still no data, fallback to MolView (The "Singular Solution" for missing 3D data)
        if (!sdfData) {
            console.warn(`⚠️ No 3D data from PubChem for CID ${cid}. Falling back to MolView...`);
            fallbackToMolView(cid, true);
            return;
        }

        console.log(`📦 Using ${fetchMethod} for visualization`);


        // Add molecule to viewer
        const model = viewer.addModel(sdfData, 'sdf');

        // Count atoms for protein detection
        const atoms = model.selectedAtoms({});
        const atomCount = atoms.length;
        console.log(`Molecule has ${atomCount} atoms`);

        // Auto-detect large molecules (proteins) and switch to cartoon mode
        let finalStyle = style;
        const styles = style.split(':');
        const isStickOrBallStick = styles.includes('stick') || (styles.includes('stick') && styles.includes('sphere'));

        if (atomCount > 100 && isStickOrBallStick) {
            console.log(`⚙️ Large molecule detected (${atomCount} atoms), switching to cartoon mode for better visualization`);
            finalStyle = 'cartoon';
        }

        // Apply style based on user settings or auto-detection
        const styleArr = finalStyle.split(':');
        const styleConfig = {};

        // Calculate scaling factor for user's sphere radius preference
        // Base VdW radius for Carbon (1.70 Å), normalized to user's preference
        const radiusScale = sphereRadius / 1.70;

        styleArr.forEach(s => {
            if (s === 'stick') {
                styleConfig.stick = {
                    radius: stickRadius,
                    colorscheme: 'Jmol'
                };
            } else if (s === 'sphere') {
                // Use Van der Waals radii for accurate atomic sizes
                // Apply per-atom styling with correct radii ratios
                styleConfig.sphere = {
                    colorscheme: 'Jmol',
                    // 3Dmol will use vdwRadii internally for proper sizing
                };
            } else if (s === 'line') {
                styleConfig.line = { colorscheme: 'Jmol' };
            } else if (s === 'cross') {
                styleConfig.cross = { colorscheme: 'Jmol' };
            } else if (s === 'cartoon') {
                styleConfig.cartoon = { color: 'spectrum' };
            }
        });

        // Default to ball-and-stick if no style specified
        if (Object.keys(styleConfig).length === 0) {
            styleConfig.stick = { radius: stickRadius, colorscheme: 'Jmol' };
            styleConfig.sphere = { colorscheme: 'Jmol' };
        }

        // Apply styles with proper atomic radii
        // For each atom, set the radius based on Van der Waals radii
        const modelAtoms = model.selectedAtoms({});
        modelAtoms.forEach(atom => {
            const element = atom.elem;
            const vdwRadius = vanDerWaalsRadii[element] || 1.70; // Default to Carbon if unknown
            const scaledRadius = (vdwRadius / 1.70) * sphereRadius; // Scale relative to Carbon

            // Apply style with correct radius for this specific atom
            const atomStyle = JSON.parse(JSON.stringify(styleConfig)); // Deep copy
            if (atomStyle.sphere) {
                atomStyle.sphere.radius = scaledRadius;
            }
            viewer.setStyle({ serial: atom.serial }, atomStyle);
        });

        // If no sphere style, just apply the config normally
        if (!styleArr.includes('sphere')) {
            viewer.setStyle({}, styleConfig);
        }

        // Center and zoom
        viewer.zoomTo();
        viewer.zoom(1.2);
        viewer.render();

        // Add custom zoom control for very precise zooming
        const viewerElement = document.getElementById('viewer');
        viewerElement.addEventListener('wheel', (event) => {
            event.preventDefault();
            event.stopPropagation();

            // Apply very slow zoom factor for precise control
            // 1% change per scroll tick for smooth, controlled zooming
            const zoomFactor = event.deltaY > 0 ? 0.99 : 1.01;

            // Apply the new zoom
            viewer.zoom(zoomFactor);
            viewer.render();
        }, { passive: false });

        // Enable rotation if user wants it
        if (autoRotate) {
            viewer.spin(true);
        }

        // Hide loading
        document.getElementById('loading').style.display = 'none';

        console.log('✅ Molecule loaded successfully:', name || cid);
    } catch (error) {
        console.error('Error loading molecule:', error);
        showError(error.message);
    }
}
