// Get CID from URL parameters
const urlParams = new URLSearchParams(window.location.search);
const cid = urlParams.get('cid');
const pdbid = urlParams.get('pdbid');  // NEW: Direct PDB ID support
const codid = urlParams.get('codid');
const name = urlParams.get('name');
const style = urlParams.get('style') || 'stick:sphere';
const stickRadius = parseFloat(urlParams.get('stickRadius') || '0.15');
const sphereRadius = parseFloat(urlParams.get('sphereRadius') || '0.3');
const autoRotate = urlParams.get('autoRotate') !== 'false';
const bgColor = decodeURIComponent(urlParams.get('bgColor') || '#1a1a2e');

// Set body background color from URL parameter
document.documentElement.style.setProperty('--bg-color', bgColor);

const MOLVIEW_API = 'http://localhost:8000';
const SEARCH_API = 'http://localhost:8001';

// =============================================================================
// CIR (Chemical Identifier Resolver) - Generates 3D coordinates from SMILES
// This is the secret sauce that makes MolView work for complex molecules!
// =============================================================================

/**
 * Fetch SMILES string from PubChem for a given CID
 * @param {number} cid - PubChem Compound ID
 * @returns {Promise<string|null>} - SMILES string or null if not found
 */
async function getSMILESFromPubChem(cid) {
    console.log(`🧪 Fetching SMILES from PubChem for CID ${cid}...`);
    try {
        const url = `https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/${cid}/property/IsomericSMILES/JSON`;
        const response = await fetch(url);
        if (response.ok) {
            const data = await response.json();
            const smiles = data?.PropertyTable?.Properties?.[0]?.IsomericSMILES;
            if (smiles) {
                console.log(`✅ Got SMILES: ${smiles.substring(0, 50)}${smiles.length > 50 ? '...' : ''}`);
                return smiles;
            }
        }
    } catch (e) {
        console.warn(`⚠️ Failed to get SMILES from PubChem: ${e.message}`);
    }
    return null;
}

/**
 * Generate 3D SDF from SMILES using CIR (Chemical Identifier Resolver)
 * This uses NCI's CADD server to compute 3D coordinates via force field optimization
 * @param {string} smiles - SMILES string
 * @returns {Promise<string|null>} - 3D SDF data or null if failed
 */
async function generate3DFromCIR(smiles) {
    console.log(`🔬 Generating 3D structure from SMILES via CIR...`);
    try {
        // Encode SMILES properly - CIR needs special handling for # and \
        const encodedSmiles = encodeURIComponent(smiles);

        // The magic parameter: get3d=True makes CIR compute 3D coordinates!
        const url = `https://cactus.nci.nih.gov/chemical/structure/${encodedSmiles}/file?format=sdf&get3d=True`;

        console.log(`📡 CIR URL: ${url.substring(0, 100)}...`);

        const response = await fetch(url);
        if (response.ok) {
            const sdfData = await response.text();

            // Validate SDF data - must have proper structure
            // CIR returns HTML error pages on failure, not valid SDF
            if (!sdfData ||
                sdfData.includes('<!DOCTYPE') ||
                sdfData.includes('<html') ||
                sdfData.includes('<h1>Page not found') ||
                sdfData.includes('Error') ||
                !sdfData.includes('V2000') && !sdfData.includes('V3000')) {
                console.warn(`⚠️ CIR returned invalid SDF data`);
                return null;
            }

            // Additional validation: check for atom block
            const lines = sdfData.split('\n');
            if (lines.length < 5) {
                console.warn(`⚠️ CIR returned too few lines (${lines.length})`);
                return null;
            }

            console.log(`✅ CIR generated valid 3D structure (${sdfData.length} bytes, ${lines.length} lines)`);
            return sdfData;
        } else {
            console.warn(`⚠️ CIR request failed with status ${response.status}`);
        }
    } catch (e) {
        console.warn(`⚠️ CIR 3D generation failed: ${e.message}`);
    }
    return null;
}

// =============================================================================
// RCSB PDB - Direct loading of protein/biomolecule structures
// =============================================================================

/**
 * Fetch PDB file directly from RCSB
 * @param {string} pdbId - PDB ID (e.g., "4INS", "4RHV")
 * @returns {Promise<string|null>} - PDB data or null if failed
 */
async function getPDBFromRCSB(pdbId) {
    console.log(`🧬 Fetching PDB ${pdbId} from RCSB...`);
    try {
        const url = `https://files.rcsb.org/view/${pdbId.toUpperCase()}.pdb`;
        const response = await fetch(url);
        if (response.ok) {
            const pdbData = await response.text();
            console.log(`✅ Got PDB data (${pdbData.length} bytes)`);
            return pdbData;
        } else {
            console.warn(`⚠️ RCSB returned status ${response.status}`);
        }
    } catch (e) {
        console.warn(`⚠️ Failed to fetch PDB from RCSB: ${e.message}`);
    }
    return null;
}

/**
 * Load and render a biomolecule from PDB ID using 3Dmol.js directly
 * @param {string} pdbId - PDB ID
 */
async function loadBiomolecule(pdbId) {
    console.log(`%c🧬 Loading biomolecule: ${pdbId}`, 'background: #9C27B0; color: white; font-weight: bold; padding: 4px;');

    try {
        const pdbData = await getPDBFromRCSB(pdbId);

        if (!pdbData) {
            console.warn(`⚠️ Failed to get PDB data for ${pdbId}`);
            showError(`Failed to load PDB: ${pdbId}`);
            return;
        }

        // Create 3Dmol viewer
        const effectiveBgColor = bgColor === 'transparent' ?
            getComputedStyle(document.body).backgroundColor : bgColor;

        let viewer = $3Dmol.createViewer('viewer', {
            backgroundColor: effectiveBgColor,
            antialias: true
        });

        // Add PDB model
        const model = viewer.addModel(pdbData, 'pdb');

        // Count atoms to determine visualization style
        const atoms = model.selectedAtoms({});
        const atomCount = atoms.length;
        console.log(`Biomolecule has ${atomCount} atoms`);

        // For proteins/large biomolecules, use cartoon representation
        if (atomCount > 500) {
            console.log(`📊 Large biomolecule - using cartoon representation`);
            viewer.setStyle({}, { cartoon: { color: 'spectrum' } });
        } else if (atomCount > 100) {
            console.log(`📊 Medium biomolecule - using cartoon + line`);
            viewer.setStyle({}, {
                cartoon: { color: 'spectrum' },
                line: { colorscheme: 'Jmol' }
            });
        } else {
            // Small peptide - use ball and stick
            console.log(`📊 Small biomolecule - using ball and stick`);
            viewer.setStyle({}, {
                stick: { radius: 0.15, colorscheme: 'Jmol' },
                sphere: { radius: 0.3, colorscheme: 'Jmol' }
            });
        }

        viewer.zoomTo();
        viewer.zoom(1.2);
        viewer.render();

        // Enable rotation
        if (autoRotate) {
            viewer.spin(true);
        }

        // Hide loading
        document.getElementById('loading').style.display = 'none';

        console.log(`✅ Biomolecule ${pdbId} loaded successfully!`);
    } catch (error) {
        console.error('Error loading biomolecule:', error);
        showError(error.message);
    }
}

// =============================================================================
// COD (Crystallography Open Database) - Direct loading of mineral structures
// =============================================================================

/**
 * Fetch CIF file directly from COD
 * @param {string} codId - COD ID (e.g., "1010928" for calcite)
 * @returns {Promise<string|null>} - CIF data or null if failed
 */
async function getCIFFromCOD(codId) {
    console.log(`💎 Fetching CIF ${codId} from COD...`);
    try {
        const url = `https://www.crystallography.net/cod/${codId}.cif`;
        const response = await fetch(url);
        if (response.ok) {
            const cifData = await response.text();
            console.log(`✅ Got CIF data (${cifData.length} bytes)`);
            return cifData;
        } else {
            console.warn(`⚠️ COD returned status ${response.status}`);
        }
    } catch (e) {
        console.warn(`⚠️ Failed to fetch CIF from COD: ${e.message}`);
    }
    return null;
}

/**
 * Load and render a mineral/crystal structure from COD ID using 3Dmol.js
 * @param {string} codId - COD ID
 */
async function loadMineral(codId) {
    console.log(`%c💎 Loading mineral: ${codId}`, 'background: #FF9800; color: white; font-weight: bold; padding: 4px;');

    try {
        const cifData = await getCIFFromCOD(codId);

        if (!cifData) {
            console.warn(`⚠️ Failed to get CIF data for ${codId}, falling back to MolView`);
            fallbackToMolView(codId, 'codid');
            return;
        }

        // Create 3Dmol viewer
        const effectiveBgColor = bgColor === 'transparent' ?
            getComputedStyle(document.body).backgroundColor : bgColor;

        let viewer = $3Dmol.createViewer('viewer', {
            backgroundColor: effectiveBgColor,
            antialias: true
        });

        // Add CIF model - 3Dmol.js supports CIF format!
        // Use options to expand unit cell for better visualization
        const model = viewer.addModel(cifData, 'cif', { doAssembly: true, normalizeAssembly: true });

        // Count atoms
        const atoms = model.selectedAtoms({});
        const atomCount = atoms.length;
        console.log(`Mineral has ${atomCount} atoms`);

        // Crystal structures look best with sphere style to show the lattice
        viewer.setStyle({}, {
            sphere: {
                colorscheme: 'Jmol',
                radius: 0.4
            }
        });

        // Add unit cell box if available
        viewer.addUnitCell();

        viewer.zoomTo();
        viewer.zoom(1.2);
        viewer.render();

        // Enable rotation
        if (autoRotate) {
            viewer.spin(true);
        }

        // Hide loading
        document.getElementById('loading').style.display = 'none';

        console.log(`✅ Mineral ${codId} loaded successfully!`);
    } catch (error) {
        console.error('Error loading mineral:', error);
        // CIF parsing can fail - fall back to MolView
        console.warn(`⚠️ CIF parsing failed, falling back to MolView`);
        fallbackToMolView(codId, 'codid');
    }
}

/**
 * Fetch 2D SDF from PubChem (fallback for flat molecule display)
 * @param {number} cid - PubChem Compound ID
 * @returns {Promise<string|null>} - 2D SDF data or null if failed
 */
async function get2DFromPubChem(cid) {
    console.log(`📄 Fetching 2D structure from PubChem for CID ${cid}...`);
    try {
        const url = `https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/${cid}/SDF?record_type=2d`;
        const response = await fetch(url);
        if (response.ok) {
            const sdfData = await response.text();
            console.log(`✅ Got 2D structure (${sdfData.length} bytes)`);
            return sdfData;
        }
    } catch (e) {
        console.warn(`⚠️ Failed to get 2D from PubChem: ${e.message}`);
    }
    return null;
}

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

// =============================================================================
// MAIN ROUTING - Determine what type of molecule to load
// =============================================================================
if (cid) {
    loadMolecule(cid);
} else if (pdbid) {
    // Direct PDB ID - load biomolecule from RCSB
    loadBiomolecule(pdbid);
} else if (codid) {
    // Crystal structure - load from COD directly with 3Dmol.js
    loadMineral(codid);
} else if (name) {
    resolveAndLoad(name);
} else {
    showError('No molecule ID or Name provided');
}

async function resolveAndLoad(name) {
    console.log(`🔍 Resolving molecule: ${name}`);

    // STEP 1: Query Search API (port 8001) for intelligent molecule detection
    try {
        console.log(`🌐 Querying Search API: ${SEARCH_API}/search?q=${encodeURIComponent(name)}`);
        const response = await fetch(`${SEARCH_API}/search?q=${encodeURIComponent(name)}`);

        if (response.ok) {
            const data = await response.json();
            console.log('📦 Search API response:', data);

            // Check both 'category' and 'primary_type' fields (API uses primary_type)
            const category = (data.category || data.primary_type || '').toLowerCase();

            if (data.found !== false && category) {
                console.log(`📋 Detected type: ${category}`);

                // CASE 1: Compound → Use CID with full fallback chain
                if (category === 'compound') {
                    console.log(`✅ Found compound: ${name}`);

                    // If we have a CID, use the full loadMolecule with CIR fallback
                    if (data.cid) {
                        loadMolecule(data.cid);
                        return;
                    } else if (data.smiles || data.canonical_smiles || data.isomeric_smiles) {
                        // Try to get CID from SMILES
                        const smiles = data.smiles || data.canonical_smiles || data.isomeric_smiles;
                        console.log('🧪 Using SMILES:', smiles);
                        const cid = await getCIDFromSMILES(smiles);
                        if (cid) {
                            loadMolecule(cid);
                            return;
                        }
                        // Fallback: try CIR directly with the SMILES
                        const sdfData = await generate3DFromCIR(smiles);
                        if (sdfData) {
                            loadMoleculeFromSDF(sdfData, name);
                            return;
                        }
                    }
                }

                // CASE 2: Mineral → Load from COD directly with 3Dmol.js
                else if (category === 'mineral') {
                    console.log(`💎 Found mineral: ${name}`);
                    const mineralId = data.codid || data.cod_id || data.id;
                    if (mineralId) {
                        console.log(`%c💎 Loading mineral directly from COD: ${mineralId}`, 'background: #FF9800; color: white; font-weight: bold; padding: 4px;');
                        loadMineral(mineralId);
                        return;
                    }
                }

                // CASE 3: Protein/Biomolecule → Use direct RCSB PDB loading!
                else if (category === 'protein' || category === 'biomolecule' || category === 'pdb') {
                    console.log(`🧬 Found ${category}: ${name}`);

                    // Get PDB ID from search response
                    const pdbId = data.pdbid || data.pdb_id || data.id;
                    if (pdbId) {
                        console.log(`%c🧬 Loading biomolecule directly from RCSB: ${pdbId}`, 'background: #9C27B0; color: white; font-weight: bold; padding: 4px;');
                        loadBiomolecule(pdbId);
                        return;
                    }
                }
            }
        }
    } catch (e) {
        console.warn(`⚠️ Search API failed: ${e.message}`);
    }

    // STEP 2: Fallback to PubChem for small molecules
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

    // STEP 3: Try RCSB search for biomolecule names
    console.log('🔄 Trying RCSB search...');
    try {
        // Simple RCSB search - try common protein names
        const pdbSearchUrl = `https://search.rcsb.org/rcsbsearch/v2/query?json={"query":{"type":"terminal","service":"text","parameters":{"value":"${encodeURIComponent(name)}"}},"return_type":"entry","request_options":{"results_content_type":["experimental"],"paginate":{"start":0,"rows":1}}}`;
        const response = await fetch(pdbSearchUrl);
        if (response.ok) {
            const data = await response.json();
            if (data.result_set && data.result_set.length > 0) {
                const pdbId = data.result_set[0].identifier;
                console.log(`✅ Found PDB ID "${pdbId}" for "${name}" via RCSB search`);
                loadBiomolecule(pdbId);
                return;
            }
        }
    } catch (e) {
        console.warn(`⚠️ RCSB search failed: ${e.message}`);
    }

    // STEP 4: Final fallback to MolView
    console.warn(`⚠️ Could not resolve "${name}". Using MolView fallback.`);
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

    // Create MolView iframe pointing to LOCAL MolView instance
    const iframe = document.createElement('iframe');

    // Use local MolView /embed/v2/ endpoint with query parameter
    let param;
    if (type === 'cid' || type === true) {
        param = `cid=${query}`;
    } else if (type === 'codid') {
        param = `codid=${query}`;
    } else {
        param = `q=${encodeURIComponent(query)}`;
    }
    iframe.src = `${MOLVIEW_API}/embed/v2/?${param}`;

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

        // =============================================================================
        // FULL 3D FALLBACK CHAIN (same as MolView V2!)
        // Priority 1: PubChem 3D conformer (pre-computed)
        // Priority 2: PubChem computed 3D
        // Priority 3: CIR 3D generation from SMILES (the magic for cardiolipin!)
        // Priority 4: 2D as 3D (flat molecule - last resort)
        // Priority 5: MolView fallback (emergency only)
        // =============================================================================

        let sdfData = null;
        let fetchMethod = '';

        // STEP 1: Try PubChem 3D conformer first
        console.log(`🔍 [Step 1/4] Trying PubChem 3D conformer for CID ${cid}...`);
        try {
            const url3d = `https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/${cid}/record/SDF?record_type=3d`;
            const response3d = await fetch(url3d);

            if (response3d.ok) {
                sdfData = await response3d.text();
                fetchMethod = 'PubChem 3D conformer';
                console.log(`✅ Successfully fetched 3D conformer for CID ${cid}`);
            } else {
                console.warn(`⚠️ 3D conformer not available (status ${response3d.status})`);
            }
        } catch (e) {
            console.warn(`⚠️ 3D conformer fetch failed: ${e.message}`);
        }

        // STEP 2: Try PubChem computed 3D
        if (!sdfData) {
            console.log(`🔍 [Step 2/4] Trying PubChem computed 3D for CID ${cid}...`);
            try {
                const urlComputed = `https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/${cid}/SDF?record_type=3d`;
                const responseComputed = await fetch(urlComputed);

                if (responseComputed.ok) {
                    sdfData = await responseComputed.text();
                    fetchMethod = 'PubChem computed 3D';
                    console.log(`✅ Successfully fetched computed 3D structure for CID ${cid}`);
                }
            } catch (e) {
                console.warn(`⚠️ Computed 3D fetch failed: ${e.message}`);
            }
        }

        // STEP 3: CIR 3D generation from SMILES (THIS IS THE KEY FOR CARDIOLIPIN!)
        if (!sdfData) {
            console.log(`🔍 [Step 3/4] PubChem has no 3D - trying CIR 3D generation...`);
            console.log(`%c🔬 Using CIR to generate 3D coordinates!`, 'background: #4CAF50; color: white; font-weight: bold; padding: 4px;');

            const smiles = await getSMILESFromPubChem(cid);
            if (smiles) {
                sdfData = await generate3DFromCIR(smiles);
                if (sdfData) {
                    fetchMethod = 'CIR 3D generation (computed from SMILES)';
                    console.log(`✅ CIR successfully generated 3D structure!`);
                }
            }
        }

        // STEP 4: 2D as 3D fallback (flat molecule - better than nothing)
        if (!sdfData) {
            console.log(`🔍 [Step 4/4] CIR failed - falling back to 2D structure...`);
            console.warn(`⚠️ Using flat 2D structure as 3D (molecule will appear flat)`);

            sdfData = await get2DFromPubChem(cid);
            if (sdfData) {
                fetchMethod = '2D structure (flat - no 3D available)';
            }
        }

        // STEP 5: Final fallback to MolView (emergency only)
        if (!sdfData) {
            console.warn(`⚠️ All methods failed for CID ${cid}. Falling back to MolView...`);
            fallbackToMolView(cid, true);
            return;
        }

        console.log(`%c📦 Using: ${fetchMethod}`, 'background: #2196F3; color: white; font-weight: bold; padding: 4px;');


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
