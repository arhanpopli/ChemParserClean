// Get CID from URL parameters
const urlParams = new URLSearchParams(window.location.search);
const cid = urlParams.get('cid');
const name = urlParams.get('name');
const style = urlParams.get('style') || 'stick:sphere';
const stickRadius = parseFloat(urlParams.get('stickRadius') || '0.15');
const sphereRadius = parseFloat(urlParams.get('sphereRadius') || '0.3');
const autoRotate = urlParams.get('autoRotate') !== 'false';
const bgColor = decodeURIComponent(urlParams.get('bgColor') || '#1a1a2e');

// Set body background color from URL parameter
document.documentElement.style.setProperty('--bg-color', bgColor);

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

if (!cid) {
    showError('No molecule ID provided');
} else {
    loadMolecule(cid);
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

            // Clear the viewer container
            document.body.innerHTML = '';
            document.body.style.margin = '0';
            document.body.style.padding = '0';
            document.body.style.overflow = 'hidden';

            // Create MolView iframe
            const iframe = document.createElement('iframe');
            // Use embed mode with balls and sticks, hiding the toolbar if possible
            iframe.src = `https://embed.molview.org/v1/?mode=balls&cid=${cid}`;
            iframe.style.width = '100vw';
            iframe.style.height = '100vh';
            iframe.style.border = 'none';
            iframe.style.display = 'block';

            document.body.appendChild(iframe);

            console.log(`✅ Switched to MolView for CID ${cid}`);
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
