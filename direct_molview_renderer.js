/**
 * Direct MolView Data Extractor and 3D.js Renderer
 * Extracts molecular data directly from source databases and renders with 3D.js
 */

class DirectMolViewRenderer {
    constructor() {
        this.scene = null;
        this.camera = null;
        this.renderer = null;
        this.container = null;
        this.dbAccess = null;

        // Initialize the database access if available
        if (typeof DirectChemicalDatabaseAccess !== 'undefined') {
            this.dbAccess = new DirectChemicalDatabaseAccess();
        }
    }

    /**
     * Extract molecular data directly from source databases (bypassing MolView)
     */
    async extractMolData(param, paramValue) {
        try {
            console.log(`Fetching molecular data directly from source databases: ${param}=${paramValue}`);

            // Use our direct database access class
            const molecularData = await this.dbAccess.getMolecularData(paramValue, param);

            if (!molecularData) {
                throw new Error(`Could not get molecular data for ${param}=${paramValue}`);
            }

            // Convert to standard format and enhance with protein-specific data if applicable
            let standardData = this.dbAccess.convertToStandardFormat(molecularData);

            // Enhance with protein-specific data if it's a PDB file
            if (molecularData.format === 'pdb') {
                standardData = this.dbAccess.getProteinData(molecularData);
            }

            if (!standardData || (standardData.atoms && standardData.atoms.length === 0) ||
                (!standardData.atoms && (!standardData.caAtoms || standardData.caAtoms.length === 0))) {
                throw new Error(`No valid atomic data found for ${param}=${paramValue}`);
            }

            // Return standardized molecular structure
            return {
                data: standardData,
                source: molecularData.source,
                format: molecularData.format || standardData.format,
                type: molecularData.format === 'pdb' ? 'protein' :
                      molecularData.format === 'cif' ? 'crystal' : 'molecule',
                param: param,
                paramValue: paramValue
            };
        } catch (error) {
            console.error(`Error extracting data from source databases: ${error.message}`);

            // Fallback: try to get from PubChem directly if CID
            if (param === 'cid') {
                try {
                    const pubchemUrl = `https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/${paramValue}/SDF?record_type=3d`;
                    const response = await fetch(pubchemUrl);

                    if (response.ok) {
                        const sdfData = await response.text();
                        const parsedData = this.dbAccess.parseSdf(sdfData);

                        return {
                            data: parsedData,
                            source: 'pubchem',
                            format: 'sdf',
                            type: 'molecule',
                            param: 'cid',
                            paramValue: paramValue
                        };
                    }
                } catch (fallbackError) {
                    console.warn('PubChem fallback also failed:', fallbackError.message);
                }
            }

            throw error;
        }
    }

    /**
     * Get CID from SMILES using PubChem API
     */
    async getCIDFromSmiles(smiles) {
        try {
            const url = `https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/smiles/${encodeURIComponent(smiles)}/cids/JSON`;
            const response = await fetch(url);
            
            if (response.ok) {
                const data = await response.json();
                if (data.IdentifierList && data.IdentifierList.CID && data.IdentifierList.CID.length > 0) {
                    return data.IdentifierList.CID[0].toString();
                }
            }
            return null;
        } catch (error) {
            console.error('Error getting CID from SMILES:', error);
            return null;
        }
    }

    /**
     * Parse SDF format to 3D coordinates
     */
    parseSdf(sdfData) {
        try {
            const lines = sdfData.split('\n');
            const atoms = [];
            const bonds = [];
            
            let readingAtoms = false;
            let readingBonds = false;
            let atomCount = 0;
            let bondCount = 0;
            
            for (let i = 0; i < lines.length; i++) {
                const line = lines[i].trim();
                
                // Parse counts line to get atom and bond counts
                if (/^\s*\d+\s+\d+/.test(line) && !atomCount) {
                    const parts = line.split(/\s+/).filter(p => p);
                    if (parts.length >= 2) {
                        atomCount = parseInt(parts[0]);
                        bondCount = parseInt(parts[1]);
                        readingAtoms = true;
                        continue;
                    }
                }
                
                // Parse atoms
                if (readingAtoms && atomCount > 0) {
                    if (/^\s*-?\d+(\.\d+)?\s+-?\d+(\.\d+)?\s+-?\d+(\.\d+)?/.test(line)) {
                        const parts = line.split(/\s+/).filter(p => p);
                        if (parts.length >= 4) {
                            const [x, y, z, element] = [parseFloat(parts[0]), parseFloat(parts[1]), parseFloat(parts[2]), parts[3]];
                            atoms.push({ x, y, z, element });
                            
                            if (atoms.length === atomCount) {
                                readingAtoms = false;
                                readingBonds = true;
                                continue;
                            }
                        }
                    }
                }
                
                // Parse bonds
                if (readingBonds && bondCount > 0) {
                    if (/^\s*\d+\s+\d+\s+\d/.test(line)) {
                        const parts = line.split(/\s+/).filter(p => p);
                        if (parts.length >= 3) {
                            const [from, to, order] = [parseInt(parts[0])-1, parseInt(parts[1])-1, parseInt(parts[2])]; // 0-indexed
                            bonds.push({ from, to, order });
                            
                            if (bonds.length === bondCount) {
                                break;
                            }
                        }
                    }
                }
            }
            
            return { atoms, bonds };
        } catch (error) {
            console.error('Error parsing SDF:', error);
            return { atoms: [], bonds: [] };
        }
    }

    /**
     * Parse CIF format to 3D coordinates
     */
    parseCif(cifData) {
        try {
            const lines = cifData.split('\n');
            const atoms = [];
            const cellParams = {};
            
            let inAtomsBlock = false;
            let atomLabels = [];
            
            for (let i = 0; i < lines.length; i++) {
                const line = lines[i].trim();
                
                if (line.startsWith('_cell_')) {
                    // Parse unit cell parameters
                    const parts = line.split(/\s+/);
                    if (parts.length >= 2) {
                        const param = parts[0].replace('_cell_', '').replace('.', '');
                        const value = parseFloat(parts[1]);
                        if (!isNaN(value)) {
                            cellParams[param] = value;
                        }
                    }
                }
                
                if (line.startsWith('loop_')) {
                    inAtomsBlock = false;
                }
                
                if (line.startsWith('_atom_site.') && !inAtomsBlock) {
                    atomLabels.push(line.replace('_atom_site_', '').replace('.', ''));
                }
                
                if (line === 'loop_' && lines[i+1] && lines[i+1].includes('_atom_site_label')) {
                    inAtomsBlock = true;
                    atomLabels = [];
                    continue;
                }
                
                if (inAtomsBlock && line && !line.startsWith('_') && !line.startsWith('#')) {
                    const parts = line.split(/\s+/).filter(p => p);
                    if (parts.length >= 4) { // At least label, x, y, z
                        const [label, x, y, z] = [parts[0], parseFloat(parts[2]), parseFloat(parts[3]), parseFloat(parts[4])];
                        if (!isNaN(x) && !isNaN(y) && !isNaN(z)) {
                            atoms.push({ x, y, z, element: label.replace(/\d/g, '') });
                        }
                    }
                }
            }
            
            return { atoms, cellParams };
        } catch (error) {
            console.error('Error parsing CIF:', error);
            return { atoms: [], cellParams: {} };
        }
    }

    /**
     * Create a basic 3D scene with Three.js
     */
    async create3DScene(moleculeData, containerId, options = {}) {
        // Dynamically load Three.js if not already loaded
        if (typeof THREE === 'undefined') {
            await this.loadThreeJS();
        }

        if (typeof THREE === 'undefined') {
            throw new Error('Three.js library not available');
        }

        this.container = document.getElementById(containerId);
        if (!this.container) {
            throw new Error(`Container ${containerId} not found`);
        }

        // Clear container
        this.container.innerHTML = '';

        // Create scene
        this.scene = new THREE.Scene();
        this.scene.background = new THREE.Color(options.backgroundColor || 0x1a1a2e);

        // Create camera
        this.camera = new THREE.PerspectiveCamera(
            options.fov || 75,
            this.container.clientWidth / this.container.clientHeight,
            options.near || 0.1,
            options.far || 1000
        );
        this.camera.position.set(0, 0, options.cameraDistance || 15);

        // Create renderer
        this.renderer = new THREE.WebGLRenderer({
            antialias: true,
            alpha: true
        });
        this.renderer.setSize(this.container.clientWidth, this.container.clientHeight);
        this.renderer.setPixelRatio(window.devicePixelRatio);
        this.container.appendChild(this.renderer.domElement);

        // Add lights
        const ambientLight = new THREE.AmbientLight(0x404040, options.ambientIntensity || 1.5);
        this.scene.add(ambientLight);

        const directionalLight = new THREE.DirectionalLight(0xffffff, options.directionalIntensity || 1);
        directionalLight.position.set(5, 5, 5);
        this.scene.add(directionalLight);

        const pointLight = new THREE.PointLight(0xffffff, options.pointIntensity || 1, 100);
        pointLight.position.set(-5, -5, -5);
        this.scene.add(pointLight);

        // Add molecule to scene based on options
        this.addMoleculeToScene(moleculeData, options);

        // Add controls for rotation
        this.addControls();

        // Start animation loop
        this.animate();
    }

    /**
     * Add molecule to scene with appropriate representation
     */
    addMoleculeToScene(moleculeData, options = {}) {
        if (!this.scene || !moleculeData || !moleculeData.atoms) {
            console.error('Invalid molecular data or scene');
            return;
        }

        // Determine the type of molecule based on available data
        const isProtein = (moleculeData.type === 'protein') ||
                         (moleculeData.format === 'pdb') ||
                         (moleculeData.caAtoms && moleculeData.caAtoms.length > 0) ||
                         options.moleculeType === 'protein';

        const isCrystal = (moleculeData.type === 'crystal') ||
                         (moleculeData.format === 'cif') ||
                         options.moleculeType === 'crystal';

        const atomCount = moleculeData.atoms ? moleculeData.atoms.length : 0;
        const isMacromolecule = atomCount > 100 || isProtein;

        // Create different representations based on molecule type
        if (isProtein) {
            // Use protein-specific representations
            this.addProteinToScene(moleculeData, options);
        } else if (isCrystal) {
            // Use crystal-specific representations
            this.addCrystalToScene(moleculeData, options);
        } else if (isMacromolecule) {
            // Use protein-style representations for large molecules even if not explicitly marked as protein
            this.addProteinToScene(moleculeData, options);
        } else {
            // Use standard atom/bond representation for small molecules
            this.addSmallMoleculeToScene(moleculeData, options);
        }
    }

    /**
     * Add small molecule with atom/bond representation
     */
    addSmallMoleculeToScene(moleculeData, options = {}) {
        const { atoms, bonds } = moleculeData;

        // Common atom colors based on element
        const atomColors = {
            'C': 0x1a1a1a,  // Dark gray
            'H': 0xffffff,  // White
            'O': 0xff0000,  // Red
            'N': 0x0000ff,  // Blue
            'S': 0xffff00,  // Yellow
            'P': 0xffa500,  // Orange
            'default': 0x808080  // Gray
        };

        // Create atom spheres and bond cylinders
        const atomGroup = new THREE.Group();
        const bondGroup = new THREE.Group();

        // Store atom positions for bond creation
        const atomPositions = [];

        // Create atoms
        atoms.forEach((atom, index) => {
            const color = atomColors[atom.element] || atomColors.default;
            const radius = options.atomRadius || (atom.element === 'H' ? 0.3 : 0.5);
            const geometry = new THREE.SphereGeometry(radius, 32, 32);
            const material = new THREE.MeshPhongMaterial({ color: color });
            const sphere = new THREE.Mesh(geometry, material);

            sphere.position.set(atom.x, atom.y, atom.z);
            atomGroup.add(sphere);

            // Store position for bond creation
            atomPositions.push(new THREE.Vector3(atom.x, atom.y, atom.z));
        });

        // Create bonds
        if (bonds && bonds.length > 0) {
            bonds.forEach(bond => {
                if (bond.from < atomPositions.length && bond.to < atomPositions.length) {
                    this.addBondToScene(bondGroup, atomPositions[bond.from], atomPositions[bond.to], bond.order);
                }
            });
        } else {
            // Estimate bonds based on proximity
            for (let i = 0; i < atomPositions.length; i++) {
                for (let j = i + 1; j < atomPositions.length; j++) {
                    const distance = atomPositions[i].distanceTo(atomPositions[j]);
                    // Connect atoms within bonding distance (adjustable)
                    const maxBondDistance = options.maxBondDistance || 2.0;
                    if (distance < maxBondDistance) {
                        this.addBondToScene(bondGroup, atomPositions[i], atomPositions[j], 1);
                    }
                }
            }
        }

        this.scene.add(atomGroup);
        this.scene.add(bondGroup);

        // Center the molecule
        this.centerMolecule(atomGroup, bondGroup);
    }

    /**
     * Add protein with protein-specific representations (ribbon, cartoon, etc.)
     */
    addProteinToScene(moleculeData, options = {}) {
        // Determine the representation type
        const representation = options.representation || 'ribbon'; // ribbon, cartoon, cylinders, trace

        // Ensure we have protein-specific data
        const proteinData = this.ensureProteinData(moleculeData);

        switch (representation) {
            case 'ribbon':
                this.addRibbonRepresentation(proteinData, options);
                break;
            case 'cartoon':
                this.addCartoonRepresentation(proteinData, options);
                break;
            case 'cylinders':
                this.addCylinderPlateRepresentation(proteinData, options);
                break;
            case 'trace':
                this.addCAlphaTraceRepresentation(proteinData, options);
                break;
            case 'bonds':
            default:
                this.addSmallMoleculeToScene(proteinData, options);
                break;
        }
    }

    /**
     * Add crystal structure visualization
     */
    addCrystalToScene(moleculeData, options = {}) {
        // For crystals, focus on atoms and unit cell
        const representation = options.representation || 'bonds';

        switch (representation) {
            case 'unit-cell':
                this.addUnitCellRepresentation(moleculeData, options);
                break;
            default:
                this.addSmallMoleculeToScene(moleculeData, { ...options, maxBondDistance: 3.0 }); // Crystals have longer bonds
                break;
        }
    }

    /**
     * Ensure we have proper protein data with C-alpha and secondary structure info
     */
    ensureProteinData(moleculeData) {
        if (!moleculeData.caAtoms) {
            // Find C-alpha atoms from the atom data
            const caAtoms = moleculeData.atoms?.filter(a =>
                typeof a === 'object' && a.name && a.name.toUpperCase() === 'CA'
            ) || [];

            return {
                ...moleculeData,
                caAtoms: caAtoms,
                backboneAtoms: moleculeData.atoms?.filter(a =>
                    typeof a === 'object' && a.name && ['N', 'CA', 'C', 'O'].includes(a.name.toUpperCase())
                ) || [],
                secondaryStructure: {
                    helices: moleculeData.helix || [],
                    sheets: moleculeData.sheet || []
                }
            };
        }
        return moleculeData;
    }

    /**
     * Add unit cell representation for crystals
     */
    addUnitCellRepresentation(crystalData, options = {}) {
        // Create a representation highlighting the crystal structure
        this.addSmallMoleculeToScene(crystalData, {
            ...options,
            atomRadius: options.atomRadius || 0.3,
            maxBondDistance: options.maxBondDistance || 4.0 // For crystal structures, larger distance
        });

        // If there are unit cell parameters, draw the cell boundaries
        if (crystalData.cellParams) {
            this.addUnitCellBox(crystalData.cellParams);
        }
    }

    /**
     * Add unit cell box visualization
     */
    addUnitCellBox(cellParams) {
        if (!cellParams.a || !cellParams.b || !cellParams.c) return;

        // Create wireframe box representing the unit cell
        const boxGeometry = new THREE.BoxGeometry(cellParams.a, cellParams.b, cellParams.c);
        const edges = new THREE.EdgesGeometry(boxGeometry);
        const lineMaterial = new THREE.LineBasicMaterial({ color: 0x00ffff, linewidth: 2 });
        const wireframe = new THREE.LineSegments(edges, lineMaterial);

        // Position the box
        wireframe.position.set(
            cellParams.a / 2,
            cellParams.b / 2,
            cellParams.c / 2
        );

        this.scene.add(wireframe);
    }

    /**
     * Add ribbon representation for proteins using C-alpha trace and tube
     */
    addRibbonRepresentation(proteinData, options = {}) {
        // Use the enhanced protein data that includes C-alpha atoms
        const caAtoms = proteinData.caAtoms && proteinData.caAtoms.length > 0 ?
                       proteinData.caAtoms :
                       proteinData.atoms?.filter(atom => atom.name === 'CA') || [];

        if (caAtoms.length < 2) {
            // If no C-alpha atoms found, fall back to backbone trace
            this.addBackboneTrace(proteinData, options);
            return;
        }

        // Create points from C-alpha atoms
        const points = caAtoms.map(atom => new THREE.Vector3(atom.x, atom.y, atom.z));

        if (points.length < 3) {
            // If we don't have enough points for a tube, fall back to simple trace
            this.addCAlphaTraceRepresentation(proteinData, options);
            return;
        }

        // Create a smooth curve through the C-alpha positions
        const curve = new THREE.CatmullRomCurve3(points, false, 'centripetal', 0.5);

        // Create the ribbon/tube geometry
        const tubeGeometry = new THREE.TubeGeometry(
            curve,
            Math.min(points.length * 2, 128),  // tubular segments
            options.ribbonRadius || 0.3,      // radius of tube
            8,                                // radial segments
            false                             // closed
        );

        // Create material with proper color based on scheme
        const material = new THREE.MeshPhongMaterial({
            color: this.getColorForScheme(options.colorScheme || 'secondary', 0, 1),
            side: THREE.DoubleSide,
            shininess: 100
        });

        // Create and add the ribbon mesh
        const ribbon = new THREE.Mesh(tubeGeometry, material);
        this.scene.add(ribbon);
    }

    /**
     * Add cartoon representation showing secondary structure elements (helices, sheets)
     */
    addCartoonRepresentation(proteinData, options = {}) {
        const caAtoms = proteinData.caAtoms && proteinData.caAtoms.length > 0 ?
                       proteinData.caAtoms :
                       proteinData.atoms?.filter(atom => atom.name === 'CA') || [];

        if (caAtoms.length < 2) {
            this.addRibbonRepresentation(proteinData, options);
            return;
        }

        // This is a simplified cartoon representation
        // In a full implementation, we'd use the secondary structure data
        // to distinguish between helices (cylinders), sheets (arrows), and loops

        // Group consecutive C-alpha atoms by chain
        const chains = {};
        caAtoms.forEach(atom => {
            const chainId = atom.chainId || 'A';
            if (!chains[chainId]) {
                chains[chainId] = [];
            }
            chains[chainId].push(atom);
        });

        // Create cartoon representation for each chain
        Object.entries(chains).forEach(([chainId, chainAtoms]) => {
            if (chainAtoms.length < 3) return; // Skip small chains

            // Sort atoms by residue number
            chainAtoms.sort((a, b) => a.resSeq - b.resSeq);

            // For now, simplify by creating a ribbon for the whole chain
            const points = chainAtoms.map(atom => new THREE.Vector3(atom.x, atom.y, atom.z));
            const curve = new THREE.CatmullRomCurve3(points, false, 'centripetal', 0.5);

            const tubeGeometry = new THREE.TubeGeometry(
                curve,
                Math.min(points.length * 2, 64),
                options.ribbonRadius || 0.4,
                6,
                false
            );

            const material = new THREE.MeshPhongMaterial({
                color: this.getColorForScheme(options.colorScheme || 'chain', chainId.charCodeAt(0), 255),
                side: THREE.DoubleSide,
                shininess: 100
            });

            const cartoon = new THREE.Mesh(tubeGeometry, material);
            this.scene.add(cartoon);
        });
    }

    /**
     * Add cylinder and plate representation (α-helices as cylinders, β-sheets as flat plates)
     */
    addCylinderPlateRepresentation(proteinData, options = {}) {
        // Use secondary structure information if available
        const secondaryStructure = proteinData.secondaryStructure || {};
        const helices = secondaryStructure.helices || [];
        const sheets = secondaryStructure.sheets || [];

        const caAtoms = proteinData.caAtoms && proteinData.caAtoms.length > 0 ?
                       proteinData.caAtoms :
                       proteinData.atoms?.filter(atom => atom.name === 'CA') || [];

        if (caAtoms.length < 2) {
            this.addRibbonRepresentation(proteinData, options);
            return;
        }

        // Create geometric representations for secondary structures
        const cartoonGroup = new THREE.Group();

        // For now, fall back to simplified cartoon representation
        this.addCartoonRepresentation(proteinData, options);
    }

    /**
     * Add C-alpha trace representation (simple line connecting C-alpha atoms)
     */
    addCAlphaTraceRepresentation(proteinData, options = {}) {
        const caAtoms = proteinData.caAtoms && proteinData.caAtoms.length > 0 ?
                       proteinData.caAtoms :
                       proteinData.atoms?.filter(atom => atom.name === 'CA') || [];

        if (caAtoms.length < 2) {
            // Fallback to basic representation
            this.addSmallMoleculeToScene(proteinData, options);
            return;
        }

        // Create points from C-alpha atoms
        const points = caAtoms.map(atom => new THREE.Vector3(atom.x, atom.y, atom.z));

        // Create line geometry connecting consecutive C-alpha atoms
        const geometry = new THREE.BufferGeometry().setFromPoints(points);
        const material = new THREE.LineBasicMaterial({
            color: 0xffffff,
            linewidth: options.lineWidth || 2
        });

        const line = new THREE.Line(geometry, material);
        this.scene.add(line);
    }

    /**
     * Add backbone trace (N, CA, C atoms in sequence)
     */
    addBackboneTrace(proteinData, options = {}) {
        const backboneAtoms = proteinData.backboneAtoms && proteinData.backboneAtoms.length > 0 ?
                             proteinData.backboneAtoms :
                             proteinData.atoms?.filter(atom => ['N', 'CA', 'C'].includes(atom.name)) || [];

        if (backboneAtoms.length < 3) {
            // Not enough backbone atoms, fall back to simple representation
            this.addSmallMoleculeToScene(proteinData, options);
            return;
        }

        // Sort backbone atoms by residue number to ensure proper connectivity
        const sortedAtoms = [...backboneAtoms].sort((a, b) => {
            if (a.chainId !== b.chainId) return a.chainId.localeCompare(b.chainId);
            return (a.resSeq || 0) - (b.resSeq || 0);
        });

        // Create points from sorted backbone atoms
        const points = sortedAtoms.map(atom => new THREE.Vector3(atom.x, atom.y, atom.z));

        // Create line geometry
        const geometry = new THREE.BufferGeometry().setFromPoints(points);
        const material = new THREE.LineBasicMaterial({
            color: 0xcccccc,
            linewidth: options.lineWidth || 1.5
        });

        const line = new THREE.Line(geometry, material);
        this.scene.add(line);
    }

    /**
     * Helper to add a bond between two points
     */
    addBondToScene(group, start, end, order = 1) {
        const direction = new THREE.Vector3().subVectors(end, start);
        const length = direction.length();
        const bondGeometry = new THREE.CylinderGeometry(0.1, 0.1, length, 8);

        const bondMaterial = new THREE.MeshPhongMaterial({ color: 0x808080 });
        const cylinder = new THREE.Mesh(bondGeometry, bondMaterial);

        // Position and orient the cylinder
        cylinder.position.set(
            start.x + direction.x / 2,
            start.y + direction.y / 2,
            start.z + direction.z / 2
        );

        // Rotate to align with direction
        cylinder.quaternion.setFromUnitVectors(
            new THREE.Vector3(0, 1, 0),  // Default up direction
            direction.clone().normalize()
        );

        group.add(cylinder);
    }

    /**
     * Helper to get color based on scheme
     */
    getColorForScheme(scheme, value, maxVal) {
        // Simplified color scheme implementation
        switch (scheme) {
            case 'spectrum': // Rainbow
                return new THREE.Color(value / maxVal, 0.5, 1 - (value / maxVal));
            case 'chain':
                return new THREE.Color(0.5 + (value / maxVal) * 0.5, 0.3, 0.7 - (value / maxVal) * 0.3);
            case 'residue':
                return new THREE.Color(0.3, 0.5 + (value / maxVal) * 0.3, 0.6);
            case 'polarity':
                // Red for polar, white for non-polar (simplified)
                return value > maxVal/2 ? new THREE.Color(1, 0, 0) : new THREE.Color(1, 1, 1);
            case 'secondary':
            default: // Color by secondary structure
                return new THREE.Color(0.5, 0.7, 1.0);
        }
    }

    /**
     * Center the molecule in the scene
     */
    centerMolecule(atomGroup, bondGroup) {
        const combined = new THREE.Group();
        combined.add(atomGroup);
        combined.add(bondGroup);

        // Compute bounding box to center the molecule
        const box = new THREE.Box3().setFromObject(combined);
        const center = box.getCenter(new THREE.Vector3());

        // Adjust positions to center at origin
        atomGroup.position.sub(center);
        bondGroup.position.sub(center);
    }

    /**
     * Load Three.js dynamically if not present
     */
    async loadThreeJS() {
        return new Promise((resolve, reject) => {
            if (typeof THREE !== 'undefined') {
                resolve();
                return;
            }
            
            const script = document.createElement('script');
            script.src = 'https://cdnjs.cloudflare.com/ajax/libs/three.js/r128/three.min.js';
            script.onload = () => {
                console.log('Three.js loaded');
                resolve();
            };
            script.onerror = () => {
                reject(new Error('Failed to load Three.js'));
            };
            document.head.appendChild(script);
        });
    }

    /**
     * Add molecule to Three.js scene
     */
    addMoleculeToScene(moleculeData) {
        if (!this.scene || !moleculeData || !moleculeData.atoms) {
            console.error('Invalid molecular data or scene');
            return;
        }
        
        const { atoms, bonds } = moleculeData;
        
        // Common atom colors based on element
        const atomColors = {
            'C': 0x1a1a1a,  // Dark gray
            'H': 0xffffff,  // White
            'O': 0xff0000,  // Red
            'N': 0x0000ff,  // Blue
            'S': 0xffff00,  // Yellow
            'P': 0xffa500,  // Orange
            'default': 0x808080  // Gray
        };
        
        // Create atom spheres and bond cylinders
        const atomGroup = new THREE.Group();
        const bondGroup = new THREE.Group();
        
        // Store atom positions for bond creation
        const atomPositions = [];
        
        // Create atoms
        atoms.forEach((atom, index) => {
            const color = atomColors[atom.element] || atomColors.default;
            const geometry = new THREE.SphereGeometry(0.5, 32, 32);
            const material = new THREE.MeshPhongMaterial({ color: color });
            const sphere = new THREE.Mesh(geometry, material);
            
            sphere.position.set(atom.x, atom.y, atom.z);
            atomGroup.add(sphere);
            
            // Store position for bond creation
            atomPositions.push(new THREE.Vector3(atom.x, atom.y, atom.z));
        });
        
        // Create bonds (if available)
        if (bonds && bonds.length > 0) {
            bonds.forEach(bond => {
                if (bond.from < atomPositions.length && bond.to < atomPositions.length) {
                    const start = atomPositions[bond.from];
                    const end = atomPositions[bond.to];
                    
                    // Create cylinder for bond
                    const direction = new THREE.Vector3().subVectors(end, start);
                    const length = direction.length();
                    const bondGeometry = new THREE.CylinderGeometry(0.1, 0.1, length, 8);
                    
                    const bondMaterial = new THREE.MeshPhongMaterial({ color: 0x808080 });
                    const cylinder = new THREE.Mesh(bondGeometry, bondMaterial);
                    
                    // Position and orient the cylinder
                    cylinder.position.set(
                        start.x + direction.x / 2,
                        start.y + direction.y / 2,
                        start.z + direction.z / 2
                    );
                    
                    // Rotate to align with direction
                    cylinder.quaternion.setFromUnitVectors(
                        new THREE.Vector3(0, 1, 0),  // Default up direction
                        direction.clone().normalize()
                    );
                    
                    bondGroup.add(cylinder);
                }
            });
        } else {
            // If no bonds data, estimate bonds based on atom proximity
            for (let i = 0; i < atomPositions.length; i++) {
                for (let j = i + 1; j < atomPositions.length; j++) {
                    const distance = atomPositions[i].distanceTo(atomPositions[j]);
                    // Connect atoms that are within 2.0 units (typical bond distance)
                    if (distance < 2.0) {
                        const direction = new THREE.Vector3().subVectors(atomPositions[j], atomPositions[i]);
                        const length = direction.length();
                        const bondGeometry = new THREE.CylinderGeometry(0.05, 0.05, length, 8);
                        
                        const bondMaterial = new THREE.MeshPhongMaterial({ color: 0x808080 });
                        const cylinder = new THREE.Mesh(bondGeometry, bondMaterial);
                        
                        cylinder.position.set(
                            atomPositions[i].x + direction.x / 2,
                            atomPositions[i].y + direction.y / 2,
                            atomPositions[i].z + direction.z / 2
                        );
                        
                        cylinder.quaternion.setFromUnitVectors(
                            new THREE.Vector3(0, 1, 0),
                            direction.clone().normalize()
                        );
                        
                        bondGroup.add(cylinder);
                    }
                }
            }
        }
        
        this.scene.add(atomGroup);
        this.scene.add(bondGroup);
        
        // Position the molecule at the center
        const center = new THREE.Vector3();
        let count = 0;
        
        atomGroup.children.forEach(child => {
            center.add(child.position);
            count++;
        });
        
        if (count > 0) {
            center.divideScalar(count);
            atomGroup.position.sub(center);
            bondGroup.position.sub(center);
        }
    }

    /**
     * Add controls for rotating the molecule
     */
    addControls() {
        if (typeof THREE.OrbitControls === 'undefined') {
            // Try to load OrbitControls
            this.loadOrbitControls().then(() => {
                this.setupOrbitControls();
            }).catch(() => {
                // Fallback: add basic mouse controls
                this.setupBasicControls();
            });
        } else {
            this.setupOrbitControls();
        }
    }

    /**
     * Load OrbitControls dynamically
     */
    async loadOrbitControls() {
        return new Promise((resolve, reject) => {
            const script = document.createElement('script');
            script.src = 'https://cdnjs.cloudflare.com/ajax/libs/three.js/r128/examples/js/controls/OrbitControls.js';
            script.onload = () => {
                console.log('OrbitControls loaded');
                resolve();
            };
            script.onerror = () => {
                reject(new Error('Failed to load OrbitControls'));
            };
            document.head.appendChild(script);
        });
    }

    /**
     * Setup OrbitControls
     */
    setupOrbitControls() {
        if (typeof THREE.OrbitControls !== 'undefined' && this.camera && this.renderer) {
            const controls = new THREE.OrbitControls(this.camera, this.renderer.domElement);
            controls.enableDamping = true;
            controls.dampingFactor = 0.05;
            controls.screenSpacePanning = false;
            controls.minDistance = 1;
            controls.maxDistance = 100;
        } else {
            this.setupBasicControls();
        }
    }

    /**
     * Setup basic rotation controls
     */
    setupBasicControls() {
        let isDragging = false;
        let previousMousePosition = {
            x: 0,
            y: 0
        };

        this.renderer.domElement.addEventListener('mousedown', (e) => {
            isDragging = true;
            previousMousePosition = {
                x: e.clientX,
                y: e.clientY
            };
        });

        this.renderer.domElement.addEventListener('mousemove', (e) => {
            if (isDragging) {
                const deltaMove = {
                    x: e.clientX - previousMousePosition.x,
                    y: e.clientY - previousMousePosition.y
                };

                // Rotate the whole scene
                this.scene.rotation.y += deltaMove.x * 0.01;
                this.scene.rotation.x += deltaMove.y * 0.01;

                previousMousePosition = {
                    x: e.clientX,
                    y: e.clientY
                };
            }
        });

        this.renderer.domElement.addEventListener('mouseup', () => {
            isDragging = false;
        });

        this.renderer.domElement.addEventListener('mouseleave', () => {
            isDragging = false;
        });
    }

    /**
     * Animation loop
     */
    animate() {
        requestAnimationFrame(() => this.animate());
        
        // Rotate slowly if no user interaction
        if (!this.isRotating) {
            this.scene.rotation.y += 0.001;
        }
        
        this.renderer.render(this.scene, this.camera);
    }

    /**
     * Main method to render molecule directly from source database data (bypassing MolView)
     */
    async renderDirectMolecule(param, paramValue, containerId, options = {}) {
        try {
            console.log(`Extracting molecular data for ${param}=${paramValue}`);

            // Determine if this is likely a protein/macromolecule based on param and paramValue
            let moleculeType = 'small';
            let representation = 'bonds'; // Default

            if (param === 'pdbid' || paramValue.length === 4 ||
                paramValue.toLowerCase().includes('protein') ||
                paramValue.toLowerCase().includes('enzyme') ||
                paramValue.toLowerCase().includes('antibody') ||
                paramValue.toLowerCase().includes('virus')) {
                moleculeType = 'protein';
                representation = 'ribbon'; // Default for proteins
            } else if (param === 'codid' || parseInt(paramValue) > 100000) {
                moleculeType = 'crystal';
                representation = 'bonds';
            }

            // Override with user options if provided
            representation = options.representation || representation;

            // Extract the molecular data directly from source databases
            const molData = await this.extractMolData(param, paramValue);

            // The data is already parsed by the database access layer
            const parsedData = molData.data;

            // Create 3D visualization with appropriate options
            const renderOptions = {
                ...options,
                moleculeType: moleculeType,
                representation: representation,
                backgroundColor: options.backgroundColor || 0x1a1a2e,
                ambientIntensity: options.ambientIntensity || 1.5,
                directionalIntensity: options.directionalIntensity || 1,
                pointIntensity: options.pointIntensity || 1,
                atomRadius: options.atomRadius || 0.5,
                maxBondDistance: options.maxBondDistance || 2.0,
                ribbonRadius: options.ribbonRadius || 0.5,
                colorScheme: options.colorScheme || 'secondary'
            };

            await this.create3DScene(parsedData, containerId, renderOptions);

            console.log(`Successfully rendered ${param}=${paramValue} in 3D from source database`);
            return true;

        } catch (error) {
            console.error('Error rendering molecule:', error);

            // Update container with error message
            const container = document.getElementById(containerId);
            if (container) {
                container.innerHTML = `
                    <div style="display: flex; align-items: center; justify-content: center; height: 100%; background: #1a1a2e; color: white; flex-direction: column;">
                        <h3>Error Rendering Molecule</h3>
                        <p>${error.message}</p>
                        <p>Could not fetch or parse data from source databases.</p>
                        <button onclick="location.reload()" style="padding: 10px 20px; background: #667eea; color: white; border: none; border-radius: 4px; margin-top: 15px; cursor: pointer;">
                            Retry
                        </button>
                    </div>
                `;
            }

            return false;
        }
    }
}

// Export for use in other modules
if (typeof window !== 'undefined') {
    window.DirectMolViewRenderer = DirectMolViewRenderer;
}