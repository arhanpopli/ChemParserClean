/**
 * Direct Chemical Database Access for Molecular Visualization
 * Bypasses MolView entirely by accessing source databases directly
 */

class DirectChemicalDatabaseAccess {
    constructor() {
        this.databases = {
            // Small molecule databases
            pubchem: {
                name: 'PubChem',
                baseUrl: 'https://pubchem.ncbi.nlm.nih.gov/rest/pug',
                endpoints: {
                    cid: '/compound/cid/{id}/SDF?record_type=3d',
                    name: '/compound/name/{name}/SDF?record_type=3d',
                    smiles: '/compound/smiles/{smiles}/SDF?record_type=3d'
                },
                categories: ['small molecules', 'organic compounds', 'inorganic compounds', 'drugs', 'metabolites']
            },

            chemspider: {
                name: 'ChemSpider',
                baseUrl: 'http://www.chemspider.com',
                endpoints: {
                    id: '/Chemical-Structure.{id}.sdf',
                    smiles: '/Search.asmx/SearchBySmiles?smiles={smiles}'
                },
                categories: ['small molecules', 'organic compounds']
            },

            chembl: {
                name: 'ChEMBL',
                baseUrl: 'https://www.ebi.ac.uk/chembl/api/data',
                endpoints: {
                    molecule: '/molecule/{id}.json'
                },
                categories: ['bioactive molecules', 'drug-like compounds', 'pharmaceuticals']
            },

            drugbank: {
                name: 'DrugBank',
                baseUrl: 'https://go.drugbank.com',
                endpoints: {
                    drug: '/structures/{id}.sdf'
                },
                categories: ['drugs', 'drug targets', 'pharmaceutical compounds']
            },

            // Macromolecule databases
            pdb: {
                name: 'RCSB Protein Data Bank',
                baseUrl: 'https://files.rcsb.org/download',
                endpoints: {
                    pdbid: '/{id}.pdb',
                    entry: '/{id}-assembly1.cif'  // For biological assemblies
                },
                categories: ['proteins', 'nucleic acids', 'enzymes', 'antibodies', 'viruses', 'ribosomes', 'complexes']
            },

            alphafolddb: {
                name: 'AlphaFold DB',
                baseUrl: 'https://alphafold.ebi.ac.uk/files',
                endpoints: {
                    afid: '/AF-{id}-F1-model_v{version}.cif'
                },
                categories: ['predicted proteins', 'alphafold structures']
            },

            // Crystal structure databases
            cod: {
                name: 'Crystallography Open Database',
                baseUrl: 'https://www.crystallography.net/cod',
                endpoints: {
                    codid: '/{id}.cif'
                },
                categories: ['crystal structures', 'minerals', 'solid state']
            },

            icdd: {
                name: 'ICDD (International Centre for Diffraction Data)',
                baseUrl: 'https://www.icdd.com',
                endpoints: {
                    pdf: '/pdfs/{id}.cif'
                },
                categories: ['crystallographic data', 'powder diffraction', 'minerals']
            },

            // Spectral databases
            nist: {
                name: 'NIST Chemistry WebBook',
                baseUrl: 'https://webbook.nist.gov/cgi/inchi',
                endpoints: {
                    inchi: '?ID={inchi}&Mask=4'
                },
                categories: ['spectra', 'thermophysical data', 'mass spec', 'IR spec', 'NMR']
            },

            spectrum: {
                name: 'SDBS (Spectral Database for Organic Compounds)',
                baseUrl: 'https://sdbs.db.aist.go.jp',
                endpoints: {
                    spectrum: '/sdbs/cgi-bin/direct_frame_top.cgi?sdbsno={id}'
                },
                categories: ['spectra', 'mass spec', 'IR spec', 'NMR', 'UV-Vis']
            }
        };
    }

    /**
     * Get compound data by CID from PubChem
     */
    async getFromPubChem(cid) {
        try {
            // Try 3D conformer first
            const url = `${this.databases.pubchem.baseUrl}${this.databases.pubchem.endpoints.cid.replace('{id}', cid)}`;
            console.log(`Fetching from PubChem: ${url}`);
            
            const response = await fetch(url);
            
            if (!response.ok) {
                // Try computed 3D
                const computedUrl = `https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/${cid}/SDF?record_type=3d`;
                const computedResponse = await fetch(computedUrl);
                
                if (computedResponse.ok) {
                    return {
                        data: await computedResponse.text(),
                        source: 'pubchem',
                        format: 'sdf',
                        cid: cid
                    };
                }
                
                throw new Error(`PubChem request failed: ${response.status} ${response.statusText}`);
            }
            
            return {
                data: await response.text(),
                source: 'pubchem',
                format: 'sdf',
                cid: cid
            };
        } catch (error) {
            console.error(`Error getting data from PubChem for CID ${cid}:`, error);
            return null;
        }
    }

    /**
     * Get compound by name from PubChem
     */
    async getFromPubChemByName(name) {
        try {
            // First get CID from name
            const cidUrl = `${this.databases.pubchem.baseUrl}/compound/name/${encodeURIComponent(name)}/cids/JSON`;
            console.log(`Getting CID for: ${name} from ${cidUrl}`);
            
            const cidResponse = await fetch(cidUrl);
            if (!cidResponse.ok) {
                throw new Error(`Could not find compound: ${name}`);
            }
            
            const cidData = await cidResponse.json();
            if (!cidData.IdentifierList || !cidData.IdentifierList.CID || cidData.IdentifierList.CID.length === 0) {
                throw new Error(`No CID found for: ${name}`);
            }
            
            const cid = cidData.IdentifierList.CID[0];
            console.log(`Found CID: ${cid} for name: ${name}`);
            
            // Now get the SDF data using the CID
            return await this.getFromPubChem(cid);
        } catch (error) {
            console.error(`Error getting data from PubChem for name ${name}:`, error);
            return null;
        }
    }

    /**
     * Get protein structure by PDB ID
     */
    async getFromPDB(pdbid) {
        try {
            const url = `${this.databases.pdb.baseUrl}${this.databases.pdb.endpoints.pdbid.replace('{id}', pdbid)}`;
            console.log(`Fetching from PDB: ${url}`);
            
            const response = await fetch(url);
            
            if (!response.ok) {
                throw new Error(`PDB request failed: ${response.status} ${response.statusText}`);
            }
            
            return {
                data: await response.text(),
                source: 'pdb',
                format: 'pdb',
                pdbid: pdbid
            };
        } catch (error) {
            console.error(`Error getting data from PDB for ${pdbid}:`, error);
            return null;
        }
    }

    /**
     * Get crystal structure from COD
     */
    async getFromCOD(codid) {
        try {
            const url = `${this.databases.cod.baseUrl}${this.databases.cod.endpoints.codid.replace('{id}', codid)}`;
            console.log(`Fetching from COD: ${url}`);
            
            const response = await fetch(url);
            
            if (!response.ok) {
                throw new Error(`COD request failed: ${response.status} ${response.statusText}`);
            }
            
            return {
                data: await response.text(),
                source: 'cod',
                format: 'cif',
                codid: codid
            };
        } catch (error) {
            console.error(`Error getting data from COD for ${codid}:`, error);
            return null;
        }
    }

    /**
     * Detect molecule type based on identifier or content
     */
    detectMoleculeType(identifier) {
        // Check if it's a protein-related identifier
        if (typeof identifier === 'string') {
            // PDB format: 4 alphanumeric characters, usually starts with digit
            if (identifier.length === 4 && /^[1-9][A-Za-z0-9]{3}$/.test(identifier)) {
                return { type: 'protein', database: 'pdb' };
            }

            // AlphaFold format: E.g., Q5VSL9 (UniProt-like)
            if (/^[A-NR-Z][0-9][A-Z][A-Z0-9]{2}[0-9]$/.test(identifier) ||
                /^[A-Z0-9]{1,6}$/.test(identifier)) {
                return { type: 'protein', database: 'alphafolddb' };
            }

            // Check for keywords in the identifier
            const lowerId = identifier.toLowerCase();
            if (lowerId.includes('enzyme') || lowerId.includes('kinase') ||
                lowerId.includes('antibody') || lowerId.includes('igg') ||
                lowerId.includes('virus') || lowerId.includes('protein')) {
                return { type: 'protein', database: 'pdb' };
            }

            // Check against known protein databases
            if (lowerId.startsWith('af-')) {  // AlphaFold DB format
                return { type: 'protein', database: 'alphafolddb' };
            }
        }

        // For numbers, determine based on range
        if (typeof identifier === 'number' || /^\d+$/.test(identifier)) {
            const numId = parseInt(identifier);
            if (numId > 1000000) { // Large numbers might be other DB IDs
                // Could be from various databases, check for context
                return { type: 'molecule', database: 'pubchem' }; // Default to small molecules
            } else if (numId > 100000) { // Medium-large numbers
                return { type: 'crystal', database: 'cod' };
            } else { // Smaller numbers are likely PubChem CIDs
                return { type: 'molecule', database: 'pubchem' };
            }
        }

        // Default to molecule type for names and other identifiers
        return { type: 'molecule', database: 'pubchem' };
    }

    /**
     * Universal method to get molecular data based on identifier type
     */
    async getMolecularData(identifier, type = 'auto') {
        try {
            let detectedType = { type: 'molecule', database: 'pubchem' };

            if (type === 'auto') {
                detectedType = this.detectMoleculeType(identifier);
            } else {
                // Use provided type
                if (type === 'cid' || type === 'name' || type === 'smiles') {
                    detectedType = { type: 'molecule', database: 'pubchem' };
                } else if (type === 'pdbid') {
                    detectedType = { type: 'protein', database: 'pdb' };
                } else if (type === 'codid') {
                    detectedType = { type: 'crystal', database: 'cod' };
                } else if (type === 'afid') {
                    detectedType = { type: 'protein', database: 'alphafolddb' };
                }
            }

            // Route to appropriate database based on detected type
            switch (detectedType.database) {
                case 'pdb':
                    return await this.getFromPDB(identifier);
                case 'cod':
                    return await this.getFromCOD(identifier);
                case 'alphafolddb':
                    return await this.getFromAlphaFoldDB(identifier);
                case 'pubchem':
                    // Determine which PubChem method to use
                    if (typeof identifier === 'number' || /^\d+$/.test(identifier)) {
                        return await this.getFromPubChem(identifier);
                    } else {
                        return await this.getFromPubChemByName(identifier);
                    }
                case 'chemspider':
                    return await this.getFromChemSpider(identifier);
                case 'chembl':
                    return await this.getFromChEMBL(identifier);
                case 'drugbank':
                    return await this.getFromDrugBank(identifier);
                default:
                    // Default to PubChem for small molecules
                    if (typeof identifier === 'number' || /^\d+$/.test(identifier)) {
                        return await this.getFromPubChem(identifier);
                    } else {
                        return await this.getFromPubChemByName(identifier);
                    }
            }
        } catch (error) {
            console.error('Error getting molecular data:', error);
            return null;
        }
    }

    /**
     * Get from ChemSpider
     */
    async getFromChemSpider(identifier) {
        try {
            let url;
            if (typeof identifier === 'number' || /^\d+$/.test(identifier)) {
                url = `${this.databases.chemspider.baseUrl}${this.databases.chemspider.endpoints.id.replace('{id}', identifier)}`;
            } else {
                // Search by SMILES or name
                url = `${this.databases.chemspider.baseUrl}${this.databases.chemspider.endpoints.smiles.replace('{smiles}', encodeURIComponent(identifier))}`;
            }

            const response = await fetch(url);
            if (!response.ok) {
                throw new Error(`ChemSpider request failed: ${response.status}`);
            }

            return {
                data: await response.text(),
                source: 'chemspider',
                format: 'sdf',
                identifier: identifier
            };
        } catch (error) {
            console.error(`Error getting data from ChemSpider for ${identifier}:`, error);
            return null;
        }
    }

    /**
     * Get from ChEMBL
     */
    async getFromChEMBL(identifier) {
        try {
            const url = `${this.databases.chembl.baseUrl}${this.databases.chembl.endpoints.molecule.replace('{id}', identifier)}`;
            const response = await fetch(url);

            if (!response.ok) {
                throw new Error(`ChEMBL request failed: ${response.status}`);
            }

            const data = await response.json();

            // Extract molecular data if available
            if (data.molecule && data.molecule.molfile) {
                return {
                    data: data.molecule.molfile,
                    source: 'chembl',
                    format: 'mol',
                    identifier: identifier
                };
            }

            throw new Error(`No molecular data found for ${identifier} in ChEMBL`);
        } catch (error) {
            console.error(`Error getting data from ChEMBL for ${identifier}:`, error);
            return null;
        }
    }

    /**
     * Get from DrugBank
     */
    async getFromDrugBank(identifier) {
        try {
            const url = `${this.databases.drugbank.baseUrl}${this.databases.drugbank.endpoints.drug.replace('{id}', identifier)}`;
            const response = await fetch(url);

            if (!response.ok) {
                throw new Error(`DrugBank request failed: ${response.status}`);
            }

            return {
                data: await response.text(),
                source: 'drugbank',
                format: 'sdf',
                identifier: identifier
            };
        } catch (error) {
            console.error(`Error getting data from DrugBank for ${identifier}:`, error);
            return null;
        }
    }

    /**
     * Get from AlphaFold DB
     */
    async getFromAlphaFoldDB(identifier) {
        try {
            // Try with version 4 (most common)
            const url = `${this.databases.alphafolddb.baseUrl}${this.databases.alphafolddb.endpoints.afid.replace('{id}', identifier).replace('{version}', '4')}`;
            const response = await fetch(url);

            if (!response.ok) {
                // Try with version 3
                const url_v3 = `${this.databases.alphafolddb.baseUrl}${this.databases.alphafolddb.endpoints.afid.replace('{id}', identifier).replace('{version}', '3')}`;
                const response_v3 = await fetch(url_v3);

                if (!response_v3.ok) {
                    throw new Error(`AlphaFold DB request failed: ${response.status}`);
                }

                return {
                    data: await response_v3.text(),
                    source: 'alphafolddb',
                    format: 'cif',
                    identifier: identifier
                };
            }

            return {
                data: await response.text(),
                source: 'alphafolddb',
                format: 'cif',
                identifier: identifier
            };
        } catch (error) {
            console.error(`Error getting data from AlphaFold DB for ${identifier}:`, error);
            return null;
        }
    }

    /**
     * Parse SDF format to molecular structure
     */
    parseSdf(sdfData) {
        const lines = sdfData.split('\n');
        const atoms = [];
        const bonds = [];
        
        let readingAtoms = false;
        let readingBonds = false;
        let atomCount = 0;
        let bondCount = 0;
        let currentLine = 0;
        
        // Skip header lines to get to counts line
        while (currentLine < lines.length) {
            if (/^\s*\d+\s+\d+/.test(lines[currentLine])) {
                // Found the counts line
                const parts = lines[currentLine].trim().split(/\s+/).filter(p => p);
                if (parts.length >= 2) {
                    atomCount = parseInt(parts[0]);
                    bondCount = parseInt(parts[1]);
                    currentLine++;
                    break;
                }
            }
            currentLine++;
        }
        
        // Read atoms
        for (let i = 0; i < atomCount && currentLine < lines.length; i++) {
            const line = lines[currentLine].trim();
            if (line && /^\s*-?\d+(\.\d+)?\s+-?\d+(\.\d+)?\s+-?\d+(\.\d+)?/.test(line)) {
                const parts = line.split(/\s+/).filter(p => p);
                if (parts.length >= 4) {
                    const [x, y, z, element] = [parseFloat(parts[0]), parseFloat(parts[1]), parseFloat(parts[2]), parts[3]];
                    if (!isNaN(x) && !isNaN(y) && !isNaN(z)) {
                        atoms.push({ x, y, z, element: element || 'C' });
                    }
                }
            }
            currentLine++;
        }
        
        // Read bonds
        for (let i = 0; i < bondCount && currentLine < lines.length; i++) {
            const line = lines[currentLine].trim();
            if (line && /^\s*\d+\s+\d+\s+\d/.test(line)) {
                const parts = line.split(/\s+/).filter(p => p);
                if (parts.length >= 3) {
                    const [from, to, order] = [parseInt(parts[0])-1, parseInt(parts[1])-1, parseInt(parts[2])]; // 0-indexed
                    if (!isNaN(from) && !isNaN(to) && !isNaN(order)) {
                        bonds.push({ from, to, order });
                    }
                }
            }
            currentLine++;
        }
        
        return { atoms, bonds, format: 'sdf' };
    }

    /**
     * Parse PDB format to molecular structure with protein-specific features
     */
    parsePdb(pdbData) {
        const lines = pdbData.split('\n');
        const atoms = [];
        const helix = [];
        const sheet = [];
        const residues = [];
        const bonds = [];

        let currentResidue = null;
        let serialToIndex = {}; // Map PDB serial numbers to array indices

        for (const line of lines) {
            if (line.startsWith('ATOM') || line.startsWith('HETATM')) {
                // Parse atom record with more detail
                const serial = parseInt(line.substring(6, 11).trim()); // PDB serial number
                const name = line.substring(12, 16).trim(); // Atom name (e.g., N, CA, C, O)
                const resName = line.substring(17, 20).trim(); // Residue name (e.g., ALA, GLY)
                const chainId = line.substring(21, 22).trim(); // Chain ID
                const resSeq = parseInt(line.substring(22, 26).trim()); // Residue sequence number
                const x = parseFloat(line.substring(30, 38).trim());
                const y = parseFloat(line.substring(38, 46).trim());
                const z = parseFloat(line.substring(46, 54).trim());
                const element = line.substring(76, 78).trim() || name.charAt(0) || 'C';

                if (!isNaN(x) && !isNaN(y) && !isNaN(z)) {
                    const atom = {
                        x, y, z,
                        element: element.toUpperCase(),
                        name: name,         // Atom name (N, CA, C, O, etc.)
                        resName: resName,   // Residue name (ALA, GLY, etc.)
                        chainId: chainId || 'A', // Chain ID
                        resSeq: resSeq,     // Residue sequence
                        serial: serial      // PDB serial number
                    };

                    // Track mapping from serial to array index
                    serialToIndex[serial] = atoms.length;
                    atoms.push(atom);

                    // Group atoms by residue
                    if (currentResidue === null ||
                        currentResidue.chainId !== chainId ||
                        currentResidue.resSeq !== resSeq) {
                        currentResidue = {
                            num: resSeq,
                            name: resName,
                            chain: chainId,
                            atoms: [],
                            startIdx: atoms.length - 1
                        };
                        residues.push(currentResidue);
                    }
                    currentResidue.atoms.push(atoms.length - 1); // store index
                }
            } else if (line.startsWith('HELIX')) {
                // Parse helix definition
                const startChain = line.substring(19, 20).trim();
                const startNum = parseInt(line.substring(21, 25).trim());
                const endChain = line.substring(31, 32).trim();
                const endNum = parseInt(line.substring(33, 37).trim());

                if (!isNaN(startNum) && !isNaN(endNum)) {
                    helix.push({
                        start: { chain: startChain, num: startNum },
                        end: { chain: endChain, num: endNum },
                        type: parseInt(line.substring(38, 40).trim())
                    });
                }
            } else if (line.startsWith('SHEET')) {
                // Parse sheet definition
                const startChain = line.substring(21, 22).trim();
                const startNum = parseInt(line.substring(22, 26).trim());
                const endChain = line.substring(32, 33).trim();
                const endNum = parseInt(line.substring(33, 37).trim());

                if (!isNaN(startNum) && !isNaN(endNum)) {
                    sheet.push({
                        start: { chain: startChain, num: startNum },
                        end: { chain: endChain, num: endNum }
                    });
                }
            } else if (line.startsWith('CONECT')) {
                // Parse bond record
                const parts = line.substring(6).trim().split(/\s+/).filter(p => p);
                if (parts.length >= 2) {
                    const fromSerial = parseInt(parts[0]);
                    if (serialToIndex.hasOwnProperty(fromSerial)) {
                        const fromIndex = serialToIndex[fromSerial];
                        for (let i = 1; i < parts.length; i++) {
                            const toSerial = parseInt(parts[i]);
                            if (serialToIndex.hasOwnProperty(toSerial)) {
                                const toIndex = serialToIndex[toSerial];
                                bonds.push({ from: fromIndex, to: toIndex, order: 1 });
                            }
                        }
                    }
                }
            }
        }

        return {
            atoms,
            bonds,
            helix,
            sheet,
            residues,
            format: 'pdb',
            type: 'protein'
        };
    }

    /**
     * Parse CIF format to molecular structure with crystal-specific features
     */
    parseCif(cifData) {
        const lines = cifData.split('\n');
        const atoms = [];
        const cellParams = {};
        const symmetryOps = [];

        let inAtomsBlock = false;
        let atomLabels = [];
        let columnMap = {};

        for (let i = 0; i < lines.length; i++) {
            const line = lines[i].trim();

            if (line.startsWith('_cell_')) {
                // Parse unit cell parameters
                const parts = line.split(/\s+/);
                if (parts.length >= 2) {
                    const param = line.replace('_cell_', '').replace('.', '');
                    const value = parseFloat(parts[parts.length - 1]); // Take last part as value
                    if (!isNaN(value)) {
                        cellParams[param] = value;
                    }
                }
            }

            if (line.startsWith('_symmetry_') && line.includes('operation')) {
                // Parse symmetry operations
                const value = line.split(/\s+/).pop();
                if (value && value.includes('x') && !symmetryOps.includes(value)) {
                    symmetryOps.push(value);
                }
            }

            if (line.includes('_atom_site') && !inAtomsBlock) {
                atomLabels.push(line.replace('_atom_site_', '').replace('.', ''));
            }

            if (line === 'loop_' && lines[i+1]?.includes('_atom_site')) {
                inAtomsBlock = true;
                // Reset and collect atom site labels
                atomLabels = [];
                let j = i + 1;
                while (j < lines.length && lines[j].trim().startsWith('_atom_site')) {
                    atomLabels.push(lines[j].replace('_atom_site_', '').replace('.', ''));
                    j++;
                }
                // Map column positions
                columnMap = {};
                atomLabels.forEach((label, idx) => {
                    columnMap[label] = idx;
                });
                i = j - 1; // Skip processed lines
                continue;
            }

            if (inAtomsBlock && line && !line.startsWith('_') && !line.startsWith('#') && !line.startsWith('loop_')) {
                const parts = line.split(/\s+/).filter(p => p);

                if (parts.length >= Math.max(...Object.values(columnMap), 0) + 1) {
                    let element = 'C'; // default
                    let x = 0, y = 0, z = 0;
                    let label = '';

                    if ('label' in columnMap && columnMap.label < parts.length) {
                        label = parts[columnMap.label];
                        element = label.replace(/\d/g, '');
                    } else if ('type_symbol' in columnMap && columnMap.type_symbol < parts.length) {
                        element = parts[columnMap.type_symbol];
                    }

                    if ('fract_x' in columnMap && columnMap.fract_x < parts.length) {
                        x = parseFloat(parts[columnMap.fract_x]);
                    } else if ('cartn_x' in columnMap && columnMap.cartn_x < parts.length) {
                        x = parseFloat(parts[columnMap.cartn_x]);
                    }

                    if ('fract_y' in columnMap && columnMap.fract_y < parts.length) {
                        y = parseFloat(parts[columnMap.fract_y]);
                    } else if ('cartn_y' in columnMap && columnMap.cartn_y < parts.length) {
                        y = parseFloat(parts[columnMap.cartn_y]);
                    }

                    if ('fract_z' in columnMap && columnMap.fract_z < parts.length) {
                        z = parseFloat(parts[columnMap.fract_z]);
                    } else if ('cartn_z' in columnMap && columnMap.cartn_z < parts.length) {
                        z = parseFloat(parts[columnMap.cartn_z]);
                    }

                    if (!isNaN(x) && !isNaN(y) && !isNaN(z)) {
                        atoms.push({
                            x, y, z,
                            element: element.toUpperCase(),
                            label: label,
                            type: element.toUpperCase()
                        });
                    }
                }
            }
        }

        return {
            atoms,
            cellParams,
            symmetryOps,
            format: 'cif',
            type: 'crystal'
        };
    }

    /**
     * Parse CIF format to molecular structure
     */
    parseCif(cifData) {
        const lines = cifData.split('\n');
        const atoms = [];
        const bonds = [];
        
        let inAtomBlock = false;
        let atomLabels = [];
        let columnMap = {};
        
        for (let i = 0; i < lines.length; i++) {
            const line = lines[i].trim();
            
            if (line.startsWith('_atom_site_')) {
                const label = line.replace('_atom_site_', '').replace('.', '');
                atomLabels.push(label);
            }
            
            if (line === 'loop_') {
                // Check if next lines are atom site labels
                if (i + 1 < lines.length && lines[i + 1].includes('_atom_site')) {
                    inAtomBlock = true;
                    // Map the column positions to data types
                    columnMap = {};
                    for (let j = 0; j < atomLabels.length; j++) {
                        columnMap[atomLabels[j]] = j;
                    }
                    i++; // Skip the 'loop_' line
                    continue;
                }
            }
            
            if (inAtomBlock && line && !line.startsWith('_') && !line.startsWith('#') && !line.startsWith('loop_')) {
                const parts = line.split(/\s+/).filter(p => p);
                
                if (parts.length >= Object.keys(columnMap).length) {
                    let element = 'C'; // default
                    let x = 0, y = 0, z = 0;
                    
                    if ('label' in columnMap && columnMap.label < parts.length) {
                        element = parts[columnMap.label].replace(/\d/g, '');
                    } else if ('type_symbol' in columnMap && columnMap.type_symbol < parts.length) {
                        element = parts[columnMap.type_symbol];
                    }
                    
                    if ('fract_x' in columnMap && columnMap.fract_x < parts.length) {
                        x = parseFloat(parts[columnMap.fract_x]);
                    } else if ('cartn_x' in columnMap && columnMap.cartn_x < parts.length) {
                        x = parseFloat(parts[columnMap.cartn_x]);
                    }
                    
                    if ('fract_y' in columnMap && columnMap.fract_y < parts.length) {
                        y = parseFloat(parts[columnMap.fract_y]);
                    } else if ('cartn_y' in columnMap && columnMap.cartn_y < parts.length) {
                        y = parseFloat(parts[columnMap.cartn_y]);
                    }
                    
                    if ('fract_z' in columnMap && columnMap.fract_z < parts.length) {
                        z = parseFloat(parts[columnMap.fract_z]);
                    } else if ('cartn_z' in columnMap && columnMap.cartn_z < parts.length) {
                        z = parseFloat(parts[columnMap.cartn_z]);
                    }
                    
                    if (!isNaN(x) && !isNaN(y) && !isNaN(z)) {
                        atoms.push({ x, y, z, element: element.toUpperCase() });
                    }
                }
            }
            
            // Look for chemical connectivity in CIF
            if (line.startsWith('_geom_bond_atom_site_id_1')) {
                // This would require more complex processing for bond information
            }
        }
        
        return { atoms, bonds, format: 'cif' };
    }

    /**
     * Convert molecular data to standardized format
     */
    convertToStandardFormat(molecularData) {
        if (!molecularData || !molecularData.data) {
            return null;
        }

        if (molecularData.format === 'sdf') {
            return this.parseSdf(molecularData.data);
        } else if (molecularData.format === 'pdb') {
            return this.parsePdb(molecularData.data);
        } else if (molecularData.format === 'cif') {
            return this.parseCif(molecularData.data);
        }

        return null;
    }

    /**
     * Get protein-specific data for visualization
     */
    getProteinData(molecularData) {
        if (molecularData.format === 'pdb') {
            const pdbData = this.parsePdb(molecularData.data);
            return {
                ...pdbData,
                caAtoms: pdbData.atoms.filter(a => a.name === 'CA'), // C-alpha atoms for backbone
                backboneAtoms: pdbData.atoms.filter(a => ['N', 'CA', 'C', 'O'].includes(a.name)),
                secondaryStructure: {
                    helices: pdbData.helix,
                    sheets: pdbData.sheet
                }
            };
        }
        return molecularData;
    }
}

// Export for use in other modules
if (typeof window !== 'undefined') {
    window.DirectChemicalDatabaseAccess = DirectChemicalDatabaseAccess;
}