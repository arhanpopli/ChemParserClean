/**
 * MolView 3D Data Extractor and Local Renderer
 * This module extracts 3D molecular data from MolView and renders it with local 3D engines
 */

class MolView3DExtractor {
    constructor() {
        this.supportedEngines = ['GLmol', 'JSmol', 'ChemDoodle'];
        this.currentEngine = null;
    }

    /**
     * Extract 3D molecular data from MolView
     * @param {string} cid - PubChem compound ID
     * @param {string} smiles - SMILES string (alternative to CID)
     * @returns {Promise<Object>} - Molecular data in appropriate format
     */
    async extractMolViewData(cid, smiles = null) {
        try {
            let molViewParam;
            if (smiles) {
                molViewParam = `smiles=${encodeURIComponent(smiles)}`;
            } else if (cid) {
                molViewParam = `cid=${cid}`;
            } else {
                throw new Error('Either CID or SMILES must be provided');
            }

            // Try to get SDF data from PubChem first (reliable source)
            if (cid) {
                const sdfData = await this.fetchSdfFromPubChem(cid);
                if (sdfData) {
                    return { 
                        type: 'sdf', 
                        data: sdfData, 
                        source: 'pubchem',
                        cid: cid
                    };
                }
            }

            // If PubChem SDF not available, try to get it from MolView 
            // by constructing the appropriate URL and extracting data
            const molViewUrl = `https://molview.org/?${molViewParam}`;
            
            // Since we cannot directly access MolView's internal data due to CORS,
            // we'll need to work with what we have or try alternative sources
            if (smiles) {
                // If we have SMILES, we can convert it to 3D using our local tools
                return await this.convertSmilesTo3D(smiles);
            }

            throw new Error('Could not extract 3D data from MolView due to CORS restrictions');
        } catch (error) {
            console.error('Error extracting data from MolView:', error);
            throw error;
        }
    }

    /**
     * Fetch SDF data from PubChem (more reliable 3D source)
     */
    async fetchSdfFromPubChem(cid) {
        try {
            // Try 3D conformer first
            const url3d = `https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/${cid}/record/SDF?record_type=3d`;
            const response3d = await fetch(url3d);

            if (response3d.ok) {
                return await response3d.text();
            }

            // Fallback to computed 3D
            const urlComputed = `https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/${cid}/SDF?record_type=3d`;
            const responseComputed = await fetch(urlComputed);

            if (responseComputed.ok) {
                return await responseComputed.text();
            }

            // Final fallback to 2D (not ideal for 3D, but better than nothing)
            const url2d = `https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/${cid}/SDF`;
            const response2d = await fetch(url2d);

            if (response2d.ok) {
                return await response2d.text();
            }

            return null;
        } catch (error) {
            console.error('Error fetching SDF from PubChem:', error);
            return null;
        }
    }

    /**
     * Convert SMILES to 3D structure using a service or local tools
     */
    async convertSmilesTo3D(smiles) {
        try {
            // For now, we'll try to get 3D data using a service that can convert SMILES to 3D
            // In a real implementation, you might want to use RDKit.js or a local service
            
            // Let's try to find the CID from SMILES first, then get 3D from PubChem
            const cid = await this.getCIDFromSmiles(smiles);
            if (cid) {
                const sdfData = await this.fetchSdfFromPubChem(cid);
                if (sdfData) {
                    return {
                        type: 'sdf',
                        data: sdfData,
                        source: 'pubchem',
                        cid: cid
                    };
                }
            }

            // If we can't get from PubChem, use SMILES directly with 3D generation
            // This is a simplified approach - in practice, you'd use RDKit.js or a similar tool
            return {
                type: 'smiles',
                data: smiles,
                source: 'smiles'
            };
        } catch (error) {
            console.error('Error converting SMILES to 3D:', error);
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
     * Render molecular data using the selected 3D engine
     */
    async renderMolecule(data, containerId, options = {}) {
        try {
            const container = document.getElementById(containerId);
            if (!container) {
                throw new Error(`Container with id ${containerId} not found`);
            }

            // Clear the container
            container.innerHTML = '';

            // Set up a default 3D viewer
            container.style.width = options.width || '100%';
            container.style.height = options.height || '400px';
            container.style.position = 'relative';

            // Add loading indicator
            const loadingDiv = document.createElement('div');
            loadingDiv.id = 'molview-3d-loading';
            loadingDiv.style.cssText = `
                position: absolute;
                top: 50%;
                left: 50%;
                transform: translate(-50%, -50%);
                text-align: center;
                color: white;
                z-index: 10;
            `;
            loadingDiv.innerHTML = `
                <div style="font-size: 18px; margin-bottom: 10px;">Loading 3D structure...</div>
                <div style="width: 40px; height: 40px; border: 3px solid #2a2a4e; border-top-color: #667eea; border-radius: 50%; animation: spin 1s linear infinite; margin: 0 auto;"></div>
                <style>
                    @keyframes spin { to { transform: rotate(360deg); } }
                </style>
            `;
            container.appendChild(loadingDiv);

            // Add a new div for the 3D viewer
            const viewerDiv = document.createElement('div');
            viewerDiv.id = 'molview-3d-viewer';
            viewerDiv.style.cssText = `width: 100%; height: 100%;`;
            container.appendChild(viewerDiv);

            // Try different approaches based on data type
            if (data.type === 'sdf') {
                // Try to use JSmol first (it handles SDF well)
                await this.renderWithJSmol(data.data, viewerDiv, options);
            } else if (data.type === 'smiles') {
                // Use SMILES with JSmol or other engines
                await this.renderWithJSmolSmiles(data.data, viewerDiv, options);
            } else {
                // Fallback
                await this.renderWithChemDoodle(data.data, viewerDiv, options);
            }

            // Remove loading indicator
            if (loadingDiv.parentNode) {
                loadingDiv.parentNode.removeChild(loadingDiv);
            }

        } catch (error) {
            console.error('Error rendering molecule:', error);
            this.showError(containerId, `Error rendering molecule: ${error.message}`);
        }
    }

    /**
     * Render using JSmol
     */
    async renderWithJSmol(sdfData, container, options = {}) {
        try {
            const containerId = container.id || 'molview-3d-viewer';

            // First check if JSmol/JSME is available
            if (typeof Jmol !== 'undefined') {
                // JSmol is available - use it for rendering
                const jmolAppletId = 'jmolApplet';

                container.innerHTML = `
                    <div id="${jmolAppletId}_infotablediv" style="display: none"></div>
                    <div id="${jmolAppletId}_maindiv" style="width: 100%; height: 100%;"></div>
                `;

                // Jmol initialization script
                Jmol.getApplet(jmolAppletId, {
                    width: container.clientWidth || 400,
                    height: container.clientHeight || 400,
                    debug: false,
                    color: "#1a1a2e",
                    disableInitialConsole: true,
                    addSelectionScript: "select *",
                    startupScript: "set antialiasDisplay ON; set frank off; load INLINE \"\"\"" + sdfData + "\"\"\"; select all; calculate hbonds; set hydrogenBonds TRUE; set hbonds type set; set measurementUnits ANGSTROMS;"
                });

                // Add the applet to the container
                document.getElementById(`${jmolAppletId}_maindiv`).appendChild(
                    Jmol.getHtml(jmolAppletId, container.clientWidth || 400, container.clientHeight || 400)
                );

                // Execute the load script
                Jmol.script(jmolAppletId, `load INLINE """${sdfData}"""`);

            } else {
                // JSmol is not available, provide a fallback with 3Dmol.js if available
                if (typeof $3Dmol !== 'undefined') {
                    // Use 3Dmol.js as fallback
                    await this.renderWith3Dmol(sdfData, container, options);
                } else if (typeof ChemDoodleWeb !== 'undefined') {
                    // Use ChemDoodle as fallback
                    await this.renderWithChemDoodleSdf(sdfData, container, options);
                } else {
                    // Pure fallback implementation showing the SDF data
                    container.innerHTML = `
                        <div style="display: flex; align-items: center; justify-content: center; height: 100%; background: #1a1a2e; color: white; flex-direction: column;">
                            <div style="text-align: center; padding: 20px; max-width: 90%;">
                                <h3>3D Structure Preview</h3>
                                <p>Molecular structure data loaded successfully</p>
                                <div style="background: rgba(0,0,0,0.3); padding: 15px; margin: 15px 0; max-height: 200px; overflow-y: auto; text-align: left; font-family: monospace; font-size: 11px; border-radius: 8px;">
                                    <div style="color: #667eea; margin-bottom: 10px;">SDF Format Data Preview:</div>
                                    <pre style="margin: 0; white-space: pre-wrap; word-break: break-word;">${sdfData.substring(0, 500)}${sdfData.length > 500 ? '...' : ''}</pre>
                                </div>

                                <div style="margin-top: 15px;">
                                    <p>To view in 3D, install one of these 3D engines:</p>
                                    <div style="display: flex; justify-content: center; gap: 10px; flex-wrap: wrap; margin-top: 10px;">
                                        <span class="engine-tag">JSmol</span>
                                        <span class="engine-tag">3Dmol.js</span>
                                        <span class="engine-tag">ChemDoodle</span>
                                    </div>
                                </div>

                                <button onclick="load3DEngine('${encodeURIComponent(sdfData)}')" style="padding: 10px 20px; background: #667eea; color: white; border: none; border-radius: 4px; margin-top: 20px; cursor: pointer;">
                                    Try 3D Visualization
                                </button>
                            </div>
                        </div>
                        <style>
                            .engine-tag {
                                background: rgba(102, 126, 234, 0.2);
                                color: #a5b4fc;
                                padding: 4px 8px;
                                border-radius: 12px;
                                font-size: 12px;
                            }
                        </style>
                    `;

                    // Define the handler function in the global scope
                    window.load3DEngine = function(encodedData) {
                        const data = decodeURIComponent(encodedData);
                        alert('In a complete implementation with 3D engines installed, this SDF data would be rendered in 3D.\n\nFirst 500 characters of data:\n' + data.substring(0, 500) + (data.length > 500 ? '...' : ''));
                    };
                }
            }
        } catch (error) {
            console.error('Error with JSmol rendering:', error);
            container.innerHTML = `<div style="color: white; padding: 20px; text-align: center;">Error initializing 3D engine: ${error.message}<br>Make sure 3D libraries are loaded.</div>`;
        }
    }

    /**
     * Render with 3Dmol.js
     */
    async renderWith3Dmol(sdfData, container, options = {}) {
        try {
            // Set up the container for 3Dmol.js
            container.innerHTML = `<div id="3dmol-container" style="width:100%; height:100%; position:relative;"></div>`;
            const viewerDiv = document.getElementById('3dmol-container');

            // Initialize 3Dmol.js viewer
            let viewer = $3Dmol.createViewer(viewerDiv, {
                backgroundColor: '0x1a1a2e',
                antialias: true
            });

            // Add the model from SDF data
            let model = viewer.addModel(sdfData, "sdf");
            viewer.setStyle({},{stick:{},sphere:{radius:0.25}});
            viewer.zoomTo();
            viewer.render();

            console.log('3Dmol.js rendering successful');
        } catch (error) {
            console.error('Error with 3Dmol.js rendering:', error);
            container.innerHTML = `<div style="color: white; padding: 20px;">Error with 3Dmol.js: ${error.message}</div>`;
        }
    }

    /**
     * Render SDF data with ChemDoodle Web Components
     */
    async renderWithChemDoodleSdf(sdfData, container, options = {}) {
        try {
            container.innerHTML = `
                <div id="chemdoodle-viewer" style="width:100%; height:100%;"></div>
                <div id="chemdoodle-info" style="position:absolute; top:10px; left:10px; color:white; background:rgba(0,0,0,0.5); padding:5px; border-radius:3px; font-size:12px;"></div>
            `;

            // If ChemDoodle is available, use it to render the SDF
            if (typeof ChemDoodle !== 'undefined' && typeof ChemDoodle.ReadWrite.molfileToObj != 'undefined') {
                // Create a 3D viewer
                let viewer = new ChemDoodle.Web3D.ViewerCanvas3D('chemdoodle-viewer', container.clientWidth || 400, container.clientHeight || 400);

                // Set up some default styles
                viewer.specs.set3DFlags(false);
                viewer.specs.backgroundColor = 'black';
                viewer.specs.atoms_displayLabels_3D = true;

                // Load the molecule from SDF data
                try {
                    // Attempt to load the SDF data
                    let molecule = ChemDoodle.readMOL(sdfData);
                    if (molecule) {
                        viewer.loadMolecule(molecule);
                        viewer.animationManager.zoomAndRotate();
                    } else {
                        document.getElementById('chemdoodle-info').textContent = 'Could not parse SDF data';
                    }
                } catch (e) {
                    document.getElementById('chemdoodle-info').textContent = 'Error parsing SDF: ' + e.message;
                }
            } else {
                container.innerHTML = `<div style="color: white; padding: 20px;">ChemDoodle not available for SDF rendering</div>`;
            }
        } catch (error) {
            console.error('Error with ChemDoodle SDF rendering:', error);
            container.innerHTML = `<div style="color: white; padding: 20px;">Error with ChemDoodle SDF: ${error.message}</div>`;
        }
    }

    /**
     * Render SMILES with available engines
     */
    async renderWithJSmolSmiles(smiles, container, options = {}) {
        try {
            // Try different engines based on availability
            if (typeof $3Dmol !== 'undefined') {
                // Use 3Dmol.js directly with SMILES
                container.innerHTML = `<div id="3dmol-smiles-container" style="width:100%; height:100%; position:relative;"></div>`;
                const viewerDiv = document.getElementById('3dmol-smiles-container');

                let viewer = $3Dmol.createViewer(viewerDiv, {
                    backgroundColor: '0x1a1a2e',
                    antialias: true
                });

                // Load from SMILES (this needs a service to convert SMILES to 3D)
                // For now, we'll need to get 3D coordinates from a service or use a fallback
                viewer.addModel(smiles, "smiles");  // This may not work directly in all 3Dmol versions
                viewer.setStyle({stick:{},sphere:{radius:0.25}});
                viewer.zoomTo();
                viewer.render();
            } else if (typeof Jmol !== 'undefined') {
                // Fallback to JSmol approach for SMILES
                container.innerHTML = `
                    <div style="display: flex; align-items: center; justify-content: center; height: 100%; background: #1a1a2e; color: white;">
                        <div style="text-align: center;">
                            <h3>3D Structure from SMILES</h3>
                            <p>SMILES: <code>${smiles}</code></p>
                            <p>This would be converted to 3D coordinates and rendered with JSmol</p>
                            <button onclick="handleSmilesRender('${encodeURIComponent(smiles)}')" style="padding: 10px 20px; background: #667eea; color: white; border: none; border-radius: 4px; margin-top: 10px;">
                                Generate 3D Structure
                            </button>
                            <div style="margin-top: 15px; font-size: 12px; color: #aaa;">
                                In a complete implementation, this would use RDKit.js or a similar tool to convert SMILES to 3D coordinates
                            </div>
                        </div>
                    </div>
                `;

                // Define the handler function in the global scope
                window.handleSmilesRender = function(encodedSmiles) {
                    const smiles = decodeURIComponent(encodedSmiles);
                    alert('In a full implementation, this would convert the SMILES string to 3D coordinates and render with JSmol.\n\nSMILES:\n' + smiles);
                };
            } else {
                // Ultimate fallback showing the SMILES with option to get 3D
                container.innerHTML = `
                    <div style="display: flex; align-items: center; justify-content: center; height: 100%; background: #1a1a2e; color: white; flex-direction: column;">
                        <div style="text-align: center; padding: 20px; max-width: 80%;">
                            <h3>SMILES String</h3>
                            <p style="word-break: break-all; margin: 15px 0; color: #667eea; font-family: monospace;">${smiles}</p>
                            <p>To visualize this in 3D, it needs to be converted to coordinates</p>

                            <div style="margin: 20px 0;">
                                <p>Available 3D engines:</p>
                                <div style="display: flex; justify-content: center; gap: 10px; flex-wrap: wrap; margin-top: 10px;">
                                    <span class="engine-tag" title="Best for SDF files">JSmol</span>
                                    <span class="engine-tag" title="Good for various formats">3Dmol.js</span>
                                    <span class="engine-tag" title="ChemDoodle Web Components">ChemDoodle</span>
                                </div>
                            </div>

                            <p style="margin-top: 15px; font-size: 12px; color: #aaa;">
                                Install one of these libraries to enable 3D visualization of SMILES
                            </p>
                        </div>
                    </div>
                    <style>
                        .engine-tag {
                            background: rgba(102, 126, 234, 0.2);
                            color: #a5b4fc;
                            padding: 4px 8px;
                            border-radius: 12px;
                            font-size: 12px;
                        }
                    </style>
                `;
            }
        } catch (error) {
            console.error('Error with SMILES rendering:', error);
            container.innerHTML = `<div style="color: white; padding: 20px;">Error rendering SMILES: ${error.message}</div>`;
        }
    }

    /**
     * Render using ChemDoodle Web Components
     */
    async renderWithChemDoodle(molData, container, options = {}) {
        try {
            container.innerHTML = `
                <div style="display: flex; align-items: center; justify-content: center; height: 100%; background: #1a1a2e; color: white;">
                    <div style="text-align: center;">
                        <h3>3D Structure with ChemDoodle</h3>
                        <p>Structure would be rendered using ChemDoodle Web Components</p>
                        <div id="chemdoodle-container" style="width: 300px; height: 300px; margin: 0 auto;"></div>
                        <p style="margin-top: 15px;">Data format: ${typeof molData === 'string' ? molData.substring(0, 50) + '...' : 'Object'}</p>
                    </div>
                </div>
            `;

            // If ChemDoodle library is available, use it to render
            if (typeof ChemDoodleWeb !== 'undefined') {
                // This would initialize the ChemDoodle viewer in a real implementation
                console.log('ChemDoodle rendering would happen here');
            }

        } catch (error) {
            console.error('Error with ChemDoodle rendering:', error);
            container.innerHTML = `<div style="color: white; padding: 20px;">Error with ChemDoodle: ${error.message}</div>`;
        }
    }

    /**
     * Show error message in the container
     */
    showError(containerId, message) {
        const container = document.getElementById(containerId);
        if (container) {
            container.innerHTML = `
                <div style="display: flex; align-items: center; justify-content: center; height: 100%; background: #1a1a2e; color: white;">
                    <div style="text-align: center; padding: 20px;">
                        <h3 style="color: #ff6b6b;">Error</h3>
                        <p>${message}</p>
                        <button onclick="location.reload()" style="padding: 10px 20px; background: #667eea; color: white; border: none; border-radius: 4px; margin-top: 10px;">
                            Retry
                        </button>
                    </div>
                </div>
            `;
        }
    }

    /**
     * Get available 3D engines
     */
    getAvailableEngines() {
        const engines = [];
        
        // Check if each engine is available
        if (typeof Jmol !== 'undefined') {
            engines.push('Jmol');
        }
        if (typeof ChemDoodleWeb !== 'undefined') {
            engines.push('ChemDoodle');
        }
        if (typeof $3Dmol !== 'undefined') {
            engines.push('3Dmol');
        }
        if (typeof GLmol !== 'undefined') {
            engines.push('GLmol');
        }
        
        return engines;
    }
}

// Create a global instance
const molViewExtractor = new MolView3DExtractor();

// Export for use in other modules if needed
if (typeof module !== 'undefined' && module.exports) {
    module.exports = MolView3DExtractor;
}