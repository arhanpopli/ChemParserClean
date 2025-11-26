/**
 * Advanced Molecular Viewer using 3Dmol.js
 * This provides all the features that MolView has by using the proven 3Dmol.js library
 */

class AdvancedMolViewer {
    constructor() {
        this.viewer = null;
        this.container = null;
    }

    /**
     * Initialize 3Dmol.js viewer
     */
    async initializeViewer(containerId, options = {}) {
        try {
            // 3Dmol.js should already be loaded
            if (typeof $3Dmol === 'undefined') {
                throw new Error('3Dmol.js library is not loaded. Please include the library before initializing the viewer.');
            }

            this.container = document.getElementById(containerId);
            if (!this.container) {
                throw new Error(`Container ${containerId} not found`);
            }

            // Create the viewer
            this.viewer = $3Dmol.createViewer(containerId, {
                backgroundColor: options.backgroundColor || 0x1a1a2e,
                antialias: true,
                ...options.viewSettings
            });

            console.log('3Dmol.js viewer initialized');
            return true;
        } catch (error) {
            console.error('Error initializing 3Dmol.js viewer:', error);
            return false;
        }
    }

    /**
     * Load molecule from data string (SDF, PDB, CIF format)
     */
    async loadMolecule(dataString, format = 'sdf', options = {}) {
        if (!this.viewer) {
            throw new Error('Viewer not initialized. Call initializeViewer() first.');
        }

        try {
            // Add the model from data string
            const model = this.viewer.addModel(dataString, format);
            
            // Determine appropriate styling based on options and molecule type
            const style = this.determineStyle(options);
            
            // Apply the style to the model
            this.viewer.setStyle({}, style);
            
            // Set background color
            if (options.backgroundColor) {
                this.viewer.setBackgroundColor(options.backgroundColor);
            }
            
            // Center and zoom
            this.viewer.zoomTo();
            this.viewer.render();
            
            console.log(`Molecule loaded successfully with format: ${format}`);
            return true;
        } catch (error) {
            console.error('Error loading molecule:', error);
            throw error;
        }
    }

    /**
     * Load molecule from URL (handles PubChem, PDB, etc. URLs)
     */
    async loadMoleculeFromUrl(url, options = {}) {
        try {
            const response = await fetch(url);
            if (!response.ok) {
                throw new Error(`Failed to fetch molecule data: ${response.status} ${response.statusText}`);
            }
            
            const data = await response.text();
            const format = this.determineFormatFromUrl(url);
            
            return await this.loadMolecule(data, format, options);
        } catch (error) {
            console.error('Error loading molecule from URL:', error);
            throw error;
        }
    }

    /**
     * Determine format from URL
     */
    determineFormatFromUrl(url) {
        if (url.includes('.pdb') || url.includes('PDB')) return 'pdb';
        if (url.includes('.sdf') || url.includes('SDF')) return 'sdf';
        if (url.includes('.cif') || url.includes('CIF')) return 'cif';
        if (url.includes('pubchem')) return 'sdf';
        return 'sdf'; // default
    }

    /**
     * Determine appropriate style based on options and molecule type
     */
    determineStyle(options = {}) {
        const representation = options.representation || 'cartoon';
        const colorScheme = options.colorScheme || 'whiteCarbon';
        const atomRadius = options.atomRadius || 0.15;

        switch (representation) {
            case 'cartoon':
                // Cartoon representation (ribbons for proteins)
                return {
                    cartoon: { 
                        color: colorScheme === 'spectrum' ? 'spectrum' : 
                                 colorScheme === 'chain' ? 'chain' : 
                                 'white' 
                    }
                };
                
            case 'ribbon':
                // Ribbon representation (simplified cartoon)
                return {
                    cartoon: { 
                        style: 'trace',
                        color: colorScheme === 'spectrum' ? 'spectrum' : 
                               colorScheme === 'chain' ? 'chain' : 
                               'white',
                        radius: atomRadius * 2
                    }
                };
                
            case 'stick':
                // Stick model
                return {
                    stick: {
                        radius: atomRadius,
                        color: colorScheme === 'whiteCarbon' ? 'whiteCarbon' : 
                               colorScheme === 'byElement' ? 'greenCarbon' : 
                               'white'
                    }
                };
                
            case 'sphere':
                // Sphere model (van der Waals spheres)
                return {
                    sphere: {
                        radius: atomRadius * 1.5,
                        color: colorScheme === 'whiteCarbon' ? 'whiteCarbon' : 
                               colorScheme === 'byElement' ? 'greenCarbon' : 
                               'white'
                    }
                };
                
            case 'line':
                // Line model
                return {
                    line: {
                        color: colorScheme === 'whiteCarbon' ? 'whiteCarbon' : 
                               colorScheme === 'byElement' ? 'greenCarbon' : 
                               'white'
                    }
                };
                
            case 'cross':
                // Cross representation
                return {
                    cross: {
                        color: colorScheme === 'whiteCarbon' ? 'whiteCarbon' : 
                               colorScheme === 'byElement' ? 'greenCarbon' : 
                               'white'
                    }
                };
                
            case 'surface':
                // Surface representation
                return {
                    surface: {
                        opacity: 0.8,
                        color: colorScheme === 'spectrum' ? 'spectrum' : 
                               colorScheme === 'chain' ? 'chain' : 
                               'lightblue'
                    }
                };
                
            default:
                // Default to cartoon for proteins, stick/sphere for small molecules
                // But let's also check the molecule to see if it has protein features
                if (options.isProtein) {
                    return {
                        cartoon: { 
                            color: colorScheme === 'spectrum' ? 'spectrum' : 
                                     colorScheme === 'chain' ? 'chain' : 
                                     'white' 
                        }
                    };
                } else {
                    return {
                        stick: {
                            radius: atomRadius,
                            color: 'whiteCarbon'
                        }
                    };
                }
        }
    }

    /**
     * Apply specific protein styling (for helices, sheets, etc.)
     */
    applyProteinStyling(options = {}) {
        if (!this.viewer) return;

        const colorScheme = options.colorScheme || 'spectrum';

        // Apply cartoon (ribbon) styling for proteins
        this.viewer.setStyle({}, {
            cartoon: {
                color: colorScheme,
                radius: 1.0
            }
        });

        // Add additional styling for specific elements if needed
        this.viewer.setStyle({elem: 'C'}, {cartoon: {color: 'white'}});
        this.viewer.setStyle({elem: 'O'}, {cartoon: {color: 'red'}});
        this.viewer.setStyle({elem: 'N'}, {cartoon: {color: 'blue'}});

        this.viewer.render();
    }

    /**
     * Apply crystal styling
     */
    applyCrystalStyling(options = {}) {
        if (!this.viewer) return;

        // For crystals, use lines and sticks with unit cell highlighting
        this.viewer.setStyle({}, {
            stick: {
                radius: 0.1,
                color: options.colorScheme || 'gray'
            }
        });

        this.viewer.render();
    }

    /**
     * Apply surface representation
     */
    applySurface(options = {}) {
        if (!this.viewer) return;

        const surfaceType = options.surfaceType || 'VDW'; // Van der Waals surface
        const opacity = options.opacity || 0.8;
        const color = options.color || 'lightblue';

        this.viewer.addSurface($3Dmol.SurfaceType.VDW, {
            opacity: opacity,
            color: color
        }, {}, {});

        this.viewer.render();
    }

    /**
     * Rotate the molecule
     */
    rotate(angle, axis = 'y') {
        if (!this.viewer) return;
        
        this.viewer.rotate(angle, axis);
        this.viewer.render();
    }

    /**
     * Set zoom level
     */
    zoom(level) {
        if (!this.viewer) return;
        
        this.viewer.zoom(level);
        this.viewer.render();
    }

    /**
     * Toggle auto-rotation
     */
    setSpin(enabled) {
        if (!this.viewer) return;
        
        this.viewer.spin(enabled);
        this.viewer.render();
    }

    /**
     * Render the current view
     */
    render() {
        if (!this.viewer) return;
        
        this.viewer.render();
    }
}

// Export for use in other modules
if (typeof window !== 'undefined') {
    window.AdvancedMolViewer = AdvancedMolViewer;
}