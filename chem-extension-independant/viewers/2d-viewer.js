// 2D Molecule Viewer - uses SmilesDrawer library
// Parse URL parameters
const params = new URLSearchParams(window.location.search);
const smiles = params.get('smiles');
const name = params.get('name') || 'Molecule';
const theme = params.get('theme') || 'light';
const maxWidth = parseInt(params.get('maxWidth')) || parseInt(params.get('width')) || 600;
const maxHeight = parseInt(params.get('maxHeight')) || parseInt(params.get('height')) || 450;
const useSvg = params.get('svg') !== 'false'; // Default to SVG

// SmilesDrawer ACTUAL supported options (from official docs)
// terminalCarbons: Show CH3 groups explicitly (this is what "Show Methyl" maps to)
// explicitHydrogens: Show H atoms explicitly
// compactDrawing: Draw concatenated terminals and pseudo elements (default true)
// atomVisualization: 'default', 'balls', or 'none'

// NOTE: SmilesDrawer does NOT support:
// - Showing ALL carbons (only terminal CH3 via terminalCarbons)
// - Aromatic circles in rings (uses alternating double bonds or dashed lines)
// - Atom numbering

const showMethyl = params.get('showMethyl') === 'true';          // Maps to terminalCarbons
const addHydrogens = params.get('addHydrogens') === 'true';       // Maps to explicitHydrogens (branched H atoms)
const showImplicitH = params.get('showImplicitHydrogens') !== 'false'; // Default true - show H counts in labels (CH3, OH, NH2)
const showCarbons = params.get('showCarbons') === 'true';        // Display C labels
const compactDrawing = params.get('compactDrawing') === 'true';  // Compact mode (linear text strings) - default false
const atomVisualization = params.get('atomVisualization') || 'default'; // 'default', 'balls', 'none'

// Apply theme
document.body.className = theme === 'dark' ? 'dark-theme' : 'light-theme';

const loadingEl = document.getElementById('loading');
const errorEl = document.getElementById('error');
const errorMsgEl = document.getElementById('errorMessage');
const svgEl = document.getElementById('molecule-svg');
const canvasEl = document.getElementById('molecule-canvas');

function showError(message) {
    loadingEl.style.display = 'none';
    errorEl.style.display = 'block';
    errorMsgEl.textContent = message;
}

function hideLoading() {
    loadingEl.style.display = 'none';
}

// Render with SmilesDrawer
function renderWithSmilesDrawer() {
    console.log('🎨 Using SmilesDrawer renderer');

    const options = {
        width: maxWidth,
        height: maxHeight,
        bondThickness: 2.0,
        bondLength: 25,
        shortBondLength: 0.85,
        bondSpacing: 6,
        atomVisualization: atomVisualization,
        isomeric: true,
        debug: false,
        terminalCarbons: showMethyl || showCarbons,
        explicitHydrogens: addHydrogens,              // Legacy option (may not be used)
        showHydrogens: addHydrogens,                  // Show H atoms as separate nodes with bonds (branching)
        showImplicitHydrogens: showImplicitH,         // Show H counts in labels (H₂O, CH₃, OH, NH₂)
        compactDrawing: compactDrawing,           // Compact mode (linear text strings)
        overlapSensitivity: 0.42,
        overlapResolutionIterations: 1,
        fontSizeLarge: 10,
        fontSizeSmall: 8,
        padding: 20.0,
        experimental: false,
        themes: {
            dark: {
                C: '#ffffff',
                O: '#ff6b6b',
                N: '#4dabf7',
                F: '#51cf66',
                CL: '#20c997',
                BR: '#fd7e14',
                I: '#be4bdb',
                P: '#fd7e14',
                S: '#fcc419',
                B: '#f59f00',
                SI: '#f59f00',
                H: '#aaaaaa',
                BACKGROUND: '#1a1a2e'
            },
            light: {
                C: '#222222',
                O: '#e74c3c',
                N: '#3498db',
                F: '#27ae60',
                CL: '#16a085',
                BR: '#d35400',
                I: '#8e44ad',
                P: '#d35400',
                S: '#f1c40f',
                B: '#e67e22',
                SI: '#e67e22',
                H: '#666666',
                BACKGROUND: '#ffffff'
            }
        }
    };

    if (useSvg) {
        svgEl.style.display = 'block';
        canvasEl.style.display = 'none';

        const svgDrawer = new SmilesDrawer.SvgDrawer(options);

        SmilesDrawer.parse(smiles, function (tree) {
            svgDrawer.draw(tree, 'molecule-svg', theme);
            hideLoading();
            console.log('✅ SVG molecule rendered successfully');

            // Let SmilesDrawer handle sizing naturally, just report what it created
            setTimeout(() => {
                const svgElement = document.getElementById('molecule-svg');
                if (svgElement) {
                    try {
                        // Get the actual width/height that SmilesDrawer set
                        const svgWidth = parseInt(svgElement.getAttribute('width')) || maxWidth;
                        const svgHeight = parseInt(svgElement.getAttribute('height')) || maxHeight;

                        console.log('📏 SVG dimensions from SmilesDrawer:', `${svgWidth}x${svgHeight}`);

                        // Update body to match exactly
                        document.body.style.width = `${svgWidth}px`;
                        document.body.style.height = `${svgHeight}px`;

                        // Send exact size to parent
                        if (window.parent !== window) {
                            window.parent.postMessage({
                                type: '2d-viewer-size',
                                width: svgWidth,
                                height: svgHeight,
                                name: name
                            }, '*');
                        }
                    } catch (e) {
                        console.warn('Could not get SVG dimensions:', e);
                        if (window.parent !== window) {
                            window.parent.postMessage({
                                type: '2d-viewer-size',
                                width: maxWidth,
                                height: maxHeight,
                                name: name
                            }, '*');
                        }
                    }
                }
            }, 150);

            if (window.parent !== window) {
                window.parent.postMessage({ type: '2d-viewer-ready', name: name }, '*');
            }
        }, function (err) {
            console.error('❌ SMILES parse error:', err);
            showError('Failed to parse SMILES: ' + (err.message || err));
        });
    } else {
        svgEl.style.display = 'none';
        canvasEl.style.display = 'block';
        canvasEl.width = maxWidth;
        canvasEl.height = maxHeight;

        const canvasDrawer = new SmilesDrawer.Drawer(options);

        SmilesDrawer.parse(smiles, function (tree) {
            canvasDrawer.draw(tree, 'molecule-canvas', theme);
            hideLoading();
            console.log('✅ Canvas molecule rendered successfully');

            // For canvas, we'll use the full dimensions
            if (window.parent !== window) {
                window.parent.postMessage({
                    type: '2d-viewer-size',
                    width: maxWidth,
                    height: maxHeight,
                    name: name
                }, '*');
                window.parent.postMessage({ type: '2d-viewer-ready', name: name }, '*');
            }
        }, function (err) {
            console.error('❌ SMILES parse error:', err);
            showError('Failed to parse SMILES: ' + (err.message || err));
        });
    }
}

function renderMolecule() {
    if (!smiles) {
        showError('No SMILES provided');
        return;
    }

    console.log('🧪 2D Viewer: Rendering SMILES:', smiles);
    console.log('📐 Max Dimensions:', maxWidth, 'x', maxHeight);
    console.log('🎨 Theme:', theme);
    console.log('⚙️ Options:', {
        terminalCarbons: showMethyl || showCarbons,
        explicitHydrogens: addHydrogens,
        showHydrogens: showImplicitH,
        compactDrawing: !showCarbons,
        atomVisualization
    });

    try {
        renderWithSmilesDrawer();
    } catch (error) {
        console.error('❌ Render error:', error);
        showError('Rendering failed: ' + error.message);
    }
}

// Start rendering
renderMolecule();
