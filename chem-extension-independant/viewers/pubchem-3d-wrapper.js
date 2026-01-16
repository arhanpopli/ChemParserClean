// Get CID and SMILES from URL parameters
const urlParams = new URLSearchParams(window.location.search);
const cid = urlParams.get('cid');
const name = urlParams.get('name');
const smiles = urlParams.get('smiles'); // Fallback SMILES if SDF not available

if (!cid) {
    showError('No molecule ID provided');
} else {
    loadPubChemViewer(cid, smiles);
}

function showError(message) {
    document.getElementById('loading').style.display = 'none';
    document.getElementById('error').style.display = 'block';
    document.getElementById('errorMessage').textContent = message;
}

/**
 * Check if PubChem has 3D SDF data for this CID
 * Falls back to MolView with SMILES if not available
 */
async function checkSDFAvailability(cid) {
    try {
        const sdfUrl = `https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/${cid}/sdf?record_type=3d`;
        const response = await fetch(sdfUrl, { method: 'HEAD' }); // Just check if exists
        return response.ok; // Returns true if 200, false if 404
    } catch (error) {
        console.error('Error checking SDF availability:', error);
        return false;
    }
}

function loadPubChemViewer(cid, fallbackSmiles) {
    const iframe = document.getElementById('viewer-frame');

    // First check if 3D SDF data is available
    checkSDFAvailability(cid).then(hasSDFData => {
        if (!hasSDFData && fallbackSmiles) {
            // Fallback to MolView with SMILES
            console.log('⚠️ No 3D SDF data found for CID:', cid, '- Using MolView fallback with SMILES');
            document.getElementById('loading').querySelector('div:last-child').textContent = 'Loading MolView 3D Viewer...';

            // Use MolView with SMILES parameter
            const molviewUrl = `https://embed.molview.org/v1/?smiles=${encodeURIComponent(fallbackSmiles)}`;
            iframe.src = molviewUrl;

            iframe.onload = function () {
                document.getElementById('loading').style.display = 'none';
                console.log('✅ MolView loaded successfully with SMILES');
            };

            iframe.onerror = function () {
                showError('Failed to load MolView viewer');
            };
        } else {
            // Use normal PubChem 3D viewer
            const pubchemUrl = `https://pubchem.ncbi.nlm.nih.gov/compound/${cid}#section=3D-Conformer`;
            iframe.src = pubchemUrl;

            // Hide loading after iframe loads
            iframe.onload = function () {
                document.getElementById('loading').style.display = 'none';

                // Try to inject CSS to hide unwanted elements
                // Note: This will fail due to cross-origin restrictions, but we'll try anyway
                try {
                    const iframeDoc = iframe.contentDocument || iframe.contentWindow.document;

                    // Inject CSS to hide everything except the 3D viewer
                    const style = iframeDoc.createElement('style');
                    style.textContent = `
                        /* Hide header and navigation */
                        header, nav, .header, .navigation, #skip-to-main {
                            display: none !important;
                        }

                        /* Hide footer */
                        footer, .footer {
                            display: none !important;
                        }

                        /* Hide sidebar */
                        aside, .sidebar, #main-navigation {
                            display: none !important;
                        }

                        /* Hide breadcrumbs */
                        .breadcrumb, .breadcrumbs {
                            display: none !important;
                        }

                        /* Hide all sections except 3D Conformer */
                        .section:not(#section-3D-Conformer) {
                            display: none !important;
                        }

                        /* Show only 3D section */
                        #section-3D-Conformer {
                            display: block !important;
                            position: fixed !important;
                            top: 0 !important;
                            left: 0 !important;
                            width: 100vw !important;
                            height: 100vh !important;
                            margin: 0 !important;
                            padding: 0 !important;
                            z-index: 9999 !important;
                        }

                        /* Hide everything else */
                        body > *:not(#section-3D-Conformer) {
                            display: none !important;
                        }

                        /* Make the 3D viewer full screen */
                        #div3d, .div3d, canvas {
                            width: 100% !important;
                            height: 100% !important;
                        }
                    `;
                    iframeDoc.head.appendChild(style);

                    console.log('✅ Successfully injected CSS to hide PubChem UI elements');
                } catch (error) {
                    console.warn('⚠️ Could not inject CSS due to cross-origin restrictions:', error);
                    console.log('📍 PubChem page will show full UI (cross-origin restriction)');
                }
            };

            iframe.onerror = function () {
                showError('Failed to load PubChem viewer');
            };

            console.log('📍 Loading PubChem 3D Viewer for CID:', cid);
        }
    });
}

