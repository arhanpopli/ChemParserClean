// Get CID from URL parameters
const urlParams = new URLSearchParams(window.location.search);
const cid = urlParams.get('cid');
const name = urlParams.get('name');

if (!cid) {
    showError('No molecule ID provided');
} else {
    loadPubChemViewer(cid);
}

function showError(message) {
    document.getElementById('loading').style.display = 'none';
    document.getElementById('error').style.display = 'block';
    document.getElementById('errorMessage').textContent = message;
}

function loadPubChemViewer(cid) {
    const iframe = document.getElementById('viewer-frame');
    const pubchemUrl = `https://pubchem.ncbi.nlm.nih.gov/compound/${cid}#section=3D-Conformer`;

    iframe.src = pubchemUrl;

    // Hide loading after iframe loads
    iframe.onload = function() {
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

    iframe.onerror = function() {
        showError('Failed to load PubChem viewer');
    };

    console.log('📍 Loading PubChem 3D Viewer for CID:', cid);
}
