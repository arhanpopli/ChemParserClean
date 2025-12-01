// QuintessenLabs V2 Logic

// --- Tab Management ---
function switchTab(tabId) {
    // Hide all tabs
    document.querySelectorAll('.tab-pane').forEach(el => {
        el.classList.remove('active');
    });

    // Deactivate all nav items
    document.querySelectorAll('.nav-item').forEach(el => {
        el.classList.remove('active');
    });

    // Show target tab
    const target = document.getElementById(tabId);
    if (target) {
        target.classList.add('active');
    }

    // Activate nav item
    const navItem = document.querySelector(`.nav-item[onclick="switchTab('${tabId}')"]`);
    if (navItem) {
        navItem.classList.add('active');
        // Update header title
        const title = navItem.querySelector('span').nextSibling.textContent.trim();
        document.getElementById('page-title').textContent = title;
    }
}

// --- Mol2ChemFig ---
async function generateMol2ChemFig() {
    // Logic to be implemented or connected to existing API
    // For now, we use the iframe interface which handles its own logic
    console.log("Mol2ChemFig generation triggered");
}

// --- MoleculeViewer ---
async function generateMoleculeViewer() {
    const smiles = document.getElementById('mv-smiles').value;
    const use3D = document.getElementById('mv-3d-toggle').checked;
    const width = document.getElementById('mv-width').value;
    const height = document.getElementById('mv-height').value;
    const resultDiv = document.getElementById('mv-result');

    resultDiv.innerHTML = '<p>Processing...</p>';

    try {
        let finalSmiles = smiles;
        if (use3D) {
            resultDiv.innerHTML = '<p>Converting to 3D SMILES...</p>';
            const resp = await fetch(`http://localhost:5001/api/opsin?name=${encodeURIComponent(smiles)}`);
            const data = await resp.json();
            if (data.smiles_3d) finalSmiles = data.smiles_3d;
        }

        const imgUrl = `http://localhost:5000/img/smiles?smiles=${encodeURIComponent(finalSmiles)}&width=${width}&height=${height}`;
        resultDiv.innerHTML = `
            <img src="${imgUrl}" style="max-width: 100%; border-radius: 8px; box-shadow: 0 4px 12px rgba(0,0,0,0.1);">
            <div style="margin-top: 1rem; font-size: 0.9rem; color: #6B7280;">
                SMILES: <code>${finalSmiles}</code>
            </div>
        `;
    } catch (e) {
        resultDiv.innerHTML = `<p style="color: #EF4444;">Error: ${e.message}</p>`;
    }
}

function clearMoleculeViewer() {
    document.getElementById('mv-result').innerHTML = '<span style="color: #9CA3AF;">Visualization will appear here</span>';
}

// --- PubChem ---
async function fetchPubChem() {
    const name = document.getElementById('pubchem-name').value;
    const type = document.getElementById('pubchem-type').value;
    const resultDiv = document.getElementById('pubchem-image');

    resultDiv.innerHTML = '<p>Fetching...</p>';

    try {
        const imgUrl = `http://localhost:5002/img/${encodeURIComponent(name)}?size=large&type=${type}`;
        resultDiv.innerHTML = `
            <img src="${imgUrl}" style="max-width: 100%; border-radius: 8px;">
        `;
    } catch (e) {
        resultDiv.innerHTML = `<p style="color: #EF4444;">Error: ${e.message}</p>`;
    }
}

function open3DViewer() {
    const name = document.getElementById('pubchem-name').value;
    // Try local first
    fetch('http://localhost:5000/')
        .then(r => {
            if (r.ok) window.open(`http://localhost:5000/?q=${encodeURIComponent(name)}`, '_blank');
            else window.open(`http://localhost:5002/static/viewer-3d.html?name=${encodeURIComponent(name)}&embed=true`, '_blank');
        })
        .catch(() => {
            window.open(`http://localhost:5002/static/viewer-3d.html?name=${encodeURIComponent(name)}&embed=true`, '_blank');
        });
}

// --- MolView ---
function loadMolView() {
    const input = document.getElementById('molview-input').value;
    const iframe = document.getElementById('molview-iframe');
    iframe.src = `http://localhost:5003/?q=${encodeURIComponent(input)}`;
}

// --- 3D SMILES ---
async function generate3DSMILES() {
    const name = document.getElementById('opsin-name').value;
    const resultDiv = document.getElementById('opsin-result');

    resultDiv.innerHTML = '<p>Generating...</p>';

    try {
        const resp = await fetch(`http://localhost:5001/api/opsin?name=${encodeURIComponent(name)}`);
        const data = await resp.json();

        if (data.smiles_3d) {
            resultDiv.innerHTML = `
                <div style="background: white; padding: 1.5rem; border-radius: 8px; border: 1px solid #E5E7EB; text-align: left;">
                    <p><strong>3D SMILES:</strong></p>
                    <code style="display: block; background: #F3F4F6; padding: 1rem; border-radius: 6px; margin-top: 0.5rem; word-break: break-all;">${data.smiles_3d}</code>
                    <div style="margin-top: 1rem; display: flex; gap: 1rem;">
                        <button class="btn btn-secondary" onclick="navigator.clipboard.writeText('${data.smiles_3d}')">Copy</button>
                    </div>
                </div>
            `;
        } else {
            throw new Error(data.error || 'Failed');
        }
    } catch (e) {
        resultDiv.innerHTML = `<p style="color: #EF4444;">Error: ${e.message}</p>`;
    }
}

// --- Tests ---
function openTest(file) {
    window.open(file, '_blank');
}

// --- Server Status ---
async function checkAllServers() {
    const servers = [
        { id: 'moleculeviewer', url: 'http://localhost:5000/health' },
        { id: 'mol2chemfig', url: 'http://localhost:5001/health' },
        { id: 'pubchem', url: 'http://localhost:5002/health' },
        { id: 'molview', url: 'http://localhost:5003/' },
        { id: 'docker', url: 'http://localhost:8000/m2cf/reset', method: 'POST' }
    ];

    for (const s of servers) {
        const el = document.getElementById(`status-${s.id}`);
        if (!el) continue;

        el.className = 'node-status checking';
        el.textContent = 'Checking...';

        try {
            const r = await fetch(s.url, { method: s.method || 'GET', mode: 'cors' });
            if (r.ok) {
                el.className = 'node-status online';
                el.textContent = 'Online';
                el.style.background = '#D1FAE5';
                el.style.color = '#065F46';
            } else throw new Error();
        } catch {
            el.className = 'node-status offline';
            el.textContent = 'Offline';
            el.style.background = '#FEE2E2';
            el.style.color = '#991B1B';
        }
    }
}

// --- Flow Diagram Logic (Simplified for V2) ---
// Note: We reuse the SVG drawing logic from previous version but adapted for new IDs if needed.
// For brevity, assuming the SVG structure matches.

// Initialize
window.addEventListener('load', () => {
    checkAllServers();
    setInterval(checkAllServers, 10000);

    // Draw flow connections if on flow tab
    if (typeof drawFlowConnections === 'function') drawFlowConnections();
});
