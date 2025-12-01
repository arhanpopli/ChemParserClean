// Tab switching
function switchTab(tabName) {
    // Hide all tabs
    document.querySelectorAll('.tab-content').forEach(tab => {
        tab.classList.remove('active');
        tab.style.display = 'none';
        setTimeout(() => {
            if (tab.classList.contains('active')) {
                tab.style.opacity = '1';
                tab.style.transform = 'translateY(0)';
            } else {
                tab.style.opacity = '0';
                tab.style.transform = 'translateY(10px)';
            }
        }, 10);
    });

    // Remove active from all buttons
    document.querySelectorAll('.tab-button').forEach(btn => {
        btn.classList.remove('active');
    });

    // Show selected tab
    const selectedTab = document.getElementById(tabName);
    if (selectedTab) {
        selectedTab.classList.add('active');
        selectedTab.style.display = 'block';
        setTimeout(() => {
            selectedTab.style.opacity = '1';
            selectedTab.style.transform = 'translateY(0)';
        }, 10);
    }

    // Activate button
    const btn = document.querySelector(`button[onclick="switchTab('${tabName}')"]`);
    if (btn) btn.classList.add('active');
}

// Mol2ChemFig Functions
async function generateMol2ChemFig() {
    const input = document.getElementById('m2cf-input').value;
    const use3D = document.getElementById('m2cf-3d-toggle').checked;
    const resultDiv = document.getElementById('m2cf-result');

    resultDiv.innerHTML = '<p>Generating...</p>';

    try {
        let smiles = input;

        // If 3D is enabled, convert name to 3D SMILES first
        if (use3D) {
            resultDiv.innerHTML = '<p>Converting to 3D SMILES via OPSIN...</p>';
            const opsinResponse = await fetch(`http://localhost:5001/api/opsin?name=${encodeURIComponent(input)}`);
            const opsinData = await opsinResponse.json();

            if (opsinData.smiles_3d) {
                smiles = opsinData.smiles_3d;
                resultDiv.innerHTML += `<p>✓ 3D SMILES: <code>${smiles}</code></p>`;
            }
        }

        // Generate with Mol2ChemFig
        const response = await fetch('http://localhost:5001/api/generate', {
            method: 'POST',
            headers: { 'Content-Type': 'application/json' },
            body: JSON.stringify({ smiles: smiles })
        });

        const data = await response.json();

        if (data.svg) {
            resultDiv.innerHTML = `
                <h3>✓ Generated Successfully</h3>
                <p><strong>Input:</strong> ${input} ${use3D ? '(3D)' : ''}</p>
                <p><strong>SMILES:</strong> <code>${smiles}</code></p>
                <div>${data.svg}</div>
            `;
        } else {
            resultDiv.innerHTML = `<p style="color: red;">Error: ${data.error || 'Failed to generate'}</p>`;
        }
    } catch (error) {
        resultDiv.innerHTML = `<p style="color: red;">Error: ${error.message}</p>`;
    }
}

function clearMol2ChemFig() {
    document.getElementById('m2cf-result').innerHTML = 'Results will appear here...';
}

// MoleculeViewer Functions
async function generateMoleculeViewer() {
    const smiles = document.getElementById('mv-smiles').value;
    const use3D = document.getElementById('mv-3d-toggle').checked;
    const width = document.getElementById('mv-width').value;
    const height = document.getElementById('mv-height').value;
    const resultDiv = document.getElementById('mv-result');

    resultDiv.innerHTML = '<p>Generating...</p>';

    try {
        let finalSmiles = smiles;

        // If 3D is enabled, convert to 3D SMILES first
        if (use3D) {
            resultDiv.innerHTML = '<p>Converting to 3D SMILES via OPSIN...</p>';
            const opsinResponse = await fetch(`http://localhost:5001/api/opsin?name=${encodeURIComponent(smiles)}`);
            const opsinData = await opsinResponse.json();

            if (opsinData.smiles_3d) {
                finalSmiles = opsinData.smiles_3d;
                resultDiv.innerHTML += `<p>✓ 3D SMILES: <code>${finalSmiles}</code></p>`;
            }
        }

        // Generate image
        const imageUrl = `http://localhost:5000/img/smiles?smiles=${encodeURIComponent(finalSmiles)}&width=${width}&height=${height}`;

        resultDiv.innerHTML = `
            <h3>✓ Generated Successfully</h3>
            <p><strong>SMILES:</strong> <code>${finalSmiles}</code> ${use3D ? '(3D)' : ''}</p>
            <img src="${imageUrl}" alt="Molecule" style="max-width: 100%;">
            <p><a href="${imageUrl}" target="_blank">Open in new tab</a></p>
        `;
    } catch (error) {
        resultDiv.innerHTML = `<p style="color: red;">Error: ${error.message}</p>`;
    }
}

function clearMoleculeViewer() {
    document.getElementById('mv-result').innerHTML = 'Results will appear here...';
}

// PubChem Functions
async function fetchPubChem() {
    const name = document.getElementById('pubchem-name').value;
    const type = document.getElementById('pubchem-type').value;
    const size = document.getElementById('pubchem-size').value;
    const imageDiv = document.getElementById('pubchem-image');

    imageDiv.innerHTML = '<p>Fetching from PubChem...</p>';

    try {
        const imageUrl = `http://localhost:5002/img/${encodeURIComponent(name)}?size=${size}&type=${type}`;

        imageDiv.innerHTML = `
            <h3>✓ ${name} (${type.toUpperCase()})</h3>
            <img src="${imageUrl}" alt="${name}" style="max-width: 100%; background: white; padding: 20px; border-radius: 8px;">
            <p><a href="${imageUrl}" target="_blank">Open in new tab</a></p>
        `;
    } catch (error) {
        imageDiv.innerHTML = `<p style="color: red;">Error: ${error.message}</p>`;
    }
}

function open3DViewer() {
    const name = document.getElementById('pubchem-name').value;

    // Try to use the local MolView server first, fallback to the old method
    fetch('http://localhost:5000/')
        .then(response => {
            if (response.ok) {
                window.open(`http://localhost:5000/?q=${encodeURIComponent(name)}`, '_blank');
            } else {
                window.open(`http://localhost:5002/static/viewer-3d.html?name=${encodeURIComponent(name)}&embed=true`, '_blank');
            }
        })
        .catch(error => {
            console.log('Local MolView server not available, using fallback');
            window.open(`http://localhost:5002/static/viewer-3d.html?name=${encodeURIComponent(name)}&embed=true`, '_blank');
        });
}

function clearPubChem() {
    document.getElementById('pubchem-image').innerHTML = '';
}

// Check PubChem server status
async function checkPubChemStatus() {
    try {
        const response = await fetch('http://localhost:5002/health');
        if (response.ok) {
            const el = document.getElementById('pubchem-status');
            if (el) el.className = 'status-indicator status-online';
            const txt = document.getElementById('pubchem-status-text');
            if (txt) txt.textContent = 'Online';
        } else {
            throw new Error('Not responding');
        }
    } catch {
        const el = document.getElementById('pubchem-status');
        if (el) el.className = 'status-indicator status-offline';
        const txt = document.getElementById('pubchem-status-text');
        if (txt) txt.textContent = 'Offline';
    }
}

// 3D SMILES Functions
async function generate3DSMILES() {
    const name = document.getElementById('opsin-name').value;
    const resultDiv = document.getElementById('opsin-result');

    resultDiv.innerHTML = '<p>Converting to 3D SMILES via OPSIN...</p>';

    try {
        const response = await fetch(`http://localhost:5001/api/opsin?name=${encodeURIComponent(name)}`);
        const data = await response.json();

        if (data.smiles_3d) {
            resultDiv.innerHTML = `
                <h3>✓ 3D SMILES Generated</h3>
                <p><strong>Name:</strong> ${name}</p>
                <p><strong>2D SMILES:</strong> <code>${data.smiles || 'N/A'}</code></p>
                <p><strong>3D SMILES:</strong> <code>${data.smiles_3d}</code></p>
                <button class="btn" onclick="useSMILES('${data.smiles_3d}')">Use in Mol2ChemFig</button>
                <button class="btn" onclick="useSMILESMV('${data.smiles_3d}')">Use in MoleculeViewer</button>
            `;
        } else {
            resultDiv.innerHTML = `<p style="color: red;">Error: ${data.error || 'Failed to generate 3D SMILES'}</p>`;
        }
    } catch (error) {
        resultDiv.innerHTML = `<p style="color: red;">Error: ${error.message}</p>`;
    }
}

function clear3DSMILES() {
    document.getElementById('opsin-result').innerHTML = 'Results will appear here...';
}

function useSMILES(smiles) {
    switchTab('mol2chemfig');
    document.getElementById('m2cf-input').value = smiles;
}

function useSMILESMV(smiles) {
    switchTab('moleculeviewer');
    document.getElementById('mv-smiles').value = smiles;
}

// Test Functions
function openTest(filename) {
    window.open(filename, '_blank');
}

// MolView Functions
function loadMolView() {
    const input = document.getElementById('molview-input').value;
    if (!input.trim()) {
        alert('Please enter a compound name, SMILES, or InChI');
        return;
    }
    const iframe = document.getElementById('molview-iframe');
    iframe.src = `http://localhost:5003/?q=${encodeURIComponent(input)}`;
}

function resetMolView() {
    const iframe = document.getElementById('molview-iframe');
    iframe.src = 'http://localhost:5003/';
}

// ========================================
// FLOW DIAGRAM FUNCTIONS
// ========================================

const flowConnections = [
    { from: 'page', to: 'extension', label: '"histamine"', color: '#8b5cf6', flow: 'all' },
    { from: 'options', to: 'extension', label: '3D? Engine?', color: '#ec4899', flow: 'all' },
    { from: 'extension', to: 'smiles-bridge', label: 'Name + Options', color: '#667eea', flow: 'all' },
    { from: 'smiles-bridge', to: 'opsin-3d', label: 'If 3D enabled', color: '#f59e0b', flow: 'all' },
    { from: 'smiles-bridge', to: 'opsin', label: 'Try first', color: '#f59e0b', flow: 'all' },
    { from: 'smiles-bridge', to: 'pubchem-api', label: 'Fallback', color: '#f59e0b', flow: 'all' },
    { from: 'smiles-bridge', to: 'selector', label: 'SMILES string', color: '#fbbf24', flow: 'all' },
    { from: 'selector', to: 'moleculeviewer', label: 'SMILES', color: '#10b981', flow: 'moleculeviewer' },
    { from: 'selector', to: 'mol2chemfig', label: 'SMILES+opts', color: '#10b981', flow: 'mol2chemfig' },
    { from: 'selector', to: 'pubchem-local', label: 'Name/CID', color: '#10b981', flow: 'pubchem' },
    { from: 'moleculeviewer', to: 'rdkit', label: 'Generate', color: '#06b6d4', flow: 'moleculeviewer' },
    { from: 'mol2chemfig', to: 'docker', label: 'ChemFig', color: '#8b5cf6', flow: 'mol2chemfig' },
    { from: 'pubchem-local', to: 'molview', label: '3D View', color: '#f59e0b', flow: 'pubchem' },
    { from: 'rdkit', to: 'cache', label: 'SVG', color: '#14b8a6', flow: 'moleculeviewer' },
    { from: 'docker', to: 'cache', label: 'SVG', color: '#14b8a6', flow: 'mol2chemfig' },
    { from: 'cache', to: 'output', label: 'URL', color: '#10b981', flow: 'all' },
    { from: 'pubchem-local', to: 'output', label: 'PNG', color: '#10b981', flow: 'pubchem' },
    { from: 'output', to: 'darkmode', label: 'If dark', color: '#6366f1', flow: 'all' },
    { from: 'darkmode', to: 'extension', label: 'Final image', color: '#8b5cf6', flow: 'all' },
    { from: 'output', to: 'extension', label: 'Final image', color: '#10b981', flow: 'all' },
];

const nodeDetails = {
    'page': { title: '📄 Web Page', content: '<p>Contains text with chemical notation patterns.</p>' },
    'options': { title: '⚙️ User Options', content: '<p>Settings from popup.js.</p>' },
    'extension': { title: '🧩 Chrome Extension', content: '<p>Central hub for detection and coordination.</p>' },
    'smiles-bridge': { title: '🌉 SMILES Bridge', content: '<p>Converts names to SMILES strings.</p>' },
    'opsin-3d': { title: '🔮 OPSIN 3D', content: '<p>Generates 3D stereochemistry SMILES.</p>' },
    'opsin': { title: '🔬 OPSIN 2D', content: '<p>Standard IUPAC to SMILES conversion.</p>' },
    'pubchem-api': { title: '☁️ PubChem API', content: '<p>Fallback for common names.</p>' },
    'selector': { title: '🔀 Renderer Selector', content: '<p>Routes to chosen engine.</p>' },
    'moleculeviewer': { title: '🧬 MoleculeViewer', content: '<p>Node.js server with RDKit.</p>' },
    'mol2chemfig': { title: '📐 Mol2ChemFig', content: '<p>Python server with Docker backend.</p>' },
    'pubchem-local': { title: '🌐 PubChem Renderer', content: '<p>Fetches images from PubChem.</p>' },
    'docker': { title: '🐳 Docker Backend', content: '<p>Runs LaTeX and dvisvgm.</p>' },
    'rdkit': { title: '⚗️ RDKit Library', content: '<p>Python library for chemical rendering.</p>' },
    'molview': { title: '🎮 MolView.org', content: '<p>Interactive 3D viewer iframe.</p>' },
    'output': { title: '🖼️ Rendered Output', content: '<p>Final SVG or PNG image.</p>' },
    'darkmode': { title: '🌙 Dark Mode', content: '<p>Inverts colors for dark themes.</p>' },
    'cache': { title: '💾 Cache', content: '<p>Stores generated images.</p>' }
};

function drawFlowConnections() {
    const svg = document.getElementById('flow-connections');
    const diagram = document.getElementById('flow-diagram');
    if (!svg || !diagram) return;

    svg.querySelectorAll('path, text').forEach(el => el.remove());

    flowConnections.forEach((conn) => {
        const fromNode = document.getElementById(`node-${conn.from}`);
        const toNode = document.getElementById(`node-${conn.to}`);

        if (!fromNode || !toNode) return;

        const fromRect = fromNode.getBoundingClientRect();
        const toRect = toNode.getBoundingClientRect();
        const diagramRect = diagram.getBoundingClientRect();

        const x1 = fromRect.left + fromRect.width / 2 - diagramRect.left;
        const y1 = fromRect.top + fromRect.height / 2 - diagramRect.top;
        const x2 = toRect.left + toRect.width / 2 - diagramRect.left;
        const y2 = toRect.top + toRect.height / 2 - diagramRect.top;

        const midX = (x1 + x2) / 2;
        const midY = (y1 + y2) / 2;
        const dx = x2 - x1;
        const dy = y2 - y1;
        const offset = Math.min(Math.abs(dx), Math.abs(dy)) * 0.3;
        const cx = midX + (dy > 0 ? offset : -offset);
        const cy = midY + (dx > 0 ? -offset : offset);

        const path = document.createElementNS('http://www.w3.org/2000/svg', 'path');
        path.setAttribute('d', `M ${x1} ${y1} Q ${cx} ${cy} ${x2} ${y2}`);
        path.setAttribute('stroke', conn.color);
        path.setAttribute('class', `flow-line flow-${conn.flow}`);
        path.setAttribute('marker-end', `url(#arrowhead)`);
        path.style.fill = 'none';

        svg.appendChild(path);

        if (conn.label) {
            const text = document.createElementNS('http://www.w3.org/2000/svg', 'text');
            text.setAttribute('x', cx);
            text.setAttribute('y', cy - 5);
            text.setAttribute('class', `flow-line-label flow-${conn.flow}`);
            text.textContent = conn.label;
            text.style.fill = 'white';
            text.style.fontSize = '10px';
            text.style.textAnchor = 'middle';
            svg.appendChild(text);
        }
    });
}

function highlightFlow(flowType) {
    const lines = document.querySelectorAll('.flow-line');
    const labels = document.querySelectorAll('.flow-line-label');

    if (flowType === 'all') {
        lines.forEach(l => { l.classList.remove('dimmed', 'highlighted'); });
        labels.forEach(l => { l.classList.remove('dimmed'); });
    } else {
        lines.forEach(l => {
            if (l.classList.contains(`flow-${flowType}`) || l.classList.contains('flow-all')) {
                l.classList.add('highlighted');
                l.classList.remove('dimmed');
            } else {
                l.classList.add('dimmed');
                l.classList.remove('highlighted');
            }
        });
    }
}

function showNodeDetails(nodeId) {
    const details = nodeDetails[nodeId];
    if (!details) return;

    const panel = document.getElementById('node-details');
    const title = document.getElementById('details-title');
    const content = document.getElementById('details-content');

    title.textContent = details.title;
    content.innerHTML = details.content;
    panel.style.display = 'block';
    panel.scrollIntoView({ behavior: 'smooth', block: 'nearest' });
}

async function checkAllServers() {
    const servers = [
        { id: 'moleculeviewer', url: 'http://localhost:5000/health' },
        { id: 'mol2chemfig', url: 'http://localhost:5001/health' },
        { id: 'pubchem', url: 'http://localhost:5002/health' },
        { id: 'molview', url: 'http://localhost:5003/' },
        { id: 'docker', url: 'http://localhost:8000/m2cf/reset', method: 'POST' }
    ];

    for (const server of servers) {
        const statusEl = document.getElementById(`status-${server.id}`);
        if (!statusEl) continue;

        statusEl.textContent = 'Checking...';
        statusEl.className = 'node-status checking';

        try {
            const response = await fetch(server.url, {
                method: server.method || 'GET',
                mode: 'cors'
            });

            if (response.ok) {
                statusEl.textContent = 'Online';
                statusEl.className = 'node-status online';
            } else {
                throw new Error('Not OK');
            }
        } catch (e) {
            statusEl.textContent = 'Offline';
            statusEl.className = 'node-status offline';
        }
    }
}

function resetFlowView() {
    document.getElementById('flow-select').value = 'all';
    highlightFlow('all');
    document.getElementById('node-details').style.display = 'none';
}

// Drag and Drop
let draggedNode = null;
let dragOffset = { x: 0, y: 0 };
let isDragging = false;

function initDraggable() {
    const nodes = document.querySelectorAll('.flow-node.draggable');
    nodes.forEach(node => {
        node.addEventListener('mousedown', startDrag);
    });
    document.addEventListener('mousemove', drag);
    document.addEventListener('mouseup', stopDrag);
}

function startDrag(e) {
    if (e.target.closest('.node-status')) return;
    draggedNode = e.currentTarget;
    isDragging = false;
    const rect = draggedNode.getBoundingClientRect();
    dragOffset.x = e.clientX - rect.left;
    dragOffset.y = e.clientY - rect.top;
    draggedNode.classList.add('dragging');
}

function drag(e) {
    if (!draggedNode) return;
    isDragging = true;
    const diagram = document.getElementById('flow-diagram');
    const diagramRect = diagram.getBoundingClientRect();

    let newX = e.clientX - diagramRect.left - dragOffset.x;
    let newY = e.clientY - diagramRect.top - dragOffset.y;

    draggedNode.style.left = newX + 'px';
    draggedNode.style.top = newY + 'px';

    drawFlowConnections();
}

function stopDrag(e) {
    if (draggedNode) {
        draggedNode.classList.remove('dragging');
        if (!isDragging) {
            const nodeId = draggedNode.id.replace('node-', '');
            showNodeDetails(nodeId);
        }
        draggedNode = null;
    }
}

// Initialize
window.addEventListener('load', () => {
    checkPubChemStatus();
    setInterval(checkPubChemStatus, 10000);
    drawFlowConnections();
    checkAllServers();
    initDraggable();
    window.addEventListener('resize', drawFlowConnections);
});
