/**
 * Popup script for Chemistry Formula Renderer v3.0
 * Enhanced with layout mode toggle
 */

// DOM elements - with null checks
const enabledToggle = document.getElementById('enabledToggle');
const mhchemToggle = document.getElementById('mhchemToggle');
const chemfigToggle = document.getElementById('chemfigToggle');
const perfModeToggle = document.getElementById('perfModeToggle');
const devModeToggle = document.getElementById('devModeToggle');
const statusDiv = document.getElementById('status');
const engineInfo = document.getElementById('engineInfo');

// MoleculeViewer options
const moleculeViewerOptions = document.getElementById('moleculeViewerOptions');
const showCarbonsToggle = document.getElementById('showCarbonsToggle');
const showMethylsToggle = document.getElementById('showMethylsToggle');
const aromaticCirclesToggle = document.getElementById('aromaticCirclesToggle');
const fancyBondsToggle = document.getElementById('fancyBondsToggle');
const atomNumbersToggle = document.getElementById('atomNumbersToggle');
const flipHorizontalToggle = document.getElementById('flipHorizontalToggle');
const flipVerticalToggle = document.getElementById('flipVerticalToggle');
const hydrogensSelect = document.getElementById('hydrogensSelect');

// mol2chemfig options
const mol2chemfigOptions = document.getElementById('mol2chemfigOptions');
const m2cfShowCarbonsToggle = document.getElementById('m2cfShowCarbonsToggle');
const m2cfAromaticCirclesToggle = document.getElementById('m2cfAromaticCirclesToggle');
const m2cfShowMethylsToggle = document.getElementById('m2cfShowMethylsToggle');
const m2cfFancyBondsToggle = document.getElementById('m2cfFancyBondsToggle');
const m2cfAtomNumbersToggle = document.getElementById('m2cfAtomNumbersToggle');
const m2cfCompactToggle = document.getElementById('m2cfCompactToggle');
const m2cfFlipHorizontalToggle = document.getElementById('m2cfFlipHorizontalToggle');
const m2cfFlipVerticalToggle = document.getElementById('m2cfFlipVerticalToggle');
const m2cfHydrogensSelect = document.getElementById('m2cfHydrogensSelect');
const use3DSmilesToggle = document.getElementById('use3DSmilesToggle');

// PubChem options
const pubchemOptions = document.getElementById('pubchemOptions');
const pubchemImageSizeSelect = document.getElementById('pubchemImageSizeSelect');
const pubchem3DToggle = document.getElementById('pubchem3DToggle');
const pubchemRecordTypeSelect = document.getElementById('pubchemRecordTypeSelect');

// Image size control options
const saveSizePerImageToggle = document.getElementById('saveSizePerImageToggle');
const saveSizeBySMILESToggle = document.getElementById('saveSizeBySMILESToggle');

// 3D Viewer options
const enable3DViewerToggle = document.getElementById('enable3DViewerToggle');
const default3DViewSelect = document.getElementById('default3DViewSelect');

// Safe selector function
function safeGetElement(id) {
  const el = document.getElementById(id);
  if (!el) console.warn(`Element ${id} not found in DOM`);
  return el;
}

// Load current settings
chrome.storage.sync.get({
  enabled: true,
  renderMhchem: true,
  renderChemfig: true,
  performanceMode: true,
  rendererEngine: 'moleculeviewer',  // Default to MoleculeViewer
  devMode: false,
  showCarbons: false,
  showMethyls: false,
  aromaticCircles: true,
  fancyBonds: true,
  atomNumbers: false,
  flipHorizontal: false,
  flipVertical: false,
  hydrogensMode: 'keep',
  // mol2chemfig options
  m2cfShowCarbons: false,
  m2cfAromaticCircles: false,
  m2cfShowMethyls: false,
  m2cfFancyBonds: false,
  m2cfAtomNumbers: false,
  m2cfCompact: false,
  m2cfFlipHorizontal: false,
  m2cfFlipVertical: false,
  m2cfHydrogensMode: 'keep',
  use3DSmiles: false,
  // Image size controls
  saveSizePerImage: false,
  saveSizeBySMILES: true,  // FIX: Enable by default - global size saving across all pages
  // PubChem options
  pubchemImageSize: 'large',
  pubchem3DEnabled: true,
  pubchemRecordType: '2d',
  // 3D Viewer options
  enable3DViewer: false,
  default3DView: '2d'
}, (settings) => {
  enabledToggle.checked = settings.enabled;
  mhchemToggle.checked = settings.renderMhchem;
  chemfigToggle.checked = settings.renderChemfig;
  perfModeToggle.checked = settings.performanceMode;
  devModeToggle.checked = settings.devMode;
  
  // Load MoleculeViewer options
  if (showCarbonsToggle) showCarbonsToggle.checked = settings.showCarbons;
  if (showMethylsToggle) showMethylsToggle.checked = settings.showMethyls;
  if (aromaticCirclesToggle) aromaticCirclesToggle.checked = settings.aromaticCircles;
  if (fancyBondsToggle) fancyBondsToggle.checked = settings.fancyBonds;
  if (atomNumbersToggle) atomNumbersToggle.checked = settings.atomNumbers;
  if (flipHorizontalToggle) flipHorizontalToggle.checked = settings.flipHorizontal;
  if (flipVerticalToggle) flipVerticalToggle.checked = settings.flipVertical;
  if (hydrogensSelect) hydrogensSelect.value = settings.hydrogensMode;
  
  // Load mol2chemfig options
  if (m2cfShowCarbonsToggle) m2cfShowCarbonsToggle.checked = settings.m2cfShowCarbons;
  if (m2cfAromaticCirclesToggle) m2cfAromaticCirclesToggle.checked = settings.m2cfAromaticCircles;
  if (m2cfShowMethylsToggle) m2cfShowMethylsToggle.checked = settings.m2cfShowMethyls;
  if (m2cfFancyBondsToggle) m2cfFancyBondsToggle.checked = settings.m2cfFancyBonds;
  if (m2cfAtomNumbersToggle) m2cfAtomNumbersToggle.checked = settings.m2cfAtomNumbers;
  if (m2cfCompactToggle) m2cfCompactToggle.checked = settings.m2cfCompact;
  if (m2cfFlipHorizontalToggle) m2cfFlipHorizontalToggle.checked = settings.m2cfFlipHorizontal;
  if (m2cfFlipVerticalToggle) m2cfFlipVerticalToggle.checked = settings.m2cfFlipVertical;
  if (m2cfHydrogensSelect) m2cfHydrogensSelect.value = settings.m2cfHydrogensMode;
  if (use3DSmilesToggle) use3DSmilesToggle.checked = settings.use3DSmiles;

  // Load PubChem options
  if (pubchemImageSizeSelect) pubchemImageSizeSelect.value = settings.pubchemImageSize;
  if (pubchem3DToggle) pubchem3DToggle.checked = settings.pubchem3DEnabled;
  if (pubchemRecordTypeSelect) pubchemRecordTypeSelect.value = settings.pubchemRecordType;

  // Load image size control options
  if (saveSizePerImageToggle) saveSizePerImageToggle.checked = settings.saveSizePerImage;
  if (saveSizeBySMILESToggle) saveSizeBySMILESToggle.checked = settings.saveSizeBySMILES;

  // Load 3D Viewer options
  if (enable3DViewerToggle) enable3DViewerToggle.checked = settings.enable3DViewer;
  if (default3DViewSelect) default3DViewSelect.value = settings.default3DView;

  // Set rendering engine radio button
  const engineRadios = document.querySelectorAll('input[name="renderingEngine"]');
  engineRadios.forEach(radio => {
    radio.checked = (radio.value === settings.rendererEngine);
  });
  
  // Show MoleculeViewer, mol2chemfig, or PubChem options based on selected engine
  if (moleculeViewerOptions) {
    moleculeViewerOptions.style.display = (settings.rendererEngine === 'moleculeviewer') ? 'block' : 'none';
  }
  if (mol2chemfigOptions) {
    mol2chemfigOptions.style.display = (settings.rendererEngine === 'mol2chemfig') ? 'block' : 'none';
  }
  if (pubchemOptions) {
    pubchemOptions.style.display = (settings.rendererEngine === 'pubchem') ? 'block' : 'none';
  }

  // Update engine info
  updateEngineInfo(settings.rendererEngine);
});

// Save settings when changed
enabledToggle.addEventListener('change', (e) => {
  chrome.storage.sync.set({ enabled: e.target.checked }, () => {
    showStatus('Extension ' + (e.target.checked ? 'enabled' : 'disabled') + '. Reload page to apply.', 'success');
  });
});

mhchemToggle.addEventListener('change', (e) => {
  chrome.storage.sync.set({ renderMhchem: e.target.checked }, () => {
    showStatus('mhchem rendering ' + (e.target.checked ? 'enabled' : 'disabled') + '.', 'success');
  });
});

chemfigToggle.addEventListener('change', (e) => {
  chrome.storage.sync.set({ renderChemfig: e.target.checked }, () => {
    showStatus('chemfig rendering ' + (e.target.checked ? 'enabled' : 'disabled') + '.', 'success');
  });
});

perfModeToggle.addEventListener('change', (e) => {
  chrome.storage.sync.set({ performanceMode: e.target.checked }, () => {
    showStatus('Lazy-loading ' + (e.target.checked ? 'enabled' : 'disabled') + '.', 'success');
  });
});

// Dev mode toggle
devModeToggle.addEventListener('change', (e) => {
  chrome.storage.sync.set({ devMode: e.target.checked }, () => {
    showStatus('Developer mode ' + (e.target.checked ? 'enabled' : 'disabled') + '. Show raw chemfig text ' + (e.target.checked ? 'ON' : 'OFF') + '. Reload page to apply.', 'success');
  });
});

// Image size control toggles
if (saveSizePerImageToggle) {
  saveSizePerImageToggle.addEventListener('change', (e) => {
    chrome.storage.sync.set({ saveSizePerImage: e.target.checked }, () => {
      showStatus('Save size per page ' + (e.target.checked ? 'enabled' : 'disabled') + '. Size changes will be remembered for each page.', 'success');
    });
  });
}

if (saveSizeBySMILESToggle) {
  saveSizeBySMILESToggle.addEventListener('change', (e) => {
    chrome.storage.sync.set({ saveSizeBySMILES: e.target.checked }, () => {
      showStatus('Save size by SMILES ' + (e.target.checked ? 'enabled' : 'disabled') + '. Same size will be used for all molecules with the same SMILES.', 'success');
    });
  });
}

// MoleculeViewer rendering options
if (showCarbonsToggle) {
  showCarbonsToggle.addEventListener('change', (e) => {
    chrome.storage.sync.set({ showCarbons: e.target.checked }, () => {
      showStatus('Show carbons ' + (e.target.checked ? 'enabled' : 'disabled') + '. Reload page to apply.', 'success');
    });
  });
}

if (showMethylsToggle) {
  showMethylsToggle.addEventListener('change', (e) => {
    chrome.storage.sync.set({ showMethyls: e.target.checked }, () => {
      showStatus('Show methyls ' + (e.target.checked ? 'enabled' : 'disabled') + '. Reload page to apply.', 'success');
    });
  });
}

if (aromaticCirclesToggle) {
  aromaticCirclesToggle.addEventListener('change', (e) => {
    chrome.storage.sync.set({ aromaticCircles: e.target.checked }, () => {
      showStatus('Aromatic circles ' + (e.target.checked ? 'enabled' : 'disabled') + '. Reload page to apply.', 'success');
    });
  });
}

if (fancyBondsToggle) {
  fancyBondsToggle.addEventListener('change', (e) => {
    chrome.storage.sync.set({ fancyBonds: e.target.checked }, () => {
      showStatus('Fancy bonds ' + (e.target.checked ? 'enabled' : 'disabled') + '. Reload page to apply.', 'success');
    });
  });
}

if (atomNumbersToggle) {
  atomNumbersToggle.addEventListener('change', (e) => {
    chrome.storage.sync.set({ atomNumbers: e.target.checked }, () => {
      showStatus('Atom numbers ' + (e.target.checked ? 'enabled' : 'disabled') + '. Reload page to apply.', 'success');
    });
  });
}

if (flipHorizontalToggle) {
  flipHorizontalToggle.addEventListener('change', (e) => {
    chrome.storage.sync.set({ flipHorizontal: e.target.checked }, () => {
      showStatus('Horizontal flip ' + (e.target.checked ? 'enabled' : 'disabled') + '. Reload page to apply.', 'success');
    });
  });
}

if (flipVerticalToggle) {
  flipVerticalToggle.addEventListener('change', (e) => {
    chrome.storage.sync.set({ flipVertical: e.target.checked }, () => {
      showStatus('Vertical flip ' + (e.target.checked ? 'enabled' : 'disabled') + '. Reload page to apply.', 'success');
    });
  });
}

if (hydrogensSelect) {
  hydrogensSelect.addEventListener('change', (e) => {
    chrome.storage.sync.set({ hydrogensMode: e.target.value }, () => {
      showStatus('Hydrogens set to ' + e.target.value + '. Reload page to apply.', 'success');
    });
  });
}

// mol2chemfig rendering options
if (m2cfShowCarbonsToggle) {
  m2cfShowCarbonsToggle.addEventListener('change', (e) => {
    chrome.storage.sync.set({ m2cfShowCarbons: e.target.checked }, () => {
      showStatus('Show carbons ' + (e.target.checked ? 'enabled' : 'disabled') + '. Reload page to apply.', 'success');
    });
  });
}

if (m2cfAromaticCirclesToggle) {
  m2cfAromaticCirclesToggle.addEventListener('change', (e) => {
    chrome.storage.sync.set({ m2cfAromaticCircles: e.target.checked }, () => {
      showStatus('Aromatic circles ' + (e.target.checked ? 'enabled' : 'disabled') + '. Reload page to apply.', 'success');
    });
  });
}

if (m2cfShowMethylsToggle) {
  m2cfShowMethylsToggle.addEventListener('change', (e) => {
    chrome.storage.sync.set({ m2cfShowMethyls: e.target.checked }, () => {
      showStatus('Show methyls ' + (e.target.checked ? 'enabled' : 'disabled') + '. Reload page to apply.', 'success');
    });
  });
}

if (m2cfFancyBondsToggle) {
  m2cfFancyBondsToggle.addEventListener('change', (e) => {
    chrome.storage.sync.set({ m2cfFancyBonds: e.target.checked }, () => {
      showStatus('Fancy bonds ' + (e.target.checked ? 'enabled' : 'disabled') + '. Reload page to apply.', 'success');
    });
  });
}

if (m2cfAtomNumbersToggle) {
  m2cfAtomNumbersToggle.addEventListener('change', (e) => {
    chrome.storage.sync.set({ m2cfAtomNumbers: e.target.checked }, () => {
      showStatus('Atom numbers ' + (e.target.checked ? 'enabled' : 'disabled') + '. Reload page to apply.', 'success');
    });
  });
}

if (m2cfCompactToggle) {
  m2cfCompactToggle.addEventListener('change', (e) => {
    chrome.storage.sync.set({ m2cfCompact: e.target.checked }, () => {
      showStatus('Compact mode ' + (e.target.checked ? 'enabled' : 'disabled') + '. Reload page to apply.', 'success');
    });
  });
}

if (m2cfFlipHorizontalToggle) {
  m2cfFlipHorizontalToggle.addEventListener('change', (e) => {
    chrome.storage.sync.set({ m2cfFlipHorizontal: e.target.checked }, () => {
      showStatus('Horizontal flip ' + (e.target.checked ? 'enabled' : 'disabled') + '. Reload page to apply.', 'success');
    });
  });
}

if (m2cfFlipVerticalToggle) {
  m2cfFlipVerticalToggle.addEventListener('change', (e) => {
    chrome.storage.sync.set({ m2cfFlipVertical: e.target.checked }, () => {
      showStatus('Vertical flip ' + (e.target.checked ? 'enabled' : 'disabled') + '. Reload page to apply.', 'success');
    });
  });
}

if (m2cfHydrogensSelect) {
  m2cfHydrogensSelect.addEventListener('change', (e) => {
    chrome.storage.sync.set({ m2cfHydrogensMode: e.target.value }, () => {
      showStatus('Hydrogens set to ' + e.target.value + '. Reload page to apply.', 'success');
    });
  });
}

if (use3DSmilesToggle) {
  use3DSmilesToggle.addEventListener('change', (e) => {
    chrome.storage.sync.set({ use3DSmiles: e.target.checked }, () => {
      showStatus('3D SMILES (OPSIN) ' + (e.target.checked ? 'enabled' : 'disabled') + '. Reload page to apply.', 'success');
    });
  });
}

// PubChem options event listeners
if (pubchemImageSizeSelect) {
  pubchemImageSizeSelect.addEventListener('change', (e) => {
    chrome.storage.sync.set({ pubchemImageSize: e.target.value }, () => {
      showStatus('Image size set to ' + e.target.value + '. Reload page to apply.', 'success');
    });
  });
}

if (pubchem3DToggle) {
  pubchem3DToggle.addEventListener('change', (e) => {
    chrome.storage.sync.set({ pubchem3DEnabled: e.target.checked }, () => {
      showStatus('3D models ' + (e.target.checked ? 'enabled' : 'disabled') + '. Reload page to apply.', 'success');
    });
  });
}

if (pubchemRecordTypeSelect) {
  pubchemRecordTypeSelect.addEventListener('change', (e) => {
    chrome.storage.sync.set({ pubchemRecordType: e.target.value }, () => {
      showStatus('Record type set to ' + e.target.value + '. Reload page to apply.', 'success');
    });
  });
}

// 3D Viewer event listeners
if (enable3DViewerToggle) {
  enable3DViewerToggle.addEventListener('change', (e) => {
    chrome.storage.sync.set({ enable3DViewer: e.target.checked }, () => {
      showStatus('3D Viewer ' + (e.target.checked ? 'enabled' : 'disabled') + '. Reload page to apply.', 'success');
    });
  });
}

if (default3DViewSelect) {
  default3DViewSelect.addEventListener('change', (e) => {
    chrome.storage.sync.set({ default3DView: e.target.value }, () => {
      showStatus('Default view set to ' + e.target.value + '. Reload page to apply.', 'success');
    });
  });
}

// Add event listeners for rendering engine radio buttons
const engineRadios = document.querySelectorAll('input[name="renderingEngine"]');
engineRadios.forEach(radio => {
  radio.addEventListener('change', (e) => {
    if (e.target.checked) {
      const engine = e.target.value;
      chrome.storage.sync.set({ rendererEngine: engine }, () => {
        updateEngineInfo(engine);
        showStatus(`Switched to ${engine === 'moleculeviewer' ? '🧪 MoleculeViewer' : '📐 mol2chemfig'}. Reload page to apply.`, 'success');
      });
    }
  });
});

/**
 * Update engine info display
 */
function updateEngineInfo(engine) {
  if (!engineInfo) return;

  if (engine === 'moleculeviewer') {
    engineInfo.textContent = '🧪 MoleculeViewer Server (localhost:5000)';
    if (moleculeViewerOptions) moleculeViewerOptions.style.display = 'block';
    if (mol2chemfigOptions) mol2chemfigOptions.style.display = 'none';
    if (pubchemOptions) pubchemOptions.style.display = 'none';
  } else if (engine === 'mol2chemfig') {
    engineInfo.textContent = '📐 mol2chemfig Server (localhost:8000)';
    if (moleculeViewerOptions) moleculeViewerOptions.style.display = 'none';
    if (mol2chemfigOptions) mol2chemfigOptions.style.display = 'block';
    if (pubchemOptions) pubchemOptions.style.display = 'none';
  } else if (engine === 'pubchem') {
    engineInfo.textContent = '🌐 PubChem Server (localhost:5002)';
    if (moleculeViewerOptions) moleculeViewerOptions.style.display = 'none';
    if (mol2chemfigOptions) mol2chemfigOptions.style.display = 'none';
    if (pubchemOptions) pubchemOptions.style.display = 'block';
  }
}

/**
 * Show status message
 */
function showStatus(message, type) {
  statusDiv.textContent = message;
  statusDiv.className = `status ${type}`;
  
  setTimeout(() => {
    statusDiv.className = 'status';
  }, 2500);
}

/**
 * Open browser console (for developer help)
 */
function openConsole() {
  // Note: In Manifest V3, we can't inject code directly
  // Instead, just show message in popup
  if (statusDiv) {
    showStatus('Open DevTools with F12 to debug. Check console for errors.', 'info');
  }
}

/**
 * Live Code Editor
 */
const codeEditor = document.getElementById('codeEditor');
const renderCodeBtn = document.getElementById('renderCodeBtn');
const codePreview = document.getElementById('codePreview');
const codePreviewContent = document.getElementById('codePreviewContent');

if (renderCodeBtn && codeEditor && codePreview && codePreviewContent) {
  renderCodeBtn.addEventListener('click', () => {
    const code = codeEditor.value.trim();
    if (!code) {
      showStatus('⚠️ Please enter some code to render', 'warning');
      return;
    }

    // Show preview section
    codePreview.style.display = 'block';
    
    // Create a temporary rendering of the code
    // We'll inject it into the preview content div
    codePreviewContent.innerHTML = code;
    
    // Send message to content script to process the preview
    chrome.tabs.query({ active: true, currentWindow: true }, (tabs) => {
      if (tabs[0]) {
        chrome.tabs.sendMessage(tabs[0].id, {
          action: 'renderPreview',
          code: code
        }, (response) => {
          if (response && response.html) {
            codePreviewContent.innerHTML = response.html;
            showStatus('✅ Preview rendered successfully', 'success');
          } else {
            showStatus('⚠️ Could not render preview. Reload page to see changes.', 'warning');
          }
        });
      }
    });
  });
}
