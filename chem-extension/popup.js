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

// Safe selector function
function safeGetElement(id) {
  const el = document.getElementById(id);
  if (!el) console.warn(`Element ${id} not found in DOM`);
  return el;
}

// Load current settings - FORCE MoleculeViewer as the only engine
chrome.storage.sync.get({
  enabled: true,
  renderMhchem: true,
  renderChemfig: true,
  performanceMode: true,
  rendererEngine: 'molecule-viewer',  // ✅ ALWAYS MoleculeViewer
  devMode: false,
  showCarbons: false,
  showMethyls: false,
  aromaticCircles: true,
  fancyBonds: true,
  atomNumbers: false,
  flipHorizontal: false,
  flipVertical: false,
  hydrogensMode: 'keep'
}, (settings) => {
  // Force MoleculeViewer
  settings.rendererEngine = 'molecule-viewer';
  
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
  
  // Show MoleculeViewer options (always visible now)
  if (moleculeViewerOptions) {
    moleculeViewerOptions.style.display = 'block';
  }
  
  // Update engine info
  if (engineInfo) {
    engineInfo.textContent = '🧪 MoleculeViewer Server (localhost:5000)';
  }
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

/**
 * Update engine info display
 */
function updateEngineInfo(engine) {
  const info = rendererInfo[engine];
  engineInfo.textContent = `${info.name}: ${info.desc}`;
  
  // Show/hide MoleculeViewer options based on selected engine
  if (moleculeViewerOptions) {
    if (engine === 'molecule-viewer') {
      moleculeViewerOptions.style.display = 'block';
    } else {
      moleculeViewerOptions.style.display = 'none';
    }
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
