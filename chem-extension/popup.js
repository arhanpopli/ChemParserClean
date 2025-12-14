/**
 * Popup script for ChemTex v6.0
 * SmilesDrawer client-side rendering
 */

// DOM elements - with null checks
const enabledToggle = document.getElementById('enabledToggle');
const mhchemToggle = document.getElementById('mhchemToggle');
const chemfigToggle = document.getElementById('chemfigToggle');
const perfModeToggle = document.getElementById('perfModeToggle');
const devModeToggle = document.getElementById('devModeToggle');
const showTagsToggle = document.getElementById('showTagsToggle');
const reloadAllBtn = document.getElementById('reloadAllBtn');
const statusDiv = document.getElementById('status');
const engineInfo = document.getElementById('engineInfo');

// SmilesDrawer options
const smilesDrawerOptions = document.getElementById('smilesDrawerOptions');
const sdShowCarbonsToggle = document.getElementById('sdShowCarbonsToggle');
const sdAromaticRingsToggle = document.getElementById('sdAromaticRingsToggle');
const sdShowMethylsToggle = document.getElementById('sdShowMethylsToggle');
const sdAtomNumbersToggle = document.getElementById('sdAtomNumbersToggle');
const sdShowExplicitHydrogensToggle = document.getElementById('sdShowExplicitHydrogensToggle');
const sdShowImplicitHydrogensToggle = document.getElementById('sdShowImplicitHydrogensToggle');
const sdCompactDrawingToggle = document.getElementById('sdCompactDrawingToggle');
const sdFlipHorizontalToggle = document.getElementById('sdFlipHorizontalToggle');
const sdFlipVerticalToggle = document.getElementById('sdFlipVerticalToggle');
const sdThemeSelect = document.getElementById('sdThemeSelect');
const sdRotateSlider = document.getElementById('sdRotateSlider');
const sdRotateValue = document.getElementById('sdRotateValue');
const sdBondThicknessSlider = document.getElementById('sdBondThicknessSlider');
const sdBondThicknessValue = document.getElementById('sdBondThicknessValue');
const sdGradientColorsToggle = document.getElementById('sdGradientColorsToggle');
const sdScaleByWeightToggle = document.getElementById('sdScaleByWeightToggle');

// Data Source toggles (IntegratedSearch)
const searchPubChemToggle = document.getElementById('searchPubChemToggle');
const searchRCSBToggle = document.getElementById('searchRCSBToggle');
const searchCODToggle = document.getElementById('searchCODToggle');

// Compound Options
const compoundOptions = document.getElementById('compoundOptions');
const molviewRepresentationSelect = document.getElementById('molviewRepresentationSelect');
const compoundMolviewBgColorSelect = document.getElementById('compoundMolviewBgColorSelect');
const molviewEngineSelect = document.getElementById('molviewEngineSelect');

// Biomolecule Options
const proteinRemoveWhiteBgToggle = document.getElementById('proteinRemoveWhiteBgToggle');
const molviewBioAssemblyToggle = document.getElementById('molviewBioAssemblyToggle');
const molviewChainTypeSelect = document.getElementById('molviewChainTypeSelect');
const molviewChainBondsToggle = document.getElementById('molviewChainBondsToggle');
const molviewChainColorSelect = document.getElementById('molviewChainColorSelect');
const proteinMolviewBgColorSelect = document.getElementById('proteinMolviewBgColorSelect');

// Mineral Options
const mineralRepresentationSelect = document.getElementById('mineralRepresentationSelect');
const mineralMolviewBgColorSelect = document.getElementById('mineralMolviewBgColorSelect');
const mineralCrystallographySelect = document.getElementById('mineralCrystallographySelect');

// Image size control options
const saveSizePerImageToggle = document.getElementById('saveSizePerImageToggle');
const saveSizeBySMILESToggle = document.getElementById('saveSizeBySMILESToggle');


// MolView-Only Mode toggle
const molSearchModeToggle = document.getElementById('molSearchModeToggle');

// Disable Formula Fallback toggle
const disableFormulaFallbackToggle = document.getElementById('disableFormulaFallbackToggle');

// Safe selector function
function safeGetElement(id) {
  const el = document.getElementById(id);
  if (!el) console.warn(`Element ${id} not found in DOM`);
  return el;
}

// Helper function to save settings with verification
function saveSettingWithVerification(settingName, value, callback) {
  console.log(`[Popup] Saving ${settingName}:`, value);
  const settingObj = {};
  settingObj[settingName] = value;

  chrome.storage.sync.set(settingObj, () => {
    if (chrome.runtime.lastError) {
      console.error(`[Popup] Error saving ${settingName}:`, chrome.runtime.lastError);
      if (callback) callback(false, chrome.runtime.lastError.message);
      return;
    }

    // Verify the save by reading it back
    chrome.storage.sync.get([settingName], (result) => {
      if (chrome.runtime.lastError) {
        console.error(`[Popup] Error verifying ${settingName}:`, chrome.runtime.lastError);
        if (callback) callback(false, chrome.runtime.lastError.message);
        return;
      }

      const savedValue = result[settingName];
      console.log(`[Popup] Verified ${settingName} saved as:`, savedValue);

      if (savedValue === value) {
        console.log(`[Popup] Successfully saved and verified ${settingName}:`, value);
        if (callback) callback(true);
      } else {
        console.error(`[Popup] Save verification failed for ${settingName}. Expected:`, value, 'Got:', savedValue);
        if (callback) callback(false, 'Verification failed');
      }
    });
  });
}

// Load current settings
chrome.storage.sync.get(null, (storedSettings) => {
  // Define defaults that will be used ONLY if a setting is undefined
  const defaults = {
    enabled: true,
    renderMhchem: true,
    renderChemfig: true,
    performanceMode: true,
    rendererEngine: 'client-side',
    devMode: false,
    sdShowCarbons: true,
    sdAromaticRings: true,
    sdShowMethyls: true,
    sdAtomNumbers: false,
    sdShowExplicitHydrogens: false,
    sdShowImplicitHydrogens: true,
    sdCompactDrawing: false,
    sdFlipHorizontal: false,
    sdFlipVertical: false,
    sdTheme: 'light',
    sdRotate: 0,
    sdBondThickness: 1.0,
    sdGradientColors: false,
    sdScaleByWeight: false,
    saveSizePerImage: false,
    saveSizeBySMILES: true,
    searchPubChem: true,
    searchRCSB: true,
    searchCOD: true,
    useStereochemistry: true,
    compoundMolviewBgColor: 'black',
    disableFormulaFallback: false,
    enable3DViewer: false,
    default3DView: '2d',
    viewer3DSource: '3dmol',
    viewer3DStyle: 'stick:sphere',
    viewer3DAutoRotate: true,
    viewer3DSize: 'normal',
    viewer3DBgColor: '#1a1a2e',
    molviewRepresentation: 'ballAndStick',
    molviewEngine: 'glmol',

    proteinRemoveWhiteBg: false,
    molviewBioAssembly: false,
    molviewChainType: 'ribbon',
    molviewChainBonds: false,
    molviewChainColor: 'ss',
    proteinMolviewBgColor: 'black',
    mineralRepresentation: 'ballAndStick',
    mineralMolviewBgColor: 'black',
    mineralCrystallography: 'supercell_2x2x2',
    viewer3DStyleSettings: {
      'stick': { stickRadius: '0.15' },
      'line': {},
      'cross': {},
      'sphere': { sphereRadius: '0.7' },
      'stick:sphere': { stickRadius: '0.15', sphereRadius: '0.3' },
      'cartoon': {}
    },
    enableAIMolecularControl: false,
    molSearchMode: false,
    showTags: false,
    // UI settings
    uiBlur: 7,
    uiOpacity: 23,
    // Molecule animation settings
    moleculeCount: 22,
    moleculeScale: 0.4,
    cursorGravityStrength: -0.9,
    initialVelocity: 10
  };

  // Merge stored settings with defaults (stored settings take precedence)
  const settings = {};
  for (const key in defaults) {
    settings[key] = storedSettings[key] !== undefined ? storedSettings[key] : defaults[key];
  }

  console.log('[Popup] ========== SETTINGS LOADING DEBUG ==========');
  console.log('[Popup] Raw storage (storedSettings):', storedSettings);
  console.log('[Popup] Defaults object:', defaults);
  console.log('[Popup] Final merged settings:', settings);
  console.log('[Popup] SmilesDrawer settings specifically:', {
    sdShowCarbons: settings.sdShowCarbons,
    sdAromaticRings: settings.sdAromaticRings,
    sdShowMethyls: settings.sdShowMethyls,
    sdAtomNumbers: settings.sdAtomNumbers,
    sdShowExplicitHydrogens: settings.sdShowExplicitHydrogens,
    sdShowImplicitHydrogens: settings.sdShowImplicitHydrogens,
    sdCompactDrawing: settings.sdCompactDrawing,
    sdFlipHorizontal: settings.sdFlipHorizontal,
    sdFlipVertical: settings.sdFlipVertical
  });
  console.log('[Popup] ================================================');

  if (enabledToggle) enabledToggle.checked = settings.enabled;
  if (mhchemToggle) mhchemToggle.checked = settings.renderMhchem;
  if (chemfigToggle) chemfigToggle.checked = settings.renderChemfig;
  if (perfModeToggle) perfModeToggle.checked = settings.performanceMode;
  if (devModeToggle) devModeToggle.checked = settings.devMode;
  if (showTagsToggle) showTagsToggle.checked = settings.showTags;

  console.log('[Popup] ========== TOGGLE ELEMENT CHECK ==========');
  console.log('[Popup] sdShowCarbonsToggle exists?', !!sdShowCarbonsToggle);
  console.log('[Popup] sdAromaticRingsToggle exists?', !!sdAromaticRingsToggle);
  console.log('[Popup] sdShowMethylsToggle exists?', !!sdShowMethylsToggle);
  console.log('[Popup] sdAtomNumbersToggle exists?', !!sdAtomNumbersToggle);
  console.log('[Popup] sdShowExplicitHydrogensToggle exists?', !!sdShowExplicitHydrogensToggle);
  console.log('[Popup] sdShowImplicitHydrogensToggle exists?', !!sdShowImplicitHydrogensToggle);
  console.log('[Popup] sdFlipHorizontalToggle exists?', !!sdFlipHorizontalToggle);
  console.log('[Popup] sdFlipVerticalToggle exists?', !!sdFlipVerticalToggle);
  console.log('[Popup] ================================================');

  // Load SmilesDrawer options
  if (sdShowCarbonsToggle) {
    sdShowCarbonsToggle.checked = settings.sdShowCarbons;
    console.log('[Popup] Set sdShowCarbonsToggle.checked to:', settings.sdShowCarbons, '| Actual value now:', sdShowCarbonsToggle.checked);
  }
  if (sdAromaticRingsToggle) {
    sdAromaticRingsToggle.checked = settings.sdAromaticRings;
    console.log('[Popup] Set sdAromaticRingsToggle.checked to:', settings.sdAromaticRings, '| Actual value now:', sdAromaticRingsToggle.checked);
  }
  if (sdShowMethylsToggle) {
    sdShowMethylsToggle.checked = settings.sdShowMethyls;
    console.log('[Popup] Set sdShowMethylsToggle.checked to:', settings.sdShowMethyls, '| Actual value now:', sdShowMethylsToggle.checked);
  }
  if (sdAtomNumbersToggle) {
    sdAtomNumbersToggle.checked = settings.sdAtomNumbers;
    console.log('[Popup] Set sdAtomNumbersToggle.checked to:', settings.sdAtomNumbers, '| Actual value now:', sdAtomNumbersToggle.checked);
  }
  if (sdShowExplicitHydrogensToggle) {
    sdShowExplicitHydrogensToggle.checked = settings.sdShowExplicitHydrogens;
    console.log('[Popup] Set sdShowExplicitHydrogensToggle.checked to:', settings.sdShowExplicitHydrogens, '| Actual value now:', sdShowExplicitHydrogensToggle.checked);
  }
  if (sdShowImplicitHydrogensToggle) {
    sdShowImplicitHydrogensToggle.checked = settings.sdShowImplicitHydrogens;
    console.log('[Popup] Set sdShowImplicitHydrogensToggle.checked to:', settings.sdShowImplicitHydrogens, '| Actual value now:', sdShowImplicitHydrogensToggle.checked);
  }
  if (sdCompactDrawingToggle) {
    sdCompactDrawingToggle.checked = settings.sdCompactDrawing;
    console.log('[Popup] Set sdCompactDrawingToggle.checked to:', settings.sdCompactDrawing, '| Actual value now:', sdCompactDrawingToggle.checked);
  }
  if (sdFlipHorizontalToggle) {
    sdFlipHorizontalToggle.checked = settings.sdFlipHorizontal;
    console.log('[Popup] Set sdFlipHorizontalToggle.checked to:', settings.sdFlipHorizontal, '| Actual value now:', sdFlipHorizontalToggle.checked);
  }
  if (sdFlipVerticalToggle) {
    sdFlipVerticalToggle.checked = settings.sdFlipVertical;
    console.log('[Popup] Set sdFlipVerticalToggle.checked to:', settings.sdFlipVertical, '| Actual value now:', sdFlipVerticalToggle.checked);
  }
  if (sdThemeSelect) sdThemeSelect.value = settings.sdTheme || 'light';
  if (sdRotateSlider) {
    sdRotateSlider.value = settings.sdRotate || 0;
    if (sdRotateValue) sdRotateValue.textContent = (settings.sdRotate || 0) + '°';
  }
  if (sdBondThicknessSlider) {
    sdBondThicknessSlider.value = settings.sdBondThickness || 1.0;
    if (sdBondThicknessValue) sdBondThicknessValue.textContent = (settings.sdBondThickness || 1.0).toFixed(1);
  }
  if (sdGradientColorsToggle) sdGradientColorsToggle.checked = settings.sdGradientColors || false;
  if (sdScaleByWeightToggle) sdScaleByWeightToggle.checked = settings.sdScaleByWeight || false;

  // Load Data Source toggles
  if (searchPubChemToggle) searchPubChemToggle.checked = settings.searchPubChem !== false;
  if (searchRCSBToggle) searchRCSBToggle.checked = settings.searchRCSB !== false;
  if (searchCODToggle) searchCODToggle.checked = settings.searchCOD !== false;

  // Load stereochemistry preference
  const useStereochemistryToggle = document.getElementById('useStereochemistryToggle');
  if (useStereochemistryToggle) useStereochemistryToggle.checked = settings.useStereochemistry !== false;

  // Load image size control options
  if (saveSizePerImageToggle) saveSizePerImageToggle.checked = settings.saveSizePerImage;
  if (saveSizeBySMILESToggle) saveSizeBySMILESToggle.checked = settings.saveSizeBySMILES;

  // Load 3D Viewer options
  if (enable3DViewerToggle) enable3DViewerToggle.checked = settings.enable3DViewer;
  if (default3DViewSelect) default3DViewSelect.value = settings.default3DView;
  if (viewer3DSourceSelect) viewer3DSourceSelect.value = settings.viewer3DSource || '3dmol';
  if (viewer3DStyleSelect) viewer3DStyleSelect.value = settings.viewer3DStyle || 'stick:sphere';
  if (viewer3DAutoRotateToggle) viewer3DAutoRotateToggle.checked = settings.viewer3DAutoRotate !== false;
  if (viewer3DSizeSelect) viewer3DSizeSelect.value = settings.viewer3DSize || 'normal';
  if (viewer3DBgColorSelect) viewer3DBgColorSelect.value = settings.viewer3DBgColor || '#1a1a2e';

  // Load MolView specific options (compound)
  console.log('[Popup] Loading compound settings:', {
    storedBgColor: storedSettings.compoundMolviewBgColor,
    mergedBgColor: settings.compoundMolviewBgColor,
    willSetTo: settings.compoundMolviewBgColor || 'black'
  });
  if (molviewRepresentationSelect) molviewRepresentationSelect.value = settings.molviewRepresentation || 'ballAndStick';
  if (compoundMolviewBgColorSelect) {
    compoundMolviewBgColorSelect.value = settings.compoundMolviewBgColor || 'black';
    console.log('[Popup] Set compoundMolviewBgColorSelect.value to:', compoundMolviewBgColorSelect.value);
  }
  if (molviewEngineSelect) molviewEngineSelect.value = settings.molviewEngine || 'glmol';

  // Load Biomolecule Options
  if (proteinRemoveWhiteBgToggle) proteinRemoveWhiteBgToggle.checked = settings.proteinRemoveWhiteBg;
  if (molviewBioAssemblyToggle) molviewBioAssemblyToggle.checked = settings.molviewBioAssembly;
  if (molviewChainTypeSelect) molviewChainTypeSelect.value = settings.molviewChainType || 'ribbon';
  if (molviewChainBondsToggle) molviewChainBondsToggle.checked = settings.molviewChainBonds;
  if (molviewChainColorSelect) molviewChainColorSelect.value = settings.molviewChainColor || 'ss';
  if (proteinMolviewBgColorSelect) proteinMolviewBgColorSelect.value = settings.proteinMolviewBgColor || 'black';

  // Load Mineral Options
  if (mineralRepresentationSelect) mineralRepresentationSelect.value = settings.mineralRepresentation || 'ballAndStick';
  if (mineralMolviewBgColorSelect) mineralMolviewBgColorSelect.value = settings.mineralMolviewBgColor || 'black';
  if (mineralCrystallographySelect) mineralCrystallographySelect.value = settings.mineralCrystallography || 'supercell_2x2x2';

  // Always show Compound Options
  if (compoundOptions) {
    compoundOptions.style.display = 'block';
  }

  // Load MolView-Only Mode option
  if (molSearchModeToggle) molSearchModeToggle.checked = settings.molSearchMode;

  // Load Disable Formula Fallback option
  if (disableFormulaFallbackToggle) disableFormulaFallbackToggle.checked = settings.disableFormulaFallback;

  // Load per-style settings for the current style
  const currentStyle = settings.viewer3DStyle || 'stick:sphere';
  const styleSettings = settings.viewer3DStyleSettings || {};
  const currentStyleSettings = styleSettings[currentStyle] || {};

  if (viewer3DStickRadiusSelect) {
    viewer3DStickRadiusSelect.value = currentStyleSettings.stickRadius || '0.15';
  }
  if (viewer3DSphereRadiusSelect) {
    viewer3DSphereRadiusSelect.value = currentStyleSettings.sphereRadius || '0.3';
  }

  // Set rendering engine radio button
  const engineRadios = document.querySelectorAll('input[name="renderingEngine"]');
  engineRadios.forEach(radio => {
    radio.checked = (radio.value === settings.rendererEngine);
  });

  // Show SmilesDrawer options (client-side is the only engine now)
  if (smilesDrawerOptions) {
    smilesDrawerOptions.style.display = 'block';
  }
  // Show Compound Options for ALL engines
  const compoundOptionsInit = document.getElementById('compoundOptions');
  if (compoundOptionsInit) {
    compoundOptionsInit.style.display = 'block';
  }
  // Show Biomolecule Options for ALL engines (Search API detects biomolecules regardless of renderer)
  const biomoleculeOptionsInit = document.getElementById('biomoleculeOptions');
  if (biomoleculeOptionsInit) {
    biomoleculeOptionsInit.style.display = 'block';
  }
  // Show Mineral Options for ALL engines
  const mineralOptionsInit = document.getElementById('mineralOptions');
  if (mineralOptionsInit) {
    mineralOptionsInit.style.display = 'block';
  }
});

// Extension enabled toggle - CRITICAL: was missing, causing settings to not persist
if (enabledToggle) {
  enabledToggle.addEventListener('change', (e) => {
    chrome.storage.sync.set({ enabled: e.target.checked }, () => {
      showStatus('ChemTex ' + (e.target.checked ? 'enabled' : 'disabled'), 'success');
    });
  });
}

if (mhchemToggle) {
  mhchemToggle.addEventListener('change', (e) => {
    chrome.storage.sync.set({ renderMhchem: e.target.checked }, () => {
      showStatus('mhchem rendering ' + (e.target.checked ? 'enabled' : 'disabled') + '.', 'success');
    });
  });
}

if (chemfigToggle) {
  chemfigToggle.addEventListener('change', (e) => {
    chrome.storage.sync.set({ renderChemfig: e.target.checked }, () => {
      showStatus('chemfig rendering ' + (e.target.checked ? 'enabled' : 'disabled') + '.', 'success');
    });
  });
}

if (perfModeToggle) {
  perfModeToggle.addEventListener('change', (e) => {
    chrome.storage.sync.set({ performanceMode: e.target.checked }, () => {
      broadcastSettingsChange({ performanceMode: e.target.checked });
      showStatus('Lazy-loading ' + (e.target.checked ? 'enabled' : 'disabled'), 'success');
    });
  });
}

// Dev mode toggle
if (devModeToggle) {
  devModeToggle.addEventListener('change', (e) => {
    chrome.storage.sync.set({ devMode: e.target.checked }, () => {
      showStatus('Developer mode ' + (e.target.checked ? 'enabled' : 'disabled') + '. Show raw chemfig text ' + (e.target.checked ? 'ON' : 'OFF') + '. Reload page to apply.', 'success');
    });
  });
}

// Show Tags toggle - always show compound name labels
if (showTagsToggle) {
  showTagsToggle.addEventListener('change', (e) => {
    const value = e.target.checked;
    chrome.storage.sync.set({ showTags: value }, () => {
      broadcastSettingsChange({ showTags: value });
      showStatus('Tags ' + (value ? 'always visible' : 'only on hover'), 'success');
    });
  });
}

// Reload All button - re-render all molecules without page reload
if (reloadAllBtn) {
  reloadAllBtn.addEventListener('click', () => {
    chrome.runtime.sendMessage({
      type: 'RELOAD_ALL_IMAGES'
    });
    showStatus('Reloading all images...', 'success');
  });
}

// Clear Cache button - clear all cached SMILES lookups
const clearCacheBtn = document.getElementById('clearCacheBtn');
if (clearCacheBtn) {
  clearCacheBtn.addEventListener('click', () => {
    // Send message to content script to clear cache
    chrome.runtime.sendMessage({
      type: 'CLEAR_CACHE'
    });

    // Also clear any local storage cache
    chrome.storage.local.get(null, (items) => {
      const keysToRemove = Object.keys(items).filter(key =>
        key.startsWith('smiles_') ||
        key.startsWith('cid_') ||
        key.startsWith('cache_') ||
        key.startsWith('search_') ||
        key === 'chemtex_smiles_cache'  // CRITICAL: The persistent SMILES cache
      );
      if (keysToRemove.length > 0) {
        chrome.storage.local.remove(keysToRemove, () => {
          console.log('[Popup] Cleared', keysToRemove.length, 'cached items from storage:', keysToRemove);
        });
      } else {
        // Even if no prefix matches, ALWAYS try to clear the main cache
        chrome.storage.local.remove('chemtex_smiles_cache', () => {
          console.log('[Popup] Cleared chemtex_smiles_cache');
        });
      }
    });

    showStatus('Cache cleared! Reload images to fetch fresh data.', 'success');
  });
}

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

// Helper function to broadcast settings change to all tabs
function broadcastSettingsChange(changedSettings) {
  chrome.runtime.sendMessage({
    type: 'SETTINGS_CHANGED',
    settings: changedSettings
  });
}

// SmilesDrawer rendering options
if (sdShowCarbonsToggle) {
  sdShowCarbonsToggle.addEventListener('change', (e) => {
    const value = e.target.checked;
    console.log('[Popup] Setting sdShowCarbons to:', value);
    chrome.storage.sync.set({ sdShowCarbons: value }, () => {
      if (chrome.runtime.lastError) {
        console.error('[Popup] Error saving sdShowCarbons:', chrome.runtime.lastError);
        showStatus('Error saving setting: ' + chrome.runtime.lastError.message, 'error');
        return;
      }
      console.log('[Popup] Successfully saved sdShowCarbons:', value);
      broadcastSettingsChange({ sdShowCarbons: value });
      showStatus('Show carbons ' + (value ? 'enabled' : 'disabled'), 'success');
    });
  });
}

if (sdAromaticRingsToggle) {
  sdAromaticRingsToggle.addEventListener('change', (e) => {
    const value = e.target.checked;
    console.log('[Popup] Setting sdAromaticRings to:', value);
    chrome.storage.sync.set({ sdAromaticRings: value }, () => {
      if (chrome.runtime.lastError) {
        console.error('[Popup] Error saving sdAromaticRings:', chrome.runtime.lastError);
        showStatus('Error saving setting: ' + chrome.runtime.lastError.message, 'error');
        return;
      }
      console.log('[Popup] Successfully saved sdAromaticRings:', value);
      broadcastSettingsChange({ sdAromaticRings: value });
      showStatus('Aromatic ring circles ' + (value ? 'enabled' : 'disabled'), 'success');
    });
  });
}

if (sdShowMethylsToggle) {
  sdShowMethylsToggle.addEventListener('change', (e) => {
    const value = e.target.checked;
    saveSettingWithVerification('sdShowMethyls', value, (success, error) => {
      if (success) {
        broadcastSettingsChange({ sdShowMethyls: value });
        showStatus('Show methyls ' + (value ? 'enabled' : 'disabled'), 'success');
      } else {
        showStatus('Error saving setting: ' + (error || 'Unknown error'), 'error');
        // Revert the toggle to its previous state
        sdShowMethylsToggle.checked = !value;
      }
    });
  });
}

if (sdAtomNumbersToggle) {
  sdAtomNumbersToggle.addEventListener('change', (e) => {
    const value = e.target.checked;
    saveSettingWithVerification('sdAtomNumbers', value, (success, error) => {
      if (success) {
        broadcastSettingsChange({ sdAtomNumbers: value });
        showStatus('Atom numbers ' + (value ? 'enabled' : 'disabled'), 'success');
      } else {
        showStatus('Error saving setting: ' + (error || 'Unknown error'), 'error');
        sdAtomNumbersToggle.checked = !value;
      }
    });
  });
}

if (sdShowExplicitHydrogensToggle) {
  sdShowExplicitHydrogensToggle.addEventListener('change', (e) => {
    const value = e.target.checked;
    saveSettingWithVerification('sdShowExplicitHydrogens', value, (success, error) => {
      if (success) {
        broadcastSettingsChange({ sdShowExplicitHydrogens: value });
        showStatus('Explicit hydrogens ' + (value ? 'enabled' : 'disabled'), 'success');
      } else {
        showStatus('Error saving setting: ' + (error || 'Unknown error'), 'error');
        sdShowExplicitHydrogensToggle.checked = !value;
      }
    });
  });
}

if (sdShowImplicitHydrogensToggle) {
  sdShowImplicitHydrogensToggle.addEventListener('change', (e) => {
    const value = e.target.checked;
    saveSettingWithVerification('sdShowImplicitHydrogens', value, (success, error) => {
      if (success) {
        broadcastSettingsChange({ sdShowImplicitHydrogens: value });
        showStatus('Implicit hydrogens ' + (value ? 'shown' : 'hidden'), 'success');
      } else {
        showStatus('Error saving setting: ' + (error || 'Unknown error'), 'error');
        sdShowImplicitHydrogensToggle.checked = !value;
      }
    });
  });
}

if (sdCompactDrawingToggle) {
  sdCompactDrawingToggle.addEventListener('change', (e) => {
    const value = e.target.checked;
    saveSettingWithVerification('sdCompactDrawing', value, (success, error) => {
      if (success) {
        broadcastSettingsChange({ sdCompactDrawing: value });
        showStatus('Compact drawing ' + (value ? 'enabled' : 'disabled'), 'success');
      } else {
        showStatus('Error saving setting: ' + (error || 'Unknown error'), 'error');
        sdCompactDrawingToggle.checked = !value;
      }
    });
  });
}

if (sdFlipHorizontalToggle) {
  sdFlipHorizontalToggle.addEventListener('change', (e) => {
    const value = e.target.checked;
    saveSettingWithVerification('sdFlipHorizontal', value, (success, error) => {
      if (success) {
        broadcastSettingsChange({ sdFlipHorizontal: value });
        showStatus('Horizontal flip ' + (value ? 'enabled' : 'disabled'), 'success');
      } else {
        showStatus('Error saving setting: ' + (error || 'Unknown error'), 'error');
        sdFlipHorizontalToggle.checked = !value;
      }
    });
  });
}

if (sdFlipVerticalToggle) {
  sdFlipVerticalToggle.addEventListener('change', (e) => {
    const value = e.target.checked;
    saveSettingWithVerification('sdFlipVertical', value, (success, error) => {
      if (success) {
        broadcastSettingsChange({ sdFlipVertical: value });
        showStatus('Vertical flip ' + (value ? 'enabled' : 'disabled'), 'success');
      } else {
        showStatus('Error saving setting: ' + (error || 'Unknown error'), 'error');
        sdFlipVerticalToggle.checked = !value;
      }
    });
  });
}

if (sdThemeSelect) {
  sdThemeSelect.addEventListener('change', (e) => {
    const value = e.target.value;
    chrome.storage.sync.set({ sdTheme: value }, () => {
      broadcastSettingsChange({ sdTheme: value });
      showStatus('Theme set to ' + value, 'success');
    });
  });
}

if (sdRotateSlider) {
  sdRotateSlider.addEventListener('input', (e) => {
    sdRotateValue.textContent = e.target.value + '°';
  });

  sdRotateSlider.addEventListener('change', (e) => {
    const value = parseInt(e.target.value);
    chrome.storage.sync.set({ sdRotate: value }, () => {
      broadcastSettingsChange({ sdRotate: value });
      showStatus('Rotation set to ' + value + '°', 'success');
    });
  });
}

if (sdBondThicknessSlider) {
  sdBondThicknessSlider.addEventListener('input', (e) => {
    sdBondThicknessValue.textContent = parseFloat(e.target.value).toFixed(1);
  });

  sdBondThicknessSlider.addEventListener('change', (e) => {
    const value = parseFloat(e.target.value);
    chrome.storage.sync.set({ sdBondThickness: value }, () => {
      broadcastSettingsChange({ sdBondThickness: value });
      showStatus('Bond thickness set to ' + value, 'success');
    });
  });
}

if (sdGradientColorsToggle) {
  sdGradientColorsToggle.addEventListener('change', (e) => {
    const value = e.target.checked;
    chrome.storage.sync.set({ sdGradientColors: value }, () => {
      broadcastSettingsChange({ sdGradientColors: value });
      showStatus('Gradient colors ' + (value ? 'enabled' : 'disabled'), 'success');
    });
  });
}

if (sdScaleByWeightToggle) {
  sdScaleByWeightToggle.addEventListener('change', (e) => {
    chrome.storage.sync.set({ sdScaleByWeight: e.target.checked }, () => {
      showStatus('Scale by weight ' + (e.target.checked ? 'enabled' : 'disabled') + '. Reload page to apply.', 'success');
    });
  });
}

// Data Source toggle event listeners
if (searchPubChemToggle) {
  searchPubChemToggle.addEventListener('change', (e) => {
    chrome.storage.sync.set({ searchPubChem: e.target.checked }, () => {
      showStatus('PubChem search ' + (e.target.checked ? 'enabled' : 'disabled') + '. Reload page to apply.', 'success');
    });
  });
}

if (searchRCSBToggle) {
  searchRCSBToggle.addEventListener('change', (e) => {
    chrome.storage.sync.set({ searchRCSB: e.target.checked }, () => {
      showStatus('RCSB search ' + (e.target.checked ? 'enabled' : 'disabled') + '. Reload page to apply.', 'success');
    });
  });
}

if (searchCODToggle) {
  searchCODToggle.addEventListener('change', (e) => {
    chrome.storage.sync.set({ searchCOD: e.target.checked }, () => {
      showStatus('COD search ' + (e.target.checked ? 'enabled' : 'disabled') + '. Reload page to apply.', 'success');
    });
  });
}

// Stereochemistry toggle
const useStereochemistryToggleHandler = document.getElementById('useStereochemistryToggle');
if (useStereochemistryToggleHandler) {
  useStereochemistryToggleHandler.addEventListener('change', (e) => {
    chrome.storage.sync.set({ useStereochemistry: e.target.checked }, () => {
      showStatus('Stereochemistry ' + (e.target.checked ? 'enabled' : 'disabled'), 'success');
    });
  });
}

// MolView specific options event listeners (compound)
if (molviewRepresentationSelect) {
  molviewRepresentationSelect.addEventListener('change', (e) => {
    chrome.storage.sync.set({ molviewRepresentation: e.target.value }, () => {
      broadcastSettingsChange({ molviewRepresentation: e.target.value });
      showStatus('Representation set to ' + e.target.options[e.target.selectedIndex].text, 'success');
    });
  });
}

if (compoundMolviewBgColorSelect) {
  compoundMolviewBgColorSelect.addEventListener('change', (e) => {
    const value = e.target.value;
    console.log('[Popup] Attempting to save compoundMolviewBgColor:', value);
    saveSettingWithVerification('compoundMolviewBgColor', value, (success, error) => {
      if (success) {
        broadcastSettingsChange({ compoundMolviewBgColor: value });
        showStatus('Compound background set to ' + e.target.options[e.target.selectedIndex].text, 'success');
      } else {
        showStatus('Error saving background color: ' + (error || 'Unknown error'), 'error');
        // Don't revert the dropdown - let user see the issue
      }
    });
  });
}

if (molviewEngineSelect) {
  molviewEngineSelect.addEventListener('change', (e) => {
    chrome.storage.sync.set({ molviewEngine: e.target.value }, () => {
      broadcastSettingsChange({ molviewEngine: e.target.value });
      showStatus('Engine set to ' + e.target.options[e.target.selectedIndex].text, 'success');
    });
  });
}



// Biomolecule Options event listeners - ALL broadcast changes for instant updates
if (proteinRemoveWhiteBgToggle) {
  proteinRemoveWhiteBgToggle.addEventListener('change', (e) => {
    chrome.storage.sync.set({ proteinRemoveWhiteBg: e.target.checked }, () => {
      broadcastSettingsChange({ proteinRemoveWhiteBg: e.target.checked });
      showStatus('White background removal ' + (e.target.checked ? 'enabled' : 'disabled'), 'success');
    });
  });
}

if (molviewBioAssemblyToggle) {
  molviewBioAssemblyToggle.addEventListener('change', (e) => {
    chrome.storage.sync.set({ molviewBioAssembly: e.target.checked }, () => {
      broadcastSettingsChange({ molviewBioAssembly: e.target.checked });
      showStatus('Bio Assembly ' + (e.target.checked ? 'shown' : 'hidden'), 'success');
    });
  });
}

if (molviewChainTypeSelect) {
  molviewChainTypeSelect.addEventListener('change', (e) => {
    chrome.storage.sync.set({ molviewChainType: e.target.value }, () => {
      broadcastSettingsChange({ molviewChainType: e.target.value });
      showStatus('Chain Type set to ' + e.target.options[e.target.selectedIndex].text, 'success');
    });
  });
}

if (molviewChainBondsToggle) {
  molviewChainBondsToggle.addEventListener('change', (e) => {
    chrome.storage.sync.set({ molviewChainBonds: e.target.checked }, () => {
      broadcastSettingsChange({ molviewChainBonds: e.target.checked });
      showStatus('Chain Bonds ' + (e.target.checked ? 'shown' : 'hidden'), 'success');
    });
  });
}

if (molviewChainColorSelect) {
  molviewChainColorSelect.addEventListener('change', (e) => {
    chrome.storage.sync.set({ molviewChainColor: e.target.value }, () => {
      broadcastSettingsChange({ molviewChainColor: e.target.value });
      showStatus('Chain Color set to ' + e.target.options[e.target.selectedIndex].text, 'success');
    });
  });
}

if (proteinMolviewBgColorSelect) {
  proteinMolviewBgColorSelect.addEventListener('change', (e) => {
    chrome.storage.sync.set({ proteinMolviewBgColor: e.target.value }, () => {
      broadcastSettingsChange({ proteinMolviewBgColor: e.target.value });
      showStatus('Protein background color set to ' + e.target.options[e.target.selectedIndex].text, 'success');
    });
  });
}

// Mineral Options event listeners
if (mineralRepresentationSelect) {
  mineralRepresentationSelect.addEventListener('change', (e) => {
    chrome.storage.sync.set({ mineralRepresentation: e.target.value }, () => {
      broadcastSettingsChange({ mineralRepresentation: e.target.value });
      showStatus('Mineral representation set to ' + e.target.options[e.target.selectedIndex].text, 'success');
    });
  });
}

if (mineralMolviewBgColorSelect) {
  mineralMolviewBgColorSelect.addEventListener('change', (e) => {
    chrome.storage.sync.set({ mineralMolviewBgColor: e.target.value }, () => {
      broadcastSettingsChange({ mineralMolviewBgColor: e.target.value });
      showStatus('Mineral background color set to ' + e.target.options[e.target.selectedIndex].text, 'success');
    });
  });
}

if (mineralCrystallographySelect) {
  mineralCrystallographySelect.addEventListener('change', (e) => {
    chrome.storage.sync.set({ mineralCrystallography: e.target.value }, () => {
      showStatus('Mineral crystallography set to ' + e.target.options[e.target.selectedIndex].text + '. Reload page to apply.', 'success');
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
        const engineNames = {
          'moleculeviewer': 'MoleculeViewer',
          'mol2chemfig': 'mol2chemfig',
          'pubchem': 'PubChem',
          'client-side': 'Client-Side'
        };
        showStatus(`Switched to ${engineNames[engine]}. Reload page to apply.`, 'success');
      });
    }
  });
});

/**
 * Update engine info display
 */
const biomoleculeOptions = document.getElementById('biomoleculeOptions');
const mineralOptions = document.getElementById('mineralOptions');

function updateEngineInfo(engine) {
  if (engine === 'cdk-depict') {
    if (engineInfo) engineInfo.textContent = '🌐 CDK Depict - External API';
    if (smilesDrawerOptions) smilesDrawerOptions.style.display = 'none';
    if (cdkDepictOptions) cdkDepictOptions.style.display = 'block';
  } else {
    // Default: SmilesDrawer (client-side)
    if (engineInfo) engineInfo.textContent = '💻 SmilesDrawer - Client-Side (No Server!)';
    if (smilesDrawerOptions) smilesDrawerOptions.style.display = 'block';
    if (cdkDepictOptions) cdkDepictOptions.style.display = 'none';
  }

  // Always show Compound, Biomolecule, and Mineral Options (independent of renderer)
  if (compoundOptions) compoundOptions.style.display = 'block';
  if (biomoleculeOptions) biomoleculeOptions.style.display = 'block';
  if (mineralOptions) mineralOptions.style.display = 'block';
}

/**
 * Show status message
 */
function showStatus(message, type) {
  statusDiv.textContent = message;
  statusDiv.className = 'status show';

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

// ==========================================
// UI SETTINGS - Glassmorphism Controls
// ==========================================

const blurSlider = document.getElementById('blurSlider');
const opacitySlider = document.getElementById('opacitySlider');
const blurValue = document.getElementById('blurValue');
const opacityValue = document.getElementById('opacityValue');

// Apply glassmorphism styles to all sections
function applyGlassmorphism(blur, opacity) {
  var sections = document.querySelectorAll('.section');
  var header = document.querySelector('.header');
  var footer = document.querySelector('.footer');
  var radioOptions = document.querySelectorAll('.radio-option');

  var blurPx = blur + 'px';
  var bgOpacity = opacity / 100;

  sections.forEach(function (section) {
    section.style.backdropFilter = 'blur(' + blurPx + ')';
    section.style.webkitBackdropFilter = 'blur(' + blurPx + ')';
    section.style.background = 'rgba(255, 255, 255, ' + bgOpacity + ')';
  });

  if (header) {
    header.style.backdropFilter = 'blur(' + blurPx + ')';
    header.style.webkitBackdropFilter = 'blur(' + blurPx + ')';
    header.style.background = 'rgba(255, 255, 255, ' + (bgOpacity * 0.9) + ')';
  }

  if (footer) {
    footer.style.backdropFilter = 'blur(' + blurPx + ')';
    footer.style.webkitBackdropFilter = 'blur(' + blurPx + ')';
    footer.style.background = 'rgba(255, 255, 255, ' + (bgOpacity * 0.9) + ')';
  }

  radioOptions.forEach(function (opt) {
    opt.style.background = 'rgba(255, 255, 255, ' + (bgOpacity * 0.8) + ')';
  });
}

// Load UI settings
chrome.storage.sync.get({
  uiBlur: 7,
  uiOpacity: 23
}, function (settings) {
  if (blurSlider) {
    blurSlider.value = settings.uiBlur;
    blurValue.textContent = settings.uiBlur + 'px';
  }
  if (opacitySlider) {
    opacitySlider.value = settings.uiOpacity;
    opacityValue.textContent = settings.uiOpacity + '%';
  }
  applyGlassmorphism(settings.uiBlur, settings.uiOpacity);
});

// Blur slider handler
if (blurSlider) {
  blurSlider.addEventListener('input', function (e) {
    var val = parseInt(e.target.value);
    blurValue.textContent = val + 'px';
    var opacity = parseInt(opacitySlider.value);
    applyGlassmorphism(val, opacity);
  });

  blurSlider.addEventListener('change', function (e) {
    chrome.storage.sync.set({ uiBlur: parseInt(e.target.value) });
  });
}

// Opacity slider handler
if (opacitySlider) {
  opacitySlider.addEventListener('input', function (e) {
    var val = parseInt(e.target.value);
    opacityValue.textContent = val + '%';
    var blur = parseInt(blurSlider.value);
    applyGlassmorphism(blur, val);
  });

  opacitySlider.addEventListener('change', function (e) {
    chrome.storage.sync.set({ uiOpacity: parseInt(e.target.value) });
  });
}

// Molecule Count slider handler
const molCountSlider = document.getElementById('molCountSlider');
const molCountValue = document.getElementById('molCountValue');

if (molCountSlider) {
  molCountSlider.addEventListener('input', function (e) {
    molCountValue.textContent = e.target.value;
  });

  molCountSlider.addEventListener('change', function (e) {
    chrome.storage.sync.set({ moleculeCount: parseInt(e.target.value) }, () => {
      // Trigger animation update if possible, or just let the user reload
      // For now, we'll rely on the animation script listening to storage changes or reload
      if (window.updateMoleculeAnimation) {
        window.updateMoleculeAnimation();
      }
    });
  });
}

// Molecule Scale slider handler
const molScaleSlider = document.getElementById('molScaleSlider');
const molScaleValue = document.getElementById('molScaleValue');

if (molScaleSlider) {
  molScaleSlider.addEventListener('input', function (e) {
    molScaleValue.textContent = e.target.value;
  });

  molScaleSlider.addEventListener('change', function (e) {
    chrome.storage.sync.set({ moleculeScale: parseFloat(e.target.value) }, () => {
      if (window.updateMoleculeAnimation) {
        window.updateMoleculeAnimation();
      }
    });
  });
}

// Load Molecule Settings
chrome.storage.sync.get({
  moleculeCount: 22,
  moleculeScale: 0.4,
  cursorGravityStrength: -0.9,
  initialVelocity: 10
}, function (settings) {
  if (molCountSlider) {
    molCountSlider.value = settings.moleculeCount;
    molCountValue.textContent = settings.moleculeCount;
  }
  if (molScaleSlider) {
    molScaleSlider.value = settings.moleculeScale;
    molScaleValue.textContent = settings.moleculeScale;
  }

  // New physics UI controls
  if (cursorGravitySlider) {
    cursorGravitySlider.value = settings.cursorGravityStrength;
    cursorGravityValue.textContent = settings.cursorGravityStrength;
  }
  if (initialVelocitySlider) {
    initialVelocitySlider.value = settings.initialVelocity;
    initialVelocityValue.textContent = settings.initialVelocity;
  }
});

// =============================================
// Popup Physics Controls (cursor gravity / initial velocity)
// =============================================

const cursorGravitySlider = document.getElementById('cursorGravitySlider');
const cursorGravityValue = document.getElementById('cursorGravityValue');
const initialVelocitySlider = document.getElementById('initialVelocitySlider');
const initialVelocityValue = document.getElementById('initialVelocityValue');

function triggerPopupPhysicsUpdate() {
  try {
    if (window.updateMoleculeAnimation) {
      window.updateMoleculeAnimation();
    }
  } catch (e) {
    // ignore
  }
}

if (cursorGravitySlider) {
  cursorGravitySlider.addEventListener('input', function (e) {
    if (cursorGravityValue) cursorGravityValue.textContent = e.target.value;
  });
  cursorGravitySlider.addEventListener('change', function (e) {
    chrome.storage.sync.set({ cursorGravityStrength: parseFloat(e.target.value) }, () => {
      triggerPopupPhysicsUpdate();
    });
  });
}

if (initialVelocitySlider) {
  initialVelocitySlider.addEventListener('input', function (e) {
    if (initialVelocityValue) initialVelocityValue.textContent = e.target.value;
  });
  initialVelocitySlider.addEventListener('change', function (e) {
    chrome.storage.sync.set({ initialVelocity: parseFloat(e.target.value) }, () => {
      triggerPopupPhysicsUpdate();
    });
  });
}

// =============================================
// Collapsible Sections
// =============================================
(function initCollapsibleSections() {
  const sectionTitles = document.querySelectorAll('.section-title');

  // Load collapsed state from storage
  chrome.storage.sync.get({ collapsedSections: [] }, function (data) {
    const collapsed = data.collapsedSections || [];

    sectionTitles.forEach(function (title) {
      const sectionName = title.querySelector('span')?.textContent?.trim();
      const section = title.closest('.section');

      // Restore collapsed state
      if (sectionName && collapsed.includes(sectionName)) {
        section.classList.add('collapsed');
      }

      // Add click handler
      title.addEventListener('click', function () {
        section.classList.toggle('collapsed');

        // Save state
        const allCollapsed = [];
        document.querySelectorAll('.section.collapsed .section-title span').forEach(function (span) {
          allCollapsed.push(span.textContent.trim());
        });
        chrome.storage.sync.set({ collapsedSections: allCollapsed });
      });
    });
  });
})();

// MetaMask Donate Button Handler
const metamaskDonateBtn = document.getElementById('metamaskDonate');
if (metamaskDonateBtn) {
  metamaskDonateBtn.addEventListener('click', async (e) => {
    e.preventDefault();

    // Your Ethereum wallet address for donations
    const donationAddress = '0xYOUR_WALLET_ADDRESS_HERE'; // Replace with your actual wallet

    // Check if MetaMask is installed
    if (typeof window.ethereum !== 'undefined') {
      try {
        // Request account access
        const accounts = await window.ethereum.request({ method: 'eth_requestAccounts' });
        const fromAddress = accounts[0];

        // Send transaction (0.001 ETH default donation)
        const txHash = await window.ethereum.request({
          method: 'eth_sendTransaction',
          params: [{
            from: fromAddress,
            to: donationAddress,
            value: '0x38D7EA4C68000' // 0.001 ETH in hex
          }]
        });

        alert('Thank you for your donation! 🎉\nTransaction: ' + txHash);
      } catch (error) {
        if (error.code === 4001) {
          // User rejected
          console.log('User cancelled donation');
        } else {
          alert('Donation failed: ' + error.message);
        }
      }
    } else {
      // No MetaMask - show address to copy
      const copyAddress = confirm(
        'MetaMask not detected!\n\n' +
        'Would you like to copy the donation address?\n\n' +
        donationAddress
      );
      if (copyAddress) {
        navigator.clipboard.writeText(donationAddress);
        alert('Address copied to clipboard! 📋');
      }
    }
  });
}

// Icon button listeners (same functions as main buttons)
const clearCacheIconBtn = document.getElementById('clearCacheIconBtn');
if (clearCacheIconBtn) {
  clearCacheIconBtn.addEventListener('click', () => {
    chrome.runtime.sendMessage({ type: 'CLEAR_CACHE' });
    chrome.storage.local.get(null, (items) => {
      const keysToRemove = Object.keys(items).filter(key =>
        key.startsWith('smiles_') || key.startsWith('cid_') ||
        key.startsWith('cache_') || key.startsWith('search_') ||
        key === 'chemtex_smiles_cache'
      );
      if (keysToRemove.length > 0) {
        chrome.storage.local.remove(keysToRemove);
      } else {
        chrome.storage.local.remove('chemtex_smiles_cache');
      }
    });
    showStatus('Cache cleared!', 'success');
  });
}

const reloadAllIconBtn = document.getElementById('reloadAllIconBtn');
if (reloadAllIconBtn) {
  reloadAllIconBtn.addEventListener('click', () => {
    chrome.runtime.sendMessage({ type: 'RELOAD_ALL_IMAGES' });
    showStatus('Reloading...', 'success');
  });
}
