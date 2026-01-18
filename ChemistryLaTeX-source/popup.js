/**
 * Popup script for ChemistryLaTeX v6.0
 */

// DOM elements - with null checks
const enabledToggle = document.getElementById('enabledToggle');
const perfModeToggle = document.getElementById('perfModeToggle');
const devModeToggle = document.getElementById('devModeToggle');
const reloadAllBtn = document.getElementById('reloadAllBtn');
const statusDiv = document.getElementById('status');


// Rendering options
const renderingOptionsSection = document.getElementById('renderingOptions');
const sdUseStereochemistryToggle = document.getElementById('sdUseStereochemistryToggle');
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
// sdAutoAdaptToggle REMOVED - no auto-adapt anymore
const sdRotateSlider = document.getElementById('sdRotateSlider');
const sdRotateValue = document.getElementById('sdRotateValue');
const sdAverageSizeSlider = document.getElementById('sdAverageSizeSlider');
const sdAverageSizeValue = document.getElementById('sdAverageSizeValue');
const sdGradientColorsToggle = document.getElementById('sdGradientColorsToggle');
const sdScaleByWeightToggle = document.getElementById('sdScaleByWeightToggle');

// Data Source options
const renderingEngineSelect = document.getElementById('renderingEngineSelect');

// Compound Options
const compoundOptions = document.getElementById('compoundOptions');
const molviewRepresentationSelect = document.getElementById('molviewRepresentationSelect');
const compoundMolviewBgColorSelect = document.getElementById('compoundMolviewBgColorSelect');
const molviewEngineSelect = document.getElementById('molviewEngineSelect');
const compound3DSizeSlider = document.getElementById('compound3DSizeSlider');
const compound3DSizeValue = document.getElementById('compound3DSizeValue');

// DEBUG: Verify compound option elements exist
console.log('%c[Popup] 🔍 DOM Element Check:', 'color: #9C27B0; font-weight: bold;', {
  compoundOptions: !!compoundOptions,
  molviewRepresentationSelect: !!molviewRepresentationSelect,
  compoundMolviewBgColorSelect: !!compoundMolviewBgColorSelect,
  molviewEngineSelect: !!molviewEngineSelect,
  compound3DSizeSlider: !!compound3DSizeSlider
});

// Biomolecule Options
const proteinRemoveWhiteBgToggle = document.getElementById('proteinRemoveWhiteBgToggle');
const molviewBioAssemblyToggle = document.getElementById('molviewBioAssemblyToggle');
const molviewChainTypeSelect = document.getElementById('molviewChainTypeSelect');
const molviewChainBondsToggle = document.getElementById('molviewChainBondsToggle');
const molviewChainColorSelect = document.getElementById('molviewChainColorSelect');
const proteinMolviewBgColorSelect = document.getElementById('proteinMolviewBgColorSelect');

// Bio Assembly Settings
const bioAssemblySettings = document.getElementById('bioAssemblySettings');
const bioAssemblyViewerSelect = document.getElementById('bioAssemblyViewerSelect');
const bioAssemblyGraphicsSelect = document.getElementById('bioAssemblyGraphicsSelect');
const bioAssemblyPdbProviderSelect = document.getElementById('bioAssemblyPdbProviderSelect');

// Mineral Options
const mineralRepresentationSelect = document.getElementById('mineralRepresentationSelect');
const mineralMolviewBgColorSelect = document.getElementById('mineralMolviewBgColorSelect');
const mineralCrystallographySelect = document.getElementById('mineralCrystallographySelect');
const mineral3DSizeSlider = document.getElementById('mineral3DSizeSlider');
const mineral3DSizeValue = document.getElementById('mineral3DSizeValue');

// Image size control options
const saveSizePerImageToggle = document.getElementById('saveSizePerImageToggle');
const saveSizeBySMILESToggle = document.getElementById('saveSizeBySMILESToggle');

// 3D Viewer options
const enable3DViewerToggle = document.getElementById('enable3DViewerToggle');
const default3DViewSelect = document.getElementById('default3DViewSelect');
const viewer3DSourceSelect = document.getElementById('viewer3DSourceSelect');
const viewer3DStyleSelect = document.getElementById('viewer3DStyleSelect');
const viewer3DAutoRotateToggle = document.getElementById('viewer3DAutoRotateToggle');
const viewer3DSizeSelect = document.getElementById('viewer3DSizeSelect');
const viewer3DBgColorSelect = document.getElementById('viewer3DBgColorSelect');

// MolView-Only Mode toggle
const molSearchModeToggle = document.getElementById('molSearchModeToggle');

// Disable Formula Fallback toggle
const disableFormulaFallbackToggle = document.getElementById('disableFormulaFallbackToggle');

// Enable AI Flag Control toggle
const enableAIFlagControlToggle = document.getElementById('enableAIFlagControlToggle');

// Search toggles (may not exist in HTML but referenced in code)
const searchCompoundsToggle = document.getElementById('searchCompoundsToggle');
const searchBiomoleculesToggle = document.getElementById('searchBiomoleculesToggle');
const searchMineralsToggle = document.getElementById('searchMineralsToggle');

// Safe selector function
function safeGetElement(id) {
  const el = document.getElementById(id);
  if (!el) console.warn(`Element ${id} not found in DOM`);
  return el;
}

// Simple setting save with error checking (used by most settings)
// Also tracks which settings the user has explicitly modified
function saveSetting(settingObj, successMessage, callback) {
  chrome.storage.sync.set(settingObj, () => {
    if (chrome.runtime.lastError) {
      const errorMsg = chrome.runtime.lastError.message || 'Unknown error';
      console.error('[Popup] ❌ Error saving settings:', errorMsg, settingObj);
      showStatus('Error saving: ' + errorMsg, 'error');
      if (callback) callback(false);
      return;
    }
    console.log('[Popup] ✅ Saved:', settingObj);
    if (successMessage) showStatus(successMessage, 'success');

    // Track that the user has explicitly modified these settings
    // Server defaults won't override user-modified settings
    trackUserModifiedSettings(Object.keys(settingObj));

    if (callback) callback(true);
  });
}

// Helper function to save settings with verification (for critical settings)
// Also tracks which settings the user has explicitly modified
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

        // Track that the user has explicitly modified this setting
        trackUserModifiedSettings([settingName]);

        if (callback) callback(true);
      } else {
        console.error(`[Popup] Save verification failed for ${settingName}. Expected:`, value, 'Got:', savedValue);
        if (callback) callback(false, 'Verification failed');
      }
    });
  });
}

// Track which settings the user has explicitly modified
// This prevents server defaults from overriding user choices
function trackUserModifiedSettings(settingNames) {
  chrome.storage.sync.get(['userModifiedSettings'], (data) => {
    const modified = data.userModifiedSettings || {};
    settingNames.forEach(name => {
      modified[name] = true;
    });
    chrome.storage.sync.set({ userModifiedSettings: modified });
  });
}

// Patch notes button click handler
(function setupPatchNotesButton() {
  const patchNotesBtn = document.getElementById('patchNotesBtn');
  const notesModal = document.getElementById('notesModal');
  const notesContent = document.getElementById('notesContent');
  const closeModal = document.getElementById('closeNotesModal');

  if (patchNotesBtn) {
    patchNotesBtn.addEventListener('click', () => {
      if (notesModal && notesContent) {
        // Use server patch notes if available, otherwise show a fallback
        const patchNotes = window._serverPatchNotes || `Version 7.7
- Fixed biomolecule 3D button visibility
- Fixed button positioning for biomolecules  
- Improved rendering performance
- Added notification system for updates

Version 7.6
- Added server-controlled configuration
- Improved SMILES/IUPAC handling
- Fixed re-rendering issues`;

        notesContent.textContent = patchNotes;
        notesModal.style.display = 'flex';

        // Mark as seen if we have a version
        if (window._serverPatchNotesVersion) {
          chrome.storage.local.set({ lastSeenPatchNotesVersion: window._serverPatchNotesVersion });
        }

        // Clear notification indicator
        const notifDot = document.getElementById('patchNotifDot');
        if (notifDot) notifDot.style.display = 'none';

        // Remove animation classes
        patchNotesBtn.classList.remove('important', 'ringing');
      }
    });
  }

  // Close modal handler
  if (closeModal && notesModal) {
    closeModal.addEventListener('click', () => {
      notesModal.style.display = 'none';
    });
    // Also close when clicking outside the modal
    notesModal.addEventListener('click', (e) => {
      if (e.target === notesModal) {
        notesModal.style.display = 'none';
      }
    });
  }
})();


// Load current settings
chrome.storage.sync.get(null, (storedSettings) => {
  // DEBUG: Show what's actually stored in chrome.storage.sync
  console.log('%c[Popup] 📦 RAW STORED SETTINGS:', 'color: #FF6B00; font-weight: bold; font-size: 14px;');
  console.log(JSON.stringify(storedSettings, null, 2));
  console.log('%c[Popup] Key settings check:', 'color: #2196F3; font-weight: bold;', {
    compoundMolviewBgColor: storedSettings.compoundMolviewBgColor,
    molviewBioAssembly: storedSettings.molviewBioAssembly,
    proteinRemoveWhiteBg: storedSettings.proteinRemoveWhiteBg,
    compound3DSize: storedSettings.compound3DSize,
    sdAverageSize: storedSettings.sdAverageSize
  });

  // VISIBLE DEBUG: Show what value was loaded from storage (won't be stripped)
  // Remove this after debugging!
  if (storedSettings.molviewRepresentation) {
    showStatus('DEBUG: Loaded representation = ' + storedSettings.molviewRepresentation, 'success');
  } else {
    showStatus('DEBUG: No representation stored, using default', 'error');
  }

  // Define defaults that will be used ONLY if a setting is undefined
  const defaults = {
    enabled: true,
    performanceMode: true,
    rendererEngine: 'chemistrylatex',
    devMode: false,
    useStereochemistry: false,  // Use isomeric SMILES (show chirality)
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
    sdAutoAdapt: true,
    sdRotate: 0,
    sdBondThickness: 1.0,
    sdAverageSize: 100,  // Percentage scale for rendered molecules
    sdGradientColors: false,
    sdScaleByWeight: false,
    saveSizePerImage: false,
    saveSizeBySMILES: true,
    searchCompounds: true,
    searchBiomolecules: true,
    searchMinerals: true,
    renderingEngine: 'chemistrylatex',
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
    bioAssemblyViewer: 'molstar',        // 'molstar', 'molstar-me', 'molview'
    bioAssemblyGraphics: 'balanced',     // 'quality', 'balanced', 'performance'
    bioAssemblyPdbProvider: 'rcsb',      // 'rcsb', 'pdbe', 'pdbj'
    molviewChainType: 'ribbon',
    molviewChainBonds: false,
    molviewChainColor: 'ss',
    proteinMolviewBgColor: 'black',
    mineralRepresentation: 'ballAndStick',
    mineralMolviewBgColor: 'black',
    mineralCrystallography: 'supercell_2x2x2',
    compound3DSize: 100,   // Percentage scale for compound 3D viewer (keeps aspect ratio)
    mineral3DSize: 100,    // Percentage scale for mineral 3D viewer (keeps aspect ratio)
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
    enableAIFlagControl: false,  // When true, flags in chem:...+d: override settings
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
  console.log('[Popup] Rendering settings:', {
    useStereochemistry: settings.useStereochemistry,
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
  if (perfModeToggle) perfModeToggle.checked = settings.performanceMode;
  if (devModeToggle) devModeToggle.checked = settings.devMode;
  if (enableAIFlagControlToggle) enableAIFlagControlToggle.checked = settings.enableAIFlagControl === true;

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

  // Load rendering options
  if (sdUseStereochemistryToggle) {
    sdUseStereochemistryToggle.checked = settings.useStereochemistry;
    console.log('[Popup] Set sdUseStereochemistryToggle.checked to:', settings.useStereochemistry, '| Actual value now:', sdUseStereochemistryToggle.checked);
  }
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
  // sdAutoAdaptToggle REMOVED - no auto-adapt anymore
  if (sdRotateSlider) {
    sdRotateSlider.value = settings.sdRotate || 0;
    if (sdRotateValue) sdRotateValue.textContent = (settings.sdRotate || 0) + '°';
  }
  if (sdAverageSizeSlider) {
    sdAverageSizeSlider.value = settings.sdAverageSize || 100;
    if (sdAverageSizeValue) sdAverageSizeValue.textContent = (settings.sdAverageSize || 100) + '%';
  }
  if (sdGradientColorsToggle) sdGradientColorsToggle.checked = settings.sdGradientColors || false;
  if (sdScaleByWeightToggle) sdScaleByWeightToggle.checked = settings.sdScaleByWeight || false;

  // Load Data Source toggles
  if (searchCompoundsToggle) searchCompoundsToggle.checked = settings.searchCompounds !== false;
  if (searchBiomoleculesToggle) searchBiomoleculesToggle.checked = settings.searchBiomolecules !== false;
  if (searchMineralsToggle) searchMineralsToggle.checked = settings.searchMinerals !== false;

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

  // Load MolView representation with debug
  if (molviewRepresentationSelect) {
    const repValue = settings.molviewRepresentation || 'ballAndStick';
    molviewRepresentationSelect.value = repValue;
    // Debug: verify it was set (this won't be stripped since it's not console.log)
    if (molviewRepresentationSelect.value !== repValue) {
      showStatus('Warning: Could not set representation to: ' + repValue, 'error');
    }
  }
  if (compoundMolviewBgColorSelect) {
    compoundMolviewBgColorSelect.value = settings.compoundMolviewBgColor || 'black';
  }
  if (molviewEngineSelect) molviewEngineSelect.value = settings.molviewEngine || 'glmol';
  if (compound3DSizeSlider) {
    compound3DSizeSlider.value = settings.compound3DSize || 100;
    if (compound3DSizeValue) compound3DSizeValue.textContent = (settings.compound3DSize || 100) + '%';
  }

  // Load Biomolecule Options
  if (proteinRemoveWhiteBgToggle) proteinRemoveWhiteBgToggle.checked = settings.proteinRemoveWhiteBg === true;
  if (molviewBioAssemblyToggle) molviewBioAssemblyToggle.checked = settings.molviewBioAssembly;

  // Load Bio Assembly sub-settings
  if (bioAssemblyViewerSelect) bioAssemblyViewerSelect.value = settings.bioAssemblyViewer || 'molstar';
  if (bioAssemblyGraphicsSelect) bioAssemblyGraphicsSelect.value = settings.bioAssemblyGraphics || 'balanced';
  if (bioAssemblyPdbProviderSelect) bioAssemblyPdbProviderSelect.value = settings.bioAssemblyPdbProvider || 'rcsb';

  // Enable/disable bio assembly sub-settings based on toggle state
  if (bioAssemblySettings) {
    if (settings.molviewBioAssembly) {
      bioAssemblySettings.style.opacity = '1';
      bioAssemblySettings.style.pointerEvents = 'auto';
    } else {
      bioAssemblySettings.style.opacity = '0.5';
      bioAssemblySettings.style.pointerEvents = 'none';
    }
  }

  // Load chain type with debug
  if (molviewChainTypeSelect) {
    const chainValue = settings.molviewChainType || 'ribbon';
    molviewChainTypeSelect.value = chainValue;
    if (molviewChainTypeSelect.value !== chainValue) {
      showStatus('Warning: Could not set chain type to: ' + chainValue, 'error');
    }
  }
  if (molviewChainBondsToggle) molviewChainBondsToggle.checked = settings.molviewChainBonds;
  if (molviewChainColorSelect) molviewChainColorSelect.value = settings.molviewChainColor || 'ss';
  if (proteinMolviewBgColorSelect) proteinMolviewBgColorSelect.value = settings.proteinMolviewBgColor || 'black';

  // Load Mineral Options
  if (mineralRepresentationSelect) mineralRepresentationSelect.value = settings.mineralRepresentation || 'ballAndStick';
  if (mineralMolviewBgColorSelect) mineralMolviewBgColorSelect.value = settings.mineralMolviewBgColor || 'black';
  if (mineralCrystallographySelect) mineralCrystallographySelect.value = settings.mineralCrystallography || 'supercell_2x2x2';
  if (mineral3DSizeSlider) {
    mineral3DSizeSlider.value = settings.mineral3DSize || 100;
    if (mineral3DSizeValue) mineral3DSizeValue.textContent = (settings.mineral3DSize || 100) + '%';
  }

  // Always show Compound Options
  if (compoundOptions) {
    compoundOptions.style.display = 'block';
  }

  // Load MolView-Only Mode option
  if (molSearchModeToggle) molSearchModeToggle.checked = settings.molSearchMode;

  // Load Disable Formula Fallback option
  if (disableFormulaFallbackToggle) disableFormulaFallbackToggle.checked = settings.disableFormulaFallback;

  // Load Rendering Engine option
  if (renderingEngineSelect) renderingEngineSelect.value = settings.renderingEngine || 'chemistrylatex';


  // Load per-style settings (removed - UI elements no longer exist)
  // const currentStyle = settings.viewer3DStyle || 'stick:sphere';
  // const styleSettings = settings.viewer3DStyleSettings || {};
  // const currentStyleSettings = styleSettings[currentStyle] || {};


  // Set rendering engine radio button
  const engineRadios = document.querySelectorAll('input[name="renderingEngine"]');
  engineRadios.forEach(radio => {
    radio.checked = (radio.value === settings.rendererEngine);
  });


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
      showStatus('ChemistryLaTeX ' + (e.target.checked ? 'enabled' : 'disabled'), 'success');
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
      showStatus('Developer mode ' + (e.target.checked ? 'enabled' : 'disabled'), 'success');
    });
  });
}


// Enable AI Flag Control toggle - when disabled, flags in chem: tags are ignored
if (enableAIFlagControlToggle) {
  enableAIFlagControlToggle.addEventListener('change', (e) => {
    const value = e.target.checked;
    chrome.storage.sync.set({ enableAIFlagControl: value }, () => {
      broadcastSettingsChange({ enableAIFlagControl: value });
      showStatus('AI flag control ' + (value ? 'enabled' : 'disabled'), 'success');
    });
  });
}

// Reload All button - re-render all molecules without page reload
console.log('[Popup] reloadAllBtn:', reloadAllBtn);
if (reloadAllBtn) {
  reloadAllBtn.addEventListener('click', () => {
    console.log('[Popup] Reload All button clicked!');
    chrome.runtime.sendMessage({
      type: 'RELOAD_ALL_IMAGES'
    }, (response) => {
      console.log('[Popup] Reload message sent, response:', response);
    });
    showStatus('Reloading all images...', 'success');
  });
} else {
  console.warn('[Popup] reloadAllBtn not found in DOM!');
}

// Clear Cache button - clear all cached SMILES lookups
const clearCacheBtn = document.getElementById('clearCacheBtn');
console.log('[Popup] clearCacheBtn:', clearCacheBtn);
if (clearCacheBtn) {
  clearCacheBtn.addEventListener('click', () => {
    console.log('[Popup] Clear Cache button clicked!');

    // Send message to ALL tabs to clear in-memory cache
    chrome.tabs.query({}, (tabs) => {
      tabs.forEach(tab => {
        chrome.tabs.sendMessage(tab.id, { type: 'CLEAR_CACHE' }).catch(() => {
          // Ignore errors for tabs without content script
        });
      });
      console.log('[Popup] Sent CLEAR_CACHE to', tabs.length, 'tabs');
    });

    // Also send to runtime (background script)
    chrome.runtime.sendMessage({ type: 'CLEAR_CACHE' }, (response) => {
      console.log('[Popup] Clear cache message sent to runtime, response:', response);
    });

    // Also clear any local storage cache
    chrome.storage.local.get(null, (items) => {
      const keysToRemove = Object.keys(items).filter(key =>
        key.startsWith('smiles_') ||
        key.startsWith('cid_') ||
        key.startsWith('cache_') ||
        key.startsWith('search_') ||
        key === 'chemistrylatex_smiles_cache' ||  // CRITICAL: The persistent SMILES cache
        key === 'chemistrylatex_svg_cache'        // NEW: The persistent SVG cache
      );
      if (keysToRemove.length > 0) {
        chrome.storage.local.remove(keysToRemove, () => {
          console.log('[Popup] Cleared', keysToRemove.length, 'cached items from storage:', keysToRemove);
        });
      } else {
        // Even if no prefix matches, ALWAYS try to clear the main caches
        chrome.storage.local.remove(['chemistrylatex_smiles_cache', 'chemistrylatex_svg_cache'], () => {
          console.log('[Popup] Cleared chemistrylatex_smiles_cache and chemistrylatex_svg_cache');
        });
      }
    });

    showStatus('Cache cleared! Reload images to fetch fresh data.', 'success');
  });
} else {
  console.warn('[Popup] clearCacheBtn not found in DOM!');
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
  // Extract the keys that changed
  const changedKeys = Object.keys(changedSettings);
  chrome.runtime.sendMessage({
    type: 'SETTINGS_CHANGED',
    settings: changedSettings,
    changedKeys: changedKeys  // Tell content script which settings changed
  });
}

// Rendering options
// Stereochemistry toggle (use isomeric SMILES)
if (sdUseStereochemistryToggle) {
  sdUseStereochemistryToggle.addEventListener('change', (e) => {
    const value = e.target.checked;
    console.log('[Popup] Setting useStereochemistry to:', value);
    chrome.storage.sync.set({ useStereochemistry: value }, () => {
      if (chrome.runtime.lastError) {
        console.error('[Popup] Error saving useStereochemistry:', chrome.runtime.lastError);
        showStatus('Error saving setting: ' + chrome.runtime.lastError.message, 'error');
        return;
      }
      console.log('[Popup] Successfully saved useStereochemistry:', value);
      broadcastSettingsChange({ useStereochemistry: value });
      showStatus('Stereochemistry ' + (value ? 'enabled' : 'disabled'), 'success');
    });
  });
}

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

// sdAutoAdaptToggle event listener REMOVED - no auto-adapt anymore

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

if (sdAverageSizeSlider) {
  sdAverageSizeSlider.addEventListener('input', (e) => {
    sdAverageSizeValue.textContent = e.target.value + '%';
    // Apply size change in real-time without saving
    broadcastSettingsChange({ sdAverageSize: parseInt(e.target.value) });
  });

  sdAverageSizeSlider.addEventListener('change', (e) => {
    const value = parseInt(e.target.value);
    chrome.storage.sync.set({ sdAverageSize: value }, () => {
      broadcastSettingsChange({ sdAverageSize: value });
      showStatus('Average size set to ' + value + '%', 'success');
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

// Data Source toggle event listeners (DISABLED - elements removed from popup.html)
// These were causing "ReferenceError: searchCompoundsToggle is not defined" which blocked all subsequent code
/*
if (searchCompoundsToggle) {
  searchCompoundsToggle.addEventListener('change', (e) => {
    chrome.storage.sync.set({ searchCompounds: e.target.checked }, () => {
      showStatus('Compound search ' + (e.target.checked ? 'enabled' : 'disabled') + '. Reload page to apply.', 'success');
    });
  });
}

if (searchBiomoleculesToggle) {
  searchBiomoleculesToggle.addEventListener('change', (e) => {
    chrome.storage.sync.set({ searchBiomolecules: e.target.checked }, () => {
      showStatus('Biomolecule search ' + (e.target.checked ? 'enabled' : 'disabled') + '. Reload page to apply.', 'success');
    });
  });
}

if (searchMineralsToggle) {
  searchMineralsToggle.addEventListener('change', (e) => {
    chrome.storage.sync.set({ searchMinerals: e.target.checked }, () => {
      showStatus('Mineral search ' + (e.target.checked ? 'enabled' : 'disabled') + '. Reload page to apply.', 'success');
    });
  });
}
*/

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
    const value = e.target.value;
    saveSetting({ molviewRepresentation: value }, 'Representation set to ' + e.target.options[e.target.selectedIndex].text, (success) => {
      if (success) broadcastSettingsChange({ molviewRepresentation: value });
    });
  });
}

if (compoundMolviewBgColorSelect) {
  console.log('[Popup] ✅ Adding event listener to compoundMolviewBgColorSelect');
  compoundMolviewBgColorSelect.addEventListener('change', (e) => {
    const value = e.target.value;
    console.log('[Popup] 🔄 compoundMolviewBgColorSelect CHANGED to:', value);
    saveSettingWithVerification('compoundMolviewBgColor', value, (success, error) => {
      if (success) {
        broadcastSettingsChange({ compoundMolviewBgColor: value });
        showStatus('Compound background set to ' + e.target.options[e.target.selectedIndex].text, 'success');
      } else {
        console.error('[Popup] ❌ Failed to save compoundMolviewBgColor:', error);
        showStatus('Error saving background color: ' + (error || 'Unknown error'), 'error');
      }
    });
  });
} else {
  console.error('[Popup] ❌ compoundMolviewBgColorSelect NOT FOUND IN DOM!');
}

if (molviewEngineSelect) {
  molviewEngineSelect.addEventListener('change', (e) => {
    const value = e.target.value;
    saveSetting({ molviewEngine: value }, 'Engine set to ' + e.target.options[e.target.selectedIndex].text, (success) => {
      if (success) broadcastSettingsChange({ molviewEngine: value });
    });
  });
}

// Rendering Engine - instant update
if (renderingEngineSelect) {
  renderingEngineSelect.addEventListener('change', (e) => {
    const value = e.target.value;
    const displayNames = { chemistrylatex: 'ChemistryLaTeX', rdkit: 'RDKit', kekule: 'Kekule' };
    saveSetting({ renderingEngine: value }, 'Rendering engine set to ' + (displayNames[value] || value), (success) => {
      if (success) broadcastSettingsChange({ renderingEngine: value });
    });
  });
}

// Biomolecule Options event listeners - ALL broadcast changes for instant updates
if (proteinRemoveWhiteBgToggle) {
  console.log('[Popup] ✅ Adding event listener to proteinRemoveWhiteBgToggle');
  proteinRemoveWhiteBgToggle.addEventListener('change', (e) => {
    try {
      console.log('[Popup] 🔄 proteinRemoveWhiteBgToggle CHANGED to:', e.target.checked);
      const value = e.target.checked;
      chrome.storage.sync.set({ proteinRemoveWhiteBg: value }, () => {
        if (chrome.runtime.lastError) {
          console.error('[Popup] ❌ ERROR saving proteinRemoveWhiteBg:', chrome.runtime.lastError);
          return;
        }
        console.log('[Popup] ✅ SAVED proteinRemoveWhiteBg:', value);
        broadcastSettingsChange({ proteinRemoveWhiteBg: value });
        showStatus('White background removal ' + (value ? 'enabled' : 'disabled'), 'success');
      });
    } catch (err) {
      console.error('[Popup] ❌ EXCEPTION in proteinRemoveWhiteBgToggle handler:', err);
    }
  });
} else {
  console.error('[Popup] ❌ proteinRemoveWhiteBgToggle NOT FOUND IN DOM!');
}

if (molviewBioAssemblyToggle) {
  molviewBioAssemblyToggle.addEventListener('change', (e) => {
    const enabled = e.target.checked;
    chrome.storage.sync.set({ molviewBioAssembly: enabled }, () => {
      broadcastSettingsChange({ molviewBioAssembly: enabled });
      showStatus('Bio Assembly ' + (enabled ? 'shown' : 'hidden'), 'success');

      // Enable/disable the sub-settings panel
      if (bioAssemblySettings) {
        if (enabled) {
          bioAssemblySettings.style.opacity = '1';
          bioAssemblySettings.style.pointerEvents = 'auto';
        } else {
          bioAssemblySettings.style.opacity = '0.5';
          bioAssemblySettings.style.pointerEvents = 'none';
        }
      }
    });
  });
}

// Bio Assembly sub-settings event listeners
if (bioAssemblyViewerSelect) {
  bioAssemblyViewerSelect.addEventListener('change', (e) => {
    const value = e.target.value;
    chrome.storage.sync.set({ bioAssemblyViewer: value }, () => {
      broadcastSettingsChange({ bioAssemblyViewer: value });
      showStatus('Bio Viewer set to ' + e.target.options[e.target.selectedIndex].text, 'success');
    });
  });
}

if (bioAssemblyGraphicsSelect) {
  bioAssemblyGraphicsSelect.addEventListener('change', (e) => {
    const value = e.target.value;
    chrome.storage.sync.set({ bioAssemblyGraphics: value }, () => {
      broadcastSettingsChange({ bioAssemblyGraphics: value });
      showStatus('Graphics mode set to ' + e.target.options[e.target.selectedIndex].text, 'success');
    });
  });
}

if (bioAssemblyPdbProviderSelect) {
  bioAssemblyPdbProviderSelect.addEventListener('change', (e) => {
    const value = e.target.value;
    chrome.storage.sync.set({ bioAssemblyPdbProvider: value }, () => {
      broadcastSettingsChange({ bioAssemblyPdbProvider: value });
      showStatus('PDB Provider set to ' + e.target.options[e.target.selectedIndex].text, 'success');
    });
  });
}

if (molviewChainTypeSelect) {
  molviewChainTypeSelect.addEventListener('change', (e) => {
    const value = e.target.value;
    saveSetting({ molviewChainType: value }, 'Chain Type set to ' + e.target.options[e.target.selectedIndex].text, (success) => {
      if (success) broadcastSettingsChange({ molviewChainType: value });
    });
  });
}

if (molviewChainBondsToggle) {
  molviewChainBondsToggle.addEventListener('change', (e) => {
    const value = e.target.checked;
    saveSetting({ molviewChainBonds: value }, 'Chain Bonds ' + (value ? 'shown' : 'hidden'), (success) => {
      if (success) broadcastSettingsChange({ molviewChainBonds: value });
    });
  });
}

if (molviewChainColorSelect) {
  molviewChainColorSelect.addEventListener('change', (e) => {
    const value = e.target.value;
    saveSetting({ molviewChainColor: value }, 'Chain Color set to ' + e.target.options[e.target.selectedIndex].text, (success) => {
      if (success) broadcastSettingsChange({ molviewChainColor: value });
    });
  });
}

if (proteinMolviewBgColorSelect) {
  proteinMolviewBgColorSelect.addEventListener('change', (e) => {
    const value = e.target.value;
    saveSetting({ proteinMolviewBgColor: value }, 'Protein background color set to ' + e.target.options[e.target.selectedIndex].text, (success) => {
      if (success) broadcastSettingsChange({ proteinMolviewBgColor: value });
    });
  });
}

// Mineral Options event listeners
if (mineralRepresentationSelect) {
  mineralRepresentationSelect.addEventListener('change', (e) => {
    const value = e.target.value;
    saveSetting({ mineralRepresentation: value }, 'Mineral representation set to ' + e.target.options[e.target.selectedIndex].text, (success) => {
      if (success) broadcastSettingsChange({ mineralRepresentation: value });
    });
  });
}

if (mineralMolviewBgColorSelect) {
  mineralMolviewBgColorSelect.addEventListener('change', (e) => {
    const value = e.target.value;
    saveSetting({ mineralMolviewBgColor: value }, 'Mineral background color set to ' + e.target.options[e.target.selectedIndex].text, (success) => {
      if (success) broadcastSettingsChange({ mineralMolviewBgColor: value });
    });
  });
}

if (mineralCrystallographySelect) {
  mineralCrystallographySelect.addEventListener('change', (e) => {
    const value = e.target.value;
    saveSetting({ mineralCrystallography: value }, 'Mineral crystallography set to ' + e.target.options[e.target.selectedIndex].text + '. Reload page to apply.');
  });
}

// Add event listeners for rendering engine radio buttons
const engineRadios = document.querySelectorAll('input[name="renderingEngine"]');
engineRadios.forEach(radio => {
  radio.addEventListener('change', (e) => {
    if (e.target.checked) {
      const engine = e.target.value;
      chrome.storage.sync.set({ rendererEngine: engine }, () => {
        updateEngineInfo();
        showStatus('Rendering engine updated. Reload page to apply.', 'success');
      });
    }
  });
});


/**
 * Update engine info display
 */
const biomoleculeOptions = document.getElementById('biomoleculeOptions');
const mineralOptions = document.getElementById('mineralOptions');

function updateEngineInfo() {
  // Always show all options (single engine now)
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

  // Check for dark mode
  var isDark = document.body.classList.contains('dark-mode');

  // Configurable base colors
  var baseR = isDark ? 30 : 255;
  var baseG = isDark ? 30 : 255;
  var baseB = isDark ? 30 : 255;

  var baseColor = 'rgba(' + baseR + ', ' + baseG + ', ' + baseB + ', ';

  sections.forEach(function (section) {
    section.style.backdropFilter = 'blur(' + blurPx + ')';
    section.style.webkitBackdropFilter = 'blur(' + blurPx + ')';
    section.style.background = baseColor + bgOpacity + ')';

    // In dark mode, ensure border is subtle
    if (isDark) {
      section.style.borderColor = 'rgba(255, 255, 255, 0.05)';
    } else {
      section.style.borderColor = 'rgba(255, 255, 255, 0.5)';
    }
  });

  if (header) {
    header.style.backdropFilter = 'blur(' + blurPx + ')';
    header.style.webkitBackdropFilter = 'blur(' + blurPx + ')';
    // Header usually needs to be a bit more opaque or distinct
    var headerOpacity = isDark ? (bgOpacity * 1.5) : (bgOpacity * 0.9);
    // Cap at 0.95
    if (headerOpacity > 0.95) headerOpacity = 0.95;

    header.style.background = baseColor + headerOpacity + ')';
  }

  if (footer) {
    footer.style.backdropFilter = 'blur(' + blurPx + ')';
    footer.style.webkitBackdropFilter = 'blur(' + blurPx + ')';
    var footerOpacity = isDark ? (bgOpacity * 1.5) : (bgOpacity * 0.9);
    if (footerOpacity > 0.95) footerOpacity = 0.95;

    footer.style.background = baseColor + footerOpacity + ')';
  }

  radioOptions.forEach(function (opt) {
    opt.style.background = baseColor + (bgOpacity * 0.8) + ')';
  });
}

// Load UI settings
// ----------------------------------------------------------------------
// MAIN INITIALIZATION: Load All Settings (Dark Mode, UI, Physics)
// ----------------------------------------------------------------------
chrome.storage.sync.get({
  // UI Defaults
  darkMode: false,
  uiBlur: 7,
  uiOpacity: 23,
  // Physics Defaults
  moleculeCount: 22,
  moleculeScale: 0.4,
  cursorGravityStrength: -0.9,
  initialVelocity: 10
}, function (settings) {

  // 1. Apply Dark Mode FIRST (Critical for base colors)
  if (settings.darkMode) {
    document.body.classList.add('dark-mode');
    updateDarkModeIcon(true);
  } else {
    updateDarkModeIcon(false);
  }

  // 2. Set UI Sliders
  if (blurSlider) {
    blurSlider.value = settings.uiBlur;
    blurValue.textContent = settings.uiBlur + 'px';
  }
  if (opacitySlider) {
    opacitySlider.value = settings.uiOpacity;
    opacityValue.textContent = settings.uiOpacity + '%';
  }

  // 3. Apply Glassmorphism (Now that dark mode class is settled)
  applyGlassmorphism(settings.uiBlur, settings.uiOpacity);

  // 4. Set Physics Sliders
  if (molCountSlider) {
    molCountSlider.value = settings.moleculeCount;
    molCountValue.textContent = settings.moleculeCount;
  }
  if (molScaleSlider) {
    molScaleSlider.value = settings.moleculeScale;
    molScaleValue.textContent = settings.moleculeScale;
  }
  if (cursorGravitySlider) {
    cursorGravitySlider.value = settings.cursorGravityStrength;
    cursorGravityValue.textContent = settings.cursorGravityStrength;
  }
  if (initialVelocitySlider) {
    initialVelocitySlider.value = settings.initialVelocity;
    initialVelocityValue.textContent = settings.initialVelocity;
  }
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
// Previous individual storage getter removed to avoid race conditions
// Settings are now loaded in the main unified block above.

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
// Copy ChatGPT Prompt button
const copyChatGPTPromptBtn = document.getElementById('copyChatGPTPrompt');
const savePromptToMemoryToggle = document.getElementById('savePromptToMemoryToggle');

// Memory instruction prefix
const MEMORY_INSTRUCTION = 'SAVE TO MEMORY WITH THE EXACT INSTRUCTIONS SPECIFIED DO NOT REMOVE OR COMPRESS, SAVE EXACTLY AS IS\n\n';

if (copyChatGPTPromptBtn) {
  copyChatGPTPromptBtn.addEventListener('click', async () => {
    try {
      let usageText;

      // Priority order for ChatGPT prompt:
      // 1. Current session server prompt (freshest)
      // 2. Cached server prompt from previous session (chrome.storage.local)
      // 3. Local USAGE.md file (last resort fallback)

      if (window._serverChatGptPrompt) {
        // Server prompt already loaded this session
        usageText = window._serverChatGptPrompt;
        console.log('[Popup] ✅ Using server ChatGPT prompt (current session)');
      } else {
        // Try to get cached prompt from previous session
        const cached = await new Promise(resolve => {
          chrome.storage.local.get(['serverChatGptPrompt'], (data) => {
            resolve(data.serverChatGptPrompt);
          });
        });

        if (cached) {
          usageText = cached;
          console.log('[Popup] ✅ Using cached server ChatGPT prompt');
        } else {
          // Last resort: use local USAGE.md
          const response = await fetch(chrome.runtime.getURL('USAGE.md'));
          usageText = await response.text();
          console.log('[Popup] ⚠️ Using local USAGE.md (server unavailable)');
        }
      }

      // Check if we should add the memory instruction
      if (savePromptToMemoryToggle && savePromptToMemoryToggle.checked) {
        usageText = MEMORY_INSTRUCTION + usageText;
      }

      // Copy to clipboard
      await navigator.clipboard.writeText(usageText);

      // Show success feedback
      const originalText = copyChatGPTPromptBtn.innerHTML;
      const originalBg = copyChatGPTPromptBtn.style.background;
      copyChatGPTPromptBtn.innerHTML = 'Copied';
      copyChatGPTPromptBtn.style.background = '#22c55e';

      // Reset after 2 seconds
      setTimeout(() => {
        copyChatGPTPromptBtn.innerHTML = originalText;
        copyChatGPTPromptBtn.style.background = originalBg;
      }, 2000);
    } catch (error) {
      console.error('Failed to copy prompt:', error);
      copyChatGPTPromptBtn.innerHTML = 'Failed';
      setTimeout(() => {
        copyChatGPTPromptBtn.innerHTML = 'Copy';
      }, 2000);
    }
  });
}

// Donation wallet display
const donationWalletBtn = document.getElementById('donationWallet');
if (donationWalletBtn) {
  donationWalletBtn.addEventListener('click', async (e) => {
    e.preventDefault();

    // Your wallet addresses (replace with your actual addresses)
    const ethAddress = '0x20E8D5f0BeA0fdffF470a6DBb89013Ff9Fd51933';
    const solAddress = '4guoDVFaQc5qJQ8geLhCgUCJxdvWqmGXFAA6Ed26wXk8';

    // Create a modal-style overlay
    const overlay = document.createElement('div');
    overlay.style.cssText = `
      position: fixed;
      top: 0;
      left: 0;
      width: 100%;
      height: 100%;
      background: rgba(0, 0, 0, 0.6);
      backdrop-filter: blur(8px);
      -webkit-backdrop-filter: blur(8px);
      display: flex;
      justify-content: center;
      align-items: center;
      z-index: 10000;
    `;

    // Check for dark mode
    const isDark = document.body.classList.contains('dark-mode');

    const modal = document.createElement('div');
    modal.style.cssText = `
      background: ${isDark ? '#1e1e1e' : 'rgba(255, 255, 255, 0.98)'};
      padding: 32px;
      border-radius: 16px;
      max-width: 420px;
      width: 90%;
      box-shadow: 0 12px 48px rgba(0, 0, 0, ${isDark ? '0.6' : '0.15'});
      color: ${isDark ? '#e0e0e0' : '#1a1a1a'};
      border: 1px solid ${isDark ? '#333' : 'rgba(0, 0, 0, 0.05)'};
      backdrop-filter: blur(10px);
    `;

    modal.innerHTML = `
      <div style="display: flex; justify-content: space-between; align-items: center; margin-bottom: 20px;">
        <h3 style="margin: 0; font-size: 18px; font-weight: 700; color: ${isDark ? '#e0e0e0' : '#1a1a1a'};">Support ChemistryLaTeX</h3>
        <button id="closeWalletModal" style="background: none; border: none; font-size: 24px; cursor: pointer; color: ${isDark ? '#aaa' : '#999'}; padding: 0; line-height: 1;">&times;</button>
      </div>
      
      <p style="margin-bottom: 20px; color: ${isDark ? '#ccc' : '#555'}; line-height: 1.5; font-size: 13px;">
        If you find this extension helpful, consider supporting development.
      </p>

      <div style="margin-bottom: 16px;">
        <label style="display: block; font-size: 11px; font-weight: 600; color: ${isDark ? '#999' : '#666'}; margin-bottom: 6px;">ETH / L2</label>
        <div style="display: flex; align-items: center; gap: 8px; background: ${isDark ? '#333' : '#f5f5f5'}; padding: 10px 12px; border-radius: 8px;">
          <code id="ethAddress" style="font-family: monospace; font-size: 11px; color: ${isDark ? '#fff' : '#333'}; flex: 1; overflow: hidden; text-overflow: ellipsis;">${ethAddress}</code>
          <button id="copyEthBtn" style="background: ${isDark ? '#444' : '#1a1a1a'}; color: #fff; border: none; padding: 5px 10px; border-radius: 6px; font-size: 11px; cursor: pointer; border: 1px solid ${isDark ? '#555' : 'transparent'};">Copy</button>
        </div>
      </div>

      <div style="margin-bottom: 16px;">
        <label style="display: block; font-size: 11px; font-weight: 600; color: ${isDark ? '#999' : '#666'}; margin-bottom: 6px;">Solana</label>
        <div style="display: flex; align-items: center; gap: 8px; background: ${isDark ? '#333' : '#f5f5f5'}; padding: 10px 12px; border-radius: 8px;">
          <code id="solAddress" style="font-family: monospace; font-size: 11px; color: ${isDark ? '#fff' : '#333'}; flex: 1; overflow: hidden; text-overflow: ellipsis;">${solAddress}</code>
          <button id="copySolBtn" style="background: ${isDark ? '#444' : '#1a1a1a'}; color: #fff; border: none; padding: 5px 10px; border-radius: 6px; font-size: 11px; cursor: pointer; border: 1px solid ${isDark ? '#555' : 'transparent'};">Copy</button>
        </div>
      </div>

      <div style="text-align: center; font-size: 11px; color: #999; margin-top: 16px;">Thank you!</div>
    `;

    overlay.appendChild(modal);
    document.body.appendChild(overlay);

    // Copy ETH address
    document.getElementById('copyEthBtn').addEventListener('click', async () => {
      await navigator.clipboard.writeText(ethAddress);
      const btn = document.getElementById('copyEthBtn');
      btn.textContent = 'Copied!';
      btn.style.background = '#27ae60';
      setTimeout(() => {
        btn.textContent = 'Copy';
        btn.style.background = '#000';
      }, 2000);
    });

    // Copy SOL address
    document.getElementById('copySolBtn').addEventListener('click', async () => {
      await navigator.clipboard.writeText(solAddress);
      const btn = document.getElementById('copySolBtn');
      btn.textContent = 'Copied!';
      btn.style.background = '#27ae60';
      setTimeout(() => {
        btn.textContent = 'Copy';
        btn.style.background = '#000';
      }, 2000);
    });

    // Close modal
    document.getElementById('closeWalletModal').addEventListener('click', () => {
      overlay.remove();
    });

    // Close on overlay click
    overlay.addEventListener('click', (e) => {
      if (e.target === overlay) {
        overlay.remove();
      }
    });
  });
}


// Compound 3D Size Slider
if (compound3DSizeSlider) {
  compound3DSizeSlider.addEventListener('input', (e) => {
    compound3DSizeValue.textContent = e.target.value + '%';
    // Broadcast for real-time update
    broadcastSettingsChange({ compound3DSize: parseInt(e.target.value) });
  });

  compound3DSizeSlider.addEventListener('change', (e) => {
    const value = parseInt(e.target.value);
    chrome.storage.sync.set({ compound3DSize: value }, () => {
      broadcastSettingsChange({ compound3DSize: value });
      showStatus('Compound 3D viewer size set to ' + value + '%', 'success');
    });
  });
}
// Mineral 3D Size Slider
if (mineral3DSizeSlider) {
  mineral3DSizeSlider.addEventListener('input', (e) => {
    mineral3DSizeValue.textContent = e.target.value + '%';
    // Broadcast for real-time update
    broadcastSettingsChange({ mineral3DSize: parseInt(e.target.value) });
  });

  mineral3DSizeSlider.addEventListener('change', (e) => {
    const value = parseInt(e.target.value);
    chrome.storage.sync.set({ mineral3DSize: value }, () => {
      broadcastSettingsChange({ mineral3DSize: value });
      showStatus('Mineral 3D viewer size set to ' + value + '%', 'success');
    });
  });
}

// === PATCH NOTES FEATURE ===
const patchNotesBtn = document.getElementById('patchNotesBtn');
const notesModal = document.getElementById('notesModal');
const closeNotesModal = document.getElementById('closeNotesModal');
const notesContent = document.getElementById('notesContent');

const NOTES_SERVER = 'https://server-chemistryrenderer.vercel.app';
let cachedNotes = null;

if (patchNotesBtn && notesModal) {
  patchNotesBtn.addEventListener('click', async () => {
    notesModal.classList.add('show');

    // Use cached notes if available, otherwise fetch
    if (cachedNotes) {
      notesContent.textContent = cachedNotes;
    } else {
      notesContent.textContent = 'Loading...';
      try {
        const response = await fetch(`${NOTES_SERVER}/api/note`);
        const data = await response.json();
        if (data.success) {
          cachedNotes = data.note;
          notesContent.textContent = data.note;
        } else {
          notesContent.textContent = 'Failed to load notes.';
        }
      } catch (err) {
        notesContent.textContent = 'Could not connect to server.';
      }
    }
  });

  closeNotesModal.addEventListener('click', () => {
    notesModal.classList.remove('show');
  });

  notesModal.addEventListener('click', (e) => {
    if (e.target === notesModal) {
      notesModal.classList.remove('show');
    }
  });
}

// === DARK MODE FEATURE ===
const darkModeBtn = document.getElementById('darkModeBtn');
const darkModeIcon = document.getElementById('darkModeIcon');

// Load dark mode setting
// Dark mode is now loaded in the main unified initialization block
// to ensure correct styling order.

function updateDarkModeIcon(isDark) {
  if (!darkModeIcon) return;

  if (isDark) {
    // Switch to sun icon
    darkModeIcon.innerHTML = '<circle cx="12" cy="12" r="5"></circle><line x1="12" y1="1" x2="12" y2="3"></line><line x1="12" y1="21" x2="12" y2="23"></line><line x1="4.22" y1="4.22" x2="5.64" y2="5.64"></line><line x1="18.36" y1="18.36" x2="19.78" y2="19.78"></line><line x1="1" y1="12" x2="3" y2="12"></line><line x1="21" y1="12" x2="23" y2="12"></line><line x1="4.22" y1="19.78" x2="5.64" y2="18.36"></line><line x1="18.36" y1="5.64" x2="19.78" y2="4.22"></line>';
  } else {
    // Switch to moon icon
    darkModeIcon.innerHTML = '<path d="M21 12.79A9 9 0 1 1 11.21 3 7 7 0 0 0 21 12.79z"></path>';
  }
}

if (darkModeBtn) {
  darkModeBtn.addEventListener('click', () => {
    const isDark = document.body.classList.toggle('dark-mode');
    updateDarkModeIcon(isDark);

    // Re-apply glassmorphism with current slider values to update base color (white <-> dark)
    if (opacitySlider && blurSlider) {
      applyGlassmorphism(parseInt(blurSlider.value), parseInt(opacitySlider.value));
    }

    // Save setting
    chrome.storage.sync.set({ darkMode: isDark });

    // Update molecule animation colors
    if (window.updateMoleculeAnimation) {
      window.updateMoleculeAnimation();
    }
  });
}

// ========================================================
// SERVER-CONTROLLED CONFIGURATION
// ========================================================
// Fetches dynamic settings from server (Discord link, default settings, announcements)
// This allows updates without going through the extension review process

const CONFIG_SERVER = 'https://server-chemistryrenderer.vercel.app';
const CONFIG_CACHE_KEY = 'chemistrylatex_server_config';
const SEEN_ANNOUNCEMENTS_KEY = 'chemistrylatex_seen_announcements';

/**
 * Fetch and apply server configuration
 * - Updates Discord link dynamically
 * - Applies default settings ONLY to users who haven't changed them
 * - Shows announcements (one-time)
 */
async function fetchServerConfig() {
  try {
    // First check the cache
    const cacheData = await new Promise(resolve => {
      chrome.storage.local.get([CONFIG_CACHE_KEY], (data) => resolve(data[CONFIG_CACHE_KEY]));
    });

    // Use cache if less than 5 minutes old
    const cacheAge = cacheData ? (Date.now() - (cacheData.fetchedAt || 0)) : Infinity;
    if (cacheData && cacheAge < 5 * 60 * 1000) {
      console.log('[Popup] 📦 Using cached server config (age:', Math.round(cacheAge / 1000), 's)');
      applyServerConfig(cacheData.config);
      return;
    }

    // Fetch fresh config from server
    console.log('[Popup] 🌐 Fetching server config...');
    const response = await fetch(`${CONFIG_SERVER}/api/config`, {
      method: 'GET',
      headers: { 'Content-Type': 'application/json' }
    });

    if (!response.ok) {
      throw new Error(`Server returned ${response.status}`);
    }

    const data = await response.json();
    if (!data.success || !data.config) {
      throw new Error('Invalid config response');
    }

    // Cache the config
    chrome.storage.local.set({
      [CONFIG_CACHE_KEY]: {
        config: data.config,
        fetchedAt: Date.now()
      }
    });

    console.log('[Popup] ✅ Server config loaded:', data.config.version);
    applyServerConfig(data.config);

  } catch (error) {
    console.warn('[Popup] ⚠️ Could not fetch server config:', error.message);
    // Silently fail - extension works fine without server config
  }
}

/**
 * Apply server configuration
 */
function applyServerConfig(config) {
  if (!config) return;

  // 1. Update Discord link dynamically
  const discordLink = document.getElementById('discordLink');
  if (discordLink && config.discordLink) {
    discordLink.href = config.discordLink;
    console.log('[Popup] 🔗 Discord link updated:', config.discordLink);
  }

  // 2. Store wallet addresses for donation button
  if (config.wallets) {
    window._serverWallets = config.wallets;
    console.log('[Popup] 💰 Wallet addresses loaded');
  }

  // 3. Store ChatGPT prompt for AI features
  if (config.chatGptPrompt) {
    window._serverChatGptPrompt = config.chatGptPrompt;  // For immediate access
    chrome.storage.local.set({ serverChatGptPrompt: config.chatGptPrompt });  // For persistence
    console.log('[Popup] 🤖 ChatGPT prompt loaded');
  }

  // 4. Handle patch notes notification
  if (config.patchNotes && config.patchNotesVersion) {
    handlePatchNotesNotification(config.patchNotes, config.patchNotesVersion);
  }

  // 5. Apply default settings ONLY to users who haven't explicitly changed them
  // We check if the user has a "userModifiedSettings" flag for each setting
  if (config.defaultSettings) {
    applyDefaultSettings(config.defaultSettings);
  }

  // 6. Apply default UI settings
  if (config.defaultUI) {
    applyDefaultUISettings(config.defaultUI);
  }

  // 7. Show announcement if not seen
  if (config.announcement) {
    showAnnouncement(config.announcement);
  }
}

/**
 * Handle patch notes notification bell
 * Shows notification dot if there's new content, rings if important
 */
function handlePatchNotesNotification(patchNotes, patchNotesVersion) {
  // Store patch notes for when user clicks the bell
  window._serverPatchNotes = patchNotes;
  window._serverPatchNotesVersion = patchNotesVersion;

  // Check if user has seen this version
  chrome.storage.local.get(['lastSeenPatchNotesVersion'], (data) => {
    const lastSeen = data.lastSeenPatchNotesVersion;
    const isNew = lastSeen !== patchNotesVersion;

    // Check for $important$ or $notif$ flags
    const isImportant = patchNotes.startsWith('$important$');
    const hasNotification = patchNotes.startsWith('$notif$') || isImportant;

    // Get notification elements
    const patchNotesBtn = document.getElementById('patchNotesBtn');
    const notifDot = document.getElementById('patchNotifDot');

    // Clean the patch notes content (remove flags)
    let cleanPatchNotes = patchNotes
      .replace('$important$', '')
      .replace('$notif$', '')
      .trim();
    window._serverPatchNotes = cleanPatchNotes;

    // Show notification if new AND has notification flag
    if (isNew && hasNotification) {
      // Show notification dot
      if (notifDot) notifDot.style.display = 'block';

      // If important, add red color and ring animation
      if (isImportant && patchNotesBtn) {
        patchNotesBtn.classList.add('important', 'ringing');
      }

      console.log('[Popup] 🔔 New patch notes available (v' + patchNotesVersion + ')');
    } else if (isNew) {
      // New content but no notification flag - still show dot but no ring
      if (notifDot) notifDot.style.display = 'block';
    }
  });
}

/**
 * Apply default settings only to users who haven't changed them
 * We track which settings the user has explicitly modified
 */
function applyDefaultSettings(defaults) {
  // Get the list of settings the user has explicitly modified
  chrome.storage.sync.get(['userModifiedSettings'], (data) => {
    const userModified = data.userModifiedSettings || {};

    // For each default setting, only apply if user hasn't modified it
    const settingsToApply = {};
    let appliedCount = 0;

    for (const [key, value] of Object.entries(defaults)) {
      if (!userModified[key]) {
        settingsToApply[key] = value;
        appliedCount++;
      }
    }

    if (appliedCount > 0) {
      chrome.storage.sync.set(settingsToApply, () => {
        console.log(`[Popup] 📝 Applied ${appliedCount} default settings from server`);
      });
    }
  });
}

/**
 * Apply default UI settings only to users who haven't changed them
 */
function applyDefaultUISettings(defaults) {
  chrome.storage.sync.get(['userModifiedUI'], (data) => {
    const userModified = data.userModifiedUI || {};

    const uiToApply = {};
    let appliedCount = 0;

    for (const [key, value] of Object.entries(defaults)) {
      if (!userModified[key]) {
        uiToApply[key] = value;
        appliedCount++;
      }
    }

    if (appliedCount > 0) {
      chrome.storage.sync.set(uiToApply, () => {
        console.log(`[Popup] 🎨 Applied ${appliedCount} default UI settings from server`);
      });
    }
  });
}

/**
 * Show announcement if not already seen
 */
function showAnnouncement(announcement) {
  if (!announcement || !announcement.id) return;

  chrome.storage.local.get([SEEN_ANNOUNCEMENTS_KEY], (data) => {
    const seenIds = data[SEEN_ANNOUNCEMENTS_KEY] || [];

    if (seenIds.includes(announcement.id)) {
      return; // Already seen
    }

    // Show the announcement
    const isDark = document.body.classList.contains('dark-mode');
    const typeColors = {
      info: '#2196F3',
      warning: '#ff9800',
      success: '#4CAF50'
    };
    const color = typeColors[announcement.type] || typeColors.info;

    // Create announcement banner
    const banner = document.createElement('div');
    banner.style.cssText = `
      position: fixed;
      top: 0;
      left: 0;
      right: 0;
      background: ${color};
      color: white;
      padding: 10px 16px;
      font-size: 13px;
      z-index: 10001;
      display: flex;
      align-items: center;
      justify-content: space-between;
      box-shadow: 0 2px 8px rgba(0,0,0,0.2);
    `;

    const content = document.createElement('div');
    content.innerHTML = `
      <strong>${announcement.title || 'Announcement'}</strong>
      <span style="margin-left: 8px;">${announcement.message || ''}</span>
      ${announcement.link ? `<a href="${announcement.link}" target="_blank" style="color: white; margin-left: 8px; text-decoration: underline;">Learn More</a>` : ''}
    `;

    const closeBtn = document.createElement('button');
    closeBtn.textContent = '×';
    closeBtn.style.cssText = `
      background: none;
      border: none;
      color: white;
      font-size: 20px;
      cursor: pointer;
      padding: 0 4px;
    `;
    closeBtn.onclick = () => {
      banner.remove();
      // Mark as seen
      seenIds.push(announcement.id);
      chrome.storage.local.set({ [SEEN_ANNOUNCEMENTS_KEY]: seenIds });
    };

    banner.appendChild(content);
    banner.appendChild(closeBtn);
    document.body.appendChild(banner);

    // Auto-dismiss after 10 seconds
    setTimeout(() => {
      if (banner.parentNode) {
        banner.remove();
        seenIds.push(announcement.id);
        chrome.storage.local.set({ [SEEN_ANNOUNCEMENTS_KEY]: seenIds });
      }
    }, 10000);
  });
}

// Fetch server config on popup load
fetchServerConfig();
