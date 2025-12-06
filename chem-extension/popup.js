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
const statusDiv = document.getElementById('status');
const engineInfo = document.getElementById('engineInfo');

// SmilesDrawer options
const smilesDrawerOptions = document.getElementById('smilesDrawerOptions');
const sdShowCarbonsToggle = document.getElementById('sdShowCarbonsToggle');
const sdAromaticRingsToggle = document.getElementById('sdAromaticRingsToggle');
const sdShowMethylsToggle = document.getElementById('sdShowMethylsToggle');
const sdAtomNumbersToggle = document.getElementById('sdAtomNumbersToggle');
const sdShowHydrogensToggle = document.getElementById('sdShowHydrogensToggle');
const sdFlipHorizontalToggle = document.getElementById('sdFlipHorizontalToggle');
const sdFlipVerticalToggle = document.getElementById('sdFlipVerticalToggle');
const sdThemeSelect = document.getElementById('sdThemeSelect');
const sdRotateSlider = document.getElementById('sdRotateSlider');
const sdRotateValue = document.getElementById('sdRotateValue');
const sdBondThicknessSlider = document.getElementById('sdBondThicknessSlider');
const sdBondThicknessValue = document.getElementById('sdBondThicknessValue');
const sdGradientColorsToggle = document.getElementById('sdGradientColorsToggle');
const sdScaleByWeightToggle = document.getElementById('sdScaleByWeightToggle');

// PubChem options
const pubchemOptions = document.getElementById('pubchemOptions');
const pubchemImageSizeSelect = document.getElementById('pubchemImageSizeSelect');
const pubchem3DToggle = document.getElementById('pubchem3DToggle');
const pubchemRecordTypeSelect = document.getElementById('pubchemRecordTypeSelect');
const pubchemDirectFetchToggle = document.getElementById('pubchemDirectFetchToggle');
const pubchemRemoveBgToggle = document.getElementById('pubchemRemoveBgToggle');
const pubchemSharpenToggle = document.getElementById('pubchemSharpenToggle');
const aiMolecularControlToggle = document.getElementById('aiMolecularControlToggle');

// CDK Depict options
const cdkDepictOptions = document.getElementById('cdkDepictOptions');
const cdkShowCarbonsToggle = document.getElementById('cdkShowCarbonsToggle');
const cdkShowMethylsToggle = document.getElementById('cdkShowMethylsToggle');
const cdkColorSchemeSelect = document.getElementById('cdkColorSchemeSelect');
const cdkHydrogenDisplaySelect = document.getElementById('cdkHydrogenDisplaySelect');
const cdkAtomNumbersToggle = document.getElementById('cdkAtomNumbersToggle');
const cdkAnnotationSelect = document.getElementById('cdkAnnotationSelect');
const cdkZoomSlider = document.getElementById('cdkZoomSlider');
const cdkZoomValue = document.getElementById('cdkZoomValue');

// 3D Viewer Settings section
const viewer3DSettings = document.getElementById('viewer3DSettings');
const viewer3DSourceSelect = document.getElementById('viewer3DSourceSelect');
const viewer3DSizeSelect = document.getElementById('viewer3DSizeSelect');
const viewer3DBgColorSelect = document.getElementById('viewer3DBgColorSelect');

// MolView specific options
const molviewOptions = document.getElementById('molviewOptions');
const molviewRepresentationSelect = document.getElementById('molviewRepresentationSelect');
const molviewEngineSelect = document.getElementById('molviewEngineSelect');
const molviewCrystallographySelect = document.getElementById('molviewCrystallographySelect');

// Protein Options
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

// 3D Viewer options
const enable3DViewerToggle = document.getElementById('enable3DViewerToggle');
const default3DViewSelect = document.getElementById('default3DViewSelect');
const viewer3DStyleSelect = document.getElementById('viewer3DStyleSelect');
const viewer3DStickRadiusSelect = document.getElementById('viewer3DStickRadiusSelect');
const viewer3DSphereRadiusSelect = document.getElementById('viewer3DSphereRadiusSelect');
const viewer3DAutoRotateToggle = document.getElementById('viewer3DAutoRotateToggle');

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

// Load current settings
chrome.storage.sync.get({
  enabled: true,
  renderMhchem: true,
  renderChemfig: true,
  performanceMode: true,
  rendererEngine: 'client-side',  // Default to SmilesDrawer client-side
  devMode: false,
  // SmilesDrawer options
  sdShowCarbons: false,
  sdAromaticRings: true,  // Show aromatic circles by default
  sdShowMethyls: false,
  sdAtomNumbers: false,
  sdShowHydrogens: false,
  sdFlipHorizontal: false,
  sdFlipVertical: false,
  sdTheme: 'light',
  sdRotate: 0,
  // Image size controls
  saveSizePerImage: false,
  saveSizeBySMILES: true,  // FIX: Enable by default - global size saving across all pages
  // PubChem options
  pubchemImageSize: 'large',
  pubchem3DEnabled: true,
  pubchemRecordType: '2d',
  pubchemDirectFetch: true,  // Default to direct fetch from PubChem (no local server needed)
  pubchemRemoveBg: false,
  pubchemSharpenImages: true,
  // CDK Depict options
  cdkShowCarbons: false,
  cdkShowMethyls: false,
  cdkColorScheme: 'coc',  // Clear/Transparent
  cdkHydrogenDisplay: 'minimal',
  cdkAtomNumbers: false,
  cdkAnnotation: 'none',
  cdkZoom: 1.5,
  // 3D Viewer options
  enable3DViewer: false,
  default3DView: '2d',
  viewer3DSource: '3dmol',
  viewer3DStyle: 'stick:sphere',
  viewer3DAutoRotate: true,
  viewer3DSize: 'normal',
  viewer3DSize: 'normal',
  viewer3DBgColor: '#1a1a2e',  // Default dark blue background
  // MolView specific options
  molviewRepresentation: 'ballAndStick',
  molviewEngine: 'glmol',
  molviewCrystallography: 'none',
  // Protein Options
  proteinRemoveWhiteBg: false,
  molviewBioAssembly: false,
  molviewChainType: 'ribbon',
  molviewChainBonds: false,
  molviewChainColor: 'ss',
  proteinMolviewBgColor: 'black',
  // Mineral Options
  mineralRepresentation: 'ballAndStick',
  mineralMolviewBgColor: 'black',
  mineralCrystallography: 'supercell_2x2x2',
  // Per-style settings for stick and sphere radius
  viewer3DStyleSettings: {
    'stick': { stickRadius: '0.15' },
    'line': {},
    'cross': {},
    'sphere': { sphereRadius: '0.7' },  // Default to larger for CPK
    'stick:sphere': { stickRadius: '0.15', sphereRadius: '0.3' },
    'cartoon': {}
  },
  // AI Molecular Control
  enableAIMolecularControl: false,
  // MolView-Only Mode
  molSearchMode: false
}, (settings) => {
  enabledToggle.checked = settings.enabled;
  mhchemToggle.checked = settings.renderMhchem;
  chemfigToggle.checked = settings.renderChemfig;
  perfModeToggle.checked = settings.performanceMode;
  devModeToggle.checked = settings.devMode;

  // Load SmilesDrawer options
  if (sdShowCarbonsToggle) sdShowCarbonsToggle.checked = settings.sdShowCarbons;
  if (sdAromaticRingsToggle) sdAromaticRingsToggle.checked = settings.sdAromaticRings;
  if (sdShowMethylsToggle) sdShowMethylsToggle.checked = settings.sdShowMethyls;
  if (sdAtomNumbersToggle) sdAtomNumbersToggle.checked = settings.sdAtomNumbers;
  if (sdShowHydrogensToggle) sdShowHydrogensToggle.checked = settings.sdShowHydrogens;
  if (sdFlipHorizontalToggle) sdFlipHorizontalToggle.checked = settings.sdFlipHorizontal;
  if (sdFlipVerticalToggle) sdFlipVerticalToggle.checked = settings.sdFlipVertical;
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

  // Load PubChem options
  if (pubchemImageSizeSelect) pubchemImageSizeSelect.value = settings.pubchemImageSize;
  if (pubchem3DToggle) pubchem3DToggle.checked = settings.pubchem3DEnabled;
  if (pubchemRecordTypeSelect) pubchemRecordTypeSelect.value = settings.pubchemRecordType;
  if (pubchemDirectFetchToggle) pubchemDirectFetchToggle.checked = settings.pubchemDirectFetch;
  if (pubchemRemoveBgToggle) pubchemRemoveBgToggle.checked = settings.pubchemRemoveBg;
  if (pubchemSharpenToggle) pubchemSharpenToggle.checked = settings.pubchemSharpenImages !== false;
  if (aiMolecularControlToggle) aiMolecularControlToggle.checked = settings.enableAIMolecularControl;

  // Load CDK Depict options
  if (cdkShowCarbonsToggle) cdkShowCarbonsToggle.checked = settings.cdkShowCarbons || false;
  if (cdkShowMethylsToggle) cdkShowMethylsToggle.checked = settings.cdkShowMethyls || false;
  if (cdkColorSchemeSelect) cdkColorSchemeSelect.value = settings.cdkColorScheme || 'coc';
  if (cdkHydrogenDisplaySelect) cdkHydrogenDisplaySelect.value = settings.cdkHydrogenDisplay || 'minimal';
  if (cdkAtomNumbersToggle) cdkAtomNumbersToggle.checked = settings.cdkAtomNumbers || false;
  if (cdkAnnotationSelect) cdkAnnotationSelect.value = settings.cdkAnnotation || 'none';
  if (cdkZoomSlider) {
    cdkZoomSlider.value = settings.cdkZoom || 1.5;
    if (cdkZoomValue) cdkZoomValue.textContent = (settings.cdkZoom || 1.5) + 'x';
  }

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

  // Load MolView specific options
  if (molviewRepresentationSelect) molviewRepresentationSelect.value = settings.molviewRepresentation || 'ballAndStick';
  if (molviewEngineSelect) molviewEngineSelect.value = settings.molviewEngine || 'glmol';
  if (molviewCrystallographySelect) molviewCrystallographySelect.value = settings.molviewCrystallography || 'none';

  // Load Protein Options
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

  // Show/hide MolView options based on current source
  if (molviewOptions) {
    molviewOptions.style.display = (settings.viewer3DSource === 'molview') ? 'block' : 'none';
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

  // Show SmilesDrawer, PubChem, or CDK Depict options based on selected engine
  if (smilesDrawerOptions) {
    smilesDrawerOptions.style.display = (settings.rendererEngine === 'client-side') ? 'block' : 'none';
  }
  if (pubchemOptions) {
    pubchemOptions.style.display = (settings.rendererEngine === 'pubchem') ? 'block' : 'none';
  }
  if (cdkDepictOptions) {
    cdkDepictOptions.style.display = (settings.rendererEngine === 'cdk-depict') ? 'block' : 'none';
  }
  // Show 3D Viewer Settings section for ALL engines (always visible)
  if (viewer3DSettings) {
    viewer3DSettings.style.display = 'block';
  }
  // Show MolView Protein Options for ALL engines (Search API detects proteins regardless of renderer)
  const molviewProteinOptionsInit = document.getElementById('molviewProteinOptions');
  if (molviewProteinOptionsInit) {
    molviewProteinOptionsInit.style.display = 'block';
  }
  // Show Mineral Options for ALL engines
  const mineralOptionsInit = document.getElementById('mineralOptions');
  if (mineralOptionsInit) {
    mineralOptionsInit.style.display = 'block';
  }
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

// SmilesDrawer rendering options
if (sdShowCarbonsToggle) {
  sdShowCarbonsToggle.addEventListener('change', (e) => {
    chrome.storage.sync.set({ sdShowCarbons: e.target.checked }, () => {
      showStatus('Show carbons ' + (e.target.checked ? 'enabled' : 'disabled') + '. Reload page to apply.', 'success');
    });
  });
}

if (sdAromaticRingsToggle) {
  sdAromaticRingsToggle.addEventListener('change', (e) => {
    chrome.storage.sync.set({ sdAromaticRings: e.target.checked }, () => {
      showStatus('Aromatic ring circles ' + (e.target.checked ? 'enabled' : 'disabled') + '. Reload page to apply.', 'success');
    });
  });
}

if (sdShowMethylsToggle) {
  sdShowMethylsToggle.addEventListener('change', (e) => {
    chrome.storage.sync.set({ sdShowMethyls: e.target.checked }, () => {
      showStatus('Show methyls ' + (e.target.checked ? 'enabled' : 'disabled') + '. Reload page to apply.', 'success');
    });
  });
}

if (sdAtomNumbersToggle) {
  sdAtomNumbersToggle.addEventListener('change', (e) => {
    chrome.storage.sync.set({ sdAtomNumbers: e.target.checked }, () => {
      showStatus('Atom numbers ' + (e.target.checked ? 'enabled' : 'disabled') + '. Reload page to apply.', 'success');
    });
  });
}

if (sdShowHydrogensToggle) {
  sdShowHydrogensToggle.addEventListener('change', (e) => {
    chrome.storage.sync.set({ sdShowHydrogens: e.target.checked }, () => {
      showStatus('Show hydrogens ' + (e.target.checked ? 'enabled' : 'disabled') + '. Reload page to apply.', 'success');
    });
  });
}

if (sdFlipHorizontalToggle) {
  sdFlipHorizontalToggle.addEventListener('change', (e) => {
    chrome.storage.sync.set({ sdFlipHorizontal: e.target.checked }, () => {
      showStatus('Horizontal flip ' + (e.target.checked ? 'enabled' : 'disabled') + '. Reload page to apply.', 'success');
    });
  });
}

if (sdFlipVerticalToggle) {
  sdFlipVerticalToggle.addEventListener('change', (e) => {
    chrome.storage.sync.set({ sdFlipVertical: e.target.checked }, () => {
      showStatus('Vertical flip ' + (e.target.checked ? 'enabled' : 'disabled') + '. Reload page to apply.', 'success');
    });
  });
}

if (sdThemeSelect) {
  sdThemeSelect.addEventListener('change', (e) => {
    chrome.storage.sync.set({ sdTheme: e.target.value }, () => {
      showStatus('Theme set to ' + e.target.value + '. Reload page to apply.', 'success');
    });
  });
}

if (sdRotateSlider) {
  sdRotateSlider.addEventListener('input', (e) => {
    sdRotateValue.textContent = e.target.value + '°';
  });

  sdRotateSlider.addEventListener('change', (e) => {
    chrome.storage.sync.set({ sdRotate: parseInt(e.target.value) }, () => {
      showStatus('Rotation set to ' + e.target.value + '°. Reload page to apply.', 'success');
    });
  });
}

if (sdBondThicknessSlider) {
  sdBondThicknessSlider.addEventListener('input', (e) => {
    sdBondThicknessValue.textContent = parseFloat(e.target.value).toFixed(1);
  });

  sdBondThicknessSlider.addEventListener('change', (e) => {
    chrome.storage.sync.set({ sdBondThickness: parseFloat(e.target.value) }, () => {
      showStatus('Bond thickness set to ' + e.target.value + '. Reload page to apply.', 'success');
    });
  });
}

if (sdGradientColorsToggle) {
  sdGradientColorsToggle.addEventListener('change', (e) => {
    chrome.storage.sync.set({ sdGradientColors: e.target.checked }, () => {
      showStatus('Gradient colors ' + (e.target.checked ? 'enabled' : 'disabled') + '. Reload page to apply.', 'success');
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

// PubChem Direct Fetch toggle
if (pubchemDirectFetchToggle) {
  pubchemDirectFetchToggle.addEventListener('change', (e) => {
    chrome.storage.sync.set({ pubchemDirectFetch: e.target.checked }, () => {
      showStatus('Direct fetch from PubChem ' + (e.target.checked ? 'enabled' : 'disabled') + '. Reload page to apply.', 'success');
    });
  });
}

// PubChem Remove Background toggle
if (pubchemRemoveBgToggle) {
  pubchemRemoveBgToggle.addEventListener('change', (e) => {
    chrome.storage.sync.set({ pubchemRemoveBg: e.target.checked }, () => {
      showStatus('Background removal ' + (e.target.checked ? 'enabled' : 'disabled') + '. Reload page to apply.', 'success');
    });
  });
}

// PubChem Sharpen Images toggle
if (pubchemSharpenToggle) {
  pubchemSharpenToggle.addEventListener('change', (e) => {
    chrome.storage.sync.set({ pubchemSharpenImages: e.target.checked }, () => {
      showStatus('Image sharpening ' + (e.target.checked ? 'enabled' : 'disabled') + '. Reload page to apply.', 'success');
    });
  });
}

// AI Molecular Control event listener
if (aiMolecularControlToggle) {
  aiMolecularControlToggle.addEventListener('change', (e) => {
    chrome.storage.sync.set({ enableAIMolecularControl: e.target.checked }, () => {
      showStatus('AI molecular control ' + (e.target.checked ? 'enabled' : 'disabled') + '. Reload page to apply.', 'success');
    });
  });
}

// CDK Depict options event listeners
if (cdkShowCarbonsToggle) {
  cdkShowCarbonsToggle.addEventListener('change', (e) => {
    chrome.storage.sync.set({ cdkShowCarbons: e.target.checked }, () => {
      showStatus('Show carbons ' + (e.target.checked ? 'enabled' : 'disabled') + '. Reload page to apply.', 'success');
    });
  });
}

if (cdkShowMethylsToggle) {
  cdkShowMethylsToggle.addEventListener('change', (e) => {
    chrome.storage.sync.set({ cdkShowMethyls: e.target.checked }, () => {
      showStatus('Show methyls ' + (e.target.checked ? 'enabled' : 'disabled') + '. Reload page to apply.', 'success');
    });
  });
}

if (cdkColorSchemeSelect) {
  cdkColorSchemeSelect.addEventListener('change', (e) => {
    chrome.storage.sync.set({ cdkColorScheme: e.target.value }, () => {
      showStatus('CDK color scheme set to ' + e.target.value + '. Reload page to apply.', 'success');
    });
  });
}

if (cdkHydrogenDisplaySelect) {
  cdkHydrogenDisplaySelect.addEventListener('change', (e) => {
    chrome.storage.sync.set({ cdkHydrogenDisplay: e.target.value }, () => {
      showStatus('Hydrogen display set to ' + e.target.value + '. Reload page to apply.', 'success');
    });
  });
}

if (cdkAtomNumbersToggle) {
  cdkAtomNumbersToggle.addEventListener('change', (e) => {
    chrome.storage.sync.set({ cdkAtomNumbers: e.target.checked }, () => {
      showStatus('Atom numbers ' + (e.target.checked ? 'enabled' : 'disabled') + '. Reload page to apply.', 'success');
    });
  });
}

if (cdkAnnotationSelect) {
  cdkAnnotationSelect.addEventListener('change', (e) => {
    chrome.storage.sync.set({ cdkAnnotation: e.target.value }, () => {
      showStatus('Annotation set to ' + e.target.value + '. Reload page to apply.', 'success');
    });
  });
}

if (cdkZoomSlider) {
  cdkZoomSlider.addEventListener('input', (e) => {
    const zoomValue = parseFloat(e.target.value);
    if (cdkZoomValue) cdkZoomValue.textContent = zoomValue + 'x';
  });

  cdkZoomSlider.addEventListener('change', (e) => {
    const zoomValue = parseFloat(e.target.value);
    chrome.storage.sync.set({ cdkZoom: zoomValue }, () => {
      showStatus('Zoom set to ' + zoomValue + 'x. Reload page to apply.', 'success');
    });
  });
}

// MolView-Only Mode event listener
if (molSearchModeToggle) {
  molSearchModeToggle.addEventListener('change', (e) => {
    chrome.storage.sync.set({ molSearchMode: e.target.checked }, () => {
      showStatus('MolView-Only Mode ' + (e.target.checked ? 'enabled' : 'disabled') + '. All data will fetch from localhost:8000. Reload page to apply.', 'success');
    });
  });
}

// Disable Formula Fallback event listener
if (disableFormulaFallbackToggle) {
  disableFormulaFallbackToggle.addEventListener('change', (e) => {
    chrome.storage.sync.set({ disableFormulaFallback: e.target.checked }, () => {
      showStatus('Formula Fallback ' + (e.target.checked ? 'disabled' : 'enabled') + '. Reload page to apply.', 'success');
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

if (viewer3DSourceSelect) {
  viewer3DSourceSelect.addEventListener('change', (e) => {
    const newValue = e.target.value;
    console.log('Saving viewer3DSource:', newValue);
    chrome.storage.sync.set({ viewer3DSource: newValue }, () => {
      if (chrome.runtime.lastError) {
        console.error('Error saving viewer3DSource:', chrome.runtime.lastError);
        showStatus('Error saving setting: ' + chrome.runtime.lastError.message, 'error');
        return;
      }
      const sourceNames = {
        '3dmol': '3Dmol.js',
        'molview': 'MolView (Local)',
        'pubchem': 'PubChem Official'
      };

      // Show/hide MolView specific options
      if (molviewOptions) {
        molviewOptions.style.display = (newValue === 'molview') ? 'block' : 'none';
      }
      console.log('Successfully saved viewer3DSource:', newValue);
      showStatus('3D viewer source set to ' + sourceNames[newValue] + '. Reload page to apply.', 'success');
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

if (viewer3DStyleSelect) {
  viewer3DStyleSelect.addEventListener('change', (e) => {
    const newStyle = e.target.value;

    // Load settings for the new style
    chrome.storage.sync.get(['viewer3DStyleSettings'], (result) => {
      const styleSettings = result.viewer3DStyleSettings || {};
      const newStyleSettings = styleSettings[newStyle] || {};

      // Update the UI with the saved settings for this style
      if (viewer3DStickRadiusSelect) {
        viewer3DStickRadiusSelect.value = newStyleSettings.stickRadius || '0.15';
      }
      if (viewer3DSphereRadiusSelect) {
        viewer3DSphereRadiusSelect.value = newStyleSettings.sphereRadius || '0.3';
      }

      // Save the new style selection
      chrome.storage.sync.set({ viewer3DStyle: newStyle }, () => {
        showStatus('3D style updated. Reload page to apply.', 'success');
      });
    });
  });
}

if (viewer3DStickRadiusSelect) {
  viewer3DStickRadiusSelect.addEventListener('change', (e) => {
    const newValue = e.target.value;

    // Get current style and update its settings
    chrome.storage.sync.get(['viewer3DStyle', 'viewer3DStyleSettings'], (result) => {
      const currentStyle = result.viewer3DStyle || 'stick:sphere';
      const styleSettings = result.viewer3DStyleSettings || {};

      // Update the settings for the current style
      if (!styleSettings[currentStyle]) {
        styleSettings[currentStyle] = {};
      }
      styleSettings[currentStyle].stickRadius = newValue;

      // Save back to storage
      chrome.storage.sync.set({ viewer3DStyleSettings: styleSettings }, () => {
        showStatus('Stick thickness updated. Reload page to apply.', 'success');
      });
    });
  });
}

if (viewer3DSphereRadiusSelect) {
  viewer3DSphereRadiusSelect.addEventListener('change', (e) => {
    const newValue = e.target.value;

    // Get current style and update its settings
    chrome.storage.sync.get(['viewer3DStyle', 'viewer3DStyleSettings'], (result) => {
      const currentStyle = result.viewer3DStyle || 'stick:sphere';
      const styleSettings = result.viewer3DStyleSettings || {};

      // Update the settings for the current style
      if (!styleSettings[currentStyle]) {
        styleSettings[currentStyle] = {};
      }
      styleSettings[currentStyle].sphereRadius = newValue;

      // Save back to storage
      chrome.storage.sync.set({ viewer3DStyleSettings: styleSettings }, () => {
        showStatus('Sphere size updated. Reload page to apply.', 'success');
      });
    });
  });
}

if (viewer3DAutoRotateToggle) {
  viewer3DAutoRotateToggle.addEventListener('change', (e) => {
    chrome.storage.sync.set({ viewer3DAutoRotate: e.target.checked }, () => {
      showStatus('Auto-rotate ' + (e.target.checked ? 'enabled' : 'disabled') + '. Reload page to apply.', 'success');
    });
  });
}

if (viewer3DSizeSelect) {
  viewer3DSizeSelect.addEventListener('change', (e) => {
    chrome.storage.sync.set({ viewer3DSize: e.target.value }, () => {
      showStatus('3D viewer size set to ' + e.target.value + '. Reload page to apply.', 'success');
    });
  });
}

if (viewer3DBgColorSelect) {
  viewer3DBgColorSelect.addEventListener('change', (e) => {
    chrome.storage.sync.set({ viewer3DBgColor: e.target.value }, () => {
      const colorName = e.target.options[e.target.selectedIndex].text;
      showStatus('3D background color set to ' + colorName + '. Reload page to apply.', 'success');
    });
  });
}

// MolView specific options event listeners
if (molviewRepresentationSelect) {
  molviewRepresentationSelect.addEventListener('change', (e) => {
    chrome.storage.sync.set({ molviewRepresentation: e.target.value }, () => {
      showStatus('Representation set to ' + e.target.options[e.target.selectedIndex].text + '. Reload page to apply.', 'success');
    });
  });
}

if (molviewEngineSelect) {
  molviewEngineSelect.addEventListener('change', (e) => {
    chrome.storage.sync.set({ molviewEngine: e.target.value }, () => {
      showStatus('Engine set to ' + e.target.options[e.target.selectedIndex].text + '. Reload page to apply.', 'success');
    });
  });
}

if (molviewCrystallographySelect) {
  molviewCrystallographySelect.addEventListener('change', (e) => {
    chrome.storage.sync.set({ molviewCrystallography: e.target.value }, () => {
      showStatus('Crystallography set to ' + e.target.options[e.target.selectedIndex].text + '. Reload page to apply.', 'success');
    });
  });
}

// Protein Options event listeners
if (proteinRemoveWhiteBgToggle) {
  proteinRemoveWhiteBgToggle.addEventListener('change', (e) => {
    chrome.storage.sync.set({ proteinRemoveWhiteBg: e.target.checked }, () => {
      showStatus('White background removal ' + (e.target.checked ? 'enabled' : 'disabled') + '. Reload page to apply.', 'success');
    });
  });
}

if (molviewBioAssemblyToggle) {
  molviewBioAssemblyToggle.addEventListener('change', (e) => {
    chrome.storage.sync.set({ molviewBioAssembly: e.target.checked }, () => {
      showStatus('Bio Assembly ' + (e.target.checked ? 'shown' : 'hidden') + '. Reload page to apply.', 'success');
    });
  });
}

if (molviewChainTypeSelect) {
  molviewChainTypeSelect.addEventListener('change', (e) => {
    chrome.storage.sync.set({ molviewChainType: e.target.value }, () => {
      showStatus('Chain Type set to ' + e.target.options[e.target.selectedIndex].text + '. Reload page to apply.', 'success');
    });
  });
}

if (molviewChainBondsToggle) {
  molviewChainBondsToggle.addEventListener('change', (e) => {
    chrome.storage.sync.set({ molviewChainBonds: e.target.checked }, () => {
      showStatus('Chain Bonds ' + (e.target.checked ? 'shown' : 'hidden') + '. Reload page to apply.', 'success');
    });
  });
}

if (molviewChainColorSelect) {
  molviewChainColorSelect.addEventListener('change', (e) => {
    chrome.storage.sync.set({ molviewChainColor: e.target.value }, () => {
      showStatus('Chain Color set to ' + e.target.options[e.target.selectedIndex].text + '. Reload page to apply.', 'success');
    });
  });
}

if (proteinMolviewBgColorSelect) {
  proteinMolviewBgColorSelect.addEventListener('change', (e) => {
    chrome.storage.sync.set({ proteinMolviewBgColor: e.target.value }, () => {
      showStatus('Protein background color set to ' + e.target.options[e.target.selectedIndex].text + '. Reload page to apply.', 'success');
    });
  });
}

// Mineral Options event listeners
if (mineralRepresentationSelect) {
  mineralRepresentationSelect.addEventListener('change', (e) => {
    chrome.storage.sync.set({ mineralRepresentation: e.target.value }, () => {
      showStatus('Mineral representation set to ' + e.target.options[e.target.selectedIndex].text + '. Reload page to apply.', 'success');
    });
  });
}

if (mineralMolviewBgColorSelect) {
  mineralMolviewBgColorSelect.addEventListener('change', (e) => {
    chrome.storage.sync.set({ mineralMolviewBgColor: e.target.value }, () => {
      showStatus('Mineral background color set to ' + e.target.options[e.target.selectedIndex].text + '. Reload page to apply.', 'success');
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
const molviewProteinOptions = document.getElementById('molviewProteinOptions');
const mineralOptions = document.getElementById('mineralOptions');

function updateEngineInfo(engine) {
  if (engine === 'client-side') {
    if (engineInfo) engineInfo.textContent = '💻 SmilesDrawer - Fast Client-Side Rendering (No Server!)';
    if (smilesDrawerOptions) smilesDrawerOptions.style.display = 'block';
    if (pubchemOptions) pubchemOptions.style.display = 'none';
    if (cdkDepictOptions) cdkDepictOptions.style.display = 'none';
  } else if (engine === 'pubchem') {
    if (engineInfo) engineInfo.textContent = '🌐 PubChem - Official Chemical Images';
    if (smilesDrawerOptions) smilesDrawerOptions.style.display = 'none';
    if (pubchemOptions) pubchemOptions.style.display = 'block';
    if (cdkDepictOptions) cdkDepictOptions.style.display = 'none';
  } else if (engine === 'cdk-depict') {
    if (engineInfo) engineInfo.textContent = '🌐 CDK Depict - Free Online API';
    if (smilesDrawerOptions) smilesDrawerOptions.style.display = 'none';
    if (pubchemOptions) pubchemOptions.style.display = 'none';
    if (cdkDepictOptions) cdkDepictOptions.style.display = 'block';
  }

  // Always show 3D Viewer Settings, Protein Options, and Mineral Options (independent of compound renderer)
  if (viewer3DSettings) viewer3DSettings.style.display = 'block';
  if (molviewProteinOptions) molviewProteinOptions.style.display = 'block';
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
  uiBlur: 10,
  uiOpacity: 45
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
  moleculeCount: 20,
  moleculeScale: 0.5
}, function (settings) {
  if (molCountSlider) {
    molCountSlider.value = settings.moleculeCount;
    molCountValue.textContent = settings.moleculeCount;
  }
  if (molScaleSlider) {
    molScaleSlider.value = settings.moleculeScale;
    molScaleValue.textContent = settings.moleculeScale;
  }
});

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
