/**
 * ChemistryLaTeX - Chemistry Formula Renderer v6.0
 * Works on ChatGPT and other strict CSP sites
 * 
 * SERVER-BASED ARCHITECTURE:
 * - SVG rendering via ChemistryLaTeX Server
 * - Biomolecule/Mineral lookups via server
 * - Extension handles 3D viewers
 * - Settings apply instantly via URL params
 */


// ============================================
// CHEMISTRYLATEX SERVER CONFIGURATION
// ============================================
// NOTE: This should match your deployed Vercel domain for the SVG renderer.
// Example: https://server-chemistryrenderer.vercel.app
const CHEMISTRYLATEX_SERVER_URL = 'https://server-chemistryrenderer.vercel.app';

// Build server URL for SVG rendering
// Cache version - increment when server SVG format changes
// v6: Marker-based theming - server returns grayscale, extension applies colors
const CHEMISTRYLATEX_CACHE_VERSION = 'v6';

// Cache bust timestamp - set when user clicks "Clear Cache"
// This forces all new requests to bypass browser cache
let cacheBustTimestamp = '';

function buildServerSvgUrl(type, value, options = {}) {
  // Types that support rendering flags (mol, smiles, iupac)
  // Minerals and biomolecules do NOT support flags
  const typesWithFlags = ['mol', 'smiles', 'iupac'];
  const supportsFlags = typesWithFlags.includes(type.toLowerCase());

  // Build flags array based on options (only for types that support them)
  let valueWithFlags = value;

  if (supportsFlags) {
    const flags = [];

    // Rendering flags (lowercase letters)
    if (options.showCarbons) flags.push('c');
    if (options.aromaticCircles) flags.push('o');
    if (options.showMethyls) flags.push('m');
    if (options.atomNumbers) flags.push('n');
    if (options.showHydrogens) flags.push('h');
    if (options.showImplicitH) flags.push('i');
    if (options.gradientColors) flags.push('g');
    if (options.useStereochemistry) flags.push('s');
    if (options.compactDrawing) flags.push('k');

    // Flip flags
    if (options.flipHorizontal) flags.push('p');
    if (options.flipVertical) flags.push('q');

    // IMPORTANT: Sort flags alphabetically for consistent cache keys
    // This prevents duplicate cache entries like mol=benzene+c+o vs mol=benzene+o+c
    flags.sort();

    // Build the value with sorted flags appended
    if (flags.length > 0) {
      valueWithFlags = `${value}+${flags.join('+')}`;
    }
  }

  // Query params - cache version and optional cache bust timestamp
  const params = new URLSearchParams();
  params.set('_v', CHEMISTRYLATEX_CACHE_VERSION);

  // Add cache bust timestamp if set (forces fresh content after clearing cache)
  if (cacheBustTimestamp) {
    params.set('_t', cacheBustTimestamp);
  }

  // Determine renderer prefix based on renderingEngine option
  // 'rdkit' -> /rdkit/mol=...
  // 'kekule' -> /kekule/mol=...
  // 'chemistrylatex' or default -> /mol=...
  const rendererPrefixes = { rdkit: '/rdkit', kekule: '/kekule' };
  const rendererPrefix = rendererPrefixes[options.renderingEngine] || '';

  const queryString = params.toString();
  const url = `${CHEMISTRYLATEX_SERVER_URL}${rendererPrefix}/${type}=${encodeURIComponent(valueWithFlags)}.svg${queryString ? '?' + queryString : ''}`;
  return url;
}


// Fetch biomolecule/mineral data from server (returns pdbid/codid for 3D viewers)
async function fetchFromChemistryLaTeXServer(type, name) {
  const url = `${CHEMISTRYLATEX_SERVER_URL}/${type}=${encodeURIComponent(name)}.json`;
  console.log('%c🌐 Fetching from ChemistryLaTeX Server:', 'color: #2196F3; font-weight: bold;', url);

  const response = await fetch(url);
  if (!response.ok) {
    throw new Error(`Server returned ${response.status}`);
  }
  return response.json();
}

// ============================================
// MARKER-BASED THEMING SYSTEM
// ============================================
// Server renders SVGs with unique gray "marker" colors for each element type.
// Extension replaces these markers with actual theme colors.
// This allows ONE cached SVG to work for ALL themes (6x cache efficiency!)

// Marker colors used by server (unique gray shades)
const ELEMENT_MARKERS = {
  C: '#101010',    // Carbon (darkest)
  H: '#181818',    // Hydrogen
  N: '#202020',    // Nitrogen
  O: '#282828',    // Oxygen
  S: '#303030',    // Sulfur
  P: '#383838',    // Phosphorus
  F: '#404040',    // Fluorine
  CL: '#484848',   // Chlorine
  BR: '#505050',   // Bromine
  I: '#585858',    // Iodine
  B: '#606060',    // Boron
  SI: '#686868',   // Silicon
  ATOM_NUMBER: '#707070'  // Atom numbering
};

// Theme color definitions (actual display colors)
const THEME_COLORS = {
  light: {
    C: '#222222', H: '#666666', N: '#3498db', O: '#e74c3c',
    S: '#f1c40f', P: '#d35400', F: '#27ae60', CL: '#16a085',
    BR: '#d35400', I: '#8e44ad', B: '#e67e22', SI: '#e67e22',
    ATOM_NUMBER: '#2196F3'
  },
  dark: {
    C: '#ffffff', H: '#aaaaaa', N: '#4dabf7', O: '#ff6b6b',
    S: '#fcc419', P: '#fd7e14', F: '#51cf66', CL: '#20c997',
    BR: '#fd7e14', I: '#be4bdb', B: '#f59f00', SI: '#f59f00',
    ATOM_NUMBER: '#ffd700'
  },
  classic_black: {
    // Pure black - classic chemistry textbook style (for white backgrounds)
    C: '#000000', H: '#000000', N: '#000000', O: '#000000',
    S: '#000000', P: '#000000', F: '#000000', CL: '#000000',
    BR: '#000000', I: '#000000', B: '#000000', SI: '#000000',
    ATOM_NUMBER: '#000000'
  },
  classic_white: {
    // Pure white - for dark backgrounds
    C: '#ffffff', H: '#ffffff', N: '#ffffff', O: '#ffffff',
    S: '#ffffff', P: '#ffffff', F: '#ffffff', CL: '#ffffff',
    BR: '#ffffff', I: '#ffffff', B: '#ffffff', SI: '#ffffff',
    ATOM_NUMBER: '#ffffff'
  },
  warm_paper: {
    C: '#586e75', H: '#657b83', N: '#268bd2', O: '#dc322f',
    S: '#b58900', P: '#d33682', F: '#859900', CL: '#16a085',
    BR: '#cb4b16', I: '#6c71c4', B: '#2aa198', SI: '#2aa198',
    ATOM_NUMBER: '#2196F3'
  },
  terminal_green: {
    // Refined Matrix theme - professional neon green on dark
    C: '#00FF41', H: '#008F11', N: '#003B00', O: '#0D0208',
    S: '#00FF41', P: '#00FF41', F: '#00FF41', CL: '#00FF41',
    BR: '#00FF41', I: '#00FF41', B: '#00FF41', SI: '#00FF41',
    ATOM_NUMBER: '#00FF41'
  },
  modern_dark: {
    C: '#c9d1d9', H: '#8b949e', N: '#58a6ff', O: '#f85149',
    S: '#d29922', P: '#db6d28', F: '#3fb950', CL: '#3fb950',
    BR: '#db6d28', I: '#a371f7', B: '#db6d28', SI: '#db6d28',
    ATOM_NUMBER: '#58a6ff'
  },
  carbon: {
    C: '#f4f4f4', H: '#a8a8a8', N: '#4589ff', O: '#fa4d56',
    S: '#f1c21b', P: '#ff832b', F: '#42be65', CL: '#42be65',
    BR: '#ff832b', I: '#a56eff', B: '#ff832b', SI: '#ff832b',
    ATOM_NUMBER: '#4589ff'
  }
};

/**
 * Apply theme colors to an SVG by replacing marker colors
 * @param {string} svgText - Raw SVG text from server
 * @param {string} themeName - Theme to apply (light, dark, matrix, etc.)
 * @returns {string} - SVG with theme colors applied
 */
function applyThemeColors(svgText, themeName) {
  // NO AUTO-ADAPT - user manually selects theme
  const themeColors = THEME_COLORS[themeName] || THEME_COLORS.light;
  let coloredSvg = svgText;

  // Replace each marker color with the theme color
  for (const [element, markerColor] of Object.entries(ELEMENT_MARKERS)) {
    const themeColor = themeColors[element] || markerColor;
    // Replace in fill, stroke, and stop-color attributes
    // Case-insensitive replacement for the hex color
    const markerRegex = new RegExp(markerColor, 'gi');
    coloredSvg = coloredSvg.replace(markerRegex, themeColor);
  }

  console.log('%c🎨 Applied theme colors:', 'color: #9C27B0;', themeName);
  return coloredSvg;
}

// ============================================
// GLOBAL ERROR HANDLING - Prevent crashes
// ============================================
// These handlers catch unhandled errors that could cause "Aw, Snap!" crashes

window.addEventListener('error', function (event) {
  // Check if error is from our extension
  if (event.filename && event.filename.includes('chrome-extension://')) {
    console.error('[ChemRenderer] 🛡️ Caught unhandled error:', event.message);
    // Prevent the error from propagating and crashing the page
    event.preventDefault();
    return true;
  }
});

window.addEventListener('unhandledrejection', function (event) {
  // Catch unhandled promise rejections
  const reason = event.reason;
  if (reason && reason.stack && reason.stack.includes('ChemRenderer')) {
    console.error('[ChemRenderer] 🛡️ Caught unhandled promise rejection:', reason);
    // Prevent the rejection from propagating
    event.preventDefault();
    return true;
  }
  // Also catch extension context invalidated errors
  if (reason && reason.message && reason.message.includes('Extension context invalidated')) {
    console.warn('[ChemRenderer] Extension context invalidated, ignoring...');
    event.preventDefault();
    return true;
  }
});

// ============================================
// LOGGING UTILITIES
// ============================================
const LOG_PREFIX = '[ChemRenderer]';
const log = {
  info: (msg, data) => console.log(`${LOG_PREFIX} [INFO] ${msg}`, data || ''),
  success: (msg, data) => console.log(`%c${LOG_PREFIX} [SUCCESS] ${msg}`, 'color: green; font-weight: bold;', data || ''),
  warn: (msg, data) => console.warn(`${LOG_PREFIX} [WARN] ${msg}`, data || ''),
  error: (msg, data) => console.error(`${LOG_PREFIX} [ERROR] ${msg}`, data || ''),
  debug: (msg, data) => console.debug(`${LOG_PREFIX} [DEBUG] ${msg}`, data || ''),
  inject: (msg, data) => console.log(`%c${LOG_PREFIX} [INJECT] ${msg}`, 'color: blue; font-weight: bold;', data || '')
};

let logHistory = [];
let extensionContextInvalid = false;


// ============================================
// 🧠 SMART CACHING SYSTEM
// ============================================
// Two-layer caching for performance:
// 1. SMILES Cache (persistent): compound name → { canonicalSmiles, isomericSmiles, cid, pdbid, codid, compoundType }
//    - canonicalSmiles: SMILES without stereochemistry (used when stereochemistry is OFF)
//    - isomericSmiles: SMILES with stereochemistry markers (used when stereochemistry is ON)
// 2. Image Cache (session): smiles + options hash → SVG/image data

/**
 * In-memory rendered image cache
 * Key: generateImageCacheKey(smiles, options)
 * Value: SVG string or image data URL
 */
const renderedImageCache = new Map();
const MAX_IMAGE_CACHE_SIZE = 100; // Limit to prevent memory bloat

/**
 * Generate a unique cache key for rendered images based on SMILES and rendering options
 */
function generateImageCacheKey(smiles, options = {}) {
  const optHash = [
    options.showCarbons ? 'c1' : 'c0',
    options.showMethyls ? 'm1' : 'm0',
    options.aromaticCircles ? 'a1' : 'a0',
    options.showHydrogens ? 'h1' : 'h0',
    options.atomNumbers ? 'n1' : 'n0',
    options.theme || 'light',
    options.width || 300,
    options.height || 240
  ].join('_');
  return `${smiles}__${optHash}`;
}

// Note: getCachedImage/setCachedImage removed - server handles all caching now

/**
 * SMILES Cache: Persistent storage for compound name → SMILES mappings
 * Stored in chrome.storage.local for persistence across sessions
 */
const SMILES_CACHE_KEY = 'chemistrylatex_smiles_cache';
const SMILES_CACHE_MAX_SIZE = 500; // Max entries

// In-memory mirror of SMILES cache for fast access
let smilesCache = {};
let smilesCacheLoaded = false;

/**
 * Load SMILES cache from storage into memory
 */
async function loadSmilesCache() {
  if (smilesCacheLoaded) return smilesCache;

  try {
    const result = await chrome.storage.local.get(SMILES_CACHE_KEY);
    smilesCache = result[SMILES_CACHE_KEY] || {};
    smilesCacheLoaded = true;
    log.info(`📦 Loaded SMILES cache: ${Object.keys(smilesCache).length} entries`);
  } catch (e) {
    log.warn('Failed to load SMILES cache:', e);
    smilesCache = {};
  }
  return smilesCache;
}

/**
 * Save SMILES cache to storage
 */
async function saveSmilesCache() {
  try {
    // Enforce size limit
    const keys = Object.keys(smilesCache);
    if (keys.length > SMILES_CACHE_MAX_SIZE) {
      // Remove oldest entries (first N keys to get under limit)
      const toRemove = keys.length - SMILES_CACHE_MAX_SIZE;
      keys.slice(0, toRemove).forEach(k => delete smilesCache[k]);
    }
    await chrome.storage.local.set({ [SMILES_CACHE_KEY]: smilesCache });
  } catch (e) {
    log.warn('Failed to save SMILES cache:', e);
  }
}
/**
 * Get cached SMILES data for a compound name
 * @param {string} name - Compound name (case-insensitive)
 * @returns {object|null} - { canonicalSmiles, isomericSmiles, cid, pdbid, codid, compoundType } or null
 */
async function getCachedSmiles(name) {
  await loadSmilesCache();
  const key = name.toLowerCase().trim();
  const cached = smilesCache[key];
  if (cached) {
    log.debug(`🎯 SMILES cache HIT: ${name}`);
    // Migrate old cache format (single 'smiles' field) to new format
    if (cached.smiles && !cached.canonicalSmiles && !cached.isomericSmiles) {
      log.debug(`🔄 Migrating old cache format for: ${name}`);
      // Old format had isomeric SMILES stored in 'smiles' field
      cached.isomericSmiles = cached.smiles;
      cached.canonicalSmiles = stripStereochemistry(cached.smiles);
      delete cached.smiles;
      // Save the migrated data
      smilesCache[key] = cached;
      clearTimeout(setCachedSmiles._saveTimeout);
      setCachedSmiles._saveTimeout = setTimeout(() => saveSmilesCache(), 1000);
    }
    return cached;
  }

  // Fallback: Static cache files have been moved to the server
  // The server now checks MineralNames, CompoundNames, BiomoleculeNames
  // before calling external APIs

  return null;
}

/**
 * Store SMILES data in cache
 * @param {string} name - Compound name
 * @param {object} data - { canonicalSmiles, isomericSmiles, cid, pdbid, codid, compoundType }
 */
async function setCachedSmiles(name, data) {
  await loadSmilesCache();
  const key = name.toLowerCase().trim();
  smilesCache[key] = {
    ...data,
    cachedAt: Date.now()
  };
  // Debounce save to avoid excessive writes
  clearTimeout(setCachedSmiles._saveTimeout);
  setCachedSmiles._saveTimeout = setTimeout(() => saveSmilesCache(), 1000);
  const displaySmiles = data.isomericSmiles || data.canonicalSmiles || 'N/A';
  log.debug(`📦 Cached SMILES: ${name} → canonical: ${data.canonicalSmiles?.substring(0, 25) || 'N/A'}, isomeric: ${data.isomericSmiles?.substring(0, 25) || 'N/A'}`);

  // Static cache files (MineralNames, CompoundNames, BiomoleculeNames) are now on the server
  // The server handles caching - no client-side static cache updates needed
}

/**
 * In-memory map of pending searches to deduplicate simultaneous requests
 * If "benzene" appears twice on page, both will wait for single API call
 * Key: compound name (lowercase), Value: Promise that resolves to search result
 */
const pendingSearches = new Map();
const MAX_PENDING_SEARCHES = 50;  // Safety limit to prevent memory bloat

/**
 * Get or create a pending search for a compound name
 * Ensures only ONE API call is made per unique compound name
 * @param {string} name - Compound name
 * @param {function} searchFn - Async function to call if no pending search exists
 * @returns {Promise} - Promise that resolves to search result
 */
async function getOrCreatePendingSearch(name, searchFn) {
  const key = name.toLowerCase().trim();

  // If search is already in progress, return existing promise
  if (pendingSearches.has(key)) {
    log.debug(`🔄 Reusing pending search for: ${name}`);
    return pendingSearches.get(key);
  }

  // SAFETY: Limit pending searches to prevent memory issues
  if (pendingSearches.size >= MAX_PENDING_SEARCHES) {
    log.warn(`⚠️ Too many pending searches (${pendingSearches.size}), clearing old ones...`);
    // Clear the oldest entries (first half)
    const keysToRemove = Array.from(pendingSearches.keys()).slice(0, Math.floor(MAX_PENDING_SEARCHES / 2));
    keysToRemove.forEach(k => pendingSearches.delete(k));
  }

  // Create new search promise
  const searchPromise = searchFn().finally(() => {
    // Clean up after search completes (success or failure)
    pendingSearches.delete(key);
  });

  pendingSearches.set(key, searchPromise);
  log.debug(`🆕 Started new search for: ${name}`);
  return searchPromise;
}

/**
 * Clear all caches and force fresh content
 * Sets cache bust timestamp so all new SVG requests bypass browser cache
 */
function clearAllCaches() {
  // Clear in-memory caches (legacy - not really used anymore)
  renderedImageCache.clear();
  smilesCache = {};
  smilesCacheLoaded = false;
  pendingSearches.clear();

  // Clear chrome storage
  chrome.storage.local.remove(SMILES_CACHE_KEY);

  // SET CACHE BUST TIMESTAMP - this is the key!
  // All new SVG URLs will include this timestamp, bypassing browser cache
  cacheBustTimestamp = Date.now().toString();

  log.info(`🗑️ Cache cleared! New cache bust timestamp: ${cacheBustTimestamp}`);
  console.log('%c🗑️ CACHE CLEARED - All new SVG requests will bypass browser cache', 'background: #F44336; color: white; font-weight: bold; padding: 8px;');
}

/**
 * Strip stereochemistry markers from SMILES
 * Removes @ and @@ (chiral centers) and / \ (double bond geometry)
 * @param {string} smiles - The SMILES string (may contain stereochemistry)
 * @returns {string} - SMILES without stereochemistry markers
 */
function stripStereochemistry(smiles) {
  if (!smiles) return smiles;
  // Remove chiral center markers: @ and @@
  // Remove double bond geometry: / and \
  return smiles.replace(/@@|@|\/|\\/g, '');
}


// ============================================
// 🔄 INSTANT SETTINGS APPLICATION
// ============================================
// Listen for settings changes from popup and re-render all molecules

/**
 * Track all rendered molecule containers for re-rendering
 */
const renderedMolecules = new Map(); // Map<element, moleculeData>

/**
 * Determine which image types are affected by changed settings
 * @param {Object} changedSettings - Settings that were changed
 * @returns {Set<string>} - Set of affected types: 'compound', 'biomolecule', 'mineral'
 */
function getAffectedImageTypes(changedSettings) {
  const affectedTypes = new Set();

  // Rendering settings only affect compounds (small molecules)
  const renderSettings = [
    'sdShowCarbons', 'sdAromaticRings', 'sdShowMethyls', 'sdAtomNumbers',
    'sdShowExplicitHydrogens', 'sdShowImplicitHydrogens', 'sdFlipHorizontal', 'sdFlipVertical', 'sdTheme',
    'sdAutoAdapt',  // Auto-adapt affects theme selection based on page dark/light mode
    'sdRotate', 'sdBondThickness', 'sdGradientColors', 'sdScaleByWeight',
    'm2cfShowCarbons', 'm2cfAromaticCircles', 'm2cfShowMethyls', 'm2cfAtomNumbers',
    'm2cfAddH2', 'm2cfShowImplicitH', 'm2cfFlipHorizontal', 'm2cfFlipVertical'
  ];

  // Compound MolView 3D settings (affects 3D views of compounds)
  const compoundMolviewSettings = [
    'molviewRepresentation', 'compoundMolviewBgColor', 'molviewEngine'
  ];

  // Protein/biomolecule specific settings
  const proteinSettings = [
    'proteinRemoveWhiteBg', 'viewer3DStyle', 'viewer3DAutoRotate',
    'viewer3DBgColor', 'viewer3DSize',
    // MolView-specific protein settings
    'molviewChainType', 'molviewChainBonds', 'molviewChainColor',
    'proteinMolviewBgColor', 'molviewBioAssembly'
  ];

  // Mineral-specific settings
  const mineralSettings = [
    'mineralRepresentation', 'mineralMolviewBgColor', 'mineralCrystallography'
  ];

  // Settings that affect ALL types
  const universalSettings = ['performanceMode', 'enableAIFlagControl'];

  // Check what's changed
  const changedKeys = Object.keys(changedSettings);

  for (const key of changedKeys) {
    if (renderSettings.includes(key) || compoundMolviewSettings.includes(key)) {
      affectedTypes.add('compound');
    }
    if (proteinSettings.includes(key)) {
      affectedTypes.add('biomolecule');
    }
    if (mineralSettings.includes(key)) {
      affectedTypes.add('mineral');
    }
    if (universalSettings.includes(key)) {
      affectedTypes.add('compound');
      affectedTypes.add('biomolecule');
      affectedTypes.add('mineral');
    }
  }

  // If no specific types were affected, assume all types (safety fallback)
  if (affectedTypes.size === 0) {
    affectedTypes.add('compound');
    affectedTypes.add('biomolecule');
    affectedTypes.add('mineral');
  }

  return affectedTypes;
}

/**
 * Re-render all molecules with new settings (without page reload)
 * Now delegates to lazyReRenderMolecules to always use lazy loading.
 * Old images stay visible with a loading spinner until the new one is ready.
 */
async function reRenderAllMolecules(newSettings) {
  log.info('🔄 Re-rendering all molecules with lazy loading...');
  // Always use lazy re-rendering to keep old images visible
  return lazyReRenderMolecules(newSettings);
}

/**
 * Mark molecules that need re-rendering (for lazy re-render)
 */
let pendingReRenderSettings = null;
let pendingUpdateIndicator = null;

// Flag to prevent mutation observer from re-scanning during image updates
// This is CRITICAL - without it, the observer would trigger scanAndRender()
// which would recreate images and cause them to disappear
let isUpdatingImages = false;

/**
 * Add a small red loading indicator to an image that needs re-rendering
 */
function addLoadingIndicator(img) {
  // Don't add if already exists
  if (img.dataset.hasLoadingIndicator === 'true') return;

  // Create a small red spinner indicator (left-middle)
  const indicator = document.createElement('div');
  indicator.className = 'chemistrylatex-image-loading-indicator';
  indicator.setAttribute('aria-label', 'Updating image');
  indicator.style.cssText = `
    position: absolute;
    top: 50%;
    left: 10px;
    transform: translateY(-50%);
    width: 14px;
    height: 14px;
    border-radius: 50%;
    border: 2px solid rgba(255, 71, 87, 0.25);
    border-top-color: #ff4757;
    z-index: 10;
    animation: chemtex-spin 0.9s linear infinite;
    pointer-events: none;
  `;

  // Position the image container relatively if needed
  const container = img.parentElement;
  if (container) {
    const computedStyle = window.getComputedStyle(container);
    if (computedStyle.position === 'static') {
      container.style.position = 'relative';
    }
    container.appendChild(indicator);
  }

  // Add animations if not already added
  if (!document.getElementById('chemtex-loading-indicator-style')) {
    const style = document.createElement('style');
    style.id = 'chemtex-loading-indicator-style';
    style.textContent = `
      @keyframes chemtex-spin {
        from { transform: translateY(-50%) rotate(0deg); }
        to { transform: translateY(-50%) rotate(360deg); }
      }
      @keyframes chemtex-fade-in {
        from { opacity: 0.65; }
        to { opacity: 1; }
      }
    `;
    document.head.appendChild(style);
  }

  img.dataset.hasLoadingIndicator = 'true';
  img.dataset.loadingIndicatorElement = 'true';
}

/**
 * Remove the loading indicator from an image
 */
function removeLoadingIndicator(img) {
  if (img.dataset.hasLoadingIndicator !== 'true') return;

  const container = img.parentElement;
  if (container) {
    const indicator = container.querySelector('.chemistrylatex-image-loading-indicator');
    if (indicator) {
      indicator.remove();
    }
  }

  delete img.dataset.hasLoadingIndicator;
  delete img.dataset.loadingIndicatorElement;
}


/**
 * Apply average size scaling to all 2D molecule images
 * This uses CSS transform to scale images in real-time without re-rendering
 * @param {number} sizePercent - Scale percentage (50-200, 100 = normal)
 */
function applyAverageSizeScaling(sizePercent) {
  const scaleFactor = sizePercent / 100;
  console.log(`%c📐 Applying average size scaling: ${sizePercent}% (factor: ${scaleFactor})`, 'color: #9C27B0; font-weight: bold;');

  // Find all molecule images
  const moleculeImages = document.querySelectorAll('img.molecule-diagram, img.molecule-viewer, img[data-molecule-viewer]');

  let scaledCount = 0;
  moleculeImages.forEach(img => {
    // Skip 3D viewer iframes
    if (img.closest('.molecule-3d-viewer')) return;
    if (img.closest('.molecule-viewer-container')) return;
    if (img.closest('.biomolecule-viewer')) return;

    // ONLY apply to SVG images (compounds and minerals)
    // Skip biomolecules, structure images, and any non-SVG sources
    const compoundType = img.dataset.compoundType || '';
    if (compoundType === 'biomolecule') return;

    // Also check if it's an SVG by looking at the src
    const src = img.src || '';
    const isSvg = src.includes('data:image/svg') || src.endsWith('.svg');
    if (!isSvg) return;

    // Find the wrapper container (could be molecule-container or chem-image-container)
    const wrapper = img.closest('.molecule-container, .chem-image-container');

    if (wrapper) {
      // Get original dimensions from the image
      const originalWidth = img.naturalWidth || parseInt(img.style.width) || img.offsetWidth || 300;
      const originalHeight = img.naturalHeight || parseInt(img.style.height) || img.offsetHeight || 240;
      const scaledWidth = Math.round(originalWidth * scaleFactor);
      const scaledHeight = Math.round(originalHeight * scaleFactor);

      // Set BOTH the wrapper size AND image size to scaled dimensions
      // This makes the container expand (affecting document flow/bubbles)
      // AND keeps buttons positioned correctly since wrapper size matches content
      wrapper.style.width = scaledWidth + 'px';
      wrapper.style.height = scaledHeight + 'px';

      // Scale the image to fill the wrapper
      img.style.width = '100%';
      img.style.height = '100%';
      img.style.transform = '';  // Clear any previous transforms
    } else {
      // No wrapper - set image size directly
      const originalWidth = img.naturalWidth || parseInt(img.style.width) || img.offsetWidth || 300;
      const originalHeight = img.naturalHeight || parseInt(img.style.height) || img.offsetHeight || 240;
      img.style.width = Math.round(originalWidth * scaleFactor) + 'px';
      img.style.height = Math.round(originalHeight * scaleFactor) + 'px';
    }

    scaledCount++;
  });

  console.log(`%c✅ Applied scaling to ${scaledCount} SVG molecule images (skipped ${moleculeImages.length - scaledCount} non-SVG/biomolecules)`, 'color: #4CAF50;');
}


/**
 * Re-render ALL molecules with new settings using the proven rendering pipeline.
 * Instead of modifying img.src directly (which fails on React pages),
 * we reset the image to trigger the original loadMoleculeImage flow.
 */
async function lazyReRenderMolecules(newSettings) {
  console.warn('[ChemRenderer] 🔄 RE-RENDERING with new settings using loadMoleculeImage pipeline');

  // CRITICAL: Block mutation observer from re-scanning while we update images
  isUpdatingImages = true;

  // Update local settings object
  Object.assign(settings, newSettings);

  // Map popup settings to internal options
  if (newSettings.sdShowCarbons !== undefined) settings.m2cfShowCarbons = newSettings.sdShowCarbons;
  if (newSettings.sdAromaticRings !== undefined) settings.m2cfAromaticCircles = newSettings.sdAromaticRings;
  if (newSettings.sdShowMethyls !== undefined) settings.m2cfShowMethyls = newSettings.sdShowMethyls;
  if (newSettings.sdAtomNumbers !== undefined) settings.m2cfAtomNumbers = newSettings.sdAtomNumbers;
  if (newSettings.sdShowExplicitHydrogens !== undefined) settings.m2cfAddH2 = newSettings.sdShowExplicitHydrogens;
  if (newSettings.sdShowImplicitHydrogens !== undefined) settings.m2cfShowImplicitH = newSettings.sdShowImplicitHydrogens;
  if (newSettings.sdCompactDrawing !== undefined) settings.m2cfCompactDrawing = newSettings.sdCompactDrawing;
  if (newSettings.sdFlipHorizontal !== undefined) settings.m2cfFlipHorizontal = newSettings.sdFlipHorizontal;
  if (newSettings.sdFlipVertical !== undefined) settings.m2cfFlipVertical = newSettings.sdFlipVertical;
  if (newSettings.sdGradientColors !== undefined) settings.sdGradientColors = newSettings.sdGradientColors;
  // Theme and bond thickness
  if (newSettings.sdTheme !== undefined) settings.sdTheme = newSettings.sdTheme;
  if (newSettings.sdAutoAdapt !== undefined) {
    settings.sdAutoAdapt = newSettings.sdAutoAdapt;
    console.log('[ChemRenderer] sdAutoAdapt setting updated:', settings.sdAutoAdapt);
  }
  if (newSettings.sdBondThickness !== undefined) settings.sdBondThickness = newSettings.sdBondThickness;
  if (newSettings.sdRotate !== undefined) settings.m2cfRotate = newSettings.sdRotate;
  // Stereochemistry toggle - affects which SMILES type is fetched
  if (newSettings.useStereochemistry !== undefined) {
    settings.m2cfUse3DSmiles = newSettings.useStereochemistry;
    console.log('[ChemRenderer] Stereochemistry setting updated:', settings.m2cfUse3DSmiles);
  }


  // Handle sdAverageSize setting - apply CSS transform scaling to all molecules
  if (newSettings.sdAverageSize !== undefined) {
    settings.sdAverageSize = newSettings.sdAverageSize;
    applyAverageSizeScaling(newSettings.sdAverageSize);
    // Don't re-render, just scale existing images
    isUpdatingImages = false;
    return;
  }

  // Handle compound3DSize and mineral3DSize settings
  if (newSettings.compound3DSize !== undefined) {
    settings.compound3DSize = newSettings.compound3DSize;
  }
  if (newSettings.mineral3DSize !== undefined) {
    settings.mineral3DSize = newSettings.mineral3DSize;
  }

  // Find all molecule images (compounds, biomolecules, and minerals)
  const moleculeImages = document.querySelectorAll('img.molecule-diagram, img.molecule-viewer, img[data-molecule-viewer]');
  console.log('[ChemRenderer] 🔍 Found molecule images to re-render:', moleculeImages.length);

  if (moleculeImages.length === 0) {
    console.log('[ChemRenderer] No images found');
    isUpdatingImages = false;
    return;
  }

  let successCount = 0;

  for (const img of moleculeImages) {
    try {
      const smiles = img.dataset.smiles || '';  // Optional now - server can look up by name
      // CRITICAL: Use original nomenclature first, then dataset.nomenclature, NEVER use alt (it may contain error messages)
      const originalNomenclature = img.dataset.originalNomenclature;
      const datasetNomenclature = img.dataset.nomenclature;

      // Validate that nomenclature doesn't contain error messages
      const isValidNomenclature = (name) => {
        if (!name || name === 'molecule') return false;
        // Reject if it looks like an error message
        if (name.includes('Error loading') || name.includes('Server returned') || name.includes('Failed to')) return false;
        return true;
      };

      // Prioritize original nomenclature, then dataset, skip alt entirely
      let nomenclature = originalNomenclature || datasetNomenclature || '';

      // If nomenclature is invalid (error message), skip this image
      if (!isValidNomenclature(nomenclature)) {
        console.log(`[ChemRenderer] ⏭️ Skipping image - invalid/error nomenclature: "${nomenclature?.substring?.(0, 50) || 'N/A'}..."`);
        continue;
      }

      const compoundType = img.dataset.compoundType || 'compound';

      // Now we just need nomenclature - server handles all lookups!
      if (!nomenclature || nomenclature === 'molecule') {
        console.log(`[ChemRenderer] ⏭️ Skipping image - no nomenclature`);
        continue;
      }

      // PRESERVE original flags (isDirectSmiles, isIUPAC, smilesValue, iupacValue, etc.)
      // to ensure proper server endpoint routing on re-render
      // NOTE: flagsLocked means "don't let AI override inline flags", NOT "never re-render"
      //       User settings from popup should ALWAYS trigger re-renders!
      let originalFlags = {};
      let originalSmiles = null;
      if (img.dataset.moleculeViewer) {
        try {
          const data = JSON.parse(atob(img.dataset.moleculeViewer));
          // CRITICAL: Preserve the original flags for proper server endpoint routing
          // Without this, chem:smiles= and chem:iupac= molecules would incorrectly use mol= (PubChem)
          originalFlags = data.flags || {};
          originalSmiles = data.smiles || null;
          console.log(`[ChemRenderer] 📦 Preserved original flags:`, {
            isDirectSmiles: originalFlags.isDirectSmiles,
            isIUPAC: originalFlags.isIUPAC,
            smilesValue: originalFlags.smilesValue,
            iupacValue: originalFlags.iupacValue,
            compoundType: originalFlags.compoundType
          });
        } catch (e) {
          // If we can't decode, just continue with re-render using defaults
        }
      }

      console.log(`[ChemRenderer] 🔄 Re-rendering: ${nomenclature} (type: ${compoundType})`);


      // Build molecule data with CURRENT settings based on compound type
      let moleculeData;

      if (compoundType === 'compound') {
        // Compound - server handles lookup by name OR renders SMILES directly if available
        // CRITICAL: Use originalSmiles if available (from chem:smiles= tags), otherwise fall back to dataset.smiles
        const effectiveSmiles = originalSmiles || smiles || null;

        moleculeData = {
          nomenclature: nomenclature,
          smiles: effectiveSmiles,  // Preserve original SMILES for direct rendering
          type: effectiveSmiles ? 'smiles' : 'nomenclature',  // Use nomenclature if no cached SMILES
          options: {
            width: parseInt(img.dataset.renderWidth) || 300,
            height: parseInt(img.dataset.renderHeight) || 240,
            aromaticCircles: settings.m2cfAromaticCircles,
            showCarbons: settings.m2cfShowCarbons,
            showMethyls: settings.m2cfShowMethyls,
            atomNumbers: settings.m2cfAtomNumbers,
            addH2: settings.m2cfAddH2,
            showImplicitH: settings.m2cfShowImplicitH,
            flipHorizontal: settings.m2cfFlipHorizontal,
            flipVertical: settings.m2cfFlipVertical,
            scale: 1,
            rotation: parseInt(img.dataset.rotation) || 0
          },
          isMoleculeViewer: true,
          flags: {
            // Merge original flags with new flags to preserve type routing
            ...originalFlags,
            useDefaults: false,
            compoundType: 'compound',
            // Explicitly preserve critical routing flags
            isDirectSmiles: originalFlags.isDirectSmiles || false,
            isIUPAC: originalFlags.isIUPAC || false,
            smilesValue: originalFlags.smilesValue || effectiveSmiles,
            iupacValue: originalFlags.iupacValue || null
          }
        };
      } else if (compoundType === 'biomolecule') {
        // Protein/biomolecule settings
        moleculeData = {
          nomenclature: nomenclature,
          type: 'nomenclature',
          options: {
            width: parseInt(img.dataset.renderWidth) || 300,
            height: parseInt(img.dataset.renderHeight) || 240,
            // Protein-specific settings
            removeWhiteBg: settings.proteinRemoveWhiteBg,
            bioAssembly: settings.molviewBioAssembly,
            chainType: settings.molviewChainType,
            chainBonds: settings.molviewChainBonds,
            chainColor: settings.molviewChainColor,
            proteinBgColor: settings.proteinMolviewBgColor
          },
          isMoleculeViewer: true,
          flags: {
            // Merge original flags to preserve routing info (biomolValue, pdbid, etc.)
            ...originalFlags,
            useDefaults: false,
            compoundType: 'biomolecule'
          }
        };
      } else if (compoundType === 'mineral') {
        // Mineral with COD settings
        moleculeData = {
          nomenclature: nomenclature,
          type: 'nomenclature',
          options: {
            width: parseInt(img.dataset.renderWidth) || 300,
            height: parseInt(img.dataset.renderHeight) || 240,
            // Mineral-specific settings
            representation: settings.mineralRepresentation,
            mineralBgColor: settings.mineralMolviewBgColor,
            crystallography: settings.mineralCrystallography
          },
          isMoleculeViewer: true,
          flags: {
            // Merge original flags to preserve routing info (mineralValue, codid, etc.)
            ...originalFlags,
            useDefaults: false,
            compoundType: 'mineral'
          }
        };
      } else {
        // Default/unknown type
        const effectiveSmiles = originalSmiles || smiles || null;
        moleculeData = {
          nomenclature: nomenclature,
          smiles: effectiveSmiles,
          type: effectiveSmiles ? 'smiles' : 'nomenclature',
          options: {
            width: parseInt(img.dataset.renderWidth) || 300,
            height: parseInt(img.dataset.renderHeight) || 240
          },
          isMoleculeViewer: true,
          flags: {
            // Merge original flags to preserve any routing info
            ...originalFlags,
            useDefaults: true,
            compoundType: compoundType
          }
        };
      }

      // Encode the updated molecule data
      const encodedData = btoa(JSON.stringify(moleculeData));

      // Create a fresh placeholder image that will be processed by loadMoleculeImage
      const newPlaceholder = document.createElement('img');
      newPlaceholder.className = 'molecule-diagram molecule-viewer';
      newPlaceholder.alt = nomenclature;
      newPlaceholder.dataset.moleculeViewer = encodedData;
      newPlaceholder.dataset.loaded = 'false';
      newPlaceholder.dataset.rotation = img.dataset.rotation || '0';
      // CRITICAL: Also preserve dataset attributes that renderClientSide uses
      newPlaceholder.dataset.nomenclature = nomenclature;
      newPlaceholder.dataset.compoundType = compoundType;
      if (moleculeData.smiles) {
        newPlaceholder.dataset.smiles = moleculeData.smiles;
      }
      newPlaceholder.style.cssText = 'display: inline-block; height: auto; margin: 0 12px 8px 0; vertical-align: middle; padding: 4px; border-radius: 4px; background-color: transparent;';

      // Find the parent - could be a wrapper div or direct parent
      const wrapper = img.closest('.molecule-container') || img.parentNode;

      if (wrapper && wrapper.parentNode) {
        // Replace the wrapper/image with the new placeholder
        wrapper.parentNode.replaceChild(newPlaceholder, wrapper);
      } else if (img.parentNode) {
        // Just replace the image directly
        img.parentNode.replaceChild(newPlaceholder, img);
      } else {
        console.error('[ChemRenderer] ❌ Image has no parent!');
        continue;
      }

      // Use the proven rendering pipeline
      if (window._loadMoleculeImage) {
        window._loadMoleculeImage(newPlaceholder);
        successCount++;
        console.log(`[ChemRenderer] ✅ Triggered loadMoleculeImage for: ${nomenclature}`);
      } else {
        console.error('[ChemRenderer] ❌ _loadMoleculeImage not available!');
      }

    } catch (e) {
      console.error('[ChemRenderer] Error re-rendering:', e);
    }
  }

  console.log(`[ChemRenderer] ✅ Re-rendered ${successCount} molecules`);

  // Re-observe any unloaded images that weren't affected by our update
  // These are images that haven't been scrolled into view yet
  const unloadedImages = document.querySelectorAll('img.molecule-diagram[data-loaded="false"], img.molecule-viewer[data-loaded="false"]');
  if (unloadedImages.length > 0 && window._lazyLoadObserver) {
    console.log(`[ChemRenderer] 👁️ Re-observing ${unloadedImages.length} unloaded images for lazy loading`);
    unloadedImages.forEach(img => {
      window._lazyLoadObserver.observe(img);
    });
  }

  // CRITICAL: Reapply average size scaling after re-render
  // Without this, the size resets to default when other settings change
  if (settings.sdAverageSize && settings.sdAverageSize !== 100) {
    // Delay slightly to allow new images to render first
    setTimeout(() => {
      console.log('[ChemRenderer] 🔄 Reapplying average size scaling after re-render:', settings.sdAverageSize + '%');
      applyAverageSizeScaling(settings.sdAverageSize);
    }, 500);
  }

  // CRITICAL: Allow mutation observer to work again after updates are done
  isUpdatingImages = false;
}


/**
 * Set up intersection observer for lazy re-rendering.
 * When an image enters the viewport:
 * 1. Show loading indicator (red spinner)
 * 2. Re-render with new settings
 * 3. Swap image only when new one is ready
 * 4. Remove loading indicator
 */
let lazyReRenderObserver = null;
const imagesBeingRendered = new Set(); // Track images currently being processed
const MAX_IMAGES_BEING_RENDERED = 20;  // Safety limit to prevent memory bloat

function setupLazyReRenderObserver() {
  // Disconnect any existing observer first to prevent memory leaks
  if (lazyReRenderObserver) {
    try {
      lazyReRenderObserver.disconnect();
    } catch (e) {
      // Ignore disconnect errors
    }
  }

  lazyReRenderObserver = new IntersectionObserver((entries) => {
    entries.forEach(entry => {
      if (entry.isIntersecting) {
        const img = entry.target;
        processImageIfNeeded(img);
      }
    });
  }, {
    rootMargin: '50px',  // Start loading slightly before fully visible
    threshold: 0.1       // Trigger when at least 10% is visible
  });

  // Observe all molecule images
  document.querySelectorAll('img[data-molecule-viewer]').forEach(img => {
    lazyReRenderObserver.observe(img);
  });

  // After setting up observer, immediately check for visible images that need re-render
  // (IntersectionObserver won't fire for images already in viewport)
  triggerVisibleImagesReRender();
}

/**
 * Process an image if it needs re-rendering and isn't already being processed
 */
function processImageIfNeeded(img) {
  const smiles = img.dataset.smiles || '(no smiles)';
  const compoundName = img.dataset.compoundName || '(no name)';

  console.log('[ChemRenderer] 🔧 processImageIfNeeded called for:', compoundName);

  // Check if this image needs re-rendering
  if (img.dataset.needsReRender !== 'true' && img.getAttribute('data-needs-re-render') !== 'true') {
    console.log('[ChemRenderer] ⏭️ Image does not need re-render (no flag)');
    return; // Doesn't need re-render
  }

  // Prevent double-processing
  if (imagesBeingRendered.has(img)) {
    console.log('[ChemRenderer] ⏭️ Image already being processed');
    log.debug(`⏭️ Image already being processed: ${compoundName}`);
    return; // Already being processed
  }

  // SAFETY: Limit concurrent rendering to prevent memory issues
  if (imagesBeingRendered.size >= MAX_IMAGES_BEING_RENDERED) {
    console.log('[ChemRenderer] ⚠️ Too many images being rendered, skipping:', compoundName);
    return;
  }

  console.log('[ChemRenderer] ✅ Image will be processed:', compoundName);
  log.info(`🔄 Processing image for lazy re-render: ${compoundName} (SMILES: ${smiles?.substring?.(0, 30) || 'N/A'}...)`);
  imagesBeingRendered.add(img);

  // NOW show the loading indicator (user is looking at this image)
  console.log('[ChemRenderer] 🔴 Adding loading indicator...');
  addLoadingIndicator(img);

  // Re-render with new settings (old image stays visible until new one is ready)
  console.log('[ChemRenderer] 🔄 Starting reRenderSingleMolecule...');
  reRenderSingleMolecule(img, pendingReRenderSettings || settings)
    .then(() => {
      // Success - clear the pending flag
      console.log('[ChemRenderer] ✅ Re-render SUCCESS for:', compoundName);
      delete img.dataset.needsReRender;
      img.removeAttribute('data-needs-re-render');
      log.debug('✅ Successfully re-rendered molecule');
    })
    .catch((err) => {
      console.log('[ChemRenderer] ❌ Re-render FAILED for:', compoundName, err);
      log.error('Failed to re-render molecule:', err);
      // Clear the flag anyway to prevent infinite retries
      delete img.dataset.needsReRender;
      img.removeAttribute('data-needs-re-render');
    })
    .finally(() => {
      // Always remove loading indicator and clear from set
      console.log('[ChemRenderer] 🔵 Removing loading indicator...');
      removeLoadingIndicator(img);
      imagesBeingRendered.delete(img);
    });
}

/**
 * Check for images that are already visible and trigger their re-render
 * This handles the case where images are in viewport when settings change
 */
function triggerVisibleImagesReRender() {
  console.log('[ChemRenderer] 👀 triggerVisibleImagesReRender called');
  const moleculeImages = document.querySelectorAll('img[data-molecule-viewer][data-needs-re-render="true"]');
  console.log('[ChemRenderer] 🔍 Found images needing re-render:', moleculeImages.length);

  let visibleCount = 0;
  moleculeImages.forEach(img => {
    const rect = img.getBoundingClientRect();
    const isVisible = rect.top < window.innerHeight && rect.bottom > 0 && rect.width > 0 && rect.height > 0;

    console.log('[ChemRenderer] 🖼️ Image check:', {
      name: img.dataset.compoundName || img.dataset.smiles?.substring(0, 20),
      isVisible,
      rect: { top: rect.top, bottom: rect.bottom, width: rect.width, height: rect.height }
    });

    if (isVisible) {
      visibleCount++;
      // This image is visible and needs re-render - process it now
      console.log('[ChemRenderer] ▶️ Processing visible image...');
      processImageIfNeeded(img);
    }
  });

  console.log('[ChemRenderer] 📊 Processed', visibleCount, 'visible images');
}

/**
 * Re-render a single molecule image using server.
 * Note: Loading indicator is handled by the caller (observer).
 * Old image stays visible until new one is ready.
 */
async function reRenderSingleMolecule(img, newSettings) {
  try {
    const smiles = img.dataset.smiles;
    const moleculeName = img.dataset.moleculeName || img.dataset.nomenclature;
    const compoundType = img.dataset.compoundType || 'compound';

    // Only re-render compounds - biomolecules and minerals use different renderers
    if (compoundType === 'biomolecule' || compoundType === 'mineral') {
      log.debug(`Skipping re-render for ${compoundType} (uses different renderer)`);
      return;
    }

    if (!smiles && !moleculeName) {
      log.debug('Skipping re-render: no SMILES or name data on image');
      return;
    }

    // CRITICAL: Decode moleculeViewer data to get original type flags
    // Without this, chem:smiles= and chem:iupac= would incorrectly use mol= (PubChem)
    let originalFlags = {};
    let originalSmiles = null;
    if (img.dataset.moleculeViewer) {
      try {
        const data = JSON.parse(atob(img.dataset.moleculeViewer));
        originalFlags = data.flags || {};
        originalSmiles = data.smiles || null;
        console.log(`[ChemRenderer] 📦 reRenderSingleMolecule - preserved flags:`, {
          isDirectSmiles: originalFlags.isDirectSmiles,
          isIUPAC: originalFlags.isIUPAC,
          smilesValue: originalFlags.smilesValue
        });
      } catch (e) {
        // If we can't decode, continue with fallback logic
      }
    }

    // Build server URL with current settings
    const renderOptions = {
      showCarbons: settings.m2cfShowCarbons,
      showMethyls: settings.m2cfShowMethyls,
      aromaticCircles: settings.m2cfAromaticCircles,
      showHydrogens: settings.m2cfAddH2,
      showImplicitH: settings.m2cfShowImplicitH,
      atomNumbers: settings.m2cfAtomNumbers,
      gradientColors: settings.sdGradientColors,
      compactDrawing: settings.m2cfCompactDrawing,
      useStereochemistry: settings.m2cfUse3DSmiles,
      flipHorizontal: settings.m2cfFlipHorizontal,
      flipVertical: settings.m2cfFlipVertical,
      renderingEngine: settings.renderingEngine
    };

    // Determine the correct server endpoint based on original type flags
    let type, value;
    const effectiveSmiles = originalSmiles || smiles;

    if (originalFlags.isDirectSmiles && (originalFlags.smilesValue || effectiveSmiles)) {
      // Direct SMILES - use smiles endpoint
      type = 'smiles';
      value = originalFlags.smilesValue || effectiveSmiles;
    } else if (originalFlags.isIUPAC && originalFlags.iupacValue) {
      // IUPAC name - use iupac endpoint (OPSIN)
      type = 'iupac';
      value = originalFlags.iupacValue;
    } else if (effectiveSmiles) {
      // Has SMILES data - use smiles endpoint
      type = 'smiles';
      value = effectiveSmiles;
    } else {
      // Fallback to mol= (PubChem lookup by name)
      type = 'mol';
      value = moleculeName;
    }

    const serverUrl = buildServerSvgUrl(type, value, renderOptions);

    log.debug(`🖼️ Re-rendering via server: ${serverUrl}`);

    // Pre-load the new image before swapping
    await new Promise((resolve, reject) => {
      const probe = new Image();
      probe.onload = () => {
        log.debug('✅ New image pre-loaded successfully');
        resolve();
      };
      probe.onerror = (err) => {
        log.error('❌ Failed to pre-load new image:', err);
        reject(err);
      };
      probe.src = serverUrl;
    });

    // Swap to new image with smooth fade-in
    log.debug('🔄 Swapping image src...');
    img.style.animation = 'chemtex-fade-in 160ms ease-out';
    img.src = serverUrl;
    log.debug('✅ Image src updated successfully');

    setTimeout(() => {
      try { img.style.animation = ''; } catch (_) { }
    }, 220);

  } catch (e) {
    log.error('Error re-rendering single molecule:', e);
    throw e;
  }
}


/**
 * Reload ALL images on the page (full re-render).
 * This is triggered by the "Reload All" button - immediate re-render for all.
 * Shows loading indicator on each image during re-render.
 */
async function reloadAllImages() {
  log.info('🔄 Reloading all images (forced refresh)...');

  // Find all molecule containers and re-render them
  const moleculeImages = document.querySelectorAll('img[data-molecule-viewer]');
  let count = 0;

  for (const img of moleculeImages) {
    // Show loading indicator during re-render
    addLoadingIndicator(img);
    try {
      await reRenderSingleMolecule(img, settings);
      count++;
    } catch (e) {
      log.error('Failed to reload image:', e);
    } finally {
      removeLoadingIndicator(img);
    }
  }

  log.success(`✅ Reloaded ${count} images`);
}

/**
 * Apply CSS-only settings changes directly to existing images without re-fetching.
 * This is used for settings that only affect visual appearance, not the actual SVG content.
 * @param {object} settings - The current settings object
 */
function applyCssOnlySettings(settings) {
  console.log('[ChemRenderer] 🎨 Applying CSS-only settings to existing images');

  // Handle 2D molecule theme and size
  const moleculeImages = document.querySelectorAll('img[data-molecule-viewer="true"]');
  moleculeImages.forEach(img => {
    // Apply theme colors if raw SVG is available
    const rawSvg = img.dataset.rawSvg;
    if (rawSvg) {
      const decodedSvg = decodeURIComponent(rawSvg);
      const coloredSvg = applyThemeColors(decodedSvg, settings.sdTheme || 'light');
      img.src = 'data:image/svg+xml;charset=utf-8,' + encodeURIComponent(coloredSvg);
      console.log('[ChemRenderer] ✅ Applied theme instantly to molecule:', img.dataset.nomenclature);
    }

    // Handle Gradient Bonds
    if (settings.sdGradientColors) {
      img.classList.add('gradient-bonds');
    } else {
      img.classList.remove('gradient-bonds');
    }
  });

  // Handle proteinRemoveWhiteBg - applies to biomolecule iframes
  const biomoleculeViewers = document.querySelectorAll('iframe.biomolecule-viewer, .molecule-3d-viewer');
  biomoleculeViewers.forEach(viewer => {
    if (settings.proteinRemoveWhiteBg) {
      viewer.style.filter = 'drop-shadow(0 0 0 transparent)';
      viewer.style.mixBlendMode = 'multiply';
      console.log('[ChemRenderer] ✅ Applied white-removal effect to biomolecule viewer');
    } else {
      viewer.style.filter = '';
      viewer.style.mixBlendMode = '';
    }
  });

  console.log(`[ChemRenderer] 🎨 Applied CSS settings to ${moleculeImages.length} molecules and ${biomoleculeViewers.length} 3D viewers`);
}

// Listen for messages from popup via background script
chrome.runtime.onMessage.addListener((request, sender, sendResponse) => {
  // SAFETY: Wrap in try-catch to prevent crashes from message handling
  try {
    console.warn('[ChemRenderer] 📩📩📩 MESSAGE RECEIVED:', request.type);
    console.log('%c[ChemRenderer] Message details:', 'background: yellow; color: black; font-weight: bold;', request);

    if (request.type === 'APPLY_SETTINGS') {
      console.warn('[ChemRenderer] ⚙️⚙️⚙️ APPLY_SETTINGS - STARTING UPDATE');
      console.log('%c[ChemRenderer] Settings to apply:', 'background: green; color: white; font-weight: bold;', request.settings);
      log.info('📬 Received settings update from popup:', request.settings);

      // Check if this is a CSS-only setting change that doesn't require re-rendering
      // NOTE: sdTheme is NOT css-only anymore - it requires re-render for marker-based colors
      const cssOnlySettings = ['proteinRemoveWhiteBg', 'sdGradientColors', 'sdAutoAdapt'];

      const changedSettings = request.changedSettings || [];
      const isCssOnlyChange = changedSettings.length > 0 &&
        changedSettings.every(key => cssOnlySettings.includes(key));

      if (isCssOnlyChange) {
        console.log('[ChemRenderer] 🎨 CSS-only change detected, applying without re-render');
        applyCssOnlySettings(request.settings);
        sendResponse({ success: true, cssOnly: true });
        return true;
      }

      // Full re-render needed for rendering settings changes
      try {
        console.warn('[ChemRenderer] 🔄🔄🔄 CALLING lazyReRenderMolecules NOW...');
        lazyReRenderMolecules(request.settings);
        console.warn('[ChemRenderer] ✅✅✅ lazyReRenderMolecules FINISHED');
      } catch (e) {
        console.error('[ChemRenderer] Error in lazyReRenderMolecules:', e);
      }

      sendResponse({ success: true });
      return true;
    }

    if (request.type === 'RELOAD_ALL_IMAGES') {
      console.warn('[ChemRenderer] 🔄 RELOAD_ALL_IMAGES received');
      log.info('📬 Received reload all images command');
      try {
        reloadAllImages();
      } catch (e) {
        console.error('[ChemRenderer] Error in reloadAllImages:', e);
      }
      sendResponse({ success: true });
      return true;
    }

    if (request.type === 'CLEAR_CACHE') {
      console.warn('[ChemRenderer] 🗑️ CLEAR_CACHE received');

      // Call the centralized cache clear function
      // This sets cacheBustTimestamp so all new SVG requests bypass browser cache
      clearAllCaches();

      // Now trigger a full re-render with fresh content
      // Get current settings and re-render all molecules
      chrome.storage.sync.get(null, (currentSettings) => {
        console.log('[ChemRenderer] 🔄 Re-rendering all molecules with fresh cache...');
        lazyReRenderMolecules(currentSettings);
      });

      sendResponse({ success: true, message: 'Cache cleared, re-rendering with fresh content' });
      return true;
    }
  } catch (e) {
    console.error('[ChemRenderer] Error handling message:', e);
    sendResponse({ success: false, error: e.message });
    return true;
  }
});

// ============================================
// DARK MODE DETECTION
// ============================================
// Check if dark mode is enabled for appropriate theming of SVGs

/**
 * Detect if the page is in dark mode
 * Returns true only if background is genuinely dark, false for light backgrounds
 */
function isDarkModeEnabled() {
  // Method 1: Check prefers-color-scheme media query
  if (window.matchMedia && window.matchMedia('(prefers-color-scheme: dark)').matches) {
    return true;
  }

  // Method 2: Check color-scheme CSS property on html/body
  try {
    const htmlStyle = window.getComputedStyle(document.documentElement);
    const bodyStyle = window.getComputedStyle(document.body);
    const colorScheme = htmlStyle.colorScheme || bodyStyle.colorScheme || '';
    if (colorScheme.includes('dark')) {
      return true;
    }
  } catch (e) { }

  // Method 3: Check html/body background color (detect dark backgrounds)
  try {
    // Check both html and body - some sites (like Wikipedia) style html element
    const elements = [document.documentElement, document.body];
    for (const el of elements) {
      if (!el) continue;
      const bg = window.getComputedStyle(el).backgroundColor;
      if (bg && bg !== 'rgba(0, 0, 0, 0)' && bg !== 'transparent') {
        const rgb = bg.match(/\d+/g);
        if (rgb && rgb.length >= 3) {
          const luminance = (0.299 * parseInt(rgb[0]) + 0.587 * parseInt(rgb[1]) + 0.114 * parseInt(rgb[2]));
          // Dark if luminance is below 100
          if (luminance < 100) {
            return true;
          }
          // Explicitly light if luminance > 180
          if (luminance > 180) {
            return false;
          }
        }
      }
    }
  } catch (e) { }

  // Method 4: Check for common dark mode classes on html or body
  const darkClasses = ['dark', 'dark-mode', 'dark-theme', 'night-mode', 'skin-minerva-dark', 'client-dark-mode'];
  for (const className of darkClasses) {
    if (document.body?.classList?.contains(className) ||
      document.documentElement?.classList?.contains(className)) {
      return true;
    }
  }

  // Method 5: Check data-theme and data-color-mode attributes (Wikipedia, GitHub, etc.)
  try {
    const dataTheme = document.documentElement.getAttribute('data-theme') ||
      document.body.getAttribute('data-theme') ||
      document.documentElement.getAttribute('data-color-mode') ||
      document.body.getAttribute('data-color-mode') || '';
    if (dataTheme.toLowerCase().includes('dark') || dataTheme.toLowerCase().includes('night')) {
      return true;
    }
  } catch (e) { }

  // Default to light mode
  return false;
}


// ============================================
// AUTO-ADAPT CODE REMOVED
// ============================================
// Users manually select their theme - no auto-adaptation to page dark/light mode



// ============================================
// BACKGROUND FETCH HELPERS (CSP bypass)
// ============================================
// These functions use the background service worker to make fetch requests,
// bypassing Content Security Policy restrictions on sites like ChatGPT

/**
 * Check if chrome.runtime is available (extension context valid)
 */
function isExtensionContextValid() {
  try {
    if (extensionContextInvalid) return false;
    if (!chrome.runtime || !chrome.runtime.id) return false;
    return true;
  } catch (e) {
    return false;
  }
}

/**
 * Fetch blob/image via background script
 * Returns { base64, type } or falls back to direct blob
 */
async function backgroundFetchBlob(url) {
  try {
    if (!isExtensionContextValid()) {
      return directFetchBlob(url);
    }
    if (chrome.runtime && chrome.runtime.sendMessage) {
      return new Promise((resolve, reject) => {
        chrome.runtime.sendMessage(
          { type: 'FETCH_BLOB', url: url },
          (response) => {
            if (chrome.runtime.lastError) {
              // Check for extension context invalidated
              if (chrome.runtime.lastError.message && chrome.runtime.lastError.message.includes('Extension context invalidated')) {
                extensionContextInvalid = true;
              }
              // console.warn('[ChemRenderer] Background blob fetch failed, trying direct:', chrome.runtime.lastError.message);
              directFetchBlob(url).then(resolve).catch(reject);
              return;
            }
            if (response && response.success) {
              resolve(response.data);
            } else {
              reject(new Error(response?.error || 'Background blob fetch failed'));
            }
          }
        );
      });
    }
  } catch (e) {
    // Check for extension context invalidated
    if (e.message && e.message.includes('Extension context invalidated')) {
      extensionContextInvalid = true;
    }
    // console.warn('[ChemRenderer] Background unavailable for blob, using direct fetch');
  }
  return directFetchBlob(url);
}

/**
 * Direct fetch for blob (used as fallback)
 */
async function directFetchBlob(url) {
  const response = await fetch(url);
  if (!response.ok) {
    throw new Error(`HTTP ${response.status}: ${response.statusText} `);
  }
  const blob = await response.blob();
  // Convert to base64 for consistency with background response
  return new Promise((resolve, reject) => {
    const reader = new FileReader();
    reader.onloadend = () => {
      const base64 = reader.result.split(',')[1];
      resolve({ base64: base64, type: blob.type });
    };
    reader.onerror = reject;
    reader.readAsDataURL(blob);
  });
}

// ============================================
// 🌉 SMILES BRIDGE - Centralized Name→SMILES Conversion
// ============================================
// All renderers should use this bridge for nomenclature→SMILES conversion
// This centralizes the conversion logic so all renderers use the same flow

/**
 * SMILES Bridge: Convert chemical nomenclature to SMILES string
 * Uses ChemTex Server search API
 *
 * @param {string} name - Chemical name/nomenclature to convert
 * @param {Object} options - Optional settings
 * @param {boolean} options.use3DSmiles - Use IsomericSMILES for stereochemistry (default: false)
 * @returns {Promise<{smiles: string, source: string, cid?: number}|null>}
 */
async function smilesBridge(name, options = {}) {
  const { use3DSmiles = false } = options;

  if (!name || typeof name !== 'string') {
    console.error('❌ SMILES Bridge: Invalid name provided');
    return null;
  }

  // Strip flags from molecule name before conversion
  const cleanName = stripFlagsFromName(name.trim());
  if (!cleanName) {
    console.error('❌ SMILES Bridge: Empty name after stripping flags');
    return null;
  }

  // Use ChemistryLaTeX Server search API
  try {
    const searchUrl = `${CHEMISTRYLATEX_SERVER_URL}/search/${encodeURIComponent(cleanName)}?type=compound`;
    console.log('%c🌉 SMILES Bridge: Using server search', 'color: #4A90D9; font-weight: bold;', searchUrl);

    const response = await fetch(searchUrl);
    if (response.ok) {
      const data = await response.json();
      if (data && data.smiles) {
        // Use isomeric SMILES if available and requested
        const smiles = use3DSmiles && data.isomericSmiles ? data.isomericSmiles : data.smiles;
        return {
          smiles: smiles,
          source: 'ChemTex-Server',
          cid: data.cid
        };
      }
    }
  } catch (e) {
    console.warn('⚠️ [Bridge] Server search failed:', e.message);
  }

  console.error('%c❌ [Bridge] Name→SMILES conversion failed', 'color: #FF0000; font-weight: bold;');
  return null;
}

// Export for global access (useful for debugging in console)
window.smilesBridge = smilesBridge;

/**
 * Get compound ID for a molecule name
 * Uses ChemTex Server search API
 * 
 * @param {string} nameOrSmiles - Chemical name or SMILES string
 * @returns {Promise<number|null>} - Compound ID or null if not found
 */
async function getCompoundID(nameOrSmiles) {
  if (!nameOrSmiles || typeof nameOrSmiles !== 'string') {
    console.error('❌ Invalid input for getCompoundID');
    return null;
  }

  const cleanInput = nameOrSmiles.trim();

  // Use smilesBridge which now returns CID from server
  const bridgeResult = await smilesBridge(cleanInput, { use3DSmiles: false });

  if (bridgeResult && bridgeResult.cid) {
    console.log('%c✅ [CID] Got ID from server:', 'color: #00FF00; font-weight: bold;', bridgeResult.cid);
    return bridgeResult.cid;
  }

  console.error('%c❌ [CID] Could not find CID for:', 'color: #FF0000; font-weight: bold;', nameOrSmiles);
  return null;
}

// Export for global access
window.getCompoundID = getCompoundID;
// Legacy alias
window.getCompoundID = getCompoundID;


// addHoverControls is defined later

// ============================================
// IMAGE SIZE CONTROLS
// These may also be defined in size-controls.js, so use var to avoid redeclaration errors
// ============================================
var SIZE_STEP = typeof SIZE_STEP !== 'undefined' ? SIZE_STEP : 20;
var MIN_SIZE = typeof MIN_SIZE !== 'undefined' ? MIN_SIZE : 100;
var MAX_SIZE = typeof MAX_SIZE !== 'undefined' ? MAX_SIZE : 800;
var DEFAULT_WIDTH = typeof DEFAULT_WIDTH !== 'undefined' ? DEFAULT_WIDTH : 400;
var DEFAULT_HEIGHT = typeof DEFAULT_HEIGHT !== 'undefined' ? DEFAULT_HEIGHT : 350;

// Non-linear (parabolic) size calculation constants
var BASE_SIZE = typeof BASE_SIZE !== 'undefined' ? BASE_SIZE : 150; // Base size for small molecules
var SIZE_SCALE_FACTOR = typeof SIZE_SCALE_FACTOR !== 'undefined' ? SIZE_SCALE_FACTOR : 0.5; // Controls the rate of size increase
var MAX_DEFAULT_SIZE = typeof MAX_DEFAULT_SIZE !== 'undefined' ? MAX_DEFAULT_SIZE : 400; // Maximum default size for very large molecules


function getImageKey(moleculeData) {
  if (!moleculeData) return null;
  if (moleculeData.smiles) return `smiles:${moleculeData.smiles}`;
  if (moleculeData.nomenclature) return `nomenclature:${moleculeData.nomenclature}`;
  return null;
}

function getPageImageKey(moleculeData, pageUrl) {
  const imageKey = getImageKey(moleculeData);
  if (!imageKey) return null;
  return `page:${pageUrl}:${imageKey}`;
}

async function loadImageSize(moleculeData, pageUrl, settings) {
  try {
    const imageKey = getImageKey(moleculeData);
    if (!imageKey) return {}; // Return empty object, let smart scaling decide
    if (!settings.saveSizePerImage && !settings.saveSizeBySMILES) {
      return {}; // Return empty object, let smart scaling decide
    }
    let storageKey;
    if (settings.saveSizePerImage && !settings.saveSizeBySMILES) {
      storageKey = getPageImageKey(moleculeData, pageUrl);
    } else {
      storageKey = imageKey;
    }
    if (!storageKey) return {}; // Return empty object, let smart scaling decide
    return new Promise((resolve) => {
      chrome.storage.local.get([storageKey], (result) => {
        if (result[storageKey]) {
          // Support both old format (width/height) and new format (scale)
          if (result[storageKey].scale) {
            resolve(result[storageKey]);
          } else {
            // Legacy: convert old width/height to scale (assume default was 400x350)
            const scale = result[storageKey].width ? result[storageKey].width / DEFAULT_WIDTH : undefined;
            resolve(scale !== undefined ? { scale } : {});
          }
        } else {
          resolve({}); // Return empty object, let smart scaling decide
        }
      });
    });
  } catch (error) {
    console.error('Error loading image size:', error);
    return {}; // Return empty object, let smart scaling decide
  }
}

async function saveImageSize(moleculeData, pageUrl, size, settings) {
  try {
    const imageKey = getImageKey(moleculeData);
    if (!imageKey) return;
    if (!settings.saveSizePerImage && !settings.saveSizeBySMILES) return;
    let storageKey;
    if (settings.saveSizePerImage && !settings.saveSizeBySMILES) {
      storageKey = getPageImageKey(moleculeData, pageUrl);
    } else {
      storageKey = imageKey;
    }
    if (!storageKey) return;
    const data = {};
    data[storageKey] = size;
    chrome.storage.local.set(data, () => {
      console.log(`Saved size for ${storageKey}:`, size);
    });
  } catch (error) {
    console.error('Error saving image size:', error);
  }
}

function createSizeControls(container, svgImg, moleculeData, settings) {
  console.log('🎛️ Creating size controls:', { settings, moleculeData });

  const controlsDiv = document.createElement('div');
  controlsDiv.className = 'chem-size-controls';
  controlsDiv.style.cssText = `
    position: absolute;
    bottom: 4px;
    left: 4px;
    display: flex !important;
    flex-direction: column;
    gap: 2px;
    opacity: 0;
    pointer-events: auto;
    transition: opacity 0.2s;
    z-index: 1000;
  `;
  const upButton = document.createElement('button');
  upButton.className = 'chem-size-btn chem-size-up';
  upButton.innerHTML = '▲';
  upButton.title = 'Increase size';
  upButton.style.cssText = `
    width: 24px;
    height: 24px;
    border: none;
    background: rgba(0, 0, 0, 0.7);
    color: white;
    border-radius: 4px;
    cursor: pointer;
    font-size: 10px;
    display: flex;
    align-items: center;
    justify-content: center;
    transition: background 0.2s;
    pointer-events: auto;
  `;
  const downButton = document.createElement('button');
  downButton.className = 'chem-size-btn chem-size-down';
  downButton.innerHTML = '▼';
  downButton.title = 'Decrease size';
  downButton.style.cssText = upButton.style.cssText;
  [upButton, downButton].forEach(btn => {
    btn.addEventListener('mouseenter', () => {
      btn.style.background = 'rgba(0, 0, 0, 0.9)';
    });
    btn.addEventListener('mouseleave', () => {
      btn.style.background = 'rgba(0, 0, 0, 0.7)';
    });
  });
  upButton.addEventListener('click', (e) => {
    e.stopPropagation();
    adjustImageSize(container, svgImg, moleculeData, SIZE_STEP, settings);
  });
  downButton.addEventListener('click', (e) => {
    e.stopPropagation();
    adjustImageSize(container, svgImg, moleculeData, -SIZE_STEP, settings);
  });
  controlsDiv.appendChild(upButton);
  controlsDiv.appendChild(downButton);
  container.addEventListener('mouseenter', () => {
    console.log('🖱️ Mouse entered container');
    controlsDiv.style.opacity = '1';
  });
  container.addEventListener('mouseleave', () => {
    controlsDiv.style.opacity = '0';
  });
  return controlsDiv;
}

function adjustImageSize(container, svgImg, moleculeData, delta, settings) {
  // Get current scale (fallback to 1.0 if not set)
  let currentScale = parseFloat(svgImg.dataset.scale) || 1.0;

  // Adjust scale by delta (convert pixel delta to scale delta)
  const scaleDelta = delta / 100; // 20px = 0.2 scale change
  let newScale = currentScale + scaleDelta;

  // Constrain scale between 0.5x and 5x
  newScale = Math.max(0.5, Math.min(5, newScale));

  // Store scale in dataset for future adjustments
  svgImg.dataset.scale = newScale.toString();

  // Apply scale using CSS transform (preserves intrinsic dimensions)
  // Don't set width/height - let SVG's intrinsic dimensions control the size
  const existingTransform = svgImg.style.transform || '';
  const withoutScale = existingTransform.replace(/scale\([^)]+\)\s*/g, '');
  svgImg.style.transform = `${withoutScale} scale(${newScale})`.trim();
  svgImg.style.transformOrigin = 'top left';

  // Save scale preference
  const pageUrl = window.location.href;
  const size = { scale: newScale };
  saveImageSize(moleculeData, pageUrl, size, settings);

  console.log(`Adjusted scale: ${newScale.toFixed(2)}x (was ${currentScale.toFixed(2)}x)`);
}

async function wrapImageWithSizeControls(svgImg, originalImg, moleculeData, settings) {
  try {
    console.log('📦 Wrapping image with size controls', { settings, moleculeData });

    const container = document.createElement('div');
    container.className = 'chem-image-container';
    container.style.cssText = `
      position: relative;
      display: inline-block;
      margin: 0 12px 8px 0;
      vertical-align: middle;
      overflow: hidden;
      max-width: 100%;
    `;
    const pageUrl = window.location.href;
    const savedSize = await loadImageSize(moleculeData, pageUrl, settings);

    // Smart scaling: Calculate default scale based on molecule size to prevent overflow
    let defaultScale = 1.5; // Default for medium molecules

    // Try to get intrinsic dimensions to determine molecule size
    const getIntrinsicDimensions = () => {
      let width = svgImg.naturalWidth || svgImg.width || 0;
      let height = svgImg.naturalHeight || svgImg.height || 0;

      // If dimensions not available yet, try to extract from SVG data URL
      if ((!width || !height) && svgImg.src && svgImg.src.startsWith('data:image/svg+xml')) {
        try {
          let svgContent;
          if (svgImg.src.includes('base64,')) {
            svgContent = atob(svgImg.src.split('base64,')[1]);
          } else {
            svgContent = decodeURIComponent(svgImg.src.split(',')[1]);
          }

          // Extract width and height
          const widthMatch = svgContent.match(/width\s*=\s*["']?(\d+(?:\.\d+)?)/i);
          const heightMatch = svgContent.match(/height\s*=\s*["']?(\d+(?:\.\d+)?)/i);

          if (widthMatch) width = parseFloat(widthMatch[1]);
          if (heightMatch) height = parseFloat(heightMatch[1]);

          // Fallback to viewBox
          if (!width || !height) {
            const viewBoxMatch = svgContent.match(/viewBox\s*=\s*["']?[\d.]+\s+[\d.]+\s+([\d.]+)\s+([\d.]+)/i);
            if (viewBoxMatch) {
              width = width || parseFloat(viewBoxMatch[1]);
              height = height || parseFloat(viewBoxMatch[2]);
            }
          }
        } catch (e) {
          console.warn('Could not extract SVG dimensions:', e);
        }
      }

      return { width, height };
    };

    const { width: intrinsicWidth, height: intrinsicHeight } = getIntrinsicDimensions();

    // Calculate smart default scale based on molecule size
    // Use non-linear scaling to prevent large molecules from overflowing
    if (intrinsicWidth > 0 && intrinsicHeight > 0) {
      // Calculate "molecule size" as the diagonal (accounts for both width and height)
      const moleculeSize = Math.sqrt(intrinsicWidth * intrinsicWidth + intrinsicHeight * intrinsicHeight);

      // Define size thresholds and corresponding scales
      // Small molecules (< 250px diagonal): scale up for visibility (1.5x - 2.0x)
      // Medium molecules (250-500px): moderate scale (1.0x - 1.5x)
      // Large molecules (500-800px): smaller scale (0.6x - 1.0x)
      // Very large molecules (> 800px): minimal scale (0.4x - 0.6x)

      // Use logarithmic scaling for smooth, non-linear decrease
      // Formula: scale = baseScale * (referenceSize / moleculeSize)^exponent
      // The exponent controls how aggressively the scale decreases

      const referenceSize = 300; // Reference size for 1.5x scale
      const baseScale = 1.5;
      const exponent = 0.5; // Square root relationship (parabolic)

      // Calculate scale with logarithmic dampening
      defaultScale = baseScale * Math.pow(referenceSize / moleculeSize, exponent);

      // Clamp to reasonable bounds
      // Min: 0.4x (very large molecules like proteins)
      // Max: 2.0x (very small molecules like water, methane)
      defaultScale = Math.max(0.4, Math.min(2.0, defaultScale));

      console.log(`%c📏 Smart Scaling: ${intrinsicWidth}x${intrinsicHeight} (diagonal: ${moleculeSize.toFixed(0)}px) → ${defaultScale.toFixed(2)}x scale`,
        'color: #9C27B0; font-weight: bold;');
    } else {
      console.log('%c📏 Using default scale (dimensions not available yet)', 'color: #FF9800;');
    }

    // Use saved scale if available, otherwise use smart default
    const scale = savedSize.scale || defaultScale;
    svgImg.dataset.scale = scale.toString();

    // NOTE: Removed 300px default width override - we now use fit-content
    // to let images wrap perfectly to their intrinsic SVG dimensions

    // Wait for image to load to get natural dimensions
    // Skip scaling for fixed-size images (protein/mineral images)
    if (!svgImg.dataset.fixedSize) {
      if (svgImg.complete) {
        applyScaleToImage(svgImg, scale);
      } else {
        svgImg.onload = () => applyScaleToImage(svgImg, scale);
        // Fallback: if onload doesn't fire (e.g. cached), force apply after timeout
        setTimeout(() => applyScaleToImage(svgImg, scale), 100);
      }
    }

    if (moleculeData) {
      svgImg.dataset.moleculeData = JSON.stringify(moleculeData);
    }
    const controls = createSizeControls(container, svgImg, moleculeData, settings);
    originalImg.parentNode.insertBefore(container, originalImg);
    container.appendChild(svgImg);
    container.appendChild(controls);
    originalImg.remove();

    // Add hover controls (name + 3D viewer button) after container is created
    // Prioritize displayName from flags, then nomenclature, then smiles
    const moleculeName = moleculeData.flags?.displayName || moleculeData.nomenclature || moleculeData.smiles || 'Unknown';
    addHoverControls(container, moleculeName, moleculeData);

    console.log('✅ Size controls wrapper created successfully with scale:', scale);
    return container;
  } catch (error) {
    console.error('❌ Error wrapping image with size controls:', error);
    // Safe fallback - check all parent nodes before replaceChild
    if (originalImg && originalImg.parentNode && svgImg) {
      try {
        originalImg.parentNode.replaceChild(svgImg, originalImg);
      } catch (e) {
        console.error('❌ Failed to replace element:', e);
      }
    }
    return svgImg;
  }
}

function applyScaleToImage(svgImg, scale) {
  // Instead of calculating a fixed pixel width, use CSS transform: scale()
  // This preserves the SVG's intrinsic dimensions (fit-content) while applying scaling

  // Get intrinsic dimensions for logging
  let intrinsicWidth = svgImg.naturalWidth || svgImg.width;
  let intrinsicHeight = svgImg.naturalHeight || svgImg.height;

  // Try to extract from SVG if not available
  if ((!intrinsicWidth || intrinsicWidth === 0) && svgImg.src && svgImg.src.startsWith('data:image/svg+xml')) {
    try {
      let svgContent;
      if (svgImg.src.includes('base64,')) {
        svgContent = atob(svgImg.src.split('base64,')[1]);
      } else {
        svgContent = decodeURIComponent(svgImg.src.split(',')[1]);
      }

      const widthMatch = svgContent.match(/width\s*=\s*["']?(\d+(?:\.\d+)?)/i);
      const heightMatch = svgContent.match(/height\s*=\s*["']?(\d+(?:\.\d+)?)/i);

      if (widthMatch) intrinsicWidth = parseFloat(widthMatch[1]);
      if (heightMatch) intrinsicHeight = parseFloat(heightMatch[1]);

      // Try viewBox as fallback
      if (!intrinsicWidth || !intrinsicHeight) {
        const viewBoxMatch = svgContent.match(/viewBox\s*=\s*["']?[\d.]+\s+[\d.]+\s+([\d.]+)\s+([\d.]+)/i);
        if (viewBoxMatch) {
          intrinsicWidth = intrinsicWidth || parseFloat(viewBoxMatch[1]);
          intrinsicHeight = intrinsicHeight || parseFloat(viewBoxMatch[2]);
        }
      }
    } catch (e) {
      console.warn('Could not extract SVG dimensions from data URL:', e);
    }
  }

  // Hard limits to prevent molecules from getting too large
  const MAX_DISPLAY_WIDTH = 600;  // Maximum width in pixels
  const MAX_DISPLAY_HEIGHT = 500; // Maximum height in pixels

  // Don't set width/height - let SVG's intrinsic dimensions control the size
  // Only set max constraints to prevent overflow
  svgImg.style.maxWidth = `${MAX_DISPLAY_WIDTH}px`;
  svgImg.style.maxHeight = `${MAX_DISPLAY_HEIGHT}px`;

  // Apply scale using CSS transform (preserves intrinsic dimensions)
  // Only apply if scale is non-default
  if (scale && scale !== 1.0) {
    // Get existing transform (may have flip/rotation applied)
    const existingTransform = svgImg.style.transform || '';
    // Remove any existing scale() and add new one
    const withoutScale = existingTransform.replace(/scale\([^)]+\)\s*/g, '');
    svgImg.style.transform = `${withoutScale} scale(${scale})`.trim();
    svgImg.style.transformOrigin = 'top left';
  }

  console.log(`Applied scale ${scale}x`);
}

// Debug interface - minimal console utilities
window.chemRendererDebug = {
  getSettings: () => settings,
  scanPage: () => scanAndRender(),
  clearCache: () => clearAllCaches(),
  async getCacheStats() {
    await loadSmilesCache();
    return { smiles: Object.keys(smilesCache).length, images: renderedImageCache.size };
  }
};

// Performance stub
window.chemRendererPerformance = { recordStructure() { }, recordLoad() { }, recordFormula() { } };

log.info('ChemTex loaded');


// Settings with all defaults
let settings = {
  enabled: true,
  performanceMode: true,
  maxVisibleSVGs: 5,
  sizePreset: 'auto',
  rendererEngine: 'chemtex',
  devMode: false,
  // Rendering options
  flipHorizontal: false,
  flipVertical: false,
  use3DSmiles: false,  // Enable 3D stereochemistry
  // Rendering options
  sdShowCarbons: true,
  sdAromaticCircles: true,
  sdShowMethyls: true,
  sdCompactDrawing: true,
  sdTheme: 'light',
  sdAutoAdapt: true,
  // 3D Viewer settings
  enable3DViewer: false,
  viewer3DSource: 'molview',  // Uses embed.molview.org
  viewer3DSize: 'normal',
  viewer3DStyle: 'stick',
  // Compound MolView options
  molviewRepresentation: 'ballAndStick',
  compoundMolviewBgColor: 'black',
  molviewEngine: 'glmol',

  // Protein/Mineral MolView options
  proteinRemoveWhiteBg: false,
  molviewChainType: 'ribbon',
  molviewChainBonds: false,
  molviewChainColor: 'ss',
  proteinMolviewBgColor: 'black',
  molviewBioAssembly: false,
  bioAssemblyViewer: 'molstar',       // 'molstar', 'molstar-me', 'molview'
  bioAssemblyGraphics: 'balanced',    // 'quality', 'balanced', 'performance'
  bioAssemblyPdbProvider: 'rcsb',     // 'rcsb', 'pdbe', 'pdbj'
  // Mineral MolView options
  mineralRepresentation: 'ballAndStick',
  mineralMolviewBgColor: 'black',
  mineralCrystallography: 'supercell_2x2x2'
};

log.info('📦 Loading settings from storage...');

// Load settings - with proper callback
chrome.storage.sync.get(null, (result) => {
  // Merge stored settings with defaults
  settings = { ...settings, ...result };

  // Force ChemTex rendering engine
  settings.rendererEngine = 'chemtex';

  // Map popup settings (sd*) to internal rendering options (m2cf*)
  // The popup uses sd* names, but rendering code uses m2cf* names
  // Default to true if undefined (matching popup.js defaults)
  settings.m2cfShowCarbons = result.sdShowCarbons !== undefined ? result.sdShowCarbons : true;
  settings.m2cfAromaticCircles = result.sdAromaticRings !== undefined ? result.sdAromaticRings : true;
  settings.m2cfShowMethyls = result.sdShowMethyls !== undefined ? result.sdShowMethyls : true;
  settings.m2cfAtomNumbers = result.sdAtomNumbers === true;
  settings.m2cfAddH2 = result.sdShowExplicitHydrogens === true;
  settings.m2cfShowImplicitH = result.sdShowImplicitHydrogens !== false; // Default to true
  settings.m2cfCompactDrawing = result.sdCompactDrawing === true; // Default to false
  settings.m2cfFlipHorizontal = result.sdFlipHorizontal === true;
  settings.m2cfFlipVertical = result.sdFlipVertical === true;
  settings.m2cfRotate = result.sdRotate || 0;
  // Gradient colors for bonds
  settings.sdGradientColors = result.sdGradientColors === true;
  // Auto-adapt theme based on page dark/light mode - default to TRUE
  settings.sdAutoAdapt = result.sdAutoAdapt !== false;
  console.log('[Content] Auto-adapt setting loaded:', settings.sdAutoAdapt);
  // Stereochemistry - map useStereochemistry to m2cfUse3DSmiles for isomeric SMILES
  // When enabled, always prefer IsomericSMILES over CanonicalSMILES
  settings.m2cfUse3DSmiles = result.useStereochemistry !== false; // Default to TRUE
  console.log('[Content] Stereochemistry setting:', {
    useStereochemistry: result.useStereochemistry,
    m2cfUse3DSmiles: settings.m2cfUse3DSmiles
  });

  console.log('[Content] Raw storage values:', {
    sdShowCarbons: result.sdShowCarbons,
    sdAromaticRings: result.sdAromaticRings,
    sdShowMethyls: result.sdShowMethyls,
    sdAtomNumbers: result.sdAtomNumbers,
    sdShowExplicitHydrogens: result.sdShowExplicitHydrogens,
    sdShowImplicitHydrogens: result.sdShowImplicitHydrogens,
    sdFlipHorizontal: result.sdFlipHorizontal,
    sdFlipVertical: result.sdFlipVertical
  });


  // Log the loaded settings for debugging
  log.info('📊 Rendering settings:', {
    showCarbons: settings.m2cfShowCarbons,
    aromaticRings: settings.m2cfAromaticCircles,
    showMethyls: settings.m2cfShowMethyls,
    atomNumbers: settings.m2cfAtomNumbers,
    explicitHydrogens: settings.m2cfAddH2,
    implicitHydrogens: settings.m2cfShowImplicitH,
    gradientColors: settings.sdGradientColors
  });

  log.info(`🔍 Server search enabled`);
  log.success('✅ Settings loaded', settings);
  log.info(`Renderer Engine: ChemTex`);
  log.info(`Performance mode: ${settings.performanceMode ? 'ON ⚡' : 'OFF'}`);

  // Pre-load SMILES cache for faster lookups
  loadSmilesCache().then(cache => {
    log.info(`📦 SMILES cache ready: ${Object.keys(cache).length} cached compounds`);
  });

  if (settings.enabled) {
    log.info('🚀 Extension enabled, initializing renderer...');
    initializeRenderer();

    // Apply saved size scaling after images render
    settings.sdAverageSize = result.sdAverageSize || 100;
    console.log('[Content] Loaded sdAverageSize from storage:', settings.sdAverageSize);

    if (settings.sdAverageSize !== 100) {
      // Apply after initial render, then again after more images may have loaded
      setTimeout(() => {
        applyAverageSizeScaling(settings.sdAverageSize);
        console.log('[Content] Applied saved size scaling (1st pass):', settings.sdAverageSize + '%');
      }, 2000);

      // Second pass for lazy-loaded images
      setTimeout(() => {
        applyAverageSizeScaling(settings.sdAverageSize);
        console.log('[Content] Applied saved size scaling (2nd pass):', settings.sdAverageSize + '%');
      }, 5000);
    }
  } else {
    log.info('⏸️  Extension disabled in settings');
  }
});

// Listen for setting changes
chrome.storage.onChanged.addListener((changes) => {
  if (changes.enabled) {
    log.info('⚙️  Settings changed, reloading...', changes);
    if (changes.enabled.newValue) {
      location.reload();
    }
  }

  // Renderer engine is now always ChemTex
  // No need to listen for changes

  // Listen for 3D viewer source changes
  if (changes.viewer3DSource) {
    log.info('🔄 3D Viewer source changed, updating...', changes.viewer3DSource);
    settings.viewer3DSource = changes.viewer3DSource.newValue;
    log.success(`✅ Switched to ${changes.viewer3DSource.newValue} 3D viewer`);
    // Update settings immediately without reload
  }

  // NOTE: Rendering options (sd*) are now handled via APPLY_SETTINGS message
  // which triggers lazyReRenderMolecules() or reRenderAllMolecules() for instant updates
  // WITHOUT a full page reload. The code below was removed to fix double-reload issue.
  // 
  // Previously this code would call location.reload() after 500ms delay, but that
  // would conflict with the APPLY_SETTINGS message handler which already handles
  // dynamic re-rendering. Removing this prevents the page from reloading.
});

// Listen for messages from popup (live preview and settings changes)
chrome.runtime.onMessage.addListener((request, sender, sendResponse) => {
  if (request.action === 'renderPreview') {
    try {
      // Process the code using wrapChemicalFormulas
      const wrappedHtml = wrapChemicalFormulas(request.code);
      sendResponse({ html: wrappedHtml });
    } catch (error) {
      console.error('Error rendering preview:', error);
      sendResponse({ html: `<span style="color: red;">Error: ${error.message}</span>` });
    }
    return true; // Keep channel open for async response
  }

  // Handle SETTINGS_CHANGED from popup (instant updates with lazy loading)
  if (request.type === 'SETTINGS_CHANGED') {
    console.log('%c⚙️ [Content] Settings changed from popup:', 'background: #2196F3; color: white; font-weight: bold; padding: 4px;', request.settings);

    // Update internal settings object
    Object.assign(settings, request.settings);

    // Check if this is a 2D rendering setting that requires re-render
    const render2DSettings = [
      'sdShowCarbons', 'sdAromaticRings', 'sdShowMethyls', 'sdAtomNumbers',
      'sdShowExplicitHydrogens', 'sdShowImplicitHydrogens', 'sdCompactDrawing',
      'sdFlipHorizontal', 'sdFlipVertical', 'sdTheme', 'sdAutoAdapt',
      'sdRotate', 'sdAverageSize', 'sdGradientColors', 'sdScaleByWeight',
      'useStereochemistry'
    ];

    // Check if any 2D settings changed
    const changed2DSettings = Object.keys(request.settings).filter(key => render2DSettings.includes(key));

    if (changed2DSettings.length > 0) {
      console.log('%c🎨 [Content] 2D settings changed - will re-render visible molecules:', 'color: #4CAF50; font-weight: bold;', changed2DSettings);
      // Only re-render 2D molecules if 2D settings changed
      // For now just update settings - images will re-render on next scan
    } else {
      console.log('%c🔧 [Content] 3D/other settings changed - no 2D re-render needed', 'color: #9C27B0;', Object.keys(request.settings));
      // 3D settings like compoundMolviewBgColor, molviewBioAssembly, etc.
      // don't require re-rendering existing 2D molecules
    }

    sendResponse({ success: true });
    return true;
  }
});

/**
 * Initialize the renderer
 */
function initializeRenderer() {
  log.inject('🔧 Initialize renderer - Starting chemical formula processing');

  // We process formulas directly using Unicode and CodeCogs
  // No MathJax injection needed - this avoids CSP issues entirely

  // Inject CSS for proper rendering FIRST
  injectStyles();

  // Setup lazy-loading BEFORE scanning (so loaders are available)
  // Performance mode ON: lazy-load images as they scroll into view
  // Performance mode OFF: load all images immediately
  if (settings.performanceMode) {
    log.inject('⚡ Performance mode ENABLED - Setting up lazy-loading for SVGs');
  } else {
    log.inject('⚡ Performance mode disabled - will load images immediately');
  }
  setupLazyLoading();

  // NOW scan and render (loaders are ready)
  log.inject('🚀 Starting initial page scan...');
  scanAndRender();

  log.inject('🔄 Setting up dynamic content observer...');
  observePageChanges();


  log.success('✅ Renderer initialized - formulas will be processed as they appear');
}

/**
 * Inject CSS styles for proper rendering
 */
function injectStyles() {
  const css = `.molecule-diagram{display:inline-block!important;position:static!important;max-width:100%;max-height:400px;margin:0 12px 8px 0!important;vertical-align:middle;padding:0!important;border:none!important;border-radius:4px;background-color:transparent!important;transition:opacity .3s ease;cursor:pointer;opacity:1}.molecule-diagram[data-loaded="error"]{background-color:transparent!important;border:1px solid rgba(239,83,80,0.3)!important}.molecule-loading{opacity:.1;background:linear-gradient(90deg,#f0f0f0 25%,#e0e0e0 50%,#f0f0f0 75%);background-size:200% 100%;animation:shimmer 2s infinite}@keyframes shimmer{0%{background-position:200% 0}100%{background-position:-200% 0}}.molecule-fadein{animation:fadeIn .4s ease-in-out}@keyframes fadeIn{from{opacity:0}to{opacity:1}}.molecule-container{display:inline-flex;flex-wrap:wrap;gap:12px;align-items:center;margin:8px 0}.molecule-container .molecule-diagram{margin:0}.molecule-container.vertical{display:flex;flex-direction:column;align-items:flex-start;gap:16px;margin:12px 0}.molecule-container.vertical .molecule-diagram{margin:0;display:block}.molecule-rotate-0{transform:rotate(0deg)}.molecule-rotate-90{transform:rotate(90deg);max-width:200px;max-height:300px}.molecule-rotate-180{transform:rotate(180deg)}.molecule-rotate-270{transform:rotate(270deg);max-width:200px;max-height:300px}.molecule-diagram:hover{background-color:rgba(100,150,255,.1)!important;box-shadow:0 0 8px rgba(100,150,255,.2)!important}.molecule-viewer-container{display:inline-block!important;margin:8px 12px 8px 0!important;vertical-align:middle;padding:0!important;background:transparent!important;border:none!important;box-shadow:none!important}.molecule-viewer-container:hover{background:transparent!important;border:none!important;box-shadow:none!important}.molecule-viewer-container .molecule-diagram{display:block;max-width:300px;max-height:200px;margin-bottom:6px;cursor:pointer;border-radius:4px}.molecule-viewer-container .molecule-diagram:hover{background-color:rgba(100,150,255,.15)!important;box-shadow:0 0 8px rgba(100,150,255,.3)!important}.molecule-download-btn:hover{background:#45a049!important;transform:scale(1.02)}.molecule-download-btn:active{transform:scale(.98)}@media(max-width:768px){.molecule-diagram{max-width:250px;max-height:150px}.molecule-container{flex-direction:column;gap:8px}.molecule-viewer-container{margin:8px 0!important}}`;

  const style = document.createElement('style');
  style.textContent = css;
  style.id = 'chem-renderer-styles';
  (document.head || document.documentElement).appendChild(style);
}


// Parse flags from chem: notation and return overrides
// NEW SYNTAX:
//   chem:smiles=CCO:           → Direct SMILES rendering
//   chem:mol=benzene:          → Compound lookup
//   chem:biomol=rhinovirus:    → Biomolecule lookup (protein)
//   chem:mineral=quartz:       → Mineral lookup (crystal)
//   chem:pbdid=4RHV:           → Direct PDB ID lookup (for biomolecules)
//   chem:cid=1234:             → Direct compound ID lookup
//   chem:codid=1234567:        → Direct mineral ID lookup
//   chem:<name>smiles=CCO:     → Custom name + SMILES (name appears in tag)
//   chem:mol=benzene+c+n-o:    → Flags without +d (completely overrides settings)
//   chem:mol=benzene+d+c+n-o:  → Flags with +d (applies to default settings)
// FLAGS:
// FLAGS:
//   c = show carbons, n = atom numbers, o = aromatic rings
//   h = show hydrogens, p = flip horizontal, q = flip vertical
//   +flag to enable, -flag to disable
//   +d = use defaults then apply flags (vs completely override)
// NOTE: pbdid does not support flags (biomolecules have different rendering)
function parseChemFlags(chemFormula, enableAIFlagControl = true) {
  const result = {
    useDefaults: false,           // +d flag: Whether to use the user's defaults as base
    flagsLocked: false,           // Whether flags are explicitly set (prevents popup from overriding)
    showCarbons: undefined,
    aromaticCircles: undefined,
    showMethyls: undefined,
    atomNumbers: undefined,
    flipHorizontal: undefined,
    flipVertical: undefined,
    addHydrogens: undefined,
    showImplicitHydrogens: undefined,
    size: undefined,
    rotation: undefined,
    isCompound: false,
    compoundType: undefined,      // 'compound', 'biomolecule', 'mineral', 'smiles', 'quantum', or 'hybrid'
    isDirectSmiles: false,
    moleculeName: null,           // The extracted molecule name/value
    isDirectID: false,            // Whether this is a direct ID lookup (pbdid, cid, codid)
    idType: undefined,            // 'pdbid', 'cid', or 'codid'
    idValue: undefined,           // The actual ID value
    smilesValue: undefined,       // The actual SMILES string (for named syntax)
    displayName: undefined,       // Custom display name for the molecule tag
    isIUPAC: false,               // Whether this is an IUPAC name (use OPSIN)
    // Quantum orbital properties
    isQuantumOrbital: false,      // Whether this is a quantum orbital visualization
    isHybridOrbital: false,       // Whether this is a hybrid orbital (sp, sp2, sp3, etc.)
    orbitalElement: undefined,    // Element symbol (H, C, O, Fe, etc.)
    orbitalType: undefined,       // Orbital type (1s, 2p, 3d, 4f, etc.)
    orbitalSubtype: undefined,    // Sub-orbital (px, py, pz, dxy, dz2, etc.)
    hybridType: undefined         // Hybrid type (sp, sp2, sp3, sp3d, sp3d2)
  };

  if (!chemFormula) return result;

  // ===== NEW SYNTAX: chem:type=value: OR chem:<name>type=value: format =====
  // Match: chem:smiles=CCO:, chem:mol=benzene+c+n:, chem:biomol=hemoglobin:, chem:mineral=quartz:
  // ALSO: chem:<name>smiles=CCO: where <name> becomes the molecule tag
  // Example: chem:ethanolsmiles=CCO: → name="ethanol", type="smiles", value="CCO"
  // NEW AI-FRIENDLY FORMATS:
  // chem:pbdid=4RHV: → Direct PDB ID lookup (for biomolecules)
  // chem:cid=1234: → Direct compound ID lookup
  // chem:codid=1234567: → Direct mineral ID lookup

  // First try the standard format: chem:type=value: (check this FIRST to avoid false named matches)
  const newSyntaxMatch = chemFormula.match(/chem:(smiles|mol|biomol|mineral|iupac|pbdid|cid|codid|quan|hybrid)=([^:]+):/i);

  // Only try named format if standard didn't match
  // Named format: chem:<name>type=value: where name can start with letter OR number
  // Examples: chem:ethanolsmiles=CCO:, chem:1-methyl-propanesmiles=CCC(C):, chem:Aspirincid=2244:
  const namedSyntaxMatch = !newSyntaxMatch ?
    chemFormula.match(/chem:([a-zA-Z0-9][a-zA-Z0-9\-,\s\(\)]{0,}?)(smiles|mol|biomol|mineral|iupac|pbdid|cid|codid|quan|hybrid)=([^:]+):/i) : null;

  // Handle standard syntax first (chem:type=value:)
  if (newSyntaxMatch) {
    const type = newSyntaxMatch[1].toLowerCase();
    const valueWithFlags = newSyntaxMatch[2];

    // Set compound type based on prefix
    switch (type) {
      case 'smiles':
        result.compoundType = 'smiles';
        result.isDirectSmiles = true;
        break;
      case 'mol':
        result.compoundType = 'compound';
        result.isCompound = true;
        break;
      case 'biomol':
        result.compoundType = 'biomolecule';
        break;
      case 'mineral':
        result.compoundType = 'mineral';
        result.isMineral = true;
        break;
      case 'iupac':
        result.compoundType = 'compound';
        result.isIUPAC = true;
        result.isCompound = true;
        break;
      case 'pbdid':
        result.compoundType = 'biomolecule';
        result.isDirectID = true;
        result.idType = 'pdbid';
        break;
      case 'cid':
        result.compoundType = 'compound';
        result.isDirectID = true;
        result.idType = 'cid';
        result.isCompound = true;
        break;
      case 'codid':
        result.compoundType = 'mineral';
        result.isDirectID = true;
        result.idType = 'codid';
        break;
      case 'quan':
        result.compoundType = 'quantum';
        result.isQuantumOrbital = true;
        break;
      case 'hybrid':
        result.compoundType = 'hybrid';
        result.isHybridOrbital = true;
        break;
    }

    // Extract molecule name/ID (everything before flags)
    // CRITICAL: Don't split on hyphens that are part of compound names like "butan-2-ol"
    // Flags are: +c, +n, -c, -o, etc. (single letter after +/-)
    // But "butan-2-ol" has "-2" and "-ol" which are NOT flags
    let moleculeName = valueWithFlags;

    // For direct ID types (pbdid, cid, codid), don't parse flags for pbdid
    // pbdid doesn't support flags (biomolecules use different rendering)
    // cid and codid support flags like normal molecules
    const isBiomoleculeID = type === 'pbdid';

    // Find the first ACTUAL flag: +letter or -letter where it's NOT part of the name
    // Look for patterns like +c, +n, -c, -o at the end, separated from the name
    // Skip flag parsing for biomolecule IDs
    const flagMatch = !isBiomoleculeID ? valueWithFlags.match(/^(.+?)([+\-][a-zA-Z](?:[+\-][a-zA-Z])*)$/) : null;
    if (flagMatch) {
      moleculeName = flagMatch[1];
      result.flagsLocked = true;
    }
    result.moleculeName = moleculeName.trim();

    // For direct SMILES, store the SMILES value
    if (type === 'smiles') {
      result.smilesValue = moleculeName.trim();
    }

    // For direct IDs, store the ID value
    if (type === 'pbdid' || type === 'cid' || type === 'codid') {
      result.idValue = moleculeName.trim();
    }

    // Parse flags (if any were found and not a biomolecule ID)
    // SKIP flag parsing if AI Flag Control is disabled
    const flagsPart = flagMatch ? flagMatch[2] : '';
    if (!enableAIFlagControl) {
      // AI Flag Control is disabled - don't parse flags, just return basic info
      console.log('%c🚫 AI Flag Control DISABLED - ignoring flags in text', 'color: #FF9800; font-weight: bold;');
      return result;
    }
    if (flagsPart.includes('+d')) result.useDefaults = true;

    const plusFlags = flagsPart.match(/\+([a-zA-Z0-9]+)/g);
    if (plusFlags) {
      plusFlags.forEach(flag => {
        const f = flag.substring(1).toLowerCase();
        if (f === 'c') result.showCarbons = true;
        if (f === 'o') result.aromaticCircles = true;
        if (f === 'm') result.showMethyls = true;
        if (f === 'n') result.atomNumbers = true;
        if (f === 'h') result.addHydrogens = true;
        if (f === 'i') result.showImplicitHydrogens = true;
        if (f === 'p') result.flipHorizontal = true;
        if (f === 'q') result.flipVertical = true;
        if (f === 'd') result.useDefaults = true;
        if (f.startsWith('s') && f.length > 1) {
          const sizeValue = parseInt(f.substring(1));
          if (!isNaN(sizeValue)) result.size = sizeValue / 100;
        }
      });
    }

    const minusFlags = flagsPart.match(/-([a-zA-Z])/g);
    if (minusFlags) {
      minusFlags.forEach(flag => {
        const f = flag.substring(1).toLowerCase();
        if (f === 'c') result.showCarbons = false;
        if (f === 'o') result.aromaticCircles = false;
        if (f === 'm') result.showMethyls = false;
        if (f === 'n') result.atomNumbers = false;
        if (f === 'h') result.addHydrogens = false;
        if (f === 'i') result.showImplicitHydrogens = false;
        if (f === 'p') result.flipHorizontal = false;
        if (f === 'q') result.flipVertical = false;
      });
    }

    console.log('%c📦 Standard syntax parsed:', 'color: #2196F3;', { type, value: valueWithFlags, result });
    return result;
  }

  // Handle named syntax (chem:<name>type=value:)
  if (namedSyntaxMatch) {
    const displayName = namedSyntaxMatch[1].trim(); // e.g., "ethanol"
    const type = namedSyntaxMatch[2].toLowerCase(); // e.g., "smiles"
    const valueWithFlags = namedSyntaxMatch[3]; // e.g., "CCO+c"

    // Set compound type based on prefix
    switch (type) {
      case 'smiles':
        result.compoundType = 'smiles';
        result.isDirectSmiles = true;
        break;
      case 'mol':
        result.compoundType = 'compound';
        result.isCompound = true;
        break;
      case 'biomol':
        result.compoundType = 'biomolecule';
        break;
      case 'mineral':
        result.compoundType = 'mineral';
        result.isMineral = true;
        break;
      case 'iupac':
        result.compoundType = 'compound';
        result.isIUPAC = true;
        result.isCompound = true;
        break;
      case 'pbdid':
        result.compoundType = 'biomolecule';
        result.isDirectID = true;
        result.idType = 'pdbid';
        break;
      case 'cid':
        result.compoundType = 'compound';
        result.isDirectID = true;
        result.idType = 'cid';
        result.isCompound = true;
        break;
      case 'codid':
        result.compoundType = 'mineral';
        result.isDirectID = true;
        result.idType = 'codid';
        result.isMineral = true;
        break;
      case 'quan':
        result.compoundType = 'quantum';
        result.isQuantumOrbital = true;
        break;
      case 'hybrid':
        result.compoundType = 'hybrid';
        result.isHybridOrbital = true;
        break;
    }


    // Use the display name for the tag (not the SMILES/value string)
    result.moleculeName = displayName;

    // For SMILES, IDs, IUPAC: store the actual value separately (before flags)
    // CRITICAL: Don't just search for + or - as they might be part of the name (e.g., "2,4,6-trinitrotoluene")
    // Instead, search for actual flag patterns: +[letter] or -[letter]
    const isBiomoleculeID = type === 'pbdid';
    let firstFlagIndex = -1;
    if (!isBiomoleculeID) {
      // Search for flag patterns: +c, +o, +m, +n, +h, +i, +p, +q, +d, +s123
      // CRITICAL: Flags MUST be followed by: end of string, another +/-, or be the last char before :
      // This prevents -m in "1-methyl-propane" from being treated as a flag (it's followed by 'e')
      // Valid flag examples: benzene+c, benzene+c+n, benzene-o
      // Invalid flag examples: 1-methyl (the -m is followed by 'e', not a flag terminator)
      const flagPattern = /[+-](c|o|m|n|h|i|p|q|d|s)(?=[+-]|$)/gi;
      const match = valueWithFlags.match(flagPattern);
      if (match) {
        // Find the position of the FIRST valid flag
        firstFlagIndex = valueWithFlags.search(flagPattern);
      }
    }
    const actualValue = firstFlagIndex !== -1 ? valueWithFlags.substring(0, firstFlagIndex).trim() : valueWithFlags.trim();

    if (type === 'smiles') {
      result.smilesValue = actualValue;
    } else if (type === 'iupac') {
      result.iupacValue = actualValue;  // Store the actual IUPAC name (e.g., "2,4,6-trinitrotoluene")
    } else if (type === 'mineral') {
      result.mineralValue = actualValue;  // Store the actual mineral name
    } else if (type === 'biomol') {
      result.biomolValue = actualValue;  // Store the actual biomolecule name
    } else if (type === 'mol') {
      result.molValue = actualValue;  // Store the actual compound name
    } else if (type === 'pbdid' || type === 'cid' || type === 'codid') {
      result.idValue = actualValue;
    }

    if (firstFlagIndex !== -1) {
      result.flagsLocked = true;
    }

    // Store the custom display name for tag rendering
    result.displayName = displayName;

    // Parse flags from the value part (skip for biomolecule IDs)
    // SKIP flag parsing if AI Flag Control is disabled
    const flagsPart = firstFlagIndex !== -1 ? valueWithFlags.substring(firstFlagIndex) : '';

    if (!enableAIFlagControl) {
      // AI Flag Control is disabled - don't parse flags, just return basic info
      console.log('%c🚫 AI Flag Control DISABLED - ignoring flags in text', 'color: #FF9800; font-weight: bold;');
      return result;
    }

    // Check for +d flag (use defaults as base)
    if (flagsPart.includes('+d')) {
      result.useDefaults = true;
    }

    // Parse + flags (enable)
    const plusFlags = flagsPart.match(/\+([a-zA-Z0-9]+)/g);
    if (plusFlags) {
      plusFlags.forEach(flag => {
        const f = flag.substring(1).toLowerCase();
        if (f === 'c') result.showCarbons = true;
        if (f === 'o') result.aromaticCircles = true;
        if (f === 'm') result.showMethyls = true;
        if (f === 'n') result.atomNumbers = true;
        if (f === 'h') result.addHydrogens = true;
        if (f === 'i') result.showImplicitHydrogens = true;
        if (f === 'p') result.flipHorizontal = true;
        if (f === 'q') result.flipVertical = true;
        if (f === 'd') result.useDefaults = true;
        // Size: +s150 = 1.5x scale
        if (f.startsWith('s') && f.length > 1) {
          const size = parseInt(f.substring(1));
          if (!isNaN(size)) result.size = size / 100;
        }
        // Rotation: +r45 = 45 degrees
        if (f.startsWith('r') && f.length > 1) {
          const rot = parseInt(f.substring(1));
          if (!isNaN(rot)) result.rotation = rot;
        }
      });
    }

    // Parse - flags (disable)
    const minusFlags = flagsPart.match(/-([a-zA-Z])/g);
    if (minusFlags) {
      minusFlags.forEach(flag => {
        const f = flag.substring(1).toLowerCase();
        if (f === 'c') result.showCarbons = false;
        if (f === 'o') result.aromaticCircles = false;
        if (f === 'm') result.showMethyls = false;
        if (f === 'n') result.atomNumbers = false;
        if (f === 'h') result.addHydrogens = false;
        if (f === 'i') result.showImplicitHydrogens = false;
        if (f === 'p') result.flipHorizontal = false;
        if (f === 'q') result.flipVertical = false;
      });
    }

    console.log('%c🏷️ Named syntax parsed:', 'color: #9C27B0;', { displayName, type, value: valueWithFlags, smilesValue: result.smilesValue, result });
    return result;
  }

  // ===== LEGACY SYNTAX: chem:name/flags: or chem:name+flags: =====
  // Example: chem:histamine/+c+o: or chem:histamine+c+o:

  // Check if it contains flags with / separator
  if (chemFormula.includes('/')) {
    const parts = chemFormula.split('/');
    if (parts.length > 1) {
      const flagsPart = parts[1];
      if (flagsPart.includes(':')) {
        const flags = flagsPart.split(':')[0].toLowerCase();

        if (flags.includes('d')) {
          result.useDefaults = true;
          const remaining = flags.replace('d', '');
          if (remaining.includes('-c')) result.showCarbons = false;
          if (remaining.includes('-o')) result.aromaticCircles = false;
          if (remaining.includes('-m')) result.showMethyls = false;
          if (remaining.includes('-n')) result.atomNumbers = false;
          if (remaining.includes('-p')) result.flipHorizontal = true;
          if (remaining.includes('-q')) result.flipVertical = true;
        }
      }
    }
  }

  // Check for + flags anywhere in the formula
  if (chemFormula.includes('+')) {
    result.flagsLocked = true;  // Has explicit flags
    const flagMatches = chemFormula.match(/\+([a-zA-Z0-9]+)/g);
    if (flagMatches) {
      flagMatches.forEach(flag => {
        const f = flag.substring(1).toLowerCase();
        if (f === 'c') result.showCarbons = true;
        if (f === 'o') result.aromaticCircles = true;
        if (f === 'm') result.showMethyls = true;
        if (f === 'n') result.atomNumbers = true;
        if (f === 'h') result.addHydrogens = true;
        if (f === 'i') result.showImplicitHydrogens = true;
        if (f === 'p') result.flipHorizontal = true;
        if (f === 'q') result.flipVertical = true;
        if (f === 'd') result.useDefaults = true;
        if (f === 'compound') result.isCompound = true;
        if (f === 'compound') result.compoundType = 'compound';
        if (f === 'biomolecule' || f === 'bio' || f === 'protein') result.compoundType = 'biomolecule';
        if (f === 'mineral' || f === 'crystal') result.compoundType = 'mineral';
        if (f === 'smiles' || f === 'smi') result.isDirectSmiles = true;
        if (f.startsWith('s') && f.length > 1) {
          const size = parseInt(f.substring(1));
          if (!isNaN(size)) result.size = size / 100;
        }
        if (f.startsWith('r') && f.length > 1) {
          const rot = parseInt(f.substring(1));
          if (!isNaN(rot)) result.rotation = rot;
        }
      });
    }
  }

  // LEGACY CODE REMOVED: The old code that looked for -[letter] anywhere in the formula
  // was incorrectly interpreting hyphens in chemical names (like "2,5-hexanedione") as flags.
  // Now flags are ONLY recognized when using explicit type=value syntax with trailing flags.
  // Example: chem:mol=benzene-c+n: has flags (-c, +n), but chem:mol=2,5-hexanedione: does NOT.

  return result;
}


// Helper function to strip flags from molecule name
// Example: "phenol/d-c" → "phenol", "histamine/+c+o" → "histamine"
function stripFlagsFromName(moleculeName) {
  if (!moleculeName) return moleculeName;

  // If there's a slash, everything before it is the molecule name
  const slashIndex = moleculeName.indexOf('/');
  if (slashIndex !== -1) {
    return moleculeName.substring(0, slashIndex);
  }

  return moleculeName;
}

function applyFlagOverrides(baseSettings, flagOverrides) {
  try {
    const newSettings = { ...baseSettings };

    // Apply overrides based on flags
    if (flagOverrides?.showCarbons !== undefined) newSettings.m2cfShowCarbons = flagOverrides.showCarbons;
    if (flagOverrides?.aromaticCircles !== undefined) newSettings.m2cfAromaticCircles = flagOverrides.aromaticCircles;
    if (flagOverrides?.showMethyls !== undefined) newSettings.m2cfShowMethyls = flagOverrides.showMethyls;
    if (flagOverrides?.atomNumbers !== undefined) newSettings.m2cfAtomNumbers = flagOverrides.atomNumbers;
    if (flagOverrides?.addHydrogens !== undefined) newSettings.m2cfAddH2 = flagOverrides.addHydrogens;
    if (flagOverrides?.flipHorizontal !== undefined) {
      newSettings.m2cfFlipHorizontal = flagOverrides.flipHorizontal;
      newSettings.flipHorizontal = flagOverrides.flipHorizontal;
    }
    if (flagOverrides?.flipVertical !== undefined) {
      newSettings.m2cfFlipVertical = flagOverrides.flipVertical;
      newSettings.flipVertical = flagOverrides.flipVertical;
    }
    // invert removed - now auto-applied based on dark mode detection

    if (flagOverrides?.rotation !== undefined) {
      newSettings.m2cfRotate = flagOverrides.rotation;
      newSettings.mvRotate = flagOverrides.rotation;
    }

    return newSettings;
  } catch (error) {
    console.error('Error applying flag overrides:', error);
    // Return original settings to not break functionality
    return baseSettings;
  }
}

/**
 * Wrap an image in a rotation container that properly handles layout
 * This prevents rotated images from overlapping content below them
 * @param {HTMLElement} img - The image element to wrap
 * @param {number} rotation - Rotation angle in degrees
 * @param {string} moleculeName - Name of the molecule for tooltip
 * @param {object} moleculeData - Full molecule data for 3D viewer
 * @returns {HTMLElement} The wrapper element
 */
function wrapImageWithRotationContainer(img, rotation, moleculeName, moleculeData) {
  // Create a MINIMAL wrapper - only for controls positioning, NO size constraints whatsoever
  const wrapper = document.createElement('div');
  wrapper.className = 'molecule-container';
  wrapper.style.cssText = `
    display: inline-block;
    position: relative;
    vertical-align: middle;
    margin: 0 12px 8px 0;
    background: none;
    background-color: transparent;
    border: none;
    box-shadow: none;
  `;
  // Store molecule data for later use when adding controls
  wrapper._moleculeName = moleculeName;
  wrapper._moleculeData = moleculeData;

  // Add the image to wrapper
  wrapper.appendChild(img);

  // NOTE: Hover controls (3D button, name tag, size arrows) are added AFTER image loads
  // to ensure proper positioning. See replaceWithLoaded() in renderClientSide().

  return wrapper;
}

/**
 * Add hover controls (molecule name + 3D viewer button) to an image
 * @param {HTMLElement} element - The element to add controls to
 * @param {string} moleculeName - Name of the molecule
 * @param {object} moleculeData - Full molecule data for 3D viewer
 */
function addHoverControls(element, moleculeName, moleculeData) {
  console.log('%c🔘 addHoverControls called', 'background: #9C27B0; color: white; font-weight: bold; padding: 4px;', {
    moleculeName,
    compoundType: moleculeData?.compoundType,
    pdbid: moleculeData?.pdbid,
    flagsCompoundType: moleculeData?.flags?.compoundType
  });

  // Skip if already has controls
  if (element.querySelector('.molecule-name-overlay')) {
    console.log('⏭️ Skipping - already has controls');
    return;
  }
  // Create name overlay (bottom-right, always visible)
  const nameOverlay = document.createElement('div');
  nameOverlay.className = 'molecule-name-overlay';
  nameOverlay.textContent = moleculeName || 'Unknown';
  nameOverlay.style.cssText = `
    position: absolute;
    bottom: 4px;
    right: 4px;
    background: rgba(0, 0, 0, 0.75);
    color: white;
    padding: 6px 10px;
    border-radius: 4px;
    font-size: 11px;
    font-family: sans-serif;
    opacity: 0.7;
    transition: opacity 0.2s;
    pointer-events: none;
    z-index: 10;
    white-space: nowrap;
    line-height: 1.4;
  `;

  // Check if this is a biomolecule (for showing toggle option)
  const isBiomolecule = moleculeData.compoundType === 'biomolecule' ||
    moleculeData.pdbid ||
    (moleculeData.flags && moleculeData.flags.compoundType === 'biomolecule');

  console.log('%c🧬 isBiomolecule:', 'color: #E91E63; font-weight: bold;', isBiomolecule);

  // SVG Icons
  const cubeIconSVG = `<svg viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="1.5" stroke-linecap="round" stroke-linejoin="round" width="16" height="16" style="display:block;">
    <path d="M21 16V8a2 2 0 0 0-1-1.73l-7-4a2 2 0 0 0-2 0l-7 4A2 2 0 0 0 3 8v8a2 2 0 0 0 1 1.73l7 4a2 2 0 0 0 2 0l7-4A2 2 0 0 0 21 16z"/>
    <polyline points="3.27 6.96 12 12.01 20.73 6.96"/>
    <line x1="12" y1="22.08" x2="12" y2="12"/>
  </svg>`;

  const asteriskIconSVG = `<svg viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="2.5" stroke-linecap="round" width="16" height="16" style="display:block;">
    <path d="M12 2v20"/>
    <path d="M2.5 7l19 10"/>
    <path d="M2.5 17l19-10"/>
  </svg>`;

  const atomIconSVG = `<svg viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="1.5" width="16" height="16" style="display:block;">
    <ellipse cx="12" cy="12" rx="10" ry="4" transform="rotate(45 12 12)"/>
    <ellipse cx="12" cy="12" rx="10" ry="4" transform="rotate(-45 12 12)"/>
    <ellipse cx="12" cy="12" rx="10" ry="4"/>
    <circle cx="12" cy="12" r="2" fill="currentColor"/>
  </svg>`;

  // Create button container (merged button group)
  const buttonGroup = document.createElement('div');
  buttonGroup.className = 'molecule-3d-button-group';
  buttonGroup.style.cssText = `
    position: absolute;
    top: 4px;
    right: 4px;
    display: flex;
    opacity: 0;
    transition: opacity 0.2s;
    z-index: 10;
    border-radius: 4px;
    overflow: hidden;
    background: rgba(0, 0, 0, 0.85);
  `;

  // Bio toggle button (asterisk/atom) - only for biomolecules
  let bioToggleBtn = null;
  // Default to MolView (useMolstar = false), UNLESS bio assembly mode is enabled
  let useMolstar = settings.molviewBioAssembly === true;

  if (isBiomolecule) {
    bioToggleBtn = document.createElement('button');
    bioToggleBtn.className = 'molecule-bio-toggle';
    // Show atom icon for MolView (default), asterisk for Mol*
    bioToggleBtn.innerHTML = useMolstar ? asteriskIconSVG : atomIconSVG;
    bioToggleBtn.title = useMolstar
      ? 'Mol* viewer (click to switch to MolView)'
      : 'MolView (click to switch to Mol*)';
    bioToggleBtn.style.cssText = `
      background: transparent;
      color: white;
      border: none;
      width: 28px;
      height: 28px;
      padding: 0;
      cursor: pointer;
      display: flex;
      align-items: center;
      justify-content: center;
      transition: background 0.2s;
      border-right: 1px solid rgba(255,255,255,0.3);
      box-sizing: border-box;
    `;

    bioToggleBtn.addEventListener('mouseenter', () => {
      bioToggleBtn.style.background = 'rgba(255,255,255,0.15)';
    });
    bioToggleBtn.addEventListener('mouseleave', () => {
      bioToggleBtn.style.background = 'transparent';
    });

    bioToggleBtn.addEventListener('click', async (e) => {
      e.stopPropagation();
      useMolstar = !useMolstar;

      // Update icon and title
      if (useMolstar) {
        bioToggleBtn.innerHTML = asteriskIconSVG;
        bioToggleBtn.title = 'Mol* viewer (click to switch to MolView)';
      } else {
        bioToggleBtn.innerHTML = atomIconSVG;
        bioToggleBtn.title = 'MolView (click to switch to Mol*)';
      }

      console.log('🔄 Bio viewer toggle:', useMolstar ? 'Mol*' : 'MolView');

      // Store preference in moleculeData
      moleculeData.preferMolstar = useMolstar;

      // Find existing 3D viewer iframe and switch it
      // The 3D viewer might be the element itself, its parent, or a sibling
      let existing3DViewer = element.closest('.molecule-3d-viewer');
      if (!existing3DViewer) {
        existing3DViewer = element.querySelector('.molecule-3d-viewer');
      }
      if (!existing3DViewer && element.classList.contains('molecule-3d-viewer')) {
        existing3DViewer = element;
      }

      // Update data attribute for persistence
      if (existing3DViewer && useMolstar) {
        existing3DViewer.setAttribute('data-viewer-type', 'molstar');
      } else if (existing3DViewer) {
        existing3DViewer.setAttribute('data-viewer-type', 'molview');
      }

      console.log('🔍 Looking for 3D viewer:', {
        element: element.className,
        found: !!existing3DViewer,
        pdbid: moleculeData.pdbid
      });

      if (existing3DViewer) {
        // Already showing 3D viewer - create NEW iframe to avoid CSP issues
        const oldIframe = existing3DViewer.querySelector('iframe');
        console.log('🔍 Found iframe:', !!oldIframe);

        if (oldIframe && moleculeData.pdbid) {
          const pdbId = moleculeData.pdbid.toLowerCase();
          let newUrl;
          if (useMolstar) {
            newUrl = `https://molstar.org/viewer/?pdb=${pdbId}&hide-controls=1`;
            console.log('%c🔄 Switching to Mol*:', 'color: #E91E63; font-weight: bold;', newUrl);
          } else {
            // Rebuild MolView URL with current settings to ensure parameters are applied
            newUrl = buildMolViewEmbedUrl(pdbId, 'pdbid', 'biomolecule');
            console.log('%c🔄 Switching to MolView:', 'color: #9C27B0; font-weight: bold;', newUrl);
          }

          // Create NEW iframe to avoid CSP (same approach as automatic fallback)
          const newIframe = document.createElement('iframe');
          newIframe.src = newUrl;
          newIframe.style.cssText = oldIframe.style.cssText;
          newIframe.setAttribute('allow', oldIframe.getAttribute('allow') || 'accelerometer; autoplay; clipboard-write; encrypted-media; gyroscope; picture-in-picture; web-share; fullscreen; xr-spatial-tracking');
          newIframe.setAttribute('allowfullscreen', 'true');
          newIframe.setAttribute('frameborder', '0');
          newIframe.title = useMolstar ? `Mol* Viewer: ${pdbId}` : `MolView: ${pdbId}`;

          // Replace old iframe with new one
          if (oldIframe.parentNode) {
            oldIframe.parentNode.replaceChild(newIframe, oldIframe);
            console.log('%c✅ Created new iframe!', 'color: green; font-weight: bold;');
          }
        } else {
          console.log('%c⚠️ Could not switch - missing iframe or pdbid', 'color: orange;');
        }
      } else {
        // Not showing 3D yet - just update the preference, next 3D click will use it
        console.log('📝 Preference saved for next 3D view');
      }
    });
  }

  // 3D Cube button (always present)
  const cubeBtn = document.createElement('button');
  cubeBtn.className = 'molecule-3d-cube-btn';
  cubeBtn.innerHTML = cubeIconSVG;
  cubeBtn.title = 'View in 3D';
  cubeBtn.style.cssText = `
    background: transparent;
    color: white;
    border: none;
    width: 28px;
    height: 28px;
    padding: 0;
    cursor: pointer;
    display: flex;
    align-items: center;
    justify-content: center;
    transition: background 0.2s;
    box-sizing: border-box;
  `;

  cubeBtn.addEventListener('mouseenter', () => {
    cubeBtn.style.background = 'rgba(255,255,255,0.15)';
  });
  cubeBtn.addEventListener('mouseleave', () => {
    cubeBtn.style.background = 'transparent';
  });

  cubeBtn.addEventListener('click', async (e) => {
    e.stopPropagation();
    console.log('🔮 3D cube button clicked!', { moleculeData, element, useMolstar });

    // Store current viewer preference
    if (isBiomolecule) {
      moleculeData.preferMolstar = useMolstar;
    }

    try {
      await show3DViewerInline(moleculeData, element);
    } catch (error) {
      console.error('❌ Error showing 3D viewer:', error);
      alert('Failed to load 3D viewer: ' + error.message);
    }
  });

  // Assemble button group
  if (bioToggleBtn) {
    buttonGroup.appendChild(bioToggleBtn);
  }
  buttonGroup.appendChild(cubeBtn);

  // Store reference to button group for hover visibility
  const viewer3DBtn = buttonGroup; // Alias for compatibility with existing hover code

  // bioViewerToggleBtn is obsolete now - set to null for compatibility
  const bioViewerToggleBtn = null;

  // Create size controls container (bottom-left, same as mol2chemfig)
  const controlsDiv = document.createElement('div');
  controlsDiv.className = 'chem-size-controls';
  controlsDiv.style.cssText = `
    position: absolute;
    bottom: 4px;
    left: 4px;
    display: flex;
    flex-direction: column;
    gap: 2px;
    opacity: 0;
    transition: opacity 0.2s;
    z-index: 1000;
  `;

  // FIX: Ensure button container positioning doesn't drift
  // The controlsDiv is absolutely positioned relative to the parent (element)
  // which might not have position: relative set if it's just an <img> tag wrapper
  if (element.style.position !== 'relative' && element.style.position !== 'absolute') {
    element.style.position = 'relative';
  }

  // FORCE REFLOW: Check position immediately to ensure browser has updated layout
  // This prevents the "slight offset" issue on initial spawn
  void element.offsetWidth;

  const upButton = document.createElement('button');
  upButton.className = 'chem-size-btn chem-size-up';
  upButton.innerHTML = '▲';
  upButton.title = 'Increase size';
  upButton.style.cssText = `
    width: 24px;
    height: 24px;
    border: none;
    background: rgba(0, 0, 0, 0.7);
    color: white;
    border-radius: 4px;
    cursor: pointer;
    font-size: 10px;
    display: flex;
    align-items: center;
    justify-content: center;
    transition: background 0.2s;
  `;

  const downButton = document.createElement('button');
  downButton.className = 'chem-size-btn chem-size-down';
  downButton.innerHTML = '▼';
  downButton.title = 'Decrease size';
  downButton.style.cssText = upButton.style.cssText;

  // Helper for resizing element
  function resizeElement(factor, minW = 100, minH = 80) {
    const img = element.querySelector('img') || element.querySelector('iframe');
    if (!img) return;
    const w = img.offsetWidth || img.width || img.naturalWidth || 300;
    const h = img.offsetHeight || img.height || img.naturalHeight || 200;
    const newW = Math.max(minW, Math.round(w * factor));
    const newH = Math.max(minH, Math.round(h * factor));
    const setSize = (el) => {
      el.style.setProperty('width', `${newW}px`, 'important');
      el.style.setProperty('height', `${newH}px`, 'important');
      el.style.setProperty('max-width', 'none', 'important');
      el.style.setProperty('max-height', 'none', 'important');
    };
    setSize(img);
    element.style.setProperty('width', `${newW}px`, 'important');
    element.style.setProperty('height', `${newH}px`, 'important');
  }

  upButton.addEventListener('click', (e) => { e.stopPropagation(); resizeElement(1.2); });
  downButton.addEventListener('click', (e) => { e.stopPropagation(); resizeElement(0.8333); });

  controlsDiv.appendChild(upButton);
  controlsDiv.appendChild(downButton);

  // Add hover effect to parent element
  element.addEventListener('mouseenter', () => {
    nameOverlay.style.opacity = '1';
    viewer3DBtn.style.opacity = '1';
    controlsDiv.style.opacity = '1';
    if (bioViewerToggleBtn) bioViewerToggleBtn.style.opacity = '1';
  });

  element.addEventListener('mouseleave', () => {
    nameOverlay.style.opacity = '0.7';
    viewer3DBtn.style.opacity = '0';
    controlsDiv.style.opacity = '0';
    if (bioViewerToggleBtn) bioViewerToggleBtn.style.opacity = '0';
  });

  // Make element positioned so absolute children work
  if (getComputedStyle(element).position === 'static') {
    element.style.position = 'relative';
  }

  // Ensure container has NO background and is transparent
  element.style.background = 'none';
  element.style.backgroundColor = 'transparent';
  element.style.border = 'none';
  element.style.boxShadow = 'none';

  // DO NOT set explicit width/height on container - let it size naturally to content
  // The absolute positioned controls will overlay the image without needing container sizing

  element.appendChild(nameOverlay);
  element.appendChild(viewer3DBtn);
  if (bioViewerToggleBtn) element.appendChild(bioViewerToggleBtn);
  element.appendChild(controlsDiv);
}

/**
 * Build RCSB structure image URL directly from PDB ID
 * Server just returns pdbid, extension builds the image URL
 * 
 * @param {string} pdbid - The PDB ID
 * @returns {string} Direct RCSB image URL
 */
function buildBioImageUrl(pdbid) {
  const pdbLower = pdbid.toLowerCase();
  // Direct RCSB CDN URL - server just provides pdbid, extension builds URL
  const imageType = settings.molviewBioAssembly ? 'assembly-1' : 'model-1';
  return `https://cdn.rcsb.org/images/structures/${pdbLower}_${imageType}.jpeg`;
}


/**
 * Build a MolView embed URL with all settings parameters
 * Uses current settings for chain type, chain color, background, etc.
 * 
 * @param {string} baseIdentifier - pdbid, codid, cid, or smiles value
 * @param {string} identifierType - 'pdbid', 'codid', 'cid', or 'smiles'
 * @param {string} moleculeType - 'biomolecule' or 'mineral'
 * @returns {string} Complete MolView embed URL with parameters
 */
function buildMolViewEmbedUrl(baseIdentifier, identifierType, moleculeType = 'compound') {
  const baseUrl = 'https://embed.molview.org/v1/';
  const params = new URLSearchParams();

  // Add the base identifier
  params.set(identifierType, baseIdentifier);

  // For proteins/biomolecules, add chain settings
  if (moleculeType === 'biomolecule') {
    // Chain representation: ribbon, cylinders, btube, ctrace, bonds
    if (settings.molviewChainType) {
      params.set('chainType', settings.molviewChainType);
    }

    // Chain bonds
    if (settings.molviewChainBonds) {
      params.set('chainBonds', 'true');
    }

    // Chain color: ss, spectrum, chain, residue, polarity, bfactor
    if (settings.molviewChainColor) {
      params.set('chainColor', settings.molviewChainColor);
    }

    // Background color for proteins
    if (settings.proteinMolviewBgColor) {
      params.set('bg', settings.proteinMolviewBgColor);
    }

    // Note: Biological assembly is not supported by embed.molview.org
    // The bioassembly setting only affects the 2D biomolecule image (assembly vs model)
  }

  // For minerals, add mineral-specific settings
  if (moleculeType === 'mineral') {
    // Representation mode: balls, stick, vdw, wireframe, line
    if (settings.mineralRepresentation) {
      // Map our values to MolView mode values
      const modeMap = {
        'ballAndStick': 'balls',
        'stick': 'stick',
        'vdwSpheres': 'vdw',
        'wireframe': 'wireframe',
        'line': 'line'
      };
      params.set('mode', modeMap[settings.mineralRepresentation] || 'balls');
    }

    // Background color for minerals
    if (settings.mineralMolviewBgColor) {
      params.set('bg', settings.mineralMolviewBgColor);
    }
  }

  // For compounds with SMILES or CID, add compound-specific settings
  if (identifierType === 'smiles' || identifierType === 'cid' || moleculeType === 'compound') {
    // Representation mode: balls, stick, vdw, wireframe, line
    if (settings.molviewRepresentation) {
      const modeMap = {
        'ballAndStick': 'balls',
        'stick': 'stick',
        'vdwSpheres': 'vdw',
        'wireframe': 'wireframe',
        'line': 'line'
      };
      params.set('mode', modeMap[settings.molviewRepresentation] || 'balls');
    } else if (!params.has('mode')) {
      params.set('mode', 'balls');
    }

    // Background color for compounds
    if (settings.compoundMolviewBgColor) {
      params.set('bg', settings.compoundMolviewBgColor);
    }
  }

  return baseUrl + '?' + params.toString();
}


/**
 * Parse orbital value from chem:quan=ELEMENT+ORBITAL: syntax
 * Examples: 
 *   "C+2p" → { element: 'C', n: 2, l: 1, m: 0, suborbital: null }
 *   "Fe+3dxy" → { element: 'Fe', n: 3, l: 2, m: 4, suborbital: 'xy' }
 *   "H+1s" → { element: 'H', n: 1, l: 0, m: 0, suborbital: null }
 * 
 * @param {string} value - The orbital value like "C+2p" or "Fe+3dxy"
 * @returns {object} Parsed orbital data
 */
function parseOrbitalValue(value) {
  // Safety check for undefined or null value
  if (!value || typeof value !== 'string') {
    console.warn('[Quantum] parseOrbitalValue called with invalid value:', value);
    return null;
  }

  // Match pattern: ELEMENT+NL or ELEMENT+NL+SUBORBITAL
  // Examples: H+1s, C+2p, C+2px, Fe+3d, Fe+3dxy, Fe+3dz2
  const match = value.match(/^([A-Za-z]{1,2})\+(\d+)([spdf])([a-z0-9]{0,3})?$/i);

  if (!match) {
    console.warn('[Quantum] Could not parse orbital value:', value);
    return null;
  }

  const element = match[1];
  const n = parseInt(match[2]);
  const orbitalLetter = match[3].toLowerCase();
  const suborbital = match[4] ? match[4].toLowerCase() : null;

  // Map orbital letters to l values
  const lMap = { s: 0, p: 1, d: 2, f: 3 };
  const l = lMap[orbitalLetter];

  // Determine m value based on suborbital
  // For Falstad, URL uses 0-indexed n (n=0 for shell 1)
  // m values in Falstad: px=1, py=2, pz=0, etc.
  let m = 0;
  if (suborbital) {
    // p orbitals: px=1, py=2, pz=0
    if (orbitalLetter === 'p') {
      const pMap = { x: 1, y: 2, z: 0 };
      m = pMap[suborbital] !== undefined ? pMap[suborbital] : 0;
    }
    // d orbitals: dz2=0, dxz=1, dyz=2, dx2-y2=3 (or dx2y2), dxy=4
    if (orbitalLetter === 'd') {
      const dMap = { 'z2': 0, 'xz': 1, 'yz': 2, 'x2y2': 3, 'x2-y2': 3, 'xy': 4 };
      m = dMap[suborbital] !== undefined ? dMap[suborbital] : 0;
    }
    // f orbitals: more complex, default to 0 for now
  }

  return { element, n, l, m, orbitalLetter, suborbital };
}


/**
 * Build Falstad quantum orbital viewer URL
 * 
 * @param {string} value - Orbital value like "C+2p" or "Fe+3dxy"
 * @param {boolean} isHybrid - Whether this is a hybrid orbital (sp, sp2, sp3)
 * @returns {string} Falstad viewer URL
 */
function buildQuantumOrbitalUrl(value, isHybrid = false) {
  console.log('[Quantum] buildQuantumOrbitalUrl called with:', { value, isHybrid });

  // Safety check for undefined or null value
  if (!value || typeof value !== 'string') {
    console.warn('[Quantum] buildQuantumOrbitalUrl called with invalid value:', value);
    // Return a default 1s orbital URL
    return 'https://falstad.com/qmatom/qmatom.html?vc=0&n=0&l=0&m=0';
  }

  const baseUrl = 'https://falstad.com/qmatom/qmatom.html';
  const params = new URLSearchParams();

  if (isHybrid) {
    // Hybrid orbital: sp, sp2, sp3, sp3d, sp3d2
    // Use Hybrid Bases mode (vc=10)
    params.set('vc', '10');

    // Parse hybrid type
    const hybridMatch = value.toLowerCase().match(/^(sp3d2|sp3d|sp3|sp2|sp)$/);
    if (hybridMatch) {
      const hybridType = hybridMatch[1];
      // For hybrid bases, we need to set n and l appropriately
      // sp, sp2, sp3 use n=2, l=1
      // sp3d uses n=3, l=2
      // sp3d2 uses n=3, l=2
      switch (hybridType) {
        case 'sp':
        case 'sp2':
        case 'sp3':
          params.set('n', '1');  // 0-indexed: n=1 → shell 2
          params.set('l', '1');  // l=1 for p involvement
          params.set('m', '0');
          break;
        case 'sp3d':
        case 'sp3d2':
          params.set('n', '2');  // 0-indexed: n=2 → shell 3
          params.set('l', '2');  // l=2 for d involvement  
          params.set('m', '0');
          break;
      }
    }
  } else {
    // Atomic orbital
    const parsed = parseOrbitalValue(value);
    if (!parsed) {
      // Fallback to 1s if parsing fails
      params.set('vc', '0');
      params.set('n', '0');
      params.set('l', '0');
      params.set('m', '0');
    } else {
      // Real orbitals mode (chemistry style)
      params.set('vc', '0');
      // Falstad uses 0-indexed n (n=0 for shell 1)
      params.set('n', String(parsed.n - 1));
      params.set('l', String(parsed.l));
      params.set('m', String(parsed.m));
    }
  }

  return baseUrl + '?' + params.toString();
}

window.parseChemFlags = parseChemFlags;
window.stripFlagsFromName = stripFlagsFromName;
window.applyFlagOverrides = applyFlagOverrides;
window.buildMolViewEmbedUrl = buildMolViewEmbedUrl;
window.buildQuantumOrbitalUrl = buildQuantumOrbitalUrl;
window.parseOrbitalValue = parseOrbitalValue;


/**
 * Setup Lazy-Loading for Molecule SVGs
 * Only renders visible SVGs, defers off-screen ones to save performance
 */
function setupLazyLoading() {
  // Disconnect any existing observer first to prevent memory leaks
  if (window._lazyLoadObserver) {
    try {
      window._lazyLoadObserver.disconnect();
    } catch (e) {
      // Ignore disconnect errors
    }
    window._lazyLoadObserver = null;
  }

  log.inject('🚀 Setting up Intersection Observer for lazy-loading with performance optimization');

  // Track loading count to prevent too many simultaneous loads
  let activeLoads = 0;
  const maxConcurrentLoads = 3;  // Reduced from 5 to prevent lag

  // Helper to decrement activeLoads and trigger queue processing
  // Use this instead of direct activeLoads-- to ensure queue is processed
  function decrementActiveLoads() {
    activeLoads--;
    // Trigger queue processing if there are waiting items
    if (typeof processLoadQueue === 'function') {
      setTimeout(processLoadQueue, 50);
    }
  }

  // Rendering handled by server
  // All rendering is now handled by the ChemTex Server via buildServerSvgUrl().
  // This reduces bundle size and prevents CSP issues on strict sites.

  // Intersection Observer for lazy-loading
  const observer = new IntersectionObserver((entries) => {
    entries.forEach(entry => {
      if (entry.isIntersecting) {
        const img = entry.target;
        // Only load if not already loaded
        if (img.dataset.loaded === 'false' && img.dataset.src && activeLoads < maxConcurrentLoads) {
          // Use requestIdleCallback to load during idle time (non-blocking)
          if (window.requestIdleCallback) {
            requestIdleCallback(() => {
              if (activeLoads < maxConcurrentLoads) {
                loadImage(img);
              }
            }, { timeout: 1000 });  // Timeout ensures it loads after 1 second even if busy
          } else {
            // Fallback for browsers without requestIdleCallback
            setTimeout(() => {
              if (activeLoads < maxConcurrentLoads) {
                loadImage(img);
              }
            }, 50);
          }
        }
      }
    });
  }, {
    rootMargin: '300px'  // Start loading 300px before entering viewport
  });

  // Helper function to load image with load tracking
  function loadImage(img) {
    activeLoads++;
    log.debug(`🖼️  Loading SVG (#${activeLoads}): ${img.dataset.src.substring(0, 60)}...`);

    // Add fade-in class for smooth appearance
    img.classList.add('molecule-fadein');
    img.classList.remove('molecule-loading');

    // Preload the image to avoid blocking the UI
    const preloadImg = new Image();
    preloadImg.onload = () => {
      img.src = img.dataset.src;
      img.dataset.loaded = 'true';
      observer.unobserve(img);
      decrementActiveLoads();
      log.debug(`✅ SVG loaded (${activeLoads} active loads remaining)`);
    };
    preloadImg.onerror = () => {
      decrementActiveLoads();
      log.debug(`❌ SVG failed to load (${activeLoads} active loads remaining)`);
    };
    preloadImg.src = img.dataset.src;
  }

  // Helper function to render biomolecule 2D image from server
  // Shows the biomolecule structure image with hover controls for 3D viewing
  async function renderBiomolecule2D(moleculeData, img) {
    console.log('%c🧬 RENDERBIOMOLECULE2D CALLED', 'background: #E91E63; color: white; font-weight: bold; padding: 4px;');

    try {
      // Create the 2D image element
      const bioImg = document.createElement('img');
      bioImg.src = moleculeData.imageUrl;
      bioImg.alt = moleculeData.nomenclature || 'Biomolecule';
      bioImg.className = 'molecule-diagram molecule-biomolecule';
      bioImg.crossOrigin = 'anonymous';

      // Get viewer size from settings
      const viewerSize = settings.viewer3DSize || 'normal';
      const sizeDimensions = {
        'small': { width: 200, height: 150 },
        'medium': { width: 300, height: 250 },
        'normal': { width: 400, height: 350 },
        'large': { width: 600, height: 450 },
        'xlarge': { width: 800, height: 600 }
      };
      const dimensions = sizeDimensions[viewerSize] || sizeDimensions['normal'];

      bioImg.style.cssText = `
        width: ${dimensions.width}px;
        height: ${dimensions.height}px;
        max-width: ${dimensions.width}px;
        max-height: ${dimensions.height}px;
        object-fit: contain;
        display: block;
        margin: 0;
        vertical-align: middle;
        cursor: pointer;
        border-radius: 8px;
      `;
      // Mark as fixed-size to prevent scaling in wrapImageWithSizeControls
      bioImg.dataset.fixedSize = 'true';

      // Apply background removal if enabled
      // Use background service worker to fetch image as base64 to bypass CORS
      if (settings.proteinRemoveWhiteBg) {
        try {
          // Fetch via background script to bypass CORS
          const imageBlob = await backgroundFetchBlob(moleculeData.imageUrl);
          if (imageBlob && imageBlob.base64) {
            // Create an image from the base64 data
            const tempImg = new Image();
            tempImg.onload = () => {
              try {
                const canvas = document.createElement('canvas');
                const ctx = canvas.getContext('2d');
                canvas.width = tempImg.naturalWidth;
                canvas.height = tempImg.naturalHeight;
                ctx.drawImage(tempImg, 0, 0);
                const imageData = ctx.getImageData(0, 0, canvas.width, canvas.height);
                const data = imageData.data;

                // Remove white background with smooth alpha gradient
                for (let i = 0; i < data.length; i += 4) {
                  const r = data[i], g = data[i + 1], b = data[i + 2];
                  // Calculate distance from pure white (255, 255, 255)
                  const dist = Math.sqrt(Math.pow(255 - r, 2) + Math.pow(255 - g, 2) + Math.pow(255 - b, 2));
                  if (dist < 15) {
                    // Pure white or near-white - make fully transparent
                    data[i + 3] = 0;
                  } else if (dist < 150) {
                    // Light colors - apply gradient transparency
                    data[i + 3] = Math.floor(255 * (dist - 15) / 135);
                  }
                  // Everything else keeps original alpha
                }
                ctx.putImageData(imageData, 0, 0);
                bioImg.src = canvas.toDataURL('image/png');
                console.log('%c✅ White background removed successfully', 'color: #4CAF50; font-weight: bold;');
              } catch (e) {
                console.warn('Canvas processing failed:', e.message);
              }
            };
            tempImg.src = `data:${imageBlob.type || 'image/jpeg'};base64,${imageBlob.base64}`;
          }
        } catch (e) {
          console.warn('Background removal failed:', e.message);
          // Fallback: just use original image (already set above)
        }
      }

      // Mark original img as loaded
      img.dataset.loaded = 'true';

      // Helper function to add controls after image is ready
      const addControlsAfterLoad = async () => {
        // Get settings and add controls
        chrome.storage.sync.get({
          saveSizePerImage: false,
          saveSizeBySMILES: true
        }, async (sizeSettings) => {
          try {
            console.log('%c🔘 Calling wrapImageWithSizeControls for biomolecule', 'background: #4CAF50; color: white; font-weight: bold;', moleculeData.nomenclature);
            await wrapImageWithSizeControls(bioImg, img, moleculeData, sizeSettings);
            console.log('%c✅ wrapImageWithSizeControls completed', 'color: #4CAF50; font-weight: bold;');
          } catch (err) {
            console.error('%c❌ wrapImageWithSizeControls FAILED:', 'background: red; color: white; font-weight: bold;', err);
          }
        });
      };

      // Wait for image to load before adding controls
      // This ensures the container has proper dimensions for button positioning
      if (bioImg.complete && bioImg.naturalWidth > 0) {
        // Image already loaded (cached or inline)
        await addControlsAfterLoad();
      } else {
        // Wait for image to load
        bioImg.onload = addControlsAfterLoad;
        bioImg.onerror = () => {
          console.error('Biomolecule image failed to load');
          // Still try to add controls even if image fails
          addControlsAfterLoad();
        };
      }


    } catch (error) {
      console.error('renderBiomolecule2D error:', error);
      img.alt = 'Failed to load biomolecule image';
      img.dataset.loaded = 'error';
    }
  }

  // NEW SIMPLIFIED renderClientSide - Server handles ALL lookups!
  // Flow:
  // - Compounds: /mol=benzene.svg → Server fetches SMILES + renders
  // - Minerals: /mineral=quartz.svg → Server fetches mineral data + renders  
  // - Biomolecules: /biomol=hemoglobin.json → Server returns pdbid → Extension shows 3D viewer
  // - Direct SMILES: /smiles=CCO.svg → Server renders directly
  // - IUPAC: /iupac=2-methylpropane.svg → Server fetches from OPSIN + renders
  async function renderClientSide(moleculeData, img) {
    activeLoads++;
    console.log('%c🎨 RENDER - Server handles ALL lookups!', 'background: #222; color: #00FF00; font-size: 16px; padding: 5px;');
    console.log('Molecule data:', moleculeData);

    try {
      // Decode molecule data if not passed directly
      if (!moleculeData && img.dataset.moleculeViewer) {
        moleculeData = JSON.parse(atob(img.dataset.moleculeViewer));
      }

      const flags = moleculeData.flags || {};

      // displayName is for UI display (e.g., "TNT")
      // The actual value for server lookup should come from specific fields
      const displayName = flags.displayName || moleculeData.nomenclature || moleculeData.name || 'molecule';

      // Determine the SERVER endpoint type based on flags
      let serverType = 'mol';  // Default: compound name lookup
      let lookupValue = displayName;  // Default: use display name for lookup

      if (flags.isDirectSmiles || moleculeData.smiles) {
        // Direct SMILES - no lookup needed
        serverType = 'smiles';
        lookupValue = flags.smilesValue || moleculeData.smiles || displayName;
      } else if (flags.compoundType === 'biomolecule' || flags.isBiomolecule) {
        // Biomolecule - server returns pdbid, we show 3D viewer
        serverType = 'biomol';
        lookupValue = flags.biomolValue || moleculeData.nomenclature || displayName;
      } else if (flags.compoundType === 'mineral' || flags.isMineral) {
        // Mineral - server fetches from COD
        serverType = 'mineral';
        lookupValue = flags.mineralValue || moleculeData.nomenclature || displayName;
      } else if (flags.compoundType === 'quantum' || flags.isQuantumOrbital) {
        // Quantum orbital - show Falstad viewer iframe
        serverType = 'quantum';
        lookupValue = moleculeData.nomenclature || displayName;
      } else if (flags.compoundType === 'hybrid' || flags.isHybridOrbital) {
        // Hybrid orbital - show Falstad viewer iframe in hybrid mode
        serverType = 'hybrid';
        lookupValue = moleculeData.nomenclature || displayName;
      } else if (flags.isIUPAC) {
        // IUPAC name - server fetches from OPSIN
        // CRITICAL: Use the actual IUPAC name, NOT the display name!
        serverType = 'iupac';
        lookupValue = flags.iupacValue || moleculeData.nomenclature || displayName;
        console.log('%c🔬 IUPAC lookup:', 'color: #FF5722; font-weight: bold;', { iupacValue: flags.iupacValue, lookupValue });
      } else {
        // Compound (default) - server fetches SMILES from PubChem
        lookupValue = flags.molValue || moleculeData.nomenclature || displayName;
      }

      console.log(`%c📍 Server endpoint: /${serverType}=`, 'color: #2196F3; font-weight: bold;', lookupValue, { displayName, lookupValue });
      console.log('%c🔍 DEBUG flags.compoundType:', 'color: orange;', flags.compoundType, 'flags.isQuantumOrbital:', flags.isQuantumOrbital);

      // ===== BIOMOLECULE: Get pdbid from server, show 3D viewer =====
      if (serverType === 'biomol') {
        try {
          const result = await fetchFromChemistryLaTeXServer('biomol', lookupValue);
          if (result.pdbid) {
            console.log(`%c🧬 Got pdbid from server: ${result.pdbid}`, 'color: #E91E63; font-weight: bold;');

            const biomoleculeData = {
              nomenclature: displayName,
              compoundType: 'biomolecule',
              embedUrl: buildMolViewEmbedUrl(result.pdbid, 'pdbid', 'biomolecule'),
              imageUrl: buildBioImageUrl(result.pdbid),
              pdbid: result.pdbid,
              is3D: moleculeData.is3D || moleculeData.show3D || flags.is3D
            };

            // Show 2D image or 3D viewer based on flags
            if (biomoleculeData.is3D) {
              await show3DViewerInline(biomoleculeData, img);
            } else {
              await renderBiomolecule2D(biomoleculeData, img);
            }
            decrementActiveLoads();
            return;
          }
        } catch (e) {
          console.error('Biomolecule lookup failed:', e.message);
          img.alt = `Biomolecule not found: ${lookupValue}`;
          decrementActiveLoads();
          return;
        }
      }

      // Note: Minerals use the same 2D SVG rendering as compounds
      // Server fetches from COD, gets SMILES, and renders 2D
      // 3D viewer is only shown if explicitly requested (is3D flag)

      // ===== QUANTUM ORBITAL: Show Falstad viewer iframe =====
      if (serverType === 'quantum' || serverType === 'hybrid') {
        console.log('%c⚛️ Quantum orbital rendering:', 'background: #9C27B0; color: white; font-weight: bold; padding: 4px;', { serverType, lookupValue });

        try {
          // Build the Falstad viewer URL
          const isHybrid = serverType === 'hybrid';
          const orbitalUrl = buildQuantumOrbitalUrl(lookupValue, isHybrid);
          console.log('%c🔗 Falstad URL:', 'color: #7C4DFF;', orbitalUrl);

          // Create iframe container - uses overflow:hidden to crop out Falstad UI
          const container = document.createElement('div');
          container.className = 'quantum-orbital-container';
          container.style.cssText = `
            display: inline-block;
            position: relative;
            margin: 8px 12px 8px 0;
            vertical-align: middle;
            background: #000;
            border-radius: 8px;
            overflow: hidden;
            box-shadow: 0 2px 8px rgba(0,0,0,0.3);
            width: 350px;
            height: 300px;
          `;

          // Create iframe wrapper - allows positioning the iframe to crop UI
          const iframeWrapper = document.createElement('div');
          iframeWrapper.style.cssText = `
            width: 100%;
            height: 100%;
            overflow: hidden;
            position: relative;
          `;

          // Create iframe for Falstad viewer
          // Make it larger than container and position to show mainly the 3D view
          const iframe = document.createElement('iframe');
          iframe.src = orbitalUrl;
          iframe.style.cssText = `
            width: 800px;
            height: 600px;
            border: none;
            display: block;
            position: absolute;
            top: -60px;
            left: -220px;
            transform: scale(0.8);
            transform-origin: top left;
          `;
          iframe.allow = 'fullscreen';
          iframe.title = displayName || 'Quantum Orbital Viewer';

          iframeWrapper.appendChild(iframe);

          // Create label overlay
          const label = document.createElement('div');
          label.className = 'quantum-orbital-label';
          label.style.cssText = `
            position: absolute;
            bottom: 8px;
            left: 8px;
            background: rgba(0,0,0,0.7);
            color: #fff;
            padding: 4px 10px;
            border-radius: 4px;
            font-size: 12px;
            font-family: system-ui, -apple-system, sans-serif;
            pointer-events: none;
            z-index: 10;
          `;
          label.textContent = displayName || lookupValue;

          // Add "Open in Falstad" button
          const openBtn = document.createElement('a');
          openBtn.href = orbitalUrl;
          openBtn.target = '_blank';
          openBtn.style.cssText = `
            position: absolute;
            top: 8px;
            right: 8px;
            background: rgba(0,0,0,0.7);
            color: #fff;
            padding: 4px 10px;
            border-radius: 4px;
            font-size: 11px;
            font-family: system-ui, -apple-system, sans-serif;
            text-decoration: none;
            z-index: 10;
            opacity: 0;
            transition: opacity 0.2s;
          `;
          openBtn.textContent = '↗ Full Viewer';

          // Show button on hover
          container.addEventListener('mouseenter', () => openBtn.style.opacity = '1');
          container.addEventListener('mouseleave', () => openBtn.style.opacity = '0');

          container.appendChild(iframeWrapper);
          container.appendChild(label);
          container.appendChild(openBtn);

          // Replace the placeholder img with our container
          if (img.parentNode) {
            img.parentNode.replaceChild(container, img);
          }

          console.log('%c✅ Quantum orbital viewer created', 'color: #4CAF50; font-weight: bold;');
          decrementActiveLoads();
          return;
        } catch (e) {
          console.error('Quantum orbital rendering failed:', e.message);
          img.alt = `Orbital not supported: ${lookupValue}`;
          decrementActiveLoads();
          return;
        }
      }

      // ===== 3D REQUEST: Redirect to 3D viewer =====
      if (moleculeData.is3D || moleculeData.show3D || flags.is3D) {
        await show3DViewerInline(moleculeData, img);
        decrementActiveLoads();
        return;
      }

      // ===== 2D RENDERING: Build server URL and fetch SVG =====

      // Get rendering options from settings
      const enableAIFlagControl = settings.enableAIFlagControl === true;

      // Determine rendering flags (from popup settings, optionally overridden by AI flags)
      let showCarbons = settings.m2cfShowCarbons === true;
      let aromaticCircles = settings.m2cfAromaticCircles === true;
      let showMethyl = settings.m2cfShowMethyls === true;
      let showAtomNumbers = settings.m2cfAtomNumbers === true;
      let flipVertical = settings.m2cfFlipVertical === true;
      let flipHorizontal = settings.m2cfFlipHorizontal === true;
      let addHydrogens = settings.m2cfAddH2 === true;
      let showImplicitH = settings.m2cfShowImplicitH !== false;
      let compactDrawing = settings.m2cfCompactDrawing === true;

      // Apply AI flags if enabled
      if (enableAIFlagControl && flags.flagsLocked) {
        if (flags.showCarbons !== undefined) showCarbons = flags.showCarbons;
        if (flags.aromaticCircles !== undefined) aromaticCircles = flags.aromaticCircles;
        if (flags.showMethyls !== undefined) showMethyl = flags.showMethyls;
        if (flags.atomNumbers !== undefined) showAtomNumbers = flags.atomNumbers;
        if (flags.addHydrogens !== undefined) addHydrogens = flags.addHydrogens;
        if (flags.showImplicitH !== undefined) showImplicitH = flags.showImplicitH;
        if (flags.flipHorizontal !== undefined) flipHorizontal = flags.flipHorizontal;
        if (flags.flipVertical !== undefined) flipVertical = flags.flipVertical;
      }

      // NO AUTO-ADAPT - use selected theme directly
      const selectedTheme = settings.sdTheme || 'light';
      const theme = selectedTheme;  // No getEffectiveTheme - manual theme selection only

      // Build server options
      const serverOptions = {
        showCarbons,
        aromaticCircles,
        showMethyls: showMethyl,
        atomNumbers: showAtomNumbers,
        showHydrogens: addHydrogens,
        showImplicitH,
        gradientColors: settings.sdGradientColors === true,
        compactDrawing,
        flipHorizontal,
        flipVertical,
        useStereochemistry: settings.m2cfUse3DSmiles === true,
        theme,
        renderingEngine: settings.renderingEngine || 'chemistrylatex'
      };

      // Get the value to send to server (use lookupValue - the actual IUPAC/mineral/etc. name)
      const valueToSend = serverType === 'smiles'
        ? (moleculeData.smiles || lookupValue)  // SMILES string
        : lookupValue;  // Name for lookup (e.g., "2,4,6-trinitrotoluene" for IUPAC)

      // Build server URL
      const serverUrl = buildServerSvgUrl(serverType, valueToSend, serverOptions);
      console.log('%c📡 Server URL:', 'color: #9C27B0;', serverUrl);

      // Fetch SVG from server
      let svgDataUrl;
      let localFlips = '';

      try {
        const response = await fetch(serverUrl);
        if (!response.ok) {
          throw new Error(`Server returned ${response.status}`);
        }

        // Read local flip header
        localFlips = response.headers.get('X-Local-Flip') || '';
        if (localFlips) {
          console.log('%c🔄 Server says apply these flips locally:', 'color: #FF9800; font-weight: bold;', localFlips);
        }

        // Read CID from server response header (for 3D viewer)
        const cidFromServer = response.headers.get('X-Compound-CID');
        if (cidFromServer) {
          moleculeData.cid = cidFromServer;
          console.log('%c🎯 Got CID from server:', 'color: #4CAF50; font-weight: bold;', cidFromServer);
        }

        // Read SMILES from server response header (for 3D viewer)
        const smilesFromServer = response.headers.get('X-Compound-SMILES');
        if (smilesFromServer) {
          moleculeData.smiles = smilesFromServer;
          console.log('%c🎯 Got SMILES from server:', 'color: #4CAF50; font-weight: bold;', smilesFromServer);
        }

        // Get SVG and apply theme colors
        const rawSvgText = await response.text();
        const coloredSvg = applyThemeColors(rawSvgText, theme);

        svgDataUrl = 'data:image/svg+xml;charset=utf-8,' + encodeURIComponent(coloredSvg);
        // Store raw SVG text for instant theme switching later
        img.dataset.rawSvg = encodeURIComponent(rawSvgText);

      } catch (fetchError) {
        console.error('%c❌ Failed to fetch SVG:', 'color: red;', fetchError);
        img.alt = `Error loading ${displayName}: ${fetchError.message}`;
        decrementActiveLoads();
        return;
      }

      // Create styled img element
      const svgImg = document.createElement('img');
      svgImg.src = svgDataUrl;
      svgImg.alt = displayName;
      svgImg.className = 'molecule-diagram molecule-viewer';  // IMPORTANT: Keep molecule-viewer for re-render selector!
      // CRITICAL: Preserve the original moleculeViewer data for re-rendering
      // Previously this was set to 'true' which lost all the original flags!
      if (img.dataset.moleculeViewer && img.dataset.moleculeViewer !== 'true') {
        svgImg.dataset.moleculeViewer = img.dataset.moleculeViewer;
      } else {
        // Re-encode from moleculeData - ALWAYS set this, never leave undefined
        svgImg.dataset.moleculeViewer = btoa(JSON.stringify(moleculeData || {}));
      }
      svgImg.dataset.nomenclature = displayName;
      svgImg.dataset.originalNomenclature = displayName;
      svgImg.dataset.compoundType = flags.compoundType || serverType;
      // Preserve SMILES if available for re-rendering
      if (moleculeData?.smiles) {
        svgImg.dataset.smiles = moleculeData.smiles;
      } else if (img.dataset.smiles) {
        svgImg.dataset.smiles = img.dataset.smiles;
      }
      // Preserve the server type used (smiles, iupac, mol, mineral, biomol)
      svgImg.dataset.serverType = serverType;
      // CRITICAL: Copy raw SVG for instant theme switching
      svgImg.dataset.rawSvg = img.dataset.rawSvg;

      svgImg.style.cssText = 'display: inline-block; margin: 0; padding: 0; border: none; width: auto; height: auto; max-width: 100%; max-height: 400px; vertical-align: middle;';

      // Apply local flips via CSS transform
      let transforms = [];
      if (localFlips.includes('q')) transforms.push('scaleY(-1)');
      if (localFlips.includes('p')) transforms.push('scaleX(-1)');
      if (flags.rotation) transforms.push(`rotate(${flags.rotation}deg)`);
      if (transforms.length > 0) {
        svgImg.style.transform = transforms.join(' ');
        console.log('%c🔄 Applied CSS transforms:', 'color: #4CAF50; font-weight: bold;', transforms.join(' '));
      }

      // CSS filters are NOT used - theme colors are applied directly in SVG
      // Removed: const cssFilter = getCSSFilterForTheme(selectedTheme, pageIsDark);
      // Removed: if (cssFilter && cssFilter !== 'none') svgImg.style.filter = cssFilter;

      // Wrap and replace
      const wrapper = wrapImageWithRotationContainer(svgImg, flags.rotation, displayName, moleculeData);

      const replaceWithLoaded = () => {
        if (img.parentNode) {
          // FIX: Enforce explicit dimensions on the wrapper to match the image
          // This prevents "inline-block" layout quirks (like extra descender space)
          // from offsetting the absolute positioned controls
          if (svgImg.naturalWidth && svgImg.naturalHeight) {
            wrapper.style.width = svgImg.naturalWidth + 'px';
            wrapper.style.height = svgImg.naturalHeight + 'px';

            // Ensure image fills the wrapper exactly
            svgImg.style.width = '100%';
            svgImg.style.height = '100%';
            svgImg.style.display = 'block'; // Removes baseline gap space
          }

          img.parentNode.replaceChild(wrapper, img);
          addHoverControls(wrapper, displayName, moleculeData);
          console.log('%c✅ Server SVG loaded successfully', 'color: green; font-weight: bold;');

          // Apply saved average size scaling to new molecule
          if (settings.sdAverageSize && settings.sdAverageSize !== 100) {
            setTimeout(() => applyAverageSizeScaling(settings.sdAverageSize), 100);
          }
        }
        decrementActiveLoads();
      };


      if (svgImg.complete && svgImg.naturalWidth > 0) {
        replaceWithLoaded();
      } else {
        svgImg.onload = replaceWithLoaded;
        svgImg.onerror = () => {
          console.error('%c❌ Failed to load SVG', 'color: red;', serverUrl);
          img.alt = `Error loading ${displayName}`;
          decrementActiveLoads();
        };
      }

    } catch (error) {
      console.error('%c❌ Rendering error:', 'color: red; font-weight: bold;', error);
      img.alt = `Error: ${error.message}`;
      decrementActiveLoads();
    }
  }


  // Helper function to show 3D viewer inline (uses embed.molview.org or 3Dmol.js)
  async function show3DViewerInline(moleculeData, targetElement) {
    // console.log('%c🔮 SHOWING 3D VIEWER INLINE', 'background: #764ba2; color: white; font-size: 14px; padding: 8px;');
    // console.log('moleculeData:', moleculeData);
    // console.log('targetElement:', targetElement);

    // Extract compound name from moleculeData, prioritizing displayName from flags
    const compoundName = moleculeData.flags?.displayName || moleculeData.nomenclature || moleculeData.smiles || '';

    // Find the container - could be the element itself or a wrapper
    let container = targetElement;

    // Check if this is already a 3D viewer container (toggle back to 2D)
    if (container.classList && container.classList.contains('molecule-3d-viewer')) {
      console.log('🔄 Toggling back to 2D view');
      const original2D = container._original2DElement;
      if (original2D && container.parentNode) {
        container.parentNode.replaceChild(original2D, container);
        console.log('✅ Restored original 2D view');
      }
      return;
    }

    // Find the img element - could be the element itself or inside a container
    // Also handle placeholder divs from client-side mode
    let img = targetElement;
    if (targetElement.tagName !== 'IMG') {
      img = targetElement.querySelector('img');
      if (!img) {
        // Might be a size control wrapper
        const sizeWrapper = targetElement.querySelector('.molecule-size-wrapper');
        if (sizeWrapper) {
          img = sizeWrapper.querySelector('img');
          container = sizeWrapper;
        }
      }
      // If still no img found, check if it's a valid placeholder/container we can use
      if (!img) {
        // For client-side mode, we may not have an img - the targetElement itself is the container
        if (targetElement.classList && (
          targetElement.classList.contains('chem-3d-placeholder') ||
          targetElement.classList.contains('chem-2d-viewer-container') ||
          targetElement.classList.contains('molecule-viewer-container')
        )) {
          console.log('%c🔄 Using container as target (no img element)', 'color: #ff9800;');
          img = null; // We'll work with the container directly
          container = targetElement;
        }
      }
    }

    // If we still don't have either img or a valid container, error out
    if (!img && !container) {
      console.error('❌ Could not find img element or valid container in:', targetElement);
      throw new Error('Could not find image element');
    }

    // console.log('Found img element:', img);
    // console.log('Container element:', container);
    // console.log('Compound name:', compoundName);

    try {
      // Determine molecule type
      const moleculeType = moleculeData.compoundType || moleculeData.flags?.compoundType || 'compound';
      const isBiomolecule = moleculeType === 'biomolecule';

      let dimensions;

      if (isBiomolecule) {
        // For biomolecules, use fixed predetermined sizes (Mol* needs larger canvas)
        const viewerSize = settings.viewer3DSize || 'normal';
        const sizeDimensions = {
          'small': { width: 200, height: 150 },
          'medium': { width: 300, height: 250 },
          'normal': { width: 400, height: 350 },
          'large': { width: 600, height: 450 },
          'xlarge': { width: 800, height: 600 }
        };
        dimensions = sizeDimensions[viewerSize] || sizeDimensions['normal'];
        console.log('%c🧬 Biomolecule 3D Viewer Size (fixed):', 'color: #E91E63; font-weight: bold;', viewerSize, dimensions);
      } else {
        // For compounds and minerals, use the 2D molecule's aspect ratio
        // Get the current 2D image dimensions
        const imgWidth = img ? (img.offsetWidth || img.naturalWidth || 300) : 300;
        const imgHeight = img ? (img.offsetHeight || img.naturalHeight || 200) : 200;
        const aspectRatio = imgWidth / imgHeight;

        // Get scale factor from settings (percentage, default 100%)
        const scaleFactor = (moleculeType === 'mineral'
          ? (settings.mineral3DSize || 100)
          : (settings.compound3DSize || 100)) / 100;

        // Base width (then calculate height from aspect ratio)
        const baseWidth = 300;  // Default base width
        const scaledWidth = Math.round(baseWidth * scaleFactor);
        const scaledHeight = Math.round(scaledWidth / aspectRatio);

        dimensions = { width: scaledWidth, height: scaledHeight };
        console.log('%c🧪 Compound/Mineral 3D Viewer Size (aspect ratio preserved):', 'color: #4CAF50; font-weight: bold;',
          { type: moleculeType, scaleFactor, aspectRatio: aspectRatio.toFixed(2), dimensions });
      }

      // Create container for 3D viewer with appropriate sizing
      const viewer3DContainer = document.createElement('div');
      viewer3DContainer.className = 'molecule-viewer-container molecule-3d-viewer';
      viewer3DContainer.style.cssText = `
        display: inline-block;
        width: ${dimensions.width}px;
        height: ${dimensions.height}px;
        margin: 0;
        padding: 0;
        vertical-align: middle;
        position: relative;
        border: none;
        outline: none;
        border-radius: 8px;
        overflow: hidden;
        background: transparent;
        box-shadow: none;
        max-width: 100%;
        box-sizing: border-box;
      `;

      // Create iframe for 3D viewer based on user's selected source
      const viewer3DIframe = document.createElement('iframe');

      let viewerUrl;

      // Check if this is a protein/mineral with a direct embed URL
      if (moleculeData.embedUrl || moleculeData.pdbid || moleculeData.codid) {
        // Dynamically rebuild the embed URL with current settings (for instant settings updates)
        const moleculeType = moleculeData.compoundType || 'compound';

        if (moleculeData.pdbid) {
          const pdbId = moleculeData.pdbid.toLowerCase();

          // Check for user toggle between Mol* and MolView (from bio toggle button)
          // Default to MolView UNLESS:
          // 1. User explicitly toggled to Mol* (preferMolstar === true)
          // 2. Bio assembly mode is enabled (settings.molviewBioAssembly)
          const useMolstar = moleculeData.preferMolstar === true || settings.molviewBioAssembly === true;

          if (useMolstar) {
            // Use Mol* viewer (better for proteins, handles more structures)
            const graphics = settings.bioAssemblyGraphics || 'balanced';
            const provider = settings.bioAssemblyPdbProvider || 'rcsb';

            if (settings.molviewBioAssembly) {
              // Bio assembly mode
              // Since we removed the dropdown, default to molstar standard (it handles assemblies well)
              const viewer = 'molstar';
              switch (viewer) {
                case 'molstar-me':
                  viewerUrl = `https://molstar.org/me/viewer/?pdb=${pdbId}&hide-controls=1&graphics-mode=${graphics}&pdb-provider=${provider}`;
                  console.log('%c🧬 Using Mol* Mesoscale Explorer:', 'color: #FF5722; font-weight: bold;', viewerUrl);
                  break;
                case 'molstar':
                default:
                  viewerUrl = `https://molstar.org/viewer/?pdb=${pdbId}&assembly-id=1&hide-controls=1`;
                  console.log('%c🧬 Using Mol* Standard Viewer:', 'color: #E91E63; font-weight: bold;', viewerUrl);
                  break;
              }
            } else {
              // Asymmetric unit mode - use Mol* (selected via toggle)
              viewerUrl = `https://molstar.org/viewer/?pdb=${pdbId}&hide-controls=1`;
              console.log('%c🧬 Using Mol* (user selected):', 'color: #E91E63; font-weight: bold;', viewerUrl);
            }
          } else {
            // DEFAULT: Use MolView
            viewerUrl = buildMolViewEmbedUrl(moleculeData.pdbid, 'pdbid', 'biomolecule');
            console.log('%c🧬 Using MolView (default):', 'color: #9C27B0; font-weight: bold;', viewerUrl);
          }
        } else if (moleculeData.codid) {
          // Mineral - use buildMolViewEmbedUrl with current settings
          viewerUrl = buildMolViewEmbedUrl(moleculeData.codid, 'codid', 'mineral');
          console.log('%c💎 Built MolView URL for mineral with settings:', 'color: #00BCD4; font-weight: bold;', viewerUrl);
        } else {
          // Fallback to stored embedUrl if no identifier available
          console.log('%c🔮 Using stored embed URL:', 'color: #9C27B0; font-weight: bold;', moleculeData.embedUrl);
          viewerUrl = moleculeData.embedUrl;
        }

        // Ensure iframe has necessary permissions for WebGL and interaction
        viewer3DIframe.setAttribute('allow', 'accelerometer; autoplay; clipboard-write; encrypted-media; gyroscope; picture-in-picture; web-share; fullscreen; xr-spatial-tracking');
      } else {
        // For compounds - use MolView embed
        console.log('%c🧪 COMPOUND 3D VIEWER PATH', 'background: #4CAF50; color: white; font-weight: bold; padding: 4px;');
        console.log('moleculeData:', moleculeData);
        console.log('compoundName:', compoundName);

        // Priority: CID > SMILES > Name query (CID is most reliable for 3D)

        if (moleculeData.cid) {
          // CID is the best option - direct lookup
          viewerUrl = buildMolViewEmbedUrl(moleculeData.cid, 'cid', 'compound');
          console.log('%c📍 MolView URL (CID):', 'color: #4CAF50; font-weight: bold;', viewerUrl);
        } else if (moleculeData.smiles) {
          // SMILES is second best
          viewerUrl = buildMolViewEmbedUrl(moleculeData.smiles, 'smiles', 'compound');
          console.log('%c📍 MolView URL (SMILES):', 'color: #0066cc; font-weight: bold;', viewerUrl);
        } else if (compoundName) {
          // No CID or SMILES - try to fetch CID from our server
          try {
            const result = await fetchFromChemistryLaTeXServer('mol', compoundName);
            if (result && result.cid) {
              viewerUrl = buildMolViewEmbedUrl(result.cid, 'cid', 'compound');
              console.log('%c📍 MolView URL (fetched CID):', 'color: #4CAF50; font-weight: bold;', viewerUrl);
            } else if (result && result.smiles) {
              viewerUrl = buildMolViewEmbedUrl(result.smiles, 'smiles', 'compound');
              console.log('%c📍 MolView URL (fetched SMILES):', 'color: #0066cc; font-weight: bold;', viewerUrl);
            } else {
              throw new Error('No CID or SMILES found');
            }
          } catch (e) {
            console.error('%c❌ Could not fetch CID for 3D view:', 'color: red;', e.message);
            viewerUrl = 'about:blank';
          }
        } else {
          console.error('%c❌ No identifier for 3D view', 'color: red;');
          viewerUrl = 'about:blank';
        }



        viewer3DIframe.setAttribute('allow', 'accelerometer; autoplay; clipboard-write; encrypted-media; gyroscope; picture-in-picture; web-share; fullscreen; xr-spatial-tracking');
        console.log('%c🎯 Final viewerUrl:', 'background: blue; color: white; font-weight: bold;', viewerUrl);
      }





      // SAFETY CHECK: Ensure viewerUrl is never undefined/null
      // If it is, the iframe would load the current page which is wrong
      if (!viewerUrl) {
        console.error('%c❌ viewerUrl is undefined - preventing iframe from loading current page', 'color: red; font-weight: bold;');
        viewerUrl = 'about:blank';
      }

      // viewer3DIframe.src = viewerUrl; // Moved to after appendChild
      // Allow necessary features for 3D rendering
      viewer3DIframe.setAttribute('allow', 'accelerometer; autoplay; clipboard-write; encrypted-media; gyroscope; picture-in-picture; web-share; fullscreen; xr-spatial-tracking');
      viewer3DIframe.setAttribute('allowfullscreen', 'true');
      viewer3DIframe.setAttribute('frameborder', '0');
      viewer3DIframe.style.cssText = `
        width: 100%;
        height: 100%;
        border: none !important;
        outline: none !important;
        margin: 0;
        padding: 0;
        display: block;
        background: transparent;
        box-shadow: none;
      `;
      viewer3DIframe.title = `3D Viewer: ${compoundName}`;

      // Create toggle button to switch back to 2D (removed - now handled by hover controls)
      // Store reference to original 2D element for toggling
      viewer3DContainer._original2DElement = container;

      // =================================================================================
      // 🚀 DIRECT LOAD IMPLEMENTATION
      // =================================================================================

      // Append the iframe
      viewer3DContainer.appendChild(viewer3DIframe);

      // Add Size Controls (Bottom Left)
      const controlsDiv = document.createElement('div');
      controlsDiv.className = 'chem-size-controls';
      controlsDiv.style.cssText = `
        position: absolute;
        bottom: 4px;
        left: 4px;
        display: flex;
        flex-direction: column;
        gap: 2px;
        z-index: 1000;
        opacity: 0;
        transition: opacity 0.2s;
      `;

      // FIX: Ensure container has positioning context so controls stay in corner
      if (viewer3DContainer.style.position !== 'relative' && viewer3DContainer.style.position !== 'absolute') {
        viewer3DContainer.style.position = 'relative';
      }

      // Show controls on hover
      viewer3DContainer.onmouseenter = () => controlsDiv.style.opacity = '1';
      viewer3DContainer.onmouseleave = () => controlsDiv.style.opacity = '0';

      const createBtn = (text, delta) => {
        const btn = document.createElement('button');
        btn.innerHTML = text;
        btn.style.cssText = `
          width: 24px;
          height: 24px;
          border: none;
          background: rgba(0, 0, 0, 0.7);
          color: white;
          border-radius: 4px;
          cursor: pointer;
          font-size: 10px;
          display: flex;
          align-items: center;
          justify-content: center;
          transition: background 0.2s;
        `;
        btn.onmouseenter = () => btn.style.background = 'rgba(0, 0, 0, 0.9)';
        btn.onmouseleave = () => btn.style.background = 'rgba(0, 0, 0, 0.7)';
        btn.onclick = (e) => {
          e.stopPropagation();
          const currentWidth = parseInt(viewer3DContainer.style.width);
          const currentHeight = parseInt(viewer3DContainer.style.height);
          // Limit min size
          if (delta < 0 && currentWidth < 200) return;

          const scale = 1 + (delta / 100); // 10% change
          viewer3DContainer.style.width = `${Math.round(currentWidth * scale)}px`;
          viewer3DContainer.style.height = `${Math.round(currentHeight * scale)}px`;
        };
        return btn;
      };


      // FIX: Ensure container has positioning context so controls stay in corner
      if (viewer3DContainer.style.position !== 'relative' && viewer3DContainer.style.position !== 'absolute') {
        viewer3DContainer.style.position = 'relative';
      }

      // FORCE REFLOW for 3D viewer container as well
      void viewer3DContainer.offsetWidth;

      controlsDiv.appendChild(createBtn('▲', 10));
      controlsDiv.appendChild(createBtn('▼', -10));
      viewer3DContainer.appendChild(controlsDiv);

      // Set src to trigger load IMMEDIATELY
      viewer3DIframe.src = viewerUrl;
      console.log('%c🚀 3D Viewer Loaded directly', 'color: #00FF00; font-weight: bold');

      // Add fallback for MolView failures on biomolecules
      // MolView sometimes can't load large structures like adenovirus
      if (moleculeData.pdbid && viewerUrl.includes('embed.molview.org')) {
        const pdbId = moleculeData.pdbid.toUpperCase();
        let fallbackTriggered = false;

        // Build fallback URL - use direct Mol* (same as bio assembly mode, but without assembly-id)
        // No assembly-id = asymmetric unit
        const fallbackUrl = `https://molstar.org/viewer/?pdb=${pdbId.toLowerCase()}&hide-controls=1`;

        // Function to replace iframe with new Mol* iframe (avoids CSP issues with src change)
        const switchToMolstar = () => {
          if (fallbackTriggered) return;
          fallbackTriggered = true;

          console.log('%c⚠️ MolView failed - switching to Mol* viewer', 'color: orange; font-weight: bold;');

          // Create a NEW iframe instead of changing src (avoids CSP block)
          const newIframe = document.createElement('iframe');
          newIframe.src = fallbackUrl;
          newIframe.style.cssText = viewer3DIframe.style.cssText;
          newIframe.setAttribute('allow', viewer3DIframe.getAttribute('allow') || '');
          newIframe.setAttribute('allowfullscreen', 'true');
          newIframe.setAttribute('frameborder', '0');
          newIframe.title = `Mol* Viewer: ${pdbId}`;

          // Replace old iframe with new one
          if (viewer3DIframe.parentNode) {
            viewer3DIframe.parentNode.replaceChild(newIframe, viewer3DIframe);
            console.log('%c🔄 Switched to Mol* viewer:', 'color: #4CAF50; font-weight: bold;', fallbackUrl);
          }
        };

        // Listen for MolView error messages (postMessage)
        // These are rare but we should still handle them
        const errorHandler = (event) => {
          // Check if message is from MolView
          if (event.origin === 'https://embed.molview.org' && !fallbackTriggered) {
            const data = event.data;

            // Check for explicit error indicators
            const isError = data && (
              data.error ||
              data.type === 'error' ||
              (typeof data === 'string' && data.toLowerCase().includes('error')) ||
              (typeof data === 'string' && data.toLowerCase().includes('failed')) ||
              (typeof data === 'string' && data.toLowerCase().includes('404')) ||
              (typeof data === 'string' && data.toLowerCase().includes('not found'))
            );

            if (isError) {
              console.log('%c❌ MolView error message detected - switching to Mol*', 'color: red; font-weight: bold;', data);
              switchToMolstar();
              window.removeEventListener('message', errorHandler);
            }
          }
        };
        window.addEventListener('message', errorHandler);

        // Clean up after 60 seconds (give MolView plenty of time to load)
        setTimeout(() => {
          window.removeEventListener('message', errorHandler);
        }, 60000);
      }

      // Replace original element with 3D viewer
      if (container.parentNode) {
        container.parentNode.replaceChild(viewer3DContainer, container);

        // Re-add hover controls to the 3D viewer container
        addHoverControls(viewer3DContainer, compoundName, moleculeData);

        console.log('%c✅ 3D viewer embedded inline', 'color: green; font-weight: bold;');
      } else {
        console.error('%c❌ Cannot replace image - no parent node', 'color: red; font-weight: bold;');
        // If no parent, insert after the img element as fallback
        if (img.nextSibling) {
          img.parentElement.insertBefore(viewer3DContainer, img.nextSibling);
        } else {
          img.parentElement.appendChild(viewer3DContainer);
        }
        img.style.display = 'none';
      }
      decrementActiveLoads();

    } catch (error) {
      console.error('%c❌ Error showing 3D viewer inline:', 'color: red; font-weight: bold;', error);
      // Fallback to server-rendered 2D image (no direct database access)
      const fallbackUrl = buildServerSvgUrl('mol', compoundName, {
        showCarbons: settings.m2cfShowCarbons,
        aromaticCircles: settings.m2cfAromaticCircles,
        renderingEngine: settings.renderingEngine
      });

      img.src = fallbackUrl;
      img.classList.add('molecule-fadein');
      img.classList.remove('molecule-loading');
      decrementActiveLoads();
    }
  }

  // Load molecule image
  function loadMoleculeImage(img) {

    console.log('%c🎨 [loadMoleculeImage] Rendering molecule', 'background: #9C27B0; color: white; font-weight: bold; padding: 4px;');

    try {
      const moleculeData = JSON.parse(atob(img.dataset.moleculeViewer || img.dataset.mol2chemfig));
      console.log('%c📦 Decoded molecule data:', 'color: #9C27B0;', moleculeData);
      renderClientSide(moleculeData, img);
    } catch (e) {
      console.error('%c❌ [loadMoleculeImage] Failed to parse molecule data:', 'color: red;', e.message);
      // Show error placeholder
      img.alt = '⚠️ Failed to render molecule';
      img.style.minWidth = '100px';
      img.style.minHeight = '50px';
      img.style.backgroundColor = '#ffebee';
      img.style.border = '1px solid #ef5350';
      img.style.borderRadius = '4px';
      img.style.padding = '10px';
      img.dataset.loaded = 'error';
    }
  }

  // Queue for molecules waiting to be loaded (when activeLoads is at max)
  const moleculeLoadQueue = [];
  let isProcessingQueue = false;

  function processLoadQueue() {
    if (isProcessingQueue || moleculeLoadQueue.length === 0) return;

    isProcessingQueue = true;

    // Process one item from queue if we have capacity
    while (moleculeLoadQueue.length > 0 && activeLoads < maxConcurrentLoads) {
      const img = moleculeLoadQueue.shift();
      // Double-check it's not already loaded (could have been loaded by another path)
      if (img.dataset.loaded !== 'true' && img.dataset.loaded !== 'loading') {
        img.dataset.loaded = 'loading';
        loadMoleculeImage(img);
      }
    }

    isProcessingQueue = false;

    // If still items in queue, schedule next check
    if (moleculeLoadQueue.length > 0) {
      setTimeout(processLoadQueue, 200);  // Check again after 200ms
    }
  }

  // Override observer callback to handle both types
  const originalObserver = observer;
  const newObserver = new IntersectionObserver((entries) => {
    entries.forEach((entry) => {
      const img = entry.target;

      // Only process if intersecting AND not already loaded/loading
      if (entry.isIntersecting && img.dataset.loaded !== 'true' && img.dataset.loaded !== 'loading') {
        // Check if this is a chemistry image (MoleculeViewer, compound, quantum, or Mol2chemfig)
        if (img.classList.contains('molecule-viewer') || img.classList.contains('molecule-compound') || img.classList.contains('molecule-legacy') || img.classList.contains('molecule-quantum')) {
          console.log(`%c🧪 Molecule image entered viewport: ${img.alt || 'unknown'}`, 'background: #0088FF; color: #FFF; padding: 4px;');

          // UNOBSERVE immediately - we've seen it, will load it
          newObserver.unobserve(img);

          if (activeLoads < maxConcurrentLoads) {
            // Load immediately
            img.dataset.loaded = 'loading';
            loadMoleculeImage(img);
          } else {
            // Queue for later - don't mark as loading yet so we can track pending
            console.log(`%c⏳ Queuing molecule (${activeLoads}/${maxConcurrentLoads} active): ${img.alt || 'unknown'}`, 'background: #FFA500; color: white; padding: 2px;');
            moleculeLoadQueue.push(img);
            // Trigger queue processing
            setTimeout(processLoadQueue, 100);
          }
        } else {
          // Standard image loading
          if (activeLoads < maxConcurrentLoads) {
            loadImage(img);
          } else {
            if (typeof requestIdleCallback !== 'undefined') {
              requestIdleCallback(() => {
                if (activeLoads < maxConcurrentLoads) {
                  loadImage(img);
                }
              });
            } else {
              setTimeout(() => {
                if (activeLoads < maxConcurrentLoads) {
                  loadImage(img);
                }
              }, 50);
            }
          }
        }
      }
    });
  }, {
    rootMargin: '300px' // Start loading 300px before image enters viewport
  });

  // Replace the old observer
  window._lazyLoadObserver = newObserver;

  // Expose loading function globally for immediate loading when performance mode is off
  window._loadMoleculeImage = loadMoleculeImage;

  // Expose queue function for chunk-based loading
  window._queueMoleculeLoad = function (img) {
    if (img.dataset.loaded === 'true' || img.dataset.loaded === 'loading') return;

    if (activeLoads < maxConcurrentLoads) {
      img.dataset.loaded = 'loading';
      loadMoleculeImage(img);
    } else {
      console.log(`%c⏳ Queuing molecule: ${img.alt || 'unknown'}`, 'background: #FFA500; color: white; padding: 2px;');
      moleculeLoadQueue.push(img);
      setTimeout(processLoadQueue, 100);
    }
  };

  // Export show3DViewerInline for hover controls
  window.show3DViewerInline = show3DViewerInline;

  // ============================================
  // CHUNK-BASED SCROLL PRELOADER
  // ============================================
  // Pre-loads molecules in chunks as you scroll (10 molecules at a time)
  // This helps load molecules above/below viewport without having to see every single one

  const CHUNK_SIZE = 10;  // Load 10 molecules at a time
  const PRELOAD_DISTANCE = 1500;  // Preload molecules within 1500px of viewport
  let lastScrollChunkTime = 0;
  const SCROLL_CHUNK_INTERVAL = 500;  // Check every 500ms during scroll

  function preloadNearbyMolecules() {
    const now = Date.now();
    if (now - lastScrollChunkTime < SCROLL_CHUNK_INTERVAL) return;
    lastScrollChunkTime = now;

    // Find all unloaded molecules
    const unloadedMolecules = document.querySelectorAll('img.molecule-diagram[data-loaded="false"], img.molecule-diagram:not([data-loaded])');
    if (unloadedMolecules.length === 0) return;

    const viewportTop = window.scrollY - PRELOAD_DISTANCE;
    const viewportBottom = window.scrollY + window.innerHeight + PRELOAD_DISTANCE;

    // Find molecules within preload distance
    const nearbyMolecules = [];
    unloadedMolecules.forEach(img => {
      const rect = img.getBoundingClientRect();
      const absoluteTop = window.scrollY + rect.top;

      if (absoluteTop >= viewportTop && absoluteTop <= viewportBottom) {
        nearbyMolecules.push({
          img,
          distance: Math.abs(absoluteTop - (window.scrollY + window.innerHeight / 2))
        });
      }
    });

    if (nearbyMolecules.length === 0) return;

    // Sort by distance and take a chunk
    nearbyMolecules.sort((a, b) => a.distance - b.distance);
    const chunk = nearbyMolecules.slice(0, CHUNK_SIZE);

    console.log(`%c📦 Preloading chunk of ${chunk.length} molecules (${nearbyMolecules.length} nearby, ${unloadedMolecules.length} total unloaded)`, 'background: #9C27B0; color: white; padding: 2px;');

    // Queue the chunk for loading
    chunk.forEach(({ img }) => {
      if (window._queueMoleculeLoad) {
        window._queueMoleculeLoad(img);
      }
    });
  }

  // Set up scroll-based chunk preloading
  let scrollPreloadTimeout = null;
  window.addEventListener('scroll', () => {
    if (scrollPreloadTimeout) clearTimeout(scrollPreloadTimeout);
    scrollPreloadTimeout = setTimeout(preloadNearbyMolecules, 150);
  }, { passive: true });

  // Also preload on initial load after a short delay
  setTimeout(preloadNearbyMolecules, 1000);

  log.success('✅ Lazy-loading observer initialized (max 3 concurrent loads)');
}


// Debounce scanAndRender to prevent excessive processing
let scanAndRenderTimeout = null;
let lastScanTime = 0;
const MIN_SCAN_INTERVAL = 1000; // Minimum 1 second between scans

function scanAndRender() {
  // Clear any pending scan
  if (scanAndRenderTimeout) {
    clearTimeout(scanAndRenderTimeout);
  }

  // Check if we scanned recently
  const now = Date.now();
  const timeSinceLastScan = now - lastScanTime;

  if (timeSinceLastScan < MIN_SCAN_INTERVAL) {
    // Schedule for later
    const delay = MIN_SCAN_INTERVAL - timeSinceLastScan;
    scanAndRenderTimeout = setTimeout(() => {
      scanAndRenderImmediate();
    }, delay);
    console.log(`[ChemRenderer] Debouncing scan, will run in ${delay}ms`);
    return;
  }

  // Run immediately
  scanAndRenderImmediate();
}

function scanAndRenderImmediate() {
  lastScanTime = Date.now();

  if (!settings.enabled) {
    log.debug('scanAndRender() called but extension is disabled');
    return;
  }

  log.inject('🔍 STARTING PAGE SCAN for chemical formulas');

  // Get all text nodes (like ReplaceR does)
  const elements = document.getElementsByTagName('*');
  let replacements = 0;
  let elementsScanned = 0;
  let textNodesFound = 0;
  let chemistryTextFound = [];

  // SAFETY: Limit total elements scanned per pass to prevent crashes on huge pages
  const MAX_ELEMENTS_PER_SCAN = 10000;

  log.debug(`Total elements in DOM: ${elements.length}`);

  for (let i = 0; i < elements.length && elementsScanned < MAX_ELEMENTS_PER_SCAN; i++) {
    const element = elements[i];
    elementsScanned++;

    // Skip scripts, styles, and already processed elements
    if (element.tagName === 'SCRIPT' ||
      element.tagName === 'STYLE' ||
      element.tagName === 'NOSCRIPT' ||
      element.classList.contains('MathJax') ||
      element.hasAttribute('data-chem-processed')) {
      continue;
    }

    // CRITICAL FIX FOR CHATGPT: ChatGPT splits text like "chem:mineral=calcite:" into separate nodes
    // But we can't use element.innerHTML because it breaks React
    // Solution: Combine all text node content, process it, then replace only the first text node

    // IMPROVED: Use innerText to get full text including spaces, then find all text nodes
    // This handles cases like "chem:mol=benzyl chloride:" where space might split nodes

    const textNodes = [];
    let combinedText = '';

    // Collect ALL text nodes recursively (not just direct children)
    // This handles deeply nested text that might be split
    const walker = document.createTreeWalker(
      element,
      NodeFilter.SHOW_TEXT,
      {
        acceptNode: (node) => {
          // Only accept text nodes that are direct children (to avoid processing nested elements twice)
          if (node.parentNode === element) {
            return NodeFilter.FILTER_ACCEPT;
          }
          return NodeFilter.FILTER_SKIP;
        }
      }
    );

    let textNode;
    while (textNode = walker.nextNode()) {
      textNodes.push(textNode);
      combinedText += textNode.nodeValue;
    }

    if (textNodes.length === 0 || !combinedText) {
      continue;
    }

    // Check if combined text contains chemistry notation
    if (combinedText.includes('chem:') || combinedText.includes('\\ce{') || combinedText.includes('ce{')) {
      textNodesFound++;
      chemistryTextFound.push(combinedText.substring(0, 100));
      log.debug(`🔬 Found potential chemistry text: "${combinedText.substring(0, 80)}..."`);

      const replacedText = wrapChemicalFormulas(combinedText);

      if (replacedText !== combinedText) {
        log.debug(`✨ Found formula in text: "${combinedText.substring(0, 50)}..."`);
        log.debug(`Wrapped result: "${replacedText.substring(0, 50)}..."`);

        // SAFE: Replace only the first text node with a span containing the result
        // Remove all other text nodes to avoid duplication
        const span = document.createElement('span');
        span.innerHTML = replacedText;
        span.setAttribute('data-chem-processed', 'true');

        // Replace first text node with span
        element.replaceChild(span, textNodes[0]);

        // Remove remaining text nodes (they're now part of the combined/replaced text)
        for (let k = 1; k < textNodes.length; k++) {
          if (textNodes[k].parentNode) {
            textNodes[k].parentNode.removeChild(textNodes[k]);
          }
        }

        // Setup loading for newly created images
        const images = span.querySelectorAll('img.molecule-diagram[data-loaded="false"]');
        console.log('%c[ChemRenderer] Found images to load:', 'color: #FF00FF; font-weight: bold;', images.length);
        if (images.length > 0) {
          console.log('%c[ChemRenderer] Performance mode:', 'color: #FF00FF;', settings.performanceMode);
          console.log('%c[ChemRenderer] _lazyLoadObserver exists:', 'color: #FF00FF;', !!window._lazyLoadObserver);
          console.log('%c[ChemRenderer] _loadMoleculeImage exists:', 'color: #FF00FF;', !!window._loadMoleculeImage);

          if (settings.performanceMode && window._lazyLoadObserver) {
            // Performance mode ON: lazy-load as images scroll into view
            images.forEach(img => {
              console.log('%c[ChemRenderer] Observing image for lazy load:', 'color: #FF00FF;', img.className);
              window._lazyLoadObserver.observe(img);
            });
            log.debug(`📊 Setup lazy-loading for ${images.length} SVGs`);
          } else if (window._loadMoleculeImage && window._queueMoleculeLoad) {
            // Performance mode OFF: Use chunk-based loading to prevent memory spikes
            // Load nearby molecules first, queue distant ones
            console.log(`%c[ChemRenderer] 🚀 Chunk-loading ${images.length} images (lazy loading disabled)`, 'color: #00FF00; font-weight: bold;');

            // Sort images by distance from current viewport
            const viewportCenter = window.scrollY + (window.innerHeight / 2);
            const imagesWithDistance = Array.from(images).map(img => {
              const rect = img.getBoundingClientRect();
              const imgCenter = window.scrollY + rect.top + (rect.height / 2);
              const distance = Math.abs(imgCenter - viewportCenter);
              return { img, distance };
            });

            // Sort by distance (closest first)
            imagesWithDistance.sort((a, b) => a.distance - b.distance);

            // Queue all images - the queue system handles throttling
            imagesWithDistance.forEach(({ img }, index) => {
              // Use the queue system which respects maxConcurrentLoads
              window._queueMoleculeLoad(img);
            });

            log.debug(`📊 Queued ${images.length} SVGs for chunk-based loading`);
          } else if (window._loadMoleculeImage) {
            // Fallback: stagger loads but respect the queue
            console.log(`%c[ChemRenderer] 🚀 Fallback loading ${images.length} images`, 'color: #FFAA00; font-weight: bold;');

            images.forEach((img, index) => {
              img.dataset.loaded = 'pending';
              setTimeout(() => {
                if (img.dataset.loaded === 'pending') {
                  img.dataset.loaded = 'loading';
                  window._loadMoleculeImage(img);
                }
              }, index * 300); // 300ms stagger as fallback
            });

            log.debug(`📊 Triggered fallback loading for ${images.length} SVGs`);
          } else {
            console.error('%c[ChemRenderer] ERROR: No loader available!', 'color: #FF0000; font-weight: bold;');
          }
        }

        replacements++;
      }
    }
  }

  log.inject(`📊 Scan results:`);
  log.inject(`  - Elements scanned: ${elementsScanned}`);
  log.inject(`  - Text nodes found: ${textNodesFound}`);
  log.inject(`  - Formulas wrapped: ${replacements}`);
  if (chemistryTextFound.length > 0) {
    log.inject(`  - Chemistry text samples: ${chemistryTextFound.length} found`);
    chemistryTextFound.forEach((t, i) => {
      log.debug(`    Sample ${i + 1}: ${t}`);
    });
  }

  if (replacements > 0) {
    log.inject(`🎨 ${replacements} formula(s) found and processed!`);

    // If MathJax is available, ask it to render any wrapped formulas
    if (window.MathJax && window.MathJax.typesetPromise) {
      log.inject('Calling MathJax.typesetPromise() for final rendering');
      MathJax.typesetPromise()
        .then(() => {
          log.success('✨ RENDERING COMPLETE! Formulas rendered with MathJax');
        })
        .catch((err) => {
          log.debug('MathJax rendering skipped (not available), using fallback');
        });
    } else {
      log.success('✨ Formulas processed - CodeCogs and Unicode rendering active');
    }
  } else {
    log.debug('No formulas found in this scan - page may not contain chemistry notation');
  }

  // Watch for dynamic content (like React/Vue apps)
  if (!window._chemRendererObserver) {
    log.inject('Setting up dynamic content observer');
    observePageChanges();
  }
}


/**
 * Wrap chemical formulas in renderable format
 * Converts plain text chemistry to Unicode with subscripts/superscripts
 */
function wrapChemicalFormulas(text) {
  if (!text || text.trim().length === 0) return text;

  // SAFETY: Limit text length to prevent excessive processing
  if (text.length > 50000) {
    console.warn('[ChemRenderer] Text too long for processing, skipping');
    return text;
  }

  let result = text;
  let patternMatches = [];

  // Pattern 0b: chem:text: (double colon format - captures everything in between)
  // Example: chem:CCO:, chem:benzene:, chem:1-chloro-benzene:, chem:CC(=O)C:, chem:histamine/+c+o+m+n:, chem:phenol/d-c:
  log.debug('🧪 Applying Pattern 0b: chem:text: → Using selected renderer engine');

  // DEBUG: Log a sample of the text being scanned
  if (result.length > 0) {
    const sample = result.substring(0, 200);
    console.log('[ChemRenderer] 🔍 DEBUG: Scanning text sample:', sample);
    console.log('[ChemRenderer] 🔍 DEBUG: Text length:', result.length);
    console.log('[ChemRenderer] 🔍 DEBUG: Looking for pattern: /\\bchem:([^:]+):/g');
  }

  // CRITICAL: Match BOTH standard and named syntax with explicit type=value format
  // Standard format: chem:smiles=CCO:, chem:mol=benzene:, chem:biomol=insulin:, chem:mineral=quartz:
  // Named format: chem:ethanolsmiles=CCO:, chem:Aspirinmol=aspirin:, chem:2,5-hexanedionesmiles=CC1=CC(=O)C=CC1=O:
  // The named format allows a display name before the type keyword
  // Invalid (no longer accepted): chem:benzene:, chem:CCO:

  // Build regex that matches BOTH formats:
  // 1. chem:type=value: (standard)
  // 2. chem:NAMEtype=value: (named - name can include numbers, letters, commas, hyphens, spaces, parentheses)
  const chemPatternRegex = /\bchem:([^=]*?)(smiles|mol|biomol|mineral|iupac|pbdid|cid|codid|quan|hybrid)=([^:]+):/gi;
  const chemMatches = result.match(chemPatternRegex);

  console.log('[ChemRenderer] 🔍 DEBUG: chemMatches result:', chemMatches);

  if (chemMatches) {
    // SAFETY: Limit number of matches to prevent excessive processing
    if (chemMatches.length > 100) {
      console.warn('[ChemRenderer] Too many chem: patterns detected, limiting to first 100');
      chemMatches.length = 100;
    }

    log.debug(`  Found ${chemMatches.length} chem:type=value: patterns`);
    console.log('[ChemRenderer] ✅ Found patterns:', chemMatches);
    log.debug(`  Renderer engine: ${settings.rendererEngine}, 3D Viewer enabled: ${settings.enable3DViewer}`);

    // CRITICAL: Match both standard and named syntax
    // chem:mol=benzene: works (name='', type='mol', value='benzene')
    // chem:ethanolsmiles=CCO: works (name='ethanol', type='smiles', value='CCO')
    result = result.replace(/\bchem:([^=]*?)(smiles|mol|biomol|mineral|iupac|pbdid|cid|codid|quan|hybrid)=([^:]+):/gi, (match, displayName, type, valueWithFlags) => {
      try {
        const value_trimmed = valueWithFlags.trim();
        const name_trimmed = displayName ? displayName.trim() : '';

        if (!value_trimmed || value_trimmed.length > 500) {
          return match; // Empty or too long, skip
        }

        // VALIDATION: Reject obviously invalid molecule names
        // Valid: benzene, CCO, 2-methylpropan-1-ol, (R)-butan-2-ol
        // Invalid: "renders, +d (default flags) should only..." (English sentences)
        // Check for spaces followed by words (indicates English text, not molecule names)
        if (/\s+[a-zA-Z]{2,}/.test(value_trimmed)) {
          console.warn('[ChemRenderer] Skipping invalid content (looks like sentence):', value_trimmed);
          return match; // Skip - looks like a sentence, not a molecule
        }
        // Reject if it contains common English words
        const invalidWords = /\b(should|only|with|that|this|the|and|for|not|use|when|from|have|are|was|were)\b/i;
        if (invalidWords.test(value_trimmed)) {
          console.warn('[ChemRenderer] Skipping invalid content (contains English words):', value_trimmed);
          return match;
        }

        // Since we now have explicit type= syntax, construct the content for parseChemFlags
        // For named syntax like chem:ethanolsmiles=CCO:, include the display name
        const content_trimmed = name_trimmed
          ? name_trimmed + type.toLowerCase() + '=' + value_trimmed
          : type.toLowerCase() + '=' + value_trimmed;

        log.debug(`  Converting molecule: ${content_trimmed} (displayName: "${name_trimmed}")`);
        window.chemRendererPerformance.recordStructure();

        // Extract compound name - use display name if provided, otherwise use the value
        let compoundName = name_trimmed || value_trimmed;


        // Process flags if any are detected (always enabled by default)
        let flagOverrides = { useDefaults: false }; // Default empty flags
        let currentSettings = settings; // Default to original settings

        // With the new explicit type=value syntax, we always have the new syntax
        const hasNewSyntax = true;

        if (hasNewSyntax) {
          try {
            // Check if AI Flag Control is enabled
            const enableAIFlagControl = settings.enableAIFlagControl === true;
            console.log('%c🔧 AI Flag Control setting:', 'color: #FF6B00; font-weight: bold;', enableAIFlagControl, 'raw value:', settings.enableAIFlagControl);

            flagOverrides = window.parseChemFlags('chem:' + content_trimmed + ':', enableAIFlagControl);
            console.log('%c🏴 Parsed flags:', 'color: #FF6B00; font-weight: bold;', flagOverrides);

            // Use moleculeName from new syntax if available, or keep display name
            if (flagOverrides.moleculeName && !name_trimmed) {
              compoundName = flagOverrides.moleculeName;
            }

            // Determine base settings:
            // - If +d flag is present: use global defaults as base, then apply flags
            // - If no +d flag AND flagsLocked: start with clean slate, only apply explicit flags
            // - If no flags at all: use global defaults
            let baseSettings;
            if (flagOverrides.useDefaults) {
              // +d flag present: use defaults and apply overrides
              baseSettings = settings;
            } else if (flagOverrides.flagsLocked) {
              // Explicit flags but no +d: start clean, only apply explicit flags
              // This prevents popup settings from overriding explicit flags
              baseSettings = {
                ...settings,
                // Clear all display options to start fresh
                m2cfShowCarbons: false,
                m2cfAromaticCircles: false,
                m2cfShowMethyls: false,
                m2cfAtomNumbers: false,
                m2cfAddH2: false,
                m2cfFlipHorizontal: false,
                m2cfFlipVertical: false,
                m2cfInvert: false
              };
            } else {
              // No explicit flags: use global settings
              baseSettings = settings;
            }

            // Apply flag overrides to base settings
            currentSettings = window.applyFlagOverrides(baseSettings, flagOverrides);
            console.log('%c⚙️ Applied settings:', 'color: #9C27B0; font-weight: bold;', {
              useDefaults: flagOverrides.useDefaults,
              flagsLocked: flagOverrides.flagsLocked,
              showCarbons: currentSettings.m2cfShowCarbons,
              aromaticCircles: currentSettings.m2cfAromaticCircles,
              showMethyls: currentSettings.m2cfShowMethyls,
              atomNumbers: currentSettings.m2cfAtomNumbers,
              addH2: currentSettings.m2cfAddH2,
              flipHorizontal: currentSettings.m2cfFlipHorizontal,
              flipVertical: currentSettings.m2cfFlipVertical,
              invert: currentSettings.m2cfInvert,
              compoundType: flagOverrides.compoundType,
              isDirectSmiles: flagOverrides.isDirectSmiles
            });

            // Extract compound name by removing flag parts (for legacy syntax)
            if (!flagOverrides.moleculeName) {
              if (content_trimmed.includes('/')) {
                compoundName = content_trimmed.split('/')[0].trim();
              } else if (content_trimmed.includes('+')) {
                // For + flags, split on first + to get compound name
                const parts = content_trimmed.split('+');
                compoundName = parts[0].trim();
              } else if (content_trimmed.includes('-')) {
                // For - flags, split on first - to get compound name
                const parts = content_trimmed.split('-');
                compoundName = parts[0].trim();
              }
            }
          } catch (e) {
            console.error('Error processing flags, using default behavior:', e);
            // Fall back to basic behavior
            compoundName = content_trimmed;
            currentSettings = settings;
          }
        }

        // =============== FLAG-BASED TYPE DETECTION ===============
        // Priority (new syntax handles this in parseChemFlags):
        // 1. chem:smiles=CCO: OR +smiles → Direct SMILES rendering
        // 2. chem:biomol=insulin: OR +biomolecule/+protein → Search biomolecule database
        // 3. chem:mineral=quartz: OR +mineral → Search mineral database
        // 4. chem:mol=benzene: OR +compound → Search for compounds
        // 5. No type specified → Auto-detection

        const isDirectSmiles = flagOverrides.isDirectSmiles;
        const specifiedType = flagOverrides.compoundType; // 'compound', 'biomolecule', 'mineral', or 'smiles'

        console.log(`%c📦 Type detection: isDirectSmiles=${isDirectSmiles}, specifiedType=${specifiedType}, compoundName=${compoundName}`, 'color: #2196F3; font-weight: bold;');

        // Build flags object for renderClientSide
        // Include flagsLocked so renderer knows not to override with popup settings
        const flagsForRender = {
          ...flagOverrides,
          compoundType: specifiedType === 'smiles' ? 'compound' : specifiedType,  // 'smiles' type uses compound rendering
          isDirectSmiles: isDirectSmiles,
          smilesValue: flagOverrides.smilesValue,  // Pass actual SMILES for named syntax (e.g., chem:ethanolsmiles=CCO:)
          displayName: flagOverrides.moleculeName  // Pass display name for tag
        };

        // For named SMILES syntax, use smilesValue as the SMILES and moleculeName for display
        const actualSmiles = flagOverrides.smilesValue || (isDirectSmiles ? compoundName : null);

        // =============== DIRECT ID LOOKUPS ===============
        // Handle chem:pbdid=4RHV:, chem:cid=1234:, chem:codid=1234567:
        if (flagOverrides.isDirectID && flagOverrides.idValue) {
          const idValue = flagOverrides.idValue;
          const displayName = flagOverrides.moleculeName || idValue;

          if (flagOverrides.idType === 'pdbid') {
            // Direct PDB ID lookup for biomolecules
            console.log('%c🧬 DIRECT PDB ID LOOKUP:', 'background: #E91E63; color: white; font-weight: bold; padding: 4px;', idValue);
            const biomolData = {
              pdbid: idValue,
              name: displayName,
              type: 'pdbid',
              compoundType: 'biomolecule',
              flags: flagsForRender
            };
            const encoded = btoa(JSON.stringify(biomolData));
            const converted = `<img src="" alt="chem" class="molecule-diagram molecule-biomolecule" data-molecule-viewer="${encoded}" data-loaded="false" style="display: inline-block; height: auto; margin: 0 12px 8px 0; vertical-align: middle; padding: 4px; border-radius: 4px; background-color: transparent;">`;
            log.debug(`  📤 Sending PDB ID to biomolecule viewer: ${idValue}`);
            return converted;
          } else if (flagOverrides.idType === 'cid') {
            // Direct Compound ID lookup for compounds
            console.log('%c🧪 DIRECT CID LOOKUP:', 'background: #4CAF50; color: white; font-weight: bold; padding: 4px;', idValue);
            const compoundData = {
              cid: idValue,
              name: displayName,
              type: 'cid',
              isCompound: true,
              directLookup: true,
              flags: {
                ...flagsForRender,
                compoundType: 'compound',
                skipIntegratedSearch: true
              }
            };
            const encoded = btoa(JSON.stringify(compoundData));
            const converted = `<img src="" alt="chem" class="molecule-diagram molecule-compound" data-molecule-viewer="${encoded}" data-loaded="false" style="display: inline-block; height: auto; margin: 0 12px 8px 0; vertical-align: middle; padding: 4px; border-radius: 4px; background-color: transparent;">`;
            log.debug(`  📤 Sending CID to server: ${idValue}`);
            return converted;
          } else if (flagOverrides.idType === 'codid') {
            // Direct COD ID lookup for minerals
            console.log('%c💎 DIRECT COD ID LOOKUP:', 'background: #9C27B0; color: white; font-weight: bold; padding: 4px;', idValue);
            const mineralData = {
              codid: idValue,
              name: displayName,
              type: 'codid',
              compoundType: 'mineral',
              flags: flagsForRender
            };
            const encoded = btoa(JSON.stringify(mineralData));
            const converted = `<img src="" alt="chem" class="molecule-diagram molecule-mineral" data-molecule-viewer="${encoded}" data-loaded="false" style="display: inline-block; height: auto; margin: 0 12px 8px 0; vertical-align: middle; padding: 4px; border-radius: 4px; background-color: transparent;">`;
            log.debug(`  📤 Sending COD ID to mineral viewer: ${idValue}`);
            return converted;
          }
        }


        // Check for special flags that override renderer selection
        // +3d flag: Use 3D viewer (compound-based)
        if (flagOverrides.is3D) {
          const compoundData = {
            // Use actualSmiles if available, otherwise use compoundName for lookup
            ...(actualSmiles ? { smiles: actualSmiles } : { nomenclature: compoundName }),
            name: compoundName, // Display name for the tag
            type: actualSmiles ? 'smiles' : 'nomenclature',
            isCompound: true,
            show3D: true,  // Signal to show 3D viewer immediately
            auto3D: true,  // Hint to viewer to auto-select best 3D mode
            flags: flagsForRender
          };

          const encoded = btoa(JSON.stringify(compoundData));
          const converted = `<img src="" alt="chem" class="molecule-diagram molecule-compound" data-molecule-viewer="${encoded}" data-loaded="false" data-show-3d="true" style="display: inline-block; height: auto; margin: 0 12px 8px 0; vertical-align: middle; padding: 4px; border-radius: 4px; background-color: transparent;">`;

          log.debug(`  📤 Sending ${actualSmiles ? 'SMILES' : 'nomenclature'} to 3D viewer: ${compoundName}`);
          return converted;
        }

        // =============== QUANTUM ORBITAL HANDLING ===============
        // Handle chem:quan=H+1s:, chem:quan=C+2p:, chem:hybrid=sp3:
        if (flagOverrides.isQuantumOrbital || flagOverrides.isHybridOrbital ||
          specifiedType === 'quantum' || specifiedType === 'hybrid') {
          const isHybrid = flagOverrides.isHybridOrbital || specifiedType === 'hybrid';
          const orbitalValue = flagOverrides.moleculeName || compoundName;
          const displayName = name_trimmed || orbitalValue;

          console.log('%c⚛️ QUANTUM ORBITAL PATH:', 'background: #9C27B0; color: white; font-weight: bold; padding: 4px;',
            { isHybrid, orbitalValue, displayName });

          const orbitalData = {
            nomenclature: orbitalValue,
            name: displayName,
            compoundType: isHybrid ? 'hybrid' : 'quantum',
            flags: {
              ...flagsForRender,
              isQuantumOrbital: !isHybrid,
              isHybridOrbital: isHybrid,
              compoundType: isHybrid ? 'hybrid' : 'quantum'
            }
          };

          const encoded = btoa(JSON.stringify(orbitalData));
          const converted = `<img src="" alt="orbital" class="molecule-diagram molecule-quantum" data-molecule-viewer="${encoded}" data-loaded="false" style="display: inline-block; height: auto; margin: 0 12px 8px 0; vertical-align: middle; padding: 4px; border-radius: 4px; background-color: transparent;">`;

          log.debug(`  📤 Sending ${isHybrid ? 'hybrid' : 'quantum'} orbital to viewer: ${orbitalValue}`);
          return converted;
        }

        // mol= syntax: Force compound lookup
        // This is triggered by chem:mol=benzene: which sets isCompound=true and compoundType='compound'
        if (flagOverrides.isCompound || specifiedType === 'compound') {
          console.log('%c🧪 DIRECT COMPOUND PATH (mol= syntax)', 'background: #4CAF50; color: white; font-weight: bold; padding: 4px;');
          const compoundData = {
            ...(actualSmiles ? { smiles: actualSmiles } : { nomenclature: compoundName }),
            name: compoundName,
            type: actualSmiles ? 'smiles' : 'nomenclature',
            isCompound: true,
            directLookup: true,  // Signal to skip IntegratedSearch
            flags: {
              ...flagsForRender,
              compoundType: 'compound',  // Ensure compound type is set
              skipIntegratedSearch: true  // Ensure IntegratedSearch is skipped
            }
          };

          const encoded = btoa(JSON.stringify(compoundData));
          const converted = `<img src="" alt="chem" class="molecule-diagram molecule-compound" data-molecule-viewer="${encoded}" data-loaded="false" style="display: inline-block; height: auto; margin: 0 12px 8px 0; vertical-align: middle; padding: 4px; border-radius: 4px; background-color: transparent;">`;

          log.debug(`  📤 Sending ${actualSmiles ? 'SMILES' : 'nomenclature'} DIRECTLY to server (mol= syntax): ${compoundName}`);
          return converted;
        }

        // Default: Send to renderer (auto-detect type)

        const rotation = currentSettings.m2cfRotate || currentSettings.mvRotate || 0;
        const scale = flagOverrides.size || 1.0;

        const moleculeViewerData = {
          // If actualSmiles available (from named syntax or direct): use smiles property
          // Otherwise: use nomenclature and let server detect the type
          ...(actualSmiles ? { smiles: actualSmiles } : { nomenclature: compoundName }),
          type: actualSmiles ? 'smiles' : 'nomenclature',
          options: {
            width: Math.round(300 * scale),
            height: Math.round(200 * scale),
            aromaticCircles: currentSettings.m2cfAromaticCircles,
            fancyBonds: currentSettings.m2cfFancyBonds,
            showCarbons: currentSettings.m2cfShowCarbons,
            showMethyls: currentSettings.m2cfShowMethyls,
            atomNumbers: currentSettings.m2cfAtomNumbers,
            hydrogensMode: currentSettings.m2cfHydrogensMode,
            addH2: currentSettings.m2cfAddH2,
            flipHorizontal: currentSettings.m2cfFlipHorizontal,
            flipVertical: currentSettings.m2cfFlipVertical,
            scale: scale,  // Pass through the scale from +s flag
            rotation: flagOverrides.rotation || rotation  // Pass through the rotation from +r flag
          },
          isMoleculeViewer: true,
          flags: flagsForRender  // Pass flags for renderClientSide to use
        };

        const encoded = btoa(JSON.stringify(moleculeViewerData));
        const converted = `<img src="" alt="chem" class="molecule-diagram molecule-viewer" data-molecule-viewer="${encoded}" data-loaded="false" data-rotation="${rotation}" style="display: inline-block; height: auto; margin: 0 12px 8px 0; vertical-align: middle; padding: 4px; border-radius: 4px; background-color: transparent;">`;

        log.debug(`  📤 Sending ${isDirectSmiles ? 'SMILES' : 'nomenclature'} to MoleculeViewer: ${compoundName} (type: ${specifiedType || 'auto-detect'})`);
        return converted;
      } catch (error) {
        console.error('[ChemRenderer] ❌ Error processing molecule:', {
          error: error.message,
          stack: error.stack,
          content: content_trimmed,
          match: match
        });
        return match; // Return original text on error
      }
    });
  }

  // Note: mhchem patterns (\ce{} and ce{}) removed - use chem: syntax instead

  if (patternMatches.length > 0) {
    log.debug(`✅ Total patterns found & converted: ${patternMatches.length}`);
  }

  return result;
}


/**
 * Observe page changes for dynamic content
 * Watches for DOM mutations and re-scans if content changes
 */
function observePageChanges() {
  // Disconnect any existing observer first to prevent memory leaks
  if (window._chemRendererObserver) {
    try {
      window._chemRendererObserver.disconnect();
      clearTimeout(window._chemRendererObserver.timer);
    } catch (e) {
      // Ignore disconnect errors
    }
    window._chemRendererObserver = null;
  }

  log.inject('🔄 Setting up mutation observer for dynamic content detection');

  let mutationCount = 0;
  let isProcessing = false;  // Prevent re-entrant scanning
  let pendingMutations = 0;  // Track pending mutations
  let throttleWarningShown = false;  // Only show warning once per throttle period
  let lastThrottleTime = 0;  // Track when we last throttled
  const MAX_PENDING_MUTATIONS = 1000;  // Increased safety limit
  const THROTTLE_COOLDOWN = 2000;  // 2 seconds between throttle resets

  const observer = new MutationObserver((mutations) => {
    // SAFETY: Skip if already processing to prevent recursive loops
    if (isProcessing) {
      return;
    }

    // SAFETY: Skip mutations from our own DOM changes (check first few only for performance)
    const samplesToCheck = Math.min(mutations.length, 5);
    for (let i = 0; i < samplesToCheck; i++) {
      const target = mutations[i].target;
      if (!target) continue;
      // Check if mutation is inside our rendered elements
      if (target.classList && (
        target.classList.contains('molecule-diagram') ||
        target.classList.contains('molecule-container') ||
        target.classList.contains('molecule-viewer') ||
        target.classList.contains('chem-size-controls') ||
        target.classList.contains('chemistrylatex-image-loading-indicator')
      )) {
        return;  // Skip our own changes
      }
      // Check parent
      if (target.closest && target.closest('.molecule-container, .chem-image-container, .molecule-viewer-container')) {
        return;  // Skip our own changes
      }
    }

    pendingMutations += mutations.length;
    mutationCount++;

    // SAFETY: If too many pending mutations, throttle aggressively
    if (pendingMutations > MAX_PENDING_MUTATIONS) {
      const now = Date.now();

      // Only show warning once per throttle period
      if (!throttleWarningShown || (now - lastThrottleTime > THROTTLE_COOLDOWN)) {
        log.warn(`⚠️ High mutation rate detected (${pendingMutations}), throttling to prevent crashes...`);
        throttleWarningShown = true;
        lastThrottleTime = now;
      }

      // Reset counter but still schedule a delayed scan
      pendingMutations = 0;

      // Clear existing timer and set a longer delay for high-mutation scenarios
      clearTimeout(observer.timer);
      observer.timer = setTimeout(() => {
        throttleWarningShown = false;
        pendingMutations = 0;

        if (settings.enabled && !isProcessing && !isUpdatingImages) {
          isProcessing = true;
          try {
            log.inject('⚡ Delayed re-scan after high mutation period');
            scanAndRender();
          } catch (e) {
            log.error('Error during page scan:', e);
          } finally {
            setTimeout(() => { isProcessing = false; }, 200);
          }
        }
      }, 1500);  // Longer delay when throttling

      return;
    }

    // Only log every 100 mutations to avoid spamming
    if (mutationCount % 100 === 0) {
      log.debug(`Detected ${mutationCount} total mutations`);
    }

    clearTimeout(observer.timer);
    // Use shorter timeout for faster detection (especially for ChatGPT streaming)
    observer.timer = setTimeout(() => {
      // Reset pending count
      pendingMutations = 0;

      // CRITICAL: Skip re-scan if we're currently updating images
      // This prevents images from disappearing when settings change
      if (isUpdatingImages) {
        log.debug('⏭️ Skipping re-scan - currently updating images');
        return;
      }

      if (settings.enabled) {
        isProcessing = true;  // Set flag to prevent re-entrant calls
        try {
          log.inject('⚡ Re-scanning page after dynamic content detected');
          scanAndRender();
        } catch (e) {
          log.error('Error during page scan:', e);
        } finally {
          // Reset flag after a delay to allow DOM to settle
          setTimeout(() => {
            isProcessing = false;
          }, 100);
        }
      }
    }, 500); // Reduced from 1000ms to 500ms for faster response
  });

  // Watch entire document for text node changes (important for ChatGPT)
  observer.observe(document.documentElement, {
    childList: true,
    subtree: true,
    characterData: false  // Don't watch text content changes, only DOM structure
  });

  window._chemRendererObserver = observer;
  log.success('✅ Dynamic content observer initialized');
}

// Cleanup function to disconnect all observers on unload
function cleanupObservers() {
  try {
    if (window._chemRendererObserver) {
      window._chemRendererObserver.disconnect();
      clearTimeout(window._chemRendererObserver.timer);
      window._chemRendererObserver = null;
    }
    if (window._lazyLoadObserver) {
      window._lazyLoadObserver.disconnect();
      window._lazyLoadObserver = null;
    }
    if (lazyReRenderObserver) {
      lazyReRenderObserver.disconnect();
      lazyReRenderObserver = null;
    }
    // Clear caches to free memory
    renderedImageCache.clear();
    pendingSearches.clear();
    imagesBeingRendered.clear();
    log.info('🧹 Cleaned up all observers and caches');
  } catch (e) {
    // Ignore cleanup errors
  }
}

// Register cleanup on page unload
window.addEventListener('beforeunload', cleanupObservers);
window.addEventListener('pagehide', cleanupObservers);

// ============================================
// DARK MODE SUPPORT
// ============================================
// Listen for dark mode changes and update molecule images
if (window.matchMedia) {
  const darkModeQuery = window.matchMedia('(prefers-color-scheme: dark)');

  function updateMoleculeColors(isDark) {
    console.log(`%c🌓 Dark mode ${isDark ? 'enabled' : 'disabled'} - updating molecule colors`, 'color: #00AAFF; font-weight: bold;');

    // Update all molecule images by re-rendering with color replacements
    const moleculeImages = document.querySelectorAll('img.molecule-diagram');
    moleculeImages.forEach(img => {
      // Get original SVG from data URI
      const src = img.src;
      if (src && src.startsWith('data:image/svg+xml;base64,')) {
        try {
          const base64 = src.replace('data:image/svg+xml;base64,', '');
          const svgContent = atob(base64);

          let updatedSvg = svgContent;
          if (isDark) {
            // Replace black with white in dark mode
            updatedSvg = updatedSvg
              .replace(/#000000/gi, '#FFFFFF')
              .replace(/#000\b/gi, '#FFF')
              .replace(/rgb\(0,\s*0,\s*0\)/gi, 'rgb(255, 255, 255)')
              .replace(/stroke="black"/gi, 'stroke="white"')
              .replace(/fill="black"/gi, 'fill="white"');
          } else {
            // In light mode, restore original (undo the replacements)
            updatedSvg = svgContent
              .replace(/#FFFFFF/gi, '#000000')
              .replace(/#FFF\b/gi, '#000')
              .replace(/rgb\(255,\s*255,\s*255\)/gi, 'rgb(0, 0, 0)')
              .replace(/stroke="white"/gi, 'stroke="black"')
              .replace(/fill="white"/gi, 'fill="black"');
          }

          // Store original src to prevent data loss on multiple toggles
          if (!img.dataset.originalSrc) {
            img.dataset.originalSrc = src;
          }

          // If switching to light mode, and we have original, use it to avoid corruption
          if (!isDark && img.dataset.originalSrc) {
            img.src = img.dataset.originalSrc;
            return;
          }

          img.src = 'data:image/svg+xml;base64,' + btoa(updatedSvg);
        } catch (e) {
          console.error('Error updating molecule colors:', e);
        }
      }
    });

    // For PNG images (like compound images), we can't easily change colors like we do with SVGs
    // But we can trigger a re-render if needed by reloading the images
    const compoundImages = document.querySelectorAll('img.compound-diagram');
    compoundImages.forEach(img => {
      // For dark mode, if background was removed, we might want to re-process the image
      // to potentially adjust how it appears in dark mode
      // For now, we'll just log that dark mode changed
      console.log(`🌓 Dark mode ${isDark ? 'enabled' : 'disabled'} - consider reprocessing compound image:`, img.src);
    });
  }

  // Listen for theme changes
  darkModeQuery.addEventListener('change', (e) => {
    updateMoleculeColors(e.matches);
  });

  log.success('✅ Dark mode listener initialized');
}

// ============================================
// END OF CONTENT SCRIPT
// ============================================

log.success('🎉 ChemRenderer loaded');

// ============================================
// CONTEXT MENU HANDLER
// ============================================


// Helper: Create and render molecule image from selection
function renderSelectionAsMolecule(text, config) {
  const selection = window.getSelection();
  if (!selection.rangeCount) return;
  const range = selection.getRangeAt(0);

  const img = document.createElement('img');
  img.className = 'molecule-diagram molecule-viewer';
  img.alt = text;
  img.title = config.title;
  img.src = '';

  const moleculeData = {
    nomenclature: config.nomenclature || text,
    type: config.type || 'nomenclature',
    options: { width: 400, height: 300, scale: 1, rotation: 0 },
    isMoleculeViewer: true,
    flags: config.flags
  };
  if (config.smiles) moleculeData.smiles = config.smiles;

  img.dataset.moleculeViewer = btoa(JSON.stringify(moleculeData));
  img.dataset.loaded = 'false';
  img.dataset.rotation = '0';
  img.dataset.compoundType = config.flags.compoundType;
  if (config.smiles) img.dataset.smiles = config.smiles;

  img.style.cssText = 'display:inline-block;vertical-align:middle;margin:0 4px;min-width:100px;min-height:80px;background:rgba(128,128,128,0.1);border-radius:4px;padding:4px';

  range.deleteContents();
  range.insertNode(img);
  selection.removeAllRanges();

  if (window._loadMoleculeImage) window._loadMoleculeImage(img);
  else { img.alt = 'Error: Renderer not ready'; }
}

chrome.runtime.onMessage.addListener((request, sender, sendResponse) => {
  const text = request.text?.trim();

  // Consolidated render handlers
  switch (request.type) {
    case 'INSPECT_MOLECULE':
      if (!text) return;
      renderSelectionAsMolecule(text, {
        title: `Searching compound database: ${text}`,
        flags: { compoundType: 'compound', isDirectSmiles: false }
      });
      break;

    case 'RENDER_SMILES':
      if (!text) return;
      renderSelectionAsMolecule(text, {
        title: `SMILES: ${text}`,
        type: 'smiles',
        smiles: text,
        flags: { compoundType: 'compound', isDirectSmiles: true }
      });
      break;

    case 'RENDER_BIOMOLECULE':
      if (!text) return;
      renderSelectionAsMolecule(text, {
        title: `Searching protein database: ${text}`,
        flags: { compoundType: 'biomolecule', isDirectSmiles: false }
      });
      break;

    case 'RENDER_MINERAL':
      if (!text) return;
      renderSelectionAsMolecule(text, {
        title: `Searching mineral database: ${text}`,
        flags: { compoundType: 'mineral', isDirectSmiles: false }
      });
      break;

    case 'RENDER_IUPAC':
      if (!text) return;
      renderSelectionAsMolecule(text, {
        title: `Converting IUPAC Name: ${text}`,
        isIUPAC: true,
        flags: { compoundType: 'compound', isDirectSmiles: false, isIUPAC: true }
      });
      break;
  }


  // Handler for "Re-render as..." (right-click on image to change rendering type)
  // Uses the SAME pipeline as chem:: syntax, but changes the compound type
  if (request.type === 'RERENDER_IMAGE') {
    const srcUrl = request.srcUrl;
    const renderAs = request.renderAs;

    console.log('%c🔄 Re-rendering image as:', 'color: #FF5722; font-weight: bold;', renderAs);
    console.log('Looking for image with src:', srcUrl?.substring(0, 100) + '...');

    // Find the image by its src URL
    const images = document.querySelectorAll('img.molecule-diagram, img.molecule-viewer, img[data-molecule-viewer]');
    let targetImg = null;

    for (const img of images) {
      if (img.src === srcUrl) {
        targetImg = img;
        break;
      }
    }

    if (!targetImg) {
      console.warn('❌ Could not find molecule image with the given src');
      return;
    }

    console.log('%c✅ Found target image:', 'color: #4CAF50;', targetImg);

    // Get the original nomenclature or SMILES from the image's data
    let originalText = targetImg.dataset.nomenclature ||
      targetImg.dataset.smiles ||
      targetImg.alt?.replace(/^(SMILES|Biomolecule|Mineral):\s*/, '') ||
      targetImg.title?.replace(/^(Searching|Rendering|PDB|COD).*?:\s*/, '');

    // Try to decode from moleculeViewer dataset
    if (!originalText && targetImg.dataset.moleculeViewer) {
      try {
        const data = JSON.parse(atob(targetImg.dataset.moleculeViewer));
        originalText = data.smiles || data.nomenclature || data.name;
      } catch (e) {
        console.warn('Could not decode moleculeViewer data:', e);
      }
    }

    if (!originalText) {
      console.warn('❌ Could not extract original text from image');
      alert('Could not determine the original molecule name/SMILES from this image.');
      return;
    }

    // Clean up the text
    originalText = originalText.trim();
    console.log('%c📝 Original text:', 'color: #2196F3;', originalText);

    // Create new molecule data based on the requested render type
    let newMoleculeData;
    let newFlags;

    switch (renderAs) {
      case 'molecule':
        // Force compound lookup for compounds (no auto-detect)
        newFlags = { compoundType: 'compound', isDirectSmiles: false };
        newMoleculeData = {
          nomenclature: originalText,
          type: 'nomenclature',
          flags: newFlags
        };
        break;


      case 'smiles':
        // Render directly as SMILES
        newFlags = { compoundType: 'compound', isDirectSmiles: true };
        newMoleculeData = {
          smiles: originalText,
          nomenclature: `SMILES: ${originalText}`,
          type: 'smiles',
          flags: newFlags
        };
        targetImg.dataset.smiles = originalText;
        break;

      case 'biomolecule':
        // Force biomolecule lookup
        newFlags = { compoundType: 'biomolecule', isDirectSmiles: false };
        newMoleculeData = {
          nomenclature: originalText,
          type: 'nomenclature',
          flags: newFlags
        };
        break;

      case 'mineral':
        // Force COD lookup
        newFlags = { compoundType: 'mineral', isDirectSmiles: false };
        newMoleculeData = {
          nomenclature: originalText,
          type: 'nomenclature',
          flags: newFlags
        };
        break;

      default:
        console.warn('Unknown render type:', renderAs);
        return;
    }

    // Add common options
    newMoleculeData.options = {
      width: 400,
      height: 300,
      scale: 1,
      rotation: 0
    };
    newMoleculeData.isMoleculeViewer = true;

    // Update the image's dataset
    targetImg.dataset.moleculeViewer = btoa(JSON.stringify(newMoleculeData));
    targetImg.dataset.loaded = 'false';
    targetImg.dataset.compoundType = newFlags.compoundType || 'auto';
    targetImg.alt = originalText;
    targetImg.title = `Re-rendering as ${renderAs}: ${originalText}`;
    targetImg.src = ''; // Clear src to trigger reload

    console.log('%c🔄 Triggering re-render via _loadMoleculeImage...', 'color: #FF5722;');

    // Re-render the image using the same pipeline
    if (window._loadMoleculeImage) {
      window._loadMoleculeImage(targetImg);
    } else {
      console.error('❌ _loadMoleculeImage not found on window object');
      targetImg.alt = 'Error: Renderer not ready';
    }
  }

  // Handler for "Edit Flags..." - opens a dialog to toggle individual rendering flags
  // Changes are temporary and will be overwritten when global settings change
  if (request.type === 'EDIT_FLAGS') {
    const srcUrl = request.srcUrl;

    console.log('%c🎛️ Opening flag editor dialog', 'color: #9C27B0; font-weight: bold;');

    // Find the image by its src URL
    const images = document.querySelectorAll('img.molecule-diagram, img.molecule-viewer, img[data-molecule-viewer]');
    let targetImg = null;

    for (const img of images) {
      if (img.src === srcUrl) {
        targetImg = img;
        break;
      }
    }

    if (!targetImg) {
      console.warn('❌ Could not find molecule image with the given src');
      return;
    }

    // Get current flags from the image's data
    let currentFlags = {};
    let moleculeName = targetImg.alt || targetImg.title || 'Unknown';

    if (targetImg.dataset.moleculeViewer) {
      try {
        const data = JSON.parse(atob(targetImg.dataset.moleculeViewer));
        currentFlags = data.flags || {};
        moleculeName = data.nomenclature || data.smiles || moleculeName;
      } catch (e) {
        console.warn('Could not decode moleculeViewer data:', e);
      }
    }

    // Remove any existing dialog
    const existingDialog = document.getElementById('chemistrylatex-flag-editor');
    if (existingDialog) existingDialog.remove();

    // Create the dialog
    const dialog = document.createElement('div');
    dialog.id = 'chemistrylatex-flag-editor';
    dialog.innerHTML = `
      <div style="position: fixed; top: 0; left: 0; right: 0; bottom: 0; background: rgba(0,0,0,0.5); z-index: 999999; display: flex; align-items: center; justify-content: center;">
        <div style="background: #fff; border-radius: 12px; padding: 24px; min-width: 350px; max-width: 450px; box-shadow: 0 8px 32px rgba(0,0,0,0.3); font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', Roboto, sans-serif;">
          <h3 style="margin: 0 0 8px 0; font-size: 18px; color: #333;">🎛️ ChemistryLaTeX Flag Editor</h3>
          <p style="margin: 0 0 16px 0; font-size: 13px; color: #666; word-break: break-all;">${moleculeName}</p>
          
          <div style="display: grid; grid-template-columns: 1fr 1fr; gap: 12px; margin-bottom: 20px;">
            <label style="display: flex; align-items: center; gap: 8px; cursor: pointer;">
              <input type="checkbox" id="flag-c" ${currentFlags.showCarbons ? 'checked' : ''} style="width: 18px; height: 18px;">
              <span style="font-size: 14px;">Show Carbons (c)</span>
            </label>
            <label style="display: flex; align-items: center; gap: 8px; cursor: pointer;">
              <input type="checkbox" id="flag-n" ${currentFlags.atomNumbers ? 'checked' : ''} style="width: 18px; height: 18px;">
              <span style="font-size: 14px;">Atom Numbers (n)</span>
            </label>
            <label style="display: flex; align-items: center; gap: 8px; cursor: pointer;">
              <input type="checkbox" id="flag-o" ${currentFlags.aromaticCircles ? 'checked' : ''} style="width: 18px; height: 18px;">
              <span style="font-size: 14px;">Aromatic Rings (o)</span>
            </label>
            <label style="display: flex; align-items: center; gap: 8px; cursor: pointer;">
              <input type="checkbox" id="flag-h" ${currentFlags.addHydrogens ? 'checked' : ''} style="width: 18px; height: 18px;">
              <span style="font-size: 14px;">Show Hydrogens (h)</span>
            </label>
            <label style="display: flex; align-items: center; gap: 8px; cursor: pointer;">
              <input type="checkbox" id="flag-m" ${currentFlags.showMethyls ? 'checked' : ''} style="width: 18px; height: 18px;">
              <span style="font-size: 14px;">Show Methyls (m)</span>
            </label>
            <label style="display: flex; align-items: center; gap: 8px; cursor: pointer;">
              <input type="checkbox" id="flag-i" ${currentFlags.showImplicitHydrogens ? 'checked' : ''} style="width: 18px; height: 18px;">
              <span style="font-size: 14px;">Implicit H (i)</span>
            </label>
            <label style="display: flex; align-items: center; gap: 8px; cursor: pointer;">
              <input type="checkbox" id="flag-p" ${currentFlags.flipHorizontal ? 'checked' : ''} style="width: 18px; height: 18px;">
              <span style="font-size: 14px;">Flip Horizontal (p)</span>
            </label>
            <label style="display: flex; align-items: center; gap: 8px; cursor: pointer;">
              <input type="checkbox" id="flag-q" ${currentFlags.flipVertical ? 'checked' : ''} style="width: 18px; height: 18px;">
              <span style="font-size: 14px;">Flip Vertical (q)</span>
            </label>
            <label style="display: flex; align-items: center; gap: 8px; cursor: pointer; grid-column: span 2;">
              <input type="checkbox" id="flag-3d" ${currentFlags.is3D ? 'checked' : ''} style="width: 18px; height: 18px;">
              <span style="font-size: 14px;">Show 3D Model First (3d)</span>
            </label>
          </div>
          
          <p style="margin: 0 0 16px 0; font-size: 11px; color: #999;">⚠️ Changes are temporary and will be reset when you change global settings in the popup.</p>
          
          <div style="display: flex; gap: 12px; justify-content: flex-end;">
            <button id="flag-cancel" style="padding: 10px 20px; border: 1px solid #ddd; background: #fff; border-radius: 6px; cursor: pointer; font-size: 14px;">Cancel</button>
            <button id="flag-apply" style="padding: 10px 20px; border: none; background: #4CAF50; color: #fff; border-radius: 6px; cursor: pointer; font-size: 14px; font-weight: 500;">Apply</button>
          </div>
        </div>
      </div>
    `;

    document.body.appendChild(dialog);

    // Handle cancel
    document.getElementById('flag-cancel').addEventListener('click', () => {
      dialog.remove();
    });

    // Handle apply
    document.getElementById('flag-apply').addEventListener('click', () => {
      // Get new flag values
      const newFlags = {
        showCarbons: document.getElementById('flag-c').checked,
        atomNumbers: document.getElementById('flag-n').checked,
        aromaticCircles: document.getElementById('flag-o').checked,
        addHydrogens: document.getElementById('flag-h').checked,
        showMethyls: document.getElementById('flag-m').checked,
        showImplicitHydrogens: document.getElementById('flag-i').checked,
        flipHorizontal: document.getElementById('flag-p').checked,
        flipVertical: document.getElementById('flag-q').checked,
        is3D: document.getElementById('flag-3d').checked,
        flagsLocked: true,  // Prevent popup from overriding
        ...currentFlags  // Preserve other flags like compoundType
      };

      console.log('%c🎛️ Applying new flags:', 'color: #9C27B0;', newFlags);

      // Update the image's molecule data
      let moleculeData = {};
      if (targetImg.dataset.moleculeViewer) {
        try {
          moleculeData = JSON.parse(atob(targetImg.dataset.moleculeViewer));
        } catch (e) { }
      }

      moleculeData.flags = newFlags;
      moleculeData.isMoleculeViewer = true;

      // Update dataset
      targetImg.dataset.moleculeViewer = btoa(JSON.stringify(moleculeData));
      targetImg.dataset.loaded = 'false';
      targetImg.src = ''; // Clear src to trigger reload

      // Re-render the image
      if (window._loadMoleculeImage) {
        window._loadMoleculeImage(targetImg);
      }

      dialog.remove();
    });

    // Close on backdrop click
    dialog.querySelector('div').addEventListener('click', (e) => {
      if (e.target === dialog.querySelector('div')) {
        dialog.remove();
      }
    });
  }
});

