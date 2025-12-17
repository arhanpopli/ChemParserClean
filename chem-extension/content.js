/**
 * Chemistry Formula Renderer v4.0 - Fully Client-Side
 * Works on ChatGPT and other strict CSP sites
 * Uses web_accessible_resources to serve scripts from extension
 *
 * NO LOCAL SERVERS NEEDED - Uses:
 * - SmilesDrawer for 2D rendering (client-side)
 * - IntegratedSearch for compound lookup (queries PubChem, RCSB, COD directly)
 * - 3Dmol.js for 3D viewing (client-side)
 * - CDK Depict as optional external API
 */

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
  info: (msg, data) => {
    try {
      console.log(`${LOG_PREFIX} [INFO] ${msg}`, data || '');
      logToStorage('info', msg, data);
    } catch (e) {
      // Silently fail if logging breaks
    }
  },
  success: (msg, data) => {
    try {
      console.log(`%c${LOG_PREFIX} [SUCCESS] ${msg}`, 'color: green; font-weight: bold;', data || '');
      logToStorage('success', msg, data);
    } catch (e) {
      // Silently fail if logging breaks
    }
  },
  warn: (msg, data) => {
    try {
      console.warn(`${LOG_PREFIX} [WARN] ${msg}`, data || '');
      logToStorage('warn', msg, data);
    } catch (e) {
      // Silently fail if logging breaks
    }
  },
  error: (msg, data) => {
    try {
      console.error(`${LOG_PREFIX} [ERROR] ${msg}`, data || '');
      logToStorage('error', msg, data);
    } catch (e) {
      // Silently fail if logging breaks
    }
  },
  debug: (msg, data) => {
    try {
      console.debug(`${LOG_PREFIX} [DEBUG] ${msg}`, data || '');
      logToStorage('debug', msg, data);
    } catch (e) {
      // Silently fail if logging breaks
    }
  },
  inject: (msg, data) => {
    try {
      console.log(`%c${LOG_PREFIX} [INJECT] ${msg}`, 'color: blue; font-weight: bold;', data || '');
      logToStorage('inject', msg, data);
    } catch (e) {
      // Silently fail if logging breaks
    }
  }
};

// Store logs for later inspection (MINIMAL - Memory optimized)
let logHistory = [];
let extensionContextInvalid = false;
const MAX_LOGS = 20;  // Reduced from 100 to save memory

function logToStorage(type, msg, data) {
  // DISABLED: Too much memory overhead
  // In production, rely on console logs instead
  // Only keep most recent logs if needed for debugging
  return;  // Skip logging to save 6GB RAM!
}

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

/**
 * Get cached rendered image
 */
function getCachedImage(smiles, options) {
  const key = generateImageCacheKey(smiles, options);
  return renderedImageCache.get(key);
}

/**
 * Store rendered image in cache
 */
function setCachedImage(smiles, options, imageData) {
  // Enforce cache size limit
  if (renderedImageCache.size >= MAX_IMAGE_CACHE_SIZE) {
    // Remove oldest entry (first key)
    const firstKey = renderedImageCache.keys().next().value;
    renderedImageCache.delete(firstKey);
  }
  const key = generateImageCacheKey(smiles, options);
  renderedImageCache.set(key, imageData);
  log.debug(`📦 Cached image: ${key.substring(0, 50)}...`);
}

/**
 * SMILES Cache: Persistent storage for compound name → SMILES mappings
 * Stored in chrome.storage.local for persistence across sessions
 */
const SMILES_CACHE_KEY = 'chemtex_smiles_cache';
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
 * Clear all caches (useful for debugging)
 */
function clearAllCaches() {
  renderedImageCache.clear();
  smilesCache = {};
  smilesCacheLoaded = false;
  pendingSearches.clear();
  chrome.storage.local.remove(SMILES_CACHE_KEY);
  log.info('🗑️ All caches cleared');
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

  // SmilesDrawer settings only affect compounds (small molecules with SMILES)
  const smilesDrawerSettings = [
    'sdShowCarbons', 'sdAromaticRings', 'sdShowMethyls', 'sdAtomNumbers',
    'sdShowExplicitHydrogens', 'sdShowImplicitHydrogens', 'sdFlipHorizontal', 'sdFlipVertical', 'sdTheme',
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
  const universalSettings = ['showTags', 'performanceMode', 'enableAIFlagControl'];

  // Check what's changed
  const changedKeys = Object.keys(changedSettings);

  for (const key of changedKeys) {
    if (smilesDrawerSettings.includes(key) || compoundMolviewSettings.includes(key)) {
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
  indicator.className = 'chemtex-image-loading-indicator';
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
    const indicator = container.querySelector('.chemtex-image-loading-indicator');
    if (indicator) {
      indicator.remove();
    }
  }

  delete img.dataset.hasLoadingIndicator;
  delete img.dataset.loadingIndicatorElement;
}

function showPendingIndicator(count) {
  // Intentionally disabled: the global top-left pending pill was a UX regression.
  // We keep stale images in place and show per-image in-box indicators instead.
  return;
}

function hidePendingIndicator() {
  // Intentionally disabled (see showPendingIndicator).
  return;
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
  if (newSettings.sdBondThickness !== undefined) settings.sdBondThickness = newSettings.sdBondThickness;
  if (newSettings.sdRotate !== undefined) settings.m2cfRotate = newSettings.sdRotate;
  // Stereochemistry toggle - affects which SMILES type is fetched
  if (newSettings.useStereochemistry !== undefined) {
    settings.m2cfUse3DSmiles = newSettings.useStereochemistry;
    console.log('[ChemRenderer] Stereochemistry setting updated:', settings.m2cfUse3DSmiles);
  }

  // Handle showTags setting
  if (newSettings.showTags !== undefined) {
    applyTagVisibility(newSettings.showTags);
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
      const smiles = img.dataset.smiles;
      const nomenclature = img.dataset.nomenclature || img.alt || 'molecule';
      const compoundType = img.dataset.compoundType || 'compound';

      // For compounds, we need SMILES to re-render
      // For biomolecules/minerals, we just need the nomenclature (they use different renderers)
      if (compoundType === 'compound' && !smiles) {
        continue;
      }

      if ((compoundType === 'biomolecule' || compoundType === 'mineral') && !nomenclature) {
        continue;
      }

      // Skip molecules with locked flags - they shouldn't be affected by global settings
      // Check if flags are locked by decoding the moleculeViewer data
      if (img.dataset.moleculeViewer) {
        try {
          const data = JSON.parse(atob(img.dataset.moleculeViewer));
          if (data.flags?.flagsLocked && !data.flags?.useDefaults) {
            console.log(`[ChemRenderer] ⏭️ Skipping ${nomenclature} - flags are locked (explicit flags set)`);
            continue;
          }
        } catch (e) {
          // If we can't decode, just continue with re-render
        }
      }

      console.log(`[ChemRenderer] 🔄 Re-rendering: ${nomenclature} (type: ${compoundType})`);


      // Build molecule data with CURRENT settings based on compound type
      let moleculeData;

      if (compoundType === 'compound') {
        // SmilesDrawer compound
        moleculeData = {
          nomenclature: nomenclature,
          smiles: smiles,
          type: 'smiles',
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
            useDefaults: false,
            compoundType: 'compound'
          }
        };
      } else if (compoundType === 'biomolecule') {
        // Protein/biomolecule with RCSB/3Dmol settings
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
            useDefaults: false,
            compoundType: 'mineral'
          }
        };
      } else {
        // Default/unknown type
        moleculeData = {
          nomenclature: nomenclature,
          smiles: smiles,
          type: smiles ? 'smiles' : 'nomenclature',
          options: {
            width: parseInt(img.dataset.renderWidth) || 300,
            height: parseInt(img.dataset.renderHeight) || 240
          },
          isMoleculeViewer: true,
          flags: {
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
  log.info(`🔄 Processing image for lazy re-render: ${compoundName} (SMILES: ${smiles.substring(0, 30)}...)`);
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
 * Re-render a single molecule image.
 * Note: Loading indicator is handled by the caller (observer).
 * Old image stays visible until new one is ready.
 */
async function reRenderSingleMolecule(img, newSettings) {
  try {
    const smiles = img.dataset.smiles;
    const compoundType = img.dataset.compoundType || 'compound';

    if (!smiles) {
      log.debug('Skipping re-render: no SMILES data on image');
      return; // No error, just nothing to do - image keeps old content
    }

    // Only apply SmilesDrawer settings to compounds (not proteins/minerals)
    if (compoundType === 'biomolecule' || compoundType === 'mineral') {
      // For proteins/minerals, they use different rendering (MolView/3Dmol)
      // Just keep the existing image - no SmilesDrawer re-render needed
      log.debug(`Skipping re-render for ${compoundType} (uses different renderer)`);
      return;
    }

    if (!window.SmilesDrawer) {
      log.warn('SmilesDrawer not available for re-render');
      throw new Error('SmilesDrawer not available');
    }

    const theme = isDarkModeEnabled() ? 'dark' : 'light';
    const renderOptions = {
      showCarbons: settings.m2cfShowCarbons,
      showMethyls: settings.m2cfShowMethyls,
      aromaticCircles: settings.m2cfAromaticCircles,
      showHydrogens: settings.m2cfAddH2,
      atomNumbers: settings.m2cfAtomNumbers,
      theme: theme,
      width: parseInt(img.dataset.renderWidth) || 300,
      height: parseInt(img.dataset.renderHeight) || 240
    };

    // Check image cache first
    let svgContent = getCachedImage(smiles, renderOptions);

    if (!svgContent) {
      const smilesDrawerOptions = {
        width: renderOptions.width,
        height: renderOptions.height,
        bondThickness: 2.0,
        bondLength: 25,
        shortBondLength: 0.85,
        bondSpacing: 6,
        atomVisualization: 'default',
        isomeric: true,
        debug: false,
        showCarbons: renderOptions.showCarbons,
        terminalCarbons: renderOptions.showMethyls,
        showHydrogens: renderOptions.showHydrogens,
        showImplicitHydrogens: settings.m2cfShowImplicitH,  // Show H counts in labels (CH₃, OH)
        showAromaticRings: renderOptions.aromaticCircles,
        atomNumbering: renderOptions.atomNumbers,
        solidBondColors: !settings.sdGradientColors,
        compactDrawing: settings.m2cfCompactDrawing === true,
        themes: {
          light: { C: '#000000', O: '#FF0000', N: '#0000FF', S: '#FFD700', P: '#FFA500', F: '#00FF00', Cl: '#00FF00', Br: '#8B0000', I: '#800080', H: '#808080', background: '#FFFFFF' },
          dark: { C: '#FFFFFF', O: '#FF6B6B', N: '#6BB5FF', S: '#FFE66D', P: '#FFB347', F: '#98FB98', Cl: '#98FB98', Br: '#FF6961', I: '#DDA0DD', H: '#AAAAAA', background: 'transparent' }
        }
      };

      const svgDrawer = new SmilesDrawer.SvgDrawer(smilesDrawerOptions);

      try {
        const parsedTree = SmilesDrawer.Parser.parse(smiles);
        const svgElement = document.createElementNS('http://www.w3.org/2000/svg', 'svg');
        svgDrawer.draw(parsedTree, svgElement, theme);
        svgContent = svgElement.outerHTML;
        setCachedImage(smiles, renderOptions, svgContent);
      } catch (parseError) {
        log.error('Failed to re-render molecule:', parseError);
        throw parseError; // Re-throw so caller knows it failed
      }
    }

    const svgBase64 = btoa(unescape(encodeURIComponent(svgContent)));
    const nextSrc = `data:image/svg+xml;base64,${svgBase64}`;

    log.debug(`🖼️ Generated new SVG (${svgContent.length} chars), preparing to swap...`);

    // Pre-load the new image before swapping (old image stays visible)
    await new Promise((resolve) => {
      const probe = new Image();
      probe.onload = () => {
        log.debug('✅ New image pre-loaded successfully');
        resolve();
      };
      probe.onerror = (err) => {
        log.error('❌ Failed to pre-load new image:', err);
        resolve(); // Still resolve to continue
      };
      probe.src = nextSrc;
    });

    // Now swap to the new image with a smooth fade-in
    log.debug('🔄 Swapping image src...');
    img.style.animation = 'chemtex-fade-in 160ms ease-out';
    img.src = nextSrc;
    log.debug('✅ Image src updated successfully');

    setTimeout(() => {
      try { img.style.animation = ''; } catch (_) { }
    }, 220);

  } catch (e) {
    log.error('Error re-rendering single molecule:', e);
    throw e; // Re-throw so caller can handle
  }
}

/**
 * Apply tag visibility setting to all molecules
 * When showTags is false: tags are hidden by default but visible on hover
 * When showTags is true: tags are always visible
 */
function applyTagVisibility(showTags) {
  const styleId = 'chemtex-tag-visibility';
  let styleEl = document.getElementById(styleId);

  if (!styleEl) {
    styleEl = document.createElement('style');
    styleEl.id = styleId;
    document.head.appendChild(styleEl);
  }

  if (showTags) {
    // Tags always visible
    styleEl.textContent = `
      .molecule-name-overlay {
        opacity: 0.7 !important;
      }
      .molecule-container:hover .molecule-name-overlay {
        opacity: 1 !important;
      }
    `;
  } else {
    // Tags hidden by default, visible on hover
    styleEl.textContent = `
      .molecule-name-overlay {
        opacity: 0 !important;
      }
      .molecule-container:hover .molecule-name-overlay {
        opacity: 0.7 !important;
      }
    `;
  }
}

/**
 * Reload ALL images on the page (full re-render).
 * This is triggered by the "Reload All" button - immediate re-render for all.
 * Shows loading indicator on each image during re-render.
 */
async function reloadAllImages() {
  log.info('🔄 Reloading all images (forced refresh)...');

  // Clear pending indicator
  hidePendingIndicator();

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

      // Always use lazy re-render: old images stay visible with loading indicator
      // Images only re-render when scrolled into view (IntersectionObserver)
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
      log.info('📬 Clearing all caches...');

      // Wrap in async IIFE to use await
      (async () => {
        try {
          // Clear the in-memory SMILES cache
          if (window.smilesCache) {
            window.smilesCache.clear();
            console.log('[ChemRenderer] Cleared window.smilesCache');
          }

          // Clear any image cache maps
          if (typeof imageCache !== 'undefined') {
            imageCache.clear();
            console.log('[ChemRenderer] Cleared imageCache');
          }

          // Clear pending search cache (to prevent stale lookups)
          if (typeof pendingSearches !== 'undefined') {
            pendingSearches.clear();
            console.log('[ChemRenderer] Cleared pendingSearches');
          }

          // Try to clear IntegratedSearch cache if available
          if (window.IntegratedSearch && window.IntegratedSearch.clearCache) {
            window.IntegratedSearch.clearCache();
            console.log('[ChemRenderer] Cleared IntegratedSearch cache');
          }

          // CRITICAL: Clear the PERSISTENT SMILES cache from chrome.storage.local
          // This is where rhinovirus with wrong SMILES was stored!
          try {
            await chrome.storage.local.remove(SMILES_CACHE_KEY);
            smilesCache = {};  // Reset in-memory cache
            smilesCacheLoaded = false;  // Force reload on next access
            console.log(`[ChemRenderer] 🗑️ Cleared PERSISTENT SMILES cache (${SMILES_CACHE_KEY})`);
          } catch (e) {
            console.error('[ChemRenderer] Failed to clear persistent cache:', e);
          }

          // CRITICAL: Remove cached data from all rendered images so fresh lookup happens
          // This is where rhinovirus with wrong SMILES was stored!
          const allMoleculeImages = document.querySelectorAll('img.molecule-diagram, img.molecule-viewer, img[data-molecule-viewer]');
          let clearedCount = 0;
          allMoleculeImages.forEach(img => {
            // Remove cached SMILES data
            if (img.dataset.smiles) {
              delete img.dataset.smiles;
              clearedCount++;
            }
            // Remove compound type so it gets re-detected
            if (img.dataset.compoundType) {
              delete img.dataset.compoundType;
            }
            // Mark as not loaded so it will be re-rendered
            img.dataset.loaded = 'false';
            // Remove any cached CID/PDBID/CODID
            if (img.dataset.cid) delete img.dataset.cid;
            if (img.dataset.pdbid) delete img.dataset.pdbid;
            if (img.dataset.codid) delete img.dataset.codid;
          });
          console.log(`[ChemRenderer] Cleared cached data from ${clearedCount} rendered images`);

          // Also remove all molecule wrappers and replace with fresh placeholders
          const moleculeWrappers = document.querySelectorAll('.molecule-container');
          moleculeWrappers.forEach(wrapper => {
            try {
              const img = wrapper.querySelector('img');
              if (img && wrapper.parentNode) {
                // Keep the original nomenclature
                const nomenclature = img.dataset.nomenclature || img.alt || 'molecule';

                // Create a fresh placeholder that will trigger a new lookup
                const newPlaceholder = document.createElement('span');
                newPlaceholder.textContent = `chem:${nomenclature}:`;
                newPlaceholder.className = 'chem-molecule-placeholder';

                // Replace wrapper with placeholder text
                wrapper.parentNode.replaceChild(newPlaceholder, wrapper);
              }
            } catch (wrapperError) {
              console.warn('[ChemRenderer] Error processing wrapper:', wrapperError);
            }
          });

          log.success('✅ All caches cleared! Click "Reload All Images" to re-render with fresh data.');
        } catch (e) {
          console.error('[ChemRenderer] Error clearing caches:', e);
        }
      })();  // End async IIFE

      sendResponse({ success: true });
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
// Check if dark mode is enabled for appropriate theming of SmilesDrawer SVGs

/**
 * Detect if the page is in dark mode
 * Returns true only if background is genuinely dark, false for light backgrounds
 */
function isDarkModeEnabled() {
  // Method 1: Check prefers-color-scheme media query
  if (window.matchMedia && window.matchMedia('(prefers-color-scheme: dark)').matches) {
    return true;
  }

  // Method 2: Check body background color (detect dark backgrounds)
  // This is critical to avoid white bonds on white backgrounds
  try {
    const bodyBg = window.getComputedStyle(document.body).backgroundColor;
    if (bodyBg) {
      const rgb = bodyBg.match(/\d+/g);
      if (rgb && rgb.length >= 3) {
        const luminance = (0.299 * parseInt(rgb[0]) + 0.587 * parseInt(rgb[1]) + 0.114 * parseInt(rgb[2]));
        // Only consider dark if luminance is below 100 (stricter threshold)
        // This prevents light gray backgrounds from being treated as dark
        if (luminance < 100) {
          return true;
        }
        // Explicitly return false for light backgrounds (luminance > 180)
        if (luminance > 180) {
          return false;
        }
      }
    }
  } catch (e) {
    // Ignore errors from accessing styles
  }

  // Method 3: Check for common dark mode classes
  const darkClasses = ['dark', 'dark-mode', 'dark-theme', 'night-mode'];
  for (const className of darkClasses) {
    if (document.body?.classList?.contains(className) ||
      document.documentElement?.classList?.contains(className)) {
      return true;
    }
  }

  // Default to light mode (false) to avoid white bonds on white backgrounds
  return false;
}



// ============================================
// BACKGROUND FETCH HELPERS (CSP bypass)
// ============================================
// These functions use the background service worker to make fetch requests,
// bypassing Content Security Policy restrictions on sites like ChatGPT

/**
 * Fetch JSON from API via background script
 * Falls back to direct fetch if background is unavailable
 */
async function backgroundFetchJSON(url, options = {}) {
  // console.log('%c[ChemRenderer] backgroundFetchJSON called', 'color: #00AAFF; font-weight: bold;', { url, options });
  try {
    // Try background script first (works on CSP-restricted sites)
    if (chrome.runtime && chrome.runtime.sendMessage) {
      // console.log('%c[ChemRenderer] Sending message to background...', 'color: #00AAFF;');
      return new Promise((resolve, reject) => {
        chrome.runtime.sendMessage(
          { type: 'FETCH_API', url: url, options: options },
          (response) => {
            // console.log('%c[ChemRenderer] Got response from background:', 'color: #00AAFF;', response);
            if (chrome.runtime.lastError) {
              // console.warn('[ChemRenderer] Background fetch failed, trying direct:', chrome.runtime.lastError.message);
              // Fall back to direct fetch
              directFetchJSON(url, options).then(resolve).catch(reject);
              return;
            }
            if (response && response.success) {
              // console.log('%c[ChemRenderer] Background fetch SUCCESS', 'color: #00FF00; font-weight: bold;');
              resolve(response.data);
            } else {
              console.error('%c[ChemRenderer] Background fetch FAILED:', 'color: #FF0000;', response?.error);
              reject(new Error(response?.error || 'Background fetch failed'));
            }
          }
        );
      });
    } else {
      // console.warn('[ChemRenderer] chrome.runtime.sendMessage not available');
    }
  } catch (e) {
    // Check for extension context invalidated
    if (e.message && e.message.includes('Extension context invalidated')) {
      extensionContextInvalid = true;
      console.warn('[ChemRenderer] Extension context invalidated, using direct fetch');
    }
    // console.warn('[ChemRenderer] Background unavailable, using direct fetch:', e.message);
  }
  // Fall back to direct fetch
  // console.log('%c[ChemRenderer] Falling back to direct fetch', 'color: #FFAA00;');
  return directFetchJSON(url, options);
}

/**
 * Check if chrome.runtime is available (extension context valid)
 */
function isExtensionContextValid() {
  try {
    // Check if extension context is still valid
    if (extensionContextInvalid) return false;
    if (!chrome.runtime || !chrome.runtime.id) return false;
    return true;
  } catch (e) {
    return false;
  }
}

/**
 * Direct fetch for JSON (used as fallback)
 */
async function directFetchJSON(url, options = {}) {
  const response = await fetch(url, options);
  if (!response.ok) {
    throw new Error(`HTTP ${response.status}: ${response.statusText} `);
  }
  return response.json();
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
// Priority: PubChem Direct → PubChem Autocomplete (for class names)
// This centralizes the conversion logic so all renderers use the same flow

/**
 * SMILES Bridge: Convert chemical nomenclature to SMILES string
 * Uses PubChem as the primary source (works for most names)
 *
 * @param {string} name - Chemical name/nomenclature to convert
 * @param {Object} options - Optional settings
 * @param {boolean} options.use3DSmiles - Use IsomericSMILES for stereochemistry (default: false)
 * @returns {Promise<{smiles: string, source: string, has_stereochemistry?: boolean}|null>}
 */
async function smilesBridge(name, options = {}) {
  const { use3DSmiles = false } = options;
  // console.log('%c🌉 SMILES BRIDGE: Converting name→SMILES', 'background: #4A90D9; color: white; font-weight: bold; padding: 2px 6px;', name, use3DSmiles ? '(3D enabled)' : '(3D disabled)');
  if (!name || typeof name !== 'string') {
    console.error('❌ SMILES Bridge: Invalid name provided');
    return null;
  }

  // Strip flags from molecule name before conversion
  // Example: "phenol/d-c" → "phenol", "histamine/+c+o" → "histamine"
  const cleanName = stripFlagsFromName(name.trim());
  if (cleanName !== name.trim()) {
    // console.log('%c🧹 Stripped flags from name:', 'color: #9c88ff;', name, '→', cleanName);
  }

  if (!cleanName) {
    console.error('❌ SMILES Bridge: Empty name after stripping flags');
    return null;
  }

  // ===== PRIORITY 0: IntegratedSearch (CLIENT-SIDE - NO SERVER NEEDED!) =====
  // Uses IntegratedSearch module which queries PubChem, RCSB, COD directly
  // Note: Main rendering in renderClientSide() uses IntegratedSearch - this function is a fallback
  // See renderClientSide() for full biomolecule/mineral handling

  // ===== PRIORITY 1: Direct PubChem API =====
  // Works for most compound names (common names, trade names, drug names)
  // console.log('%c🌐 [Bridge] Priority 1: Trying Direct PubChem API...', 'color: #0088FF; font-weight: bold;');
  try {
    const pubchemUrl = `https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/${encodeURIComponent(cleanName)}/property/CanonicalSMILES,IsomericSMILES/JSON`;
    const data = await backgroundFetchJSON(pubchemUrl);
    if (data && data.PropertyTable && data.PropertyTable.Properties && data.PropertyTable.Properties.length > 0) {
      const props = data.PropertyTable.Properties[0];
      const smiles = use3DSmiles ? (props.IsomericSMILES || props.CanonicalSMILES) : props.CanonicalSMILES;
      // console.log('%c✅ [Bridge] Direct PubChem SUCCESS:', 'color: #00FF00; font-weight: bold;', smiles);
      return { smiles: smiles, source: 'PubChem' };
    }
  } catch (ePubChem) {
    // console.warn('⚠️ [Bridge] Direct PubChem API request failed:', ePubChem.message);
  }

  // ===== PRIORITY 2: PubChem Autocomplete (for class/generic names) =====
  // Handles names like "sphingomyelin" → finds "Sphingomyelin 16:0"
  // console.log('%c🌉 [Bridge] Priority 2: Trying PubChem Autocomplete for class name...', 'color: #FF8800; font-weight: bold;');
  try {
    const autocompleteUrl = `https://pubchem.ncbi.nlm.nih.gov/rest/autocomplete/compound/${encodeURIComponent(cleanName)}/json?limit=5`;
    const autocompleteData = await backgroundFetchJSON(autocompleteUrl);

    if (autocompleteData && autocompleteData.dictionary_terms && autocompleteData.dictionary_terms.compound && autocompleteData.dictionary_terms.compound.length > 0) {
      const bestMatch = autocompleteData.dictionary_terms.compound[0];
      // console.log('%c🔗 [Bridge] Found related compound:', 'color: #FF8800;', bestMatch);

      // Fetch SMILES for the matched compound
      const matchUrl = `https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/${encodeURIComponent(bestMatch)}/property/CanonicalSMILES,IsomericSMILES/JSON`;
      const matchData = await backgroundFetchJSON(matchUrl);

      if (matchData && matchData.PropertyTable && matchData.PropertyTable.Properties && matchData.PropertyTable.Properties.length > 0) {
        const props = matchData.PropertyTable.Properties[0];
        const smiles = use3DSmiles ? (props.IsomericSMILES || props.CanonicalSMILES) : props.CanonicalSMILES;
        // console.log('%c✅ [Bridge] PubChem Autocomplete SUCCESS:', 'color: #00FF00; font-weight: bold;', smiles);
        return { smiles: smiles, source: 'PubChem-Autocomplete' };
      }
    }
  } catch (eAutocomplete) {
    // console.warn('⚠️ [Bridge] PubChem Autocomplete failed:', eAutocomplete.message);
  }

  // All conversion attempts failed
  console.error('%c❌ [Bridge] All name→SMILES conversion attempts failed', 'color: #FF0000; font-weight: bold;');
  return null;
}

// Export for global access (useful for debugging in console)
window.smilesBridge = smilesBridge;

/**
 * Get PubChem CID (Compound ID) for a molecule name or SMILES
 * Uses smilesBridge for complex molecules, direct lookup for simple ones
 * Handles autocomplete fallback for class names like "phosphatidylcholine"
 * 
 * @param {string} nameOrSmiles - Chemical name or SMILES string
 * @returns {Promise<number|null>} - PubChem CID or null if not found
 */
async function getPubChemCID(nameOrSmiles) {
  // console.log('%c🔍 Getting PubChem CID for:', 'color: #0088FF; font-weight: bold;', nameOrSmiles);

  if (!nameOrSmiles || typeof nameOrSmiles !== 'string') {
    console.error('❌ Invalid input for getPubChemCID');
    return null;
  }

  const cleanInput = nameOrSmiles.trim();

  // Helper to fetch CID from name
  const fetchCidByName = async (name) => {
    try {
      const url = `https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/${encodeURIComponent(name)}/cids/JSON`;
      const data = await backgroundFetchJSON(url);
      if (data && data.IdentifierList && data.IdentifierList.CID && data.IdentifierList.CID.length > 0) {
        return data.IdentifierList.CID[0];
      }
    } catch (e) {
      return null;
    }
    return null;
  };

  // Priority 1: Direct lookup
  // console.log('%c🌐 [CID] Priority 1: Direct PubChem CID lookup...', 'color: #0088FF; font-weight: bold;');
  let cid = await fetchCidByName(cleanInput);
  if (cid) {
    // console.log('%c✅ [CID] Direct lookup SUCCESS:', 'color: #00FF00; font-weight: bold;', cid);
    return cid;
  }

  // Priority 1.5: Try appending "_1" (common PubChem pattern for representative compounds like Phosphatidylcholine)
  // This fixes "phosphatidylcholine" -> "Phosphatidylcholine_1" (CID 10425706 which has 3D)
  if (!cleanInput.includes('_')) {
    // console.log('%c🌐 [CID] Priority 1.5: Trying with _1 suffix...', 'color: #0088FF;');
    cid = await fetchCidByName(cleanInput + '_1');
    if (cid) {
      // console.log('%c✅ [CID] Suffix lookup SUCCESS:', 'color: #00FF00; font-weight: bold;', cid);
      return cid;
    }
  }

  // Priority 2: PubChem Autocomplete (The "Singular Solution" for complex molecules)
  // If direct lookup fails, ask PubChem for suggestions and use the first one
  // This handles generic names like "sphingomyelin" -> "Sphingomyelin 16:0" -> CID
  // console.log('%c🔍 [CID] Priority 2: Trying PubChem Autocomplete...', 'color: #FF8800; font-weight: bold;');
  try {
    const autocompleteUrl = `https://pubchem.ncbi.nlm.nih.gov/rest/autocomplete/compound/${encodeURIComponent(cleanInput)}/json?limit=1`;
    const autocompleteData = await backgroundFetchJSON(autocompleteUrl);

    if (autocompleteData && autocompleteData.dictionary_terms && autocompleteData.dictionary_terms.compound && autocompleteData.dictionary_terms.compound.length > 0) {
      const bestMatch = autocompleteData.dictionary_terms.compound[0];
      // console.log('%c🔗 [CID] Autocomplete found match:', 'color: #FF8800; font-weight: bold;', bestMatch);

      // Now get CID for this best match
      cid = await fetchCidByName(bestMatch);
      if (cid) {
        console.log('%c✅ [CID] Autocomplete match SUCCESS:', 'color: #00FF00; font-weight: bold;', cid);
        return cid;
      }
    }
  } catch (e) {
    // console.warn('⚠️ [CID] Autocomplete failed:', e.message);
  }

  // Priority 2.5: Check for Protein/Biomolecule (PDB)
  // If it's a protein name like "rhinovirus", "hemoglobin", etc., PubChem might not have it or it might be a substance.
  // We can try to guess if it's a protein by checking if it fails standard lookup but looks like a biological name.
  // For now, we will rely on the 3D viewer's fallback to MolView which handles PDB search better.
  // But we can add a specific check here if needed in the future.

  // Priority 3: SMILES Bridge fallback (Last resort)
  // console.log('%c🌉 [CID] Priority 3: Using SMILES Bridge fallback...', 'color: #FF8800; font-weight: bold;');
  const bridgeResult = await smilesBridge(cleanInput, { use3DSmiles: false });

  if (bridgeResult && bridgeResult.smiles) {
    // console.log('%c✅ [CID] Got SMILES from bridge:', 'color: #00FF00;', bridgeResult.smiles, `(source: ${bridgeResult.source})`);

    // Now get CID from SMILES
    try {
      const smilesUrl = `https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/smiles/${encodeURIComponent(bridgeResult.smiles)}/cids/JSON`;
      const smilesData = await backgroundFetchJSON(smilesUrl);

      if (smilesData && smilesData.IdentifierList && smilesData.IdentifierList.CID && smilesData.IdentifierList.CID.length > 0) {
        cid = smilesData.IdentifierList.CID[0];
        // console.log('%c✅ [CID] Found CID from SMILES:', 'color: #00FF00; font-weight: bold;', cid);
        return cid;
      }
    } catch (e) {
      // console.warn('⚠️ [CID] Lookup from SMILES failed:', e.message);
    }
  }

  console.error('%c❌ [CID] Could not find CID for:', 'color: #FF0000; font-weight: bold;', nameOrSmiles);
  return null;
}

// Export for global access
window.getPubChemCID = getPubChemCID;

// addHoverControls is defined later in the file (line ~1900)

// Performance monitoring system (MEMORY OPTIMIZED)
window.chemRendererPerformance = {
  metrics: {
    totalStructuresFound: 0,
    totalSVGsRendered: 0,
    totalFormulasTransformed: 0,
    totalLoadTime: 0,
    averageLoadTime: 0,
    maxLoadTime: 0,
    minLoadTime: Infinity,
    startTime: performance.now(),
    loadTimes: [],  // FIXED: Now limited to 50 items (was 1000)
  },

  recordLoad: function (duration) {
    this.metrics.totalSVGsRendered++;
    this.metrics.totalLoadTime += duration;
    this.metrics.loadTimes.push(duration);
    this.metrics.maxLoadTime = Math.max(this.metrics.maxLoadTime, duration);
    this.metrics.minLoadTime = Math.min(this.metrics.minLoadTime, duration);
    this.metrics.averageLoadTime = this.metrics.totalLoadTime / this.metrics.totalSVGsRendered;
    // FIXED: Reduced from 1000 to 50 to save memory
    if (this.metrics.loadTimes.length > 50) this.metrics.loadTimes.shift();
  },

  recordStructure: function () {
    this.metrics.totalStructuresFound++;
  },

  recordFormula: function () {
    this.metrics.totalFormulasTransformed++;
  },

  getStats: function () {
    const uptime = performance.now() - this.metrics.startTime;
    return {
      uptime: `${(uptime / 1000).toFixed(2)}s`,
      structures: this.metrics.totalStructuresFound,
      svgsRendered: this.metrics.totalSVGsRendered,
      formulasTransformed: this.metrics.totalFormulasTransformed,
      avgLoadTime: `${this.metrics.averageLoadTime.toFixed(1)}ms`,
      maxLoadTime: `${this.metrics.maxLoadTime === -Infinity ? 'N/A' : this.metrics.maxLoadTime.toFixed(1) + 'ms'}`,
      minLoadTime: `${this.metrics.minLoadTime === Infinity ? 'N/A' : this.metrics.minLoadTime.toFixed(1) + 'ms'}`,
      totalLoadTime: `${this.metrics.totalLoadTime.toFixed(0)}ms`,
    };
  },

  logStats: function () {
    const stats = this.getStats();
    console.table(stats);
    console.log('%c📊 Performance Analysis:', 'color: #0066cc; font-weight: bold; font-size: 14px;');
    console.log(`   • ${stats.structures} chemical structures detected`);
    console.log(`   • ${stats.svgsRendered} SVGs rendered (avg ${stats.avgLoadTime})`);
    console.log(`   • ${stats.formulasTransformed} formulas transformed locally`);
    console.log(`   • Load time range: ${stats.minLoadTime} - ${stats.maxLoadTime}`);
  },
};

// ============================================
// IMAGE SIZE CONTROLS
// ============================================
const SIZE_STEP = 20;
const MIN_SIZE = 100;
const MAX_SIZE = 800;
const DEFAULT_WIDTH = 400;
const DEFAULT_HEIGHT = 350;

// Non-linear (parabolic) size calculation constants
const BASE_SIZE = 150; // Base size for small molecules
const SIZE_SCALE_FACTOR = 0.5; // Controls the rate of size increase
const MAX_DEFAULT_SIZE = 400; // Maximum default size for very large molecules


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

  console.log(`Applied scale ${scale}x: intrinsic ${intrinsicWidth || '?'}x${intrinsicHeight || '?'}px (max: ${MAX_DISPLAY_WIDTH}x${MAX_DISPLAY_HEIGHT})`);
}

// Create debug interface EARLY so it's always available
window.chemRendererDebug = {
  getLogs: function () {
    console.table(logHistory);
    return logHistory;
  },
  getSettings: function () {
    console.table(settings);
    return settings;
  },
  getPerformanceStats: function () {
    return window.chemRendererPerformance.logStats();
  },
  scanPage: function () {
    console.log('[Manual Scan] Triggering page scan...');
    scanAndRender();
  },
  testFormulas: function () {
    log.info('Testing formula detection...');
    const testCases = [
      '\\ce{H2O}',
      '\\ce{C6H12O6}',
      'ce{Na+ + Cl-}',
      '\\chemfig{-C(-[::30]H)(-[::-30]H)-}',
      'chemfig{C=C}'
    ];
    testCases.forEach(test => {
      const result = wrapChemicalFormulas(test);
      console.log(`${test} → ${result}`);
    });
  },
  getCurrentFormulas: function () {
    console.log('[getCurrentFormulas] Scanning page for chemistry text...');
    const formulas = [];
    const elements = document.getElementsByTagName('*');
    for (let i = 0; i < elements.length; i++) {
      const el = elements[i];
      for (let j = 0; j < el.childNodes.length; j++) {
        const node = el.childNodes[j];
        if (node.nodeType === 3) {
          const text = node.nodeValue;
          if (text.includes('\\') || text.includes('ce{') || text.includes('chem:')) {
            formulas.push(text.substring(0, 200));
          }
        }
      }
    }
    console.log(`Found ${formulas.length} potential chemistry texts:`);
    console.table(formulas);
    return formulas;
  },
  rotateFormulas: function (angle = 90) {
    console.log(`[rotateFormulas] Rotating all structures by ${angle}°`);
    const images = document.querySelectorAll('.molecule-diagram');
    let rotated = 0;
    images.forEach(img => {
      const currentRotation = parseInt(img.style.transform.match(/\d+/)?.[0] || 0);
      const newRotation = (currentRotation + angle) % 360;
      img.style.transform = `rotate(${newRotation}deg)`;
      rotated++;
    });
    console.log(`✅ Rotated ${rotated} structures to ${angle}°`);
    return rotated;
  },
  getRotationHelp: function () {
    console.log(`
    🔄 ROTATION CONTROLS:
    
    Rotate all structures 90° clockwise:
      window.chemRendererDebug.rotateFormulas(90)
    
    Rotate all structures 180°:
      window.chemRendererDebug.rotateFormulas(180)
    
    Rotate all structures 270° (or -90°):
      window.chemRendererDebug.rotateFormulas(270)
    
    Reset to 0°:
      window.chemRendererDebug.rotateFormulas(-270)
    
    Or tell ChatGPT:
      "Rotate the structures ↻" or "Rotate 90° clockwise ↻"
      "Flip upside down ↕" or "Rotate 180° ↕"
      "Rotate counter-clockwise ↺" or "Rotate 270° ↺"
    `);
  },
  togglePerformanceMode: function (enable) {
    if (enable !== undefined) {
      settings.performanceMode = enable;
      log.success(`⚡ Performance mode ${enable ? 'ENABLED' : 'DISABLED'}`);
    } else {
      settings.performanceMode = !settings.performanceMode;
      log.success(`⚡ Performance mode toggled to: ${settings.performanceMode ? 'ON' : 'OFF'}`);
    }
    chrome.storage.sync.set({ performanceMode: settings.performanceMode });
    location.reload();
  },
  setMaxVisibleSVGs: function (count) {
    settings.maxVisibleSVGs = count;
    log.success(`📊 Max visible SVGs set to: ${count}`);
    chrome.storage.sync.set({ maxVisibleSVGs: count });
  },
  getPerformanceStats: function () {
    const stats = {
      performanceModeEnabled: settings.performanceMode,
      maxVisibleSVGs: settings.maxVisibleSVGs,
      totalSVGsOnPage: document.querySelectorAll('img.molecule-diagram').length,
      visibleSVGs: document.querySelectorAll('img.molecule-diagram[data-loaded="true"]').length,
      pendingSVGs: document.querySelectorAll('img.molecule-diagram[data-loaded="false"]').length
    };
    console.table(stats);
    return stats;
  },
  setLayoutMode: function (mode) {
    if (mode !== 'horizontal' && mode !== 'vertical') {
      console.error('❌ Invalid layout mode. Use "horizontal" or "vertical"');
      return;
    }
    settings.layoutMode = mode;
    applyLayoutMode();
    log.success(`📐 Layout mode set to: ${mode.toUpperCase()}`);
    chrome.storage.sync.set({ layoutMode: mode });
    return mode;
  },
  toggleLayoutMode: function () {
    settings.layoutMode = settings.layoutMode === 'horizontal' ? 'vertical' : 'horizontal';
    applyLayoutMode();
    log.success(`📐 Layout mode toggled to: ${settings.layoutMode.toUpperCase()}`);
    chrome.storage.sync.set({ layoutMode: settings.layoutMode });
    return settings.layoutMode;
  },
  getLayoutSettings: function () {
    const layoutInfo = {
      currentMode: settings.layoutMode,
      containersOnPage: document.querySelectorAll('.molecule-container').length,
      structuresOnPage: document.querySelectorAll('.molecule-diagram').length
    };
    console.table(layoutInfo);
    console.log(`📐 Current Layout: ${settings.layoutMode === 'vertical' ? '↕️ VERTICAL (stacked)' : '↔️ HORIZONTAL (side-by-side)'}`);
    return layoutInfo;
  },
  // MEMORY OPTIMIZATION: Clear unnecessary caches
  clearMemory: function () {
    console.log('🧹 Clearing extension memory...');
    let clearedBytes = 0;

    // Clear logs (disabled anyway but just in case)
    const oldLogSize = JSON.stringify(logHistory).length;
    logHistory = [];
    clearedBytes += oldLogSize;

    // Clear old load times (keep only last 10)
    const oldMetricsSize = JSON.stringify(window.chemRendererPerformance.metrics.loadTimes).length;
    if (window.chemRendererPerformance.metrics.loadTimes.length > 10) {
      window.chemRendererPerformance.metrics.loadTimes =
        window.chemRendererPerformance.metrics.loadTimes.slice(-10);
    }
    const newMetricsSize = JSON.stringify(window.chemRendererPerformance.metrics.loadTimes).length;
    clearedBytes += oldMetricsSize - newMetricsSize;

    // Force garbage collection hint
    if (window.gc) {
      window.gc();
      console.log('🗑️  Garbage collection triggered');
    }

    console.log(`✅ Cleared ~${(clearedBytes / 1024 / 1024).toFixed(2)}MB of memory`);
    return {
      cleared: true,
      estimatedBytes: clearedBytes
    };
  },
  // Check current memory usage
  checkMemory: function () {
    if (performance.memory) {
      const used = (performance.memory.usedJSHeapSize / 1048576).toFixed(2);
      const limit = (performance.memory.jsHeapSizeLimit / 1048576).toFixed(2);
      console.log(`💾 Memory Usage: ${used}MB / ${limit}MB`);
      return {
        usedMB: parseFloat(used),
        limitMB: parseFloat(limit),
        percentUsed: ((used / limit) * 100).toFixed(1)
      };
    } else {
      console.log('❌ Memory API not available');
      return null;
    }
  },
  // Toggle carbon label rendering
  toggleCarbonLabels: function (enable) {
    if (enable !== undefined) {
      settings.renderCarbonsAsSticks = enable;
      log.success(`🧪 Carbon labels ${enable ? 'DISABLED (sticks only)' : 'ENABLED (show CH, CH2, CH3)'}`);
    } else {
      settings.renderCarbonsAsSticks = !settings.renderCarbonsAsSticks;
      log.success(`🧪 Carbon labels toggled to: ${settings.renderCarbonsAsSticks ? 'STICKS ONLY' : 'FULL LABELS'}`);
    }
    chrome.storage.sync.set({ renderCarbonsAsSticks: settings.renderCarbonsAsSticks });
    location.reload();
  },
  // Set size preset for responsive rendering
  setSizePreset: function (preset) {
    const validPresets = ['auto', 'small', 'medium', 'large'];
    if (!validPresets.includes(preset)) {
      console.error(`❌ Invalid preset. Use: ${validPresets.join(', ')}`);
      return;
    }
    settings.sizePreset = preset;
    log.success(`📏 Size preset set to: ${preset}`);
    chrome.storage.sync.set({ sizePreset: preset });
    location.reload();
  },

  getRenderingHelp: function () {
    console.log(`
    🧪 RENDERING OPTIONS:
    
    Toggle carbon labels (show sticks vs CH2, CH3):
      window.chemRendererDebug.toggleCarbonLabels()
      window.chemRendererDebug.toggleCarbonLabels(true)   // Sticks only
      window.chemRendererDebug.toggleCarbonLabels(false)  // Show labels
    
    Set size preset:
      window.chemRendererDebug.setSizePreset('auto')      // Auto-detect
      window.chemRendererDebug.setSizePreset('small')     // Inline (150x100)
      window.chemRendererDebug.setSizePreset('medium')    // Default (300x200)
      window.chemRendererDebug.setSizePreset('large')     // Examples (400x300)
    
    Switch rendering engine:
      window.chemRendererDebug.setRendererEngine('codecogs')      // Standard (current)
      window.chemRendererDebug.setRendererEngine('latex-online')  // Alternative
      window.chemRendererDebug.setRendererEngine('quicklatex')    // Fast mode
    
    Or tell ChatGPT:
      "Show carbon atoms in structures 🧪"  → toggles carbon labels
      "Use sticks only for carbons 🔗"      → disables carbon labels
      "Make structures smaller 📉"          → sets small preset
      "Make structures larger 📈"           → sets large preset
      "Use quicklatex rendering ⚡"         → switches to QuickLaTeX
    `);
  },

  // ============================================
  // 🧠 CACHE MANAGEMENT
  // ============================================

  // Clear all caches (SMILES + rendered images)
  clearCache: function () {
    clearAllCaches();
    console.log('✅ All caches cleared');
  },

  // View SMILES cache contents
  viewSmilesCache: async function () {
    await loadSmilesCache();
    const entries = Object.entries(smilesCache);
    console.log(`📦 SMILES Cache (${entries.length} entries):`);
    entries.forEach(([name, data]) => {
      console.log(`  ${name}: ${data.smiles?.substring(0, 40) || 'N/A'}...`);
    });
    return smilesCache;
  },

  // View rendered image cache contents
  viewImageCache: function () {
    console.log(`🖼️ Image Cache (${renderedImageCache.size} entries):`);
    renderedImageCache.forEach((value, key) => {
      console.log(`  ${key.substring(0, 60)}...`);
    });
    return renderedImageCache;
  },

  // Get cache statistics
  getCacheStats: async function () {
    await loadSmilesCache();
    const stats = {
      smilesCache: {
        entries: Object.keys(smilesCache).length,
        maxSize: SMILES_CACHE_MAX_SIZE
      },
      imageCache: {
        entries: renderedImageCache.size,
        maxSize: MAX_IMAGE_CACHE_SIZE
      }
    };
    console.log('📊 Cache Statistics:', stats);
    return stats;
  }

};

// Make logs accessible from console (legacy)
window.chemRendererLogs = () => {
  console.table(logHistory);
  return logHistory;
};

log.inject('🚀 Content script loaded - Extension context injection method');

log.info('🧪 Content script loaded!');
log.info(`Page URL: ${window.location.href}`);
log.info(`Document readyState: ${document.readyState}`);

// Settings with all defaults - Client-side rendering (no servers needed!)
let settings = {
  enabled: true,
  performanceMode: true,
  maxVisibleSVGs: 5,
  sizePreset: 'auto',
  rendererEngine: 'client-side',  // Always SmilesDrawer (no server needed!)
  devMode: false,
  // Rendering options
  flipHorizontal: false,
  flipVertical: false,
  use3DSmiles: false,  // Enable 3D stereochemistry
  // SmilesDrawer options
  sdShowCarbons: true,
  sdAromaticCircles: true,
  sdShowMethyls: true,
  sdCompactDrawing: true,
  sdTheme: 'light',
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

  // Force client-side rendering (only SmilesDrawer supported now)
  settings.rendererEngine = 'client-side';

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

  // Apply showTags setting on startup
  if (result.showTags) {
    applyTagVisibility(true);
  }

  // Log the loaded settings for debugging
  log.info('📊 SmilesDrawer settings:', {
    showCarbons: settings.m2cfShowCarbons,
    aromaticRings: settings.m2cfAromaticCircles,
    showMethyls: settings.m2cfShowMethyls,
    atomNumbers: settings.m2cfAtomNumbers,
    explicitHydrogens: settings.m2cfAddH2,
    implicitHydrogens: settings.m2cfShowImplicitH,
    gradientColors: settings.sdGradientColors
  });

  log.info(`🔍 IntegratedSearch: Queries PubChem, RCSB, COD directly (no local server needed!)`);
  log.success('✅ Settings loaded', settings);
  log.info(`Renderer Engine: 💻 SmilesDrawer (Client-Side)`);
  log.info(`Performance mode: ${settings.performanceMode ? 'ON ⚡' : 'OFF'}`);

  // Pre-load SMILES cache for faster lookups
  loadSmilesCache().then(cache => {
    log.info(`📦 SMILES cache ready: ${Object.keys(cache).length} cached compounds`);
  });

  if (settings.enabled) {
    log.info('🚀 Extension enabled, initializing renderer...');
    initializeRenderer();
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

  // Renderer engine is now always client-side (SmilesDrawer)
  // No need to listen for changes

  // Listen for 3D viewer source changes
  if (changes.viewer3DSource) {
    log.info('🔄 3D Viewer source changed, updating...', changes.viewer3DSource);
    settings.viewer3DSource = changes.viewer3DSource.newValue;
    log.success(`✅ Switched to ${changes.viewer3DSource.newValue} 3D viewer`);
    // Update settings immediately without reload
  }

  // NOTE: SmilesDrawer rendering options (sd*) are now handled via APPLY_SETTINGS message
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

    // Use lazy re-render: old images stay visible with loading indicator
    // Images only re-render when scrolled into view
    lazyReRenderMolecules(request.settings);

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

  // Apply layout mode to any containers
  applyLayoutMode();

  log.inject('🔄 Setting up dynamic content observer...');
  observePageChanges();

  log.success('✅ Renderer initialized - formulas will be processed as they appear');
}

/**
 * Inject CSS styles for proper rendering
 */
function injectStyles() {
  const css = `
    /* Molecule structure diagrams (SmilesDrawer) - inline display with consistent margins */
    .molecule-diagram {
      display: inline-block;
      max-width: 300px;
      max-height: 200px;
      margin: 0 12px 8px 0 !important;  /* Consistent professional spacing */
      vertical-align: middle;
      padding: 4px;
      border-radius: 4px;
      background-color: transparent;  /* Transparent background for dark mode compatibility */
      transition: opacity 0.3s ease;
      cursor: pointer;
      opacity: 1;
    }
    
    /* Loading placeholder - prevents layout shift */
    .molecule-loading {
      opacity: 0.1;
      background: linear-gradient(90deg, #f0f0f0 25%, #e0e0e0 50%, #f0f0f0 75%);
      background-size: 200% 100%;
      animation: shimmer 2s infinite;
    }
    
    @keyframes shimmer {
      0% { background-position: 200% 0; }
      100% { background-position: -200% 0; }
    }
    
    /* Fade-in animation when image loads */
    .molecule-fadein {
      animation: fadeIn 0.4s ease-in-out;
    }
    
    @keyframes fadeIn {
      from {
        opacity: 0;
      }
      to {
        opacity: 1;
      }
    }
    
    /* Container for multiple molecule structures - horizontal layout (DEFAULT) */
    .molecule-container {
      display: inline-flex;
      flex-wrap: wrap;
      gap: 12px;
      align-items: center;
      margin: 8px 0;
    }
    
    .molecule-container .molecule-diagram {
      margin: 0;
    }
    
    /* Vertical layout mode - stack structures top to bottom */
    .molecule-container.vertical {
      display: flex;
      flex-direction: column;
      align-items: flex-start;
      gap: 16px;
      margin: 12px 0;
    }
    
    .molecule-container.vertical .molecule-diagram {
      margin: 0;
      display: block;
    }
    
    /* Rotation utilities - for orientation control */
    .molecule-rotate-0 {
      transform: rotate(0deg);
    }
    
    .molecule-rotate-90 {
      transform: rotate(90deg);
      max-width: 200px;
      max-height: 300px;
    }
    
    .molecule-rotate-180 {
      transform: rotate(180deg);
    }
    
    .molecule-rotate-270 {
      transform: rotate(270deg);
      max-width: 200px;
      max-height: 300px;
    }
    
    /* Hover effects */
    .molecule-diagram:hover {
      background-color: rgba(100, 150, 255, 0.1);
      box-shadow: 0 0 8px rgba(100, 150, 255, 0.2);
    }
    
    /* Molecule Viewer Container - displays image + link + download button */
    .molecule-viewer-container {
      display: inline-block !important;
      margin: 8px 12px 8px 0 !important;
      vertical-align: middle;
      padding: 8px;
      background: #f9f9f9;
      border-radius: 6px;
      border: 1px solid #e0e0e0;
      box-shadow: 0 2px 4px rgba(0,0,0,0.1);
    }
    
    .molecule-viewer-container:hover {
      background: #f0f0f0;
      border-color: #0066cc;
      box-shadow: 0 2px 8px rgba(0,102,204,0.2);
    }
    
    .molecule-viewer-container .molecule-diagram {
      display: block;
      max-width: 300px;
      max-height: 200px;
      margin-bottom: 6px;
      cursor: pointer;
      border-radius: 4px;
    }
    
    .molecule-viewer-container .molecule-diagram:hover {
      background-color: rgba(100, 150, 255, 0.15);
      box-shadow: 0 0 8px rgba(100, 150, 255, 0.3);
    }
    
    .molecule-download-btn:hover {
      background: #45a049 !important;
      transform: scale(1.02);
    }
    
    .molecule-download-btn:active {
      transform: scale(0.98);
    }
    
    /* Mobile responsive */
    @media (max-width: 768px) {
      .molecule-diagram {
        max-width: 250px;
        max-height: 150px;
      }
      
      .molecule-container {
        flex-direction: column;
        gap: 8px;
      }
      
      .molecule-viewer-container {
        margin: 8px 0 !important;
      }
    }
  `;

  const style = document.createElement('style');
  style.textContent = css;
  style.id = 'chem-renderer-styles';
  (document.head || document.documentElement).appendChild(style);
  log.inject('✅ CSS styles injected for formula rendering');
}

// Parse flags from chem: notation and return overrides
// NEW SYNTAX:
//   chem:smiles=CCO:           → Direct SMILES to SmilesDrawer
//   chem:mol=benzene:          → Uses PubChem for SMILES (compound)
//   chem:biomol=rhinovirus:    → Uses RCSB for pdbid
//   chem:mineral=quartz:       → Uses COD for codid/smiles
//   chem:pbdid=4RHV:           → Direct PDB ID lookup (for biomolecules)
//   chem:cid=1234:             → Direct PubChem CID lookup (for compounds)
//   chem:codid=1234567:        → Direct COD ID lookup (for minerals)
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
    isPubchem: false,
    compoundType: undefined,      // 'compound', 'biomolecule', 'mineral', or 'smiles'
    isDirectSmiles: false,
    moleculeName: null,           // The extracted molecule name/value
    isDirectID: false,            // Whether this is a direct ID lookup (pbdid, cid, codid)
    idType: undefined,            // 'pdbid', 'cid', or 'codid'
    idValue: undefined,           // The actual ID value
    smilesValue: undefined,       // The actual SMILES string (for named syntax)
    displayName: undefined        // Custom display name for the molecule tag
  };

  if (!chemFormula) return result;

  // ===== NEW SYNTAX: chem:type=value: OR chem:<name>type=value: format =====
  // Match: chem:smiles=CCO:, chem:mol=benzene+c+n:, chem:biomol=hemoglobin:, chem:mineral=quartz:
  // ALSO: chem:<name>smiles=CCO: where <name> becomes the molecule tag
  // Example: chem:ethanolsmiles=CCO: → name="ethanol", type="smiles", value="CCO"
  // NEW AI-FRIENDLY FORMATS:
  // chem:pbdid=4RHV: → Direct PDB ID lookup (for biomolecules)
  // chem:cid=1234: → Direct PubChem CID lookup (for compounds)
  // chem:codid=1234567: → Direct COD ID lookup (for minerals)

  // First try the standard format: chem:type=value: (check this FIRST to avoid false named matches)
  const newSyntaxMatch = chemFormula.match(/chem:(smiles|mol|biomol|mineral|pbdid|cid|codid)=([^:]+):/i);

  // Only try named format if standard didn't match
  // Named format: chem:<name>type=value: where name must start with a letter and be 2+ chars
  // Examples: chem:ethanolsmiles=CCO:, chem:Aspirinmol=aspirin:, chem:Aspirincid=2244:
  const namedSyntaxMatch = !newSyntaxMatch ?
    chemFormula.match(/chem:([a-zA-Z][a-zA-Z0-9\s]{1,}?)(smiles|mol|biomol|mineral|pbdid|cid|codid)=([^:]+):/i) : null;

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
        result.isPubchem = true;
        break;
      case 'biomol':
        result.compoundType = 'biomolecule';
        break;
      case 'mineral':
        result.compoundType = 'mineral';
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
        result.isPubchem = true;
        break;
      case 'codid':
        result.compoundType = 'mineral';
        result.isDirectID = true;
        result.idType = 'codid';
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
        result.isPubchem = true;
        break;
      case 'biomol':
        result.compoundType = 'biomolecule';
        break;
      case 'mineral':
        result.compoundType = 'mineral';
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
        result.isPubchem = true;
        break;
      case 'codid':
        result.compoundType = 'mineral';
        result.isDirectID = true;
        result.idType = 'codid';
        break;
    }

    // Use the display name for the tag (not the SMILES/value string)
    result.moleculeName = displayName;

    // For SMILES, IDs: store the actual value separately (before flags)
    const isBiomoleculeID = type === 'pbdid';
    const firstFlagIndex = !isBiomoleculeID ? valueWithFlags.search(/[+-]/) : -1;
    const actualValue = firstFlagIndex !== -1 ? valueWithFlags.substring(0, firstFlagIndex).trim() : valueWithFlags.trim();
    
    if (type === 'smiles') {
      result.smilesValue = actualValue;
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
        if (f === 'pubchem') result.isPubchem = true;
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

  // Check for - flags (disable)
  const slashIndex = chemFormula.indexOf('/');
  const afterSlash = slashIndex !== -1 ? chemFormula.substring(slashIndex) : '';
  const minusFlagsMatch = chemFormula.match(/-([cmnhopqi])/gi);

  if (minusFlagsMatch) {
    result.flagsLocked = true;  // Has explicit flags
    minusFlagsMatch.forEach(flag => {
      if (!afterSlash.toLowerCase().includes(flag.toLowerCase())) {
        const f = flag.substring(1).toLowerCase();
        if (f === 'c') result.showCarbons = false;
        if (f === 'o') result.aromaticCircles = false;
        if (f === 'm') result.showMethyls = false;
        if (f === 'n') result.atomNumbers = false;
        if (f === 'h') result.addHydrogens = false;
        if (f === 'i') result.showImplicitHydrogens = false;
        if (f === 'p') result.flipHorizontal = false;
        if (f === 'q') result.flipVertical = false;
      }
    });
  }

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
  `;
  // CRITICAL: NO width, NO height, NO max-width, NO max-height, NO fit-content
  // The wrapper will naturally wrap whatever size the image is

  // Add the image to wrapper
  wrapper.appendChild(img);

  // Add hover controls to wrapper (3D button, name tag, size arrows)
  addHoverControls(wrapper, moleculeName, moleculeData);

  return wrapper;
}

/**
 * Add hover controls (molecule name + 3D viewer button) to an image
 * @param {HTMLElement} element - The element to add controls to
 * @param {string} moleculeName - Name of the molecule
 * @param {object} moleculeData - Full molecule data for 3D viewer
 */
function addHoverControls(element, moleculeName, moleculeData) {
  // Skip if already has controls
  if (element.querySelector('.molecule-name-overlay')) {
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

  // Create 3D viewer button (top-right on hover)
  const viewer3DBtn = document.createElement('button');
  viewer3DBtn.className = 'molecule-3d-btn';
  viewer3DBtn.textContent = '3D';
  viewer3DBtn.title = 'View in 3D';
  viewer3DBtn.style.cssText = `
    position: absolute;
    top: 4px;
    right: 4px;
    background: rgba(0, 120, 255, 0.9);
    color: white;
    border: none;
    padding: 6px 12px;
    border-radius: 4px;
    font-size: 11px;
    font-weight: bold;
    font-family: sans-serif;
    cursor: pointer;
    opacity: 0;
    transition: opacity 0.2s, background 0.2s;
    z-index: 10;
    line-height: 1.4;
  `;

  viewer3DBtn.addEventListener('click', async (e) => {
    e.stopPropagation();
    console.log('🔮 3D button clicked!', { moleculeData, element });
    try {
      await show3DViewerInline(moleculeData, element);
    } catch (error) {
      console.error('❌ Error showing 3D viewer:', error);
      alert('Failed to load 3D viewer: ' + error.message);
    }
  });

  viewer3DBtn.addEventListener('mouseenter', () => {
    viewer3DBtn.style.background = 'rgba(0, 100, 220, 0.9)';
  });

  viewer3DBtn.addEventListener('mouseleave', () => {
    viewer3DBtn.style.background = 'rgba(0, 120, 255, 0.9)';
  });

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

  // Hover effects for buttons
  [upButton, downButton].forEach(btn => {
    btn.addEventListener('mouseenter', () => {
      btn.style.background = 'rgba(0, 0, 0, 0.9)';
    });
    btn.addEventListener('mouseleave', () => {
    });
  });

  // Find the img OR iframe element (iframe for 3D viewer, img for 2D)
  const imgElement = element.querySelector('img');
  const iframeElement = element.querySelector('iframe');
  const targetElement = imgElement || iframeElement;

  // DEBUG: Log what we find
  console.log('🔍 Size control setup:', {
    element: element,
    elementClass: element.className,
    imgElement: imgElement,
    iframeElement: iframeElement,
    targetElement: targetElement,
    targetType: targetElement ? targetElement.tagName : 'none',
    children: Array.from(element.children).map(c => c.tagName)
  });

  upButton.addEventListener('click', (e) => {
    e.stopPropagation();
    console.log('🔼 UP BUTTON CLICKED!');

    // Re-query in case it changed (2D/3D toggle)
    const img = element.querySelector('img') || element.querySelector('iframe');
    console.log('🔍 Looking for img/iframe:', { found: !!img, type: img ? img.tagName : 'none' });

    if (img) {
      // Get current rendered width as baseline (handles SVG data URLs properly)
      const currentWidth = img.offsetWidth || img.width || img.naturalWidth || 300;
      const currentStyle = window.getComputedStyle(img);

      console.log('📐 Current state:', {
        offsetWidth: img.offsetWidth,
        width: img.width,
        naturalWidth: img.naturalWidth,
        styleWidth: img.style.width,
        computedWidth: currentStyle.width,
        using: currentWidth
      });

      const newWidth = Math.round(currentWidth * 1.2); // 20% increase

      // Get current height to maintain aspect ratio
      const currentHeight = img.offsetHeight || img.height || img.naturalHeight || 200;
      const newHeight = Math.round(currentHeight * 1.2); // Same 20% increase

      // FORCE both width AND height - don't use auto, set explicitly!
      img.style.setProperty('width', `${newWidth}px`, 'important');
      img.style.setProperty('height', `${newHeight}px`, 'important');
      img.style.setProperty('max-width', 'none', 'important');
      img.style.setProperty('max-height', 'none', 'important');
      img.style.setProperty('min-width', `${newWidth}px`, 'important');
      img.style.setProperty('min-height', `${newHeight}px`, 'important');

      // CRITICAL: Sync wrapper (element) size so controls stay in corners
      // 'element' is the wrapper that addHoverControls was called on
      element.style.setProperty('width', `${newWidth}px`, 'important');
      element.style.setProperty('height', `${newHeight}px`, 'important');
      element.style.setProperty('max-width', 'none', 'important');
      element.style.setProperty('max-height', 'none', 'important');

      console.log(`✅ Size UP: ${currentWidth}x${currentHeight} → ${newWidth}x${newHeight} (wrapper synced)`);

      // Verify it worked
      setTimeout(() => {
        console.log('✔️ After resize:', {
          offsetWidth: img.offsetWidth,
          styleWidth: img.style.width
        });
      }, 100);
    } else {
      console.error('❌ No img element found for size control');
      console.error('Element structure:', element.innerHTML.substring(0, 200));
    }
  });

  downButton.addEventListener('click', (e) => {
    e.stopPropagation();
    console.log('🔽 DOWN BUTTON CLICKED!');

    const img = element.querySelector('img') || element.querySelector('iframe');
    if (img) {
      const currentWidth = img.offsetWidth || img.width || img.naturalWidth || 300;
      const currentHeight = img.offsetHeight || img.height || img.naturalHeight || 200;

      const newWidth = Math.max(100, Math.round(currentWidth * 0.8333)); // ~20% decrease, min 100px
      const newHeight = Math.max(80, Math.round(currentHeight * 0.8333)); // ~20% decrease, min 80px

      img.style.setProperty('width', `${newWidth}px`, 'important');
      img.style.setProperty('height', `${newHeight}px`, 'important');
      img.style.setProperty('max-width', 'none', 'important');
      img.style.setProperty('max-height', 'none', 'important');
      img.style.setProperty('min-width', `${newWidth}px`, 'important');
      img.style.setProperty('min-height', `${newHeight}px`, 'important');

      // CRITICAL: Sync wrapper (element) size so controls stay in corners
      element.style.setProperty('width', `${newWidth}px`, 'important');
      element.style.setProperty('height', `${newHeight}px`, 'important');
      element.style.setProperty('max-width', 'none', 'important');
      element.style.setProperty('max-height', 'none', 'important');

      console.log(`✅ Size DOWN: ${currentWidth}x${currentHeight} → ${newWidth}x${newHeight} (wrapper synced)`);
    } else {
      console.error('❌ No img element found for size control');
    }
  });

  controlsDiv.appendChild(upButton);
  controlsDiv.appendChild(downButton);

  // Add hover effect to parent element
  element.addEventListener('mouseenter', () => {
    nameOverlay.style.opacity = '1';
    viewer3DBtn.style.opacity = '1';
    controlsDiv.style.opacity = '1';
  });

  element.addEventListener('mouseleave', () => {
    nameOverlay.style.opacity = '0.7';
    viewer3DBtn.style.opacity = '0';
    controlsDiv.style.opacity = '0';
  });

  // Make element positioned so absolute children work
  if (getComputedStyle(element).position === 'static') {
    element.style.position = 'relative';
  }

  element.appendChild(nameOverlay);
  element.appendChild(viewer3DBtn);
  element.appendChild(controlsDiv);
}

/**
 * Get invert filter for dark mode compatibility
 * @returns {string} CSS filter string for inverting colors
 */
function getInvertFilter() {
  // Placeholder - returns empty filter (no inversion)
  // TODO: Implement dark mode detection and proper filter
  return '';
}

/**
 * Apply layout mode to molecule containers
 * Handles compact/expanded view modes
 */
function applyLayoutMode() {
  // Placeholder - no layout changes applied
  // TODO: Implement layout mode switching if needed
  log.info('applyLayoutMode called (placeholder - no action taken)');
}

/**
 * Build RCSB structure image URL
 * Returns either asymmetric unit (_model-1) or biological assembly (_assembly-1) based on settings
 * 
 * @param {string} pdbid - The PDB ID
 * @returns {string} RCSB image URL
 */
function buildRCSBImageUrl(pdbid) {
  const pdbLower = pdbid.toLowerCase();
  // Use biological assembly image if setting is enabled
  if (settings.molviewBioAssembly) {
    return `https://cdn.rcsb.org/images/structures/${pdbLower}_assembly-1.jpeg`;
  }
  // Default: asymmetric unit
  return `https://cdn.rcsb.org/images/structures/${pdbLower}_model-1.jpeg`;
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
    // The bioassembly setting only affects the 2D RCSB image (assembly vs model)
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

// Export helper functions for global access
window.parseChemFlags = parseChemFlags;
window.stripFlagsFromName = stripFlagsFromName;
window.applyFlagOverrides = applyFlagOverrides;
window.getInvertFilter = getInvertFilter;
window.buildMolViewEmbedUrl = buildMolViewEmbedUrl;

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

  // ============================================
  // GLOBAL SMILESDRAWER INSTANCE REUSE
  // ============================================
  // Reuse a single SvgDrawer instance instead of creating new ones for each molecule
  // This significantly reduces memory usage and prevents crashes with many molecules
  let _globalSvgDrawerOptions = null;
  let _globalSvgDrawer = null;

  function getOrCreateSvgDrawer(options) {
    // Check if we need to create a new instance (first time or options changed significantly)
    // We compare key options that affect rendering, not dimensions (those change per molecule)
    const optionsKey = JSON.stringify({
      showCarbons: options.showCarbons,
      terminalCarbons: options.terminalCarbons,
      showHydrogens: options.showHydrogens,
      showImplicitHydrogens: options.showImplicitHydrogens,
      showAromaticRings: options.showAromaticRings,
      atomNumbering: options.atomNumbering,
      solidBondColors: options.solidBondColors,
      compactDrawing: options.compactDrawing,
      bondThickness: options.bondThickness
    });

    if (!_globalSvgDrawer || _globalSvgDrawerOptions !== optionsKey) {
      console.log('%c🔧 Creating/updating global SvgDrawer instance', 'color: #9C27B0; font-weight: bold;');
      _globalSvgDrawer = new SmilesDrawer.SvgDrawer(options);
      _globalSvgDrawerOptions = optionsKey;
    }

    // Return a drawer with updated dimensions (dimensions don't require new instance)
    // Actually, SmilesDrawer uses options at draw time, so we can reuse the instance
    return _globalSvgDrawer;
  }

  // Queue for limiting concurrent SmilesDrawer.parse() calls
  // This prevents browser overload when many molecules are visible at once
  const PARSE_QUEUE_MAX = 3;  // Max concurrent parse operations
  let parseQueueActive = 0;
  const parseQueue = [];

  function queueSmilesDrawerParse(smiles, successCallback, errorCallback) {
    return new Promise((resolve, reject) => {
      const task = { smiles, successCallback, errorCallback, resolve, reject };

      if (parseQueueActive < PARSE_QUEUE_MAX) {
        executeParseTask(task);
      } else {
        parseQueue.push(task);
      }
    });
  }

  function executeParseTask(task) {
    parseQueueActive++;

    try {
      SmilesDrawer.parse(task.smiles,
        function (tree) {
          parseQueueActive--;
          processNextParseTask();
          try {
            task.successCallback(tree);
            task.resolve(tree);
          } catch (e) {
            console.error('[ChemRenderer] Error in parse success callback:', e);
            task.reject(e);
          }
        },
        function (err) {
          parseQueueActive--;
          processNextParseTask();
          try {
            task.errorCallback(err);
            task.reject(err);
          } catch (e) {
            console.error('[ChemRenderer] Error in parse error callback:', e);
            task.reject(e);
          }
        }
      );
    } catch (e) {
      parseQueueActive--;
      processNextParseTask();
      console.error('[ChemRenderer] Error calling SmilesDrawer.parse:', e);
      task.reject(e);
    }
  }

  function processNextParseTask() {
    if (parseQueue.length > 0 && parseQueueActive < PARSE_QUEUE_MAX) {
      const nextTask = parseQueue.shift();
      executeParseTask(nextTask);
    }
  }

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

  // Helper function to render biomolecule 2D image from RCSB
  // Shows the RCSB structure image with hover controls for 3D viewing
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

      // Wrap with size controls and hover controls (including 3D button)
      chrome.storage.sync.get({
        saveSizePerImage: false,
        saveSizeBySMILES: true
      }, async (sizeSettings) => {
        await wrapImageWithSizeControls(bioImg, img, moleculeData, sizeSettings);
      });

    } catch (error) {
      console.error('renderBiomolecule2D error:', error);
      img.alt = 'Failed to load biomolecule image';
      img.dataset.loaded = 'error';
    }
  }

  // Helper function to render molecule client-side using SmilesDrawer (SVG)
  // Uses IntegratedSearch for name→SMILES conversion (no server needed!)
  async function renderClientSide(moleculeData, img) {
    activeLoads++;
    console.log('%c🎨 RENDERCLIENTSIDE CALLED (SmilesDrawer SVG)!', 'background: #222; color: #00FF00; font-size: 16px; padding: 5px;');
    console.log('Image element:', img);
    console.log('Dataset:', img?.dataset);
    log.debug(`🎨 Rendering Client-Side SVG via SmilesDrawer (#${activeLoads})`);

    try {
      // Decode molecule data if not passed directly
      if (!moleculeData && img.dataset.moleculeViewer) {
        moleculeData = JSON.parse(atob(img.dataset.moleculeViewer));
      }

      console.log('%c📦 Molecule data:', 'color: #FF00FF; font-weight: bold;', moleculeData);

      // Check for 3D request - redirect to 3D viewer
      // Pass both moleculeData and original img element (not a container)
      if (moleculeData.is3D || moleculeData.show3D || (moleculeData.options && moleculeData.options.is3D) || (moleculeData.flags && moleculeData.flags.is3D)) {
        // console.log('%c🔮 3D View requested in Client-Side mode', 'color: #FF00FF;');
        // For 3D viewer, we need to pass the actual target element
        // Create a temporary placeholder img if needed for 3D viewer to replace
        let targetForViewer = img;
        if (!img || img.tagName !== 'IMG') {
          // If img is not an actual image element, create a placeholder
          const placeholder = document.createElement('div');
          placeholder.className = 'chem-3d-placeholder';
          placeholder.style.cssText = 'display: inline-block; width: 300px; height: 200px;';
          placeholder._moleculeData = moleculeData;
          if (img && img.parentNode) {
            img.parentNode.replaceChild(placeholder, img);
          }
          targetForViewer = placeholder;
        }
        await show3DViewerInline(moleculeData, targetForViewer);
        decrementActiveLoads();
        return;
      }

      let smiles = moleculeData.smiles;
      // Track both canonical and isomeric SMILES for proper dual-cache storage
      let fetchedCanonicalSmiles = null;
      let fetchedIsomericSmiles = null;
      
      // Use displayName from flags if available (for named syntax like chem:Benzenesmiles=c1ccccc1:)
      // Otherwise fall back to nomenclature, name, or 'molecule'
      const compoundName = moleculeData.flags?.displayName || moleculeData.nomenclature || moleculeData.name || 'molecule';

      // If no SMILES, get it using Search API (MolView-Only) or PubChem
      if (!smiles && moleculeData.nomenclature) {
        const cleanName = moleculeData.nomenclature.trim();

        // Check if 3D stereochemistry is enabled
        const use3DSmiles = settings.m2cfUse3DSmiles === true;

        // ===== CHECK FOR DIRECT SMILES FLAG (+smiles) =====
        // If ChatGPT/AI specified +smiles, treat the input as SMILES directly (skip all API lookups)
        if (moleculeData.flags?.isDirectSmiles) {
          console.log('%c🧪 [Client] Direct SMILES flag detected - using input as SMILES', 'background: #00BCD4; color: white; font-weight: bold; padding: 4px;');
          smiles = cleanName;
          // Skip all API lookups, fall through to rendering
        }

        // ===== CHECK SMILES CACHE FIRST =====
        // This avoids redundant API calls for previously looked up compounds
        // Cache stores: { canonicalSmiles, isomericSmiles, pdbid, codid, compoundType }
        if (!smiles && !moleculeData.flags?.isDirectSmiles) {
          const cachedData = await getCachedSmiles(cleanName);
          if (cachedData) {
            // ===== CACHED BIOMOLECULE (has pdbid) =====
            if (cachedData.pdbid && cachedData.compoundType === 'biomolecule') {
              console.log('%c🎯 [Client] CACHE HIT - BIOMOLECULE!', 'background: #E91E63; color: white; font-weight: bold; padding: 4px;', cleanName, cachedData.pdbid);

              // Render directly using cached pdbid - no need to search RCSB again!
              const biomoleculeData = {
                nomenclature: cleanName,
                compoundType: 'biomolecule',
                embedUrl: buildMolViewEmbedUrl(cachedData.pdbid, 'pdbid', 'biomolecule'),
                imageUrl: buildRCSBImageUrl(cachedData.pdbid),
                pdbid: cachedData.pdbid,
                is3D: false
              };
              await renderBiomolecule2D(biomoleculeData, img);
              decrementActiveLoads();
              return;
            }

            // ===== CACHED MINERAL (has codid) =====
            if (cachedData.codid && cachedData.compoundType === 'mineral') {
              console.log('%c🎯 [Client] CACHE HIT - MINERAL!', 'background: #00BCD4; color: white; font-weight: bold; padding: 4px;', cleanName, cachedData.codid);

              // Render directly using cached codid - no need to search COD again!
              const mineralData = {
                nomenclature: cleanName,
                compoundType: 'mineral',
                embedUrl: `https://embed.molview.org/v1/?codid=${cachedData.codid}`,
                codid: cachedData.codid,
                is3D: true
              };
              await show3DViewerInline(mineralData, img);
              decrementActiveLoads();
              return;
            }

            // ===== CACHED COMPOUND (has canonicalSmiles or isomericSmiles) =====
            if (cachedData.canonicalSmiles || cachedData.isomericSmiles) {
              // Select the appropriate SMILES based on stereochemistry setting
              const selectedSmiles = use3DSmiles
                ? (cachedData.isomericSmiles || cachedData.canonicalSmiles)  // Prefer isomeric when stereo ON
                : (cachedData.canonicalSmiles || cachedData.isomericSmiles); // Prefer canonical when stereo OFF

              // Validate cached SMILES - if it looks like a compound name, skip it (bad cache entry)
              const looksLikeName = /\b(acid|ol|ene|ane|ine|one|ide|ate|ose|yl|methyl|ethyl|propyl|butyl|amino|hydroxy|tartaric|butanol|ethanol)\b/i.test(selectedSmiles);
              const hasNomenclaturePattern = /^\(?[RSrs],?[RSrs]?\)?-?\d*-?[a-zA-Z]{4,}/i.test(selectedSmiles);

              if (looksLikeName || hasNomenclaturePattern) {
                console.warn('%c⚠️ [Client] BAD CACHE ENTRY DETECTED - cached value looks like a name, not SMILES!', 'background: #ff9800; color: black; font-weight: bold; padding: 4px;');
                console.warn('Cached value:', selectedSmiles);
                console.warn('Will re-fetch SMILES for:', cleanName);
                // Don't use this cached value, let it fall through to API lookup
              } else {
                console.log('%c🎯 [Client] SMILES CACHE HIT!', 'background: #4CAF50; color: white; font-weight: bold; padding: 4px;', cleanName);
                smiles = selectedSmiles;
                if (use3DSmiles) {
                  console.log('%c🔬 [Client] Using cached IsomericSMILES (stereochemistry ON):', 'color: #4CAF50;', smiles);
                } else {
                  console.log('%c🔬 [Client] Using cached CanonicalSMILES (stereochemistry OFF):', 'color: #FF9800;', smiles);
                }
                // Skip all API lookups, fall through to rendering
              }
            }
          }
        }

        // ===== CHECK FOR COMPOUND TYPE FLAG (skip IntegratedSearch if specified) =====
        // If ChatGPT/AI specified the type, use that directly instead of searching
        // Skip if we already have SMILES from +smiles flag
        if (!smiles && moleculeData.flags?.compoundType) {
          const specifiedType = moleculeData.flags.compoundType;
          console.log(`%c🏷️ [Client] Compound type specified by flag: ${specifiedType}`, 'background: #FF9800; color: white; font-weight: bold; padding: 4px;');

          if (specifiedType === 'biomolecule') {
            // For biomolecules, search RCSB only
            try {
              const rcsbSearchUrl = `https://search.rcsb.org/rcsbsearch/v2/query`;
              const rcsbQuery = {
                "query": { "type": "terminal", "service": "full_text", "parameters": { "value": cleanName } },
                "return_type": "entry",
                "request_options": { "return_all_hits": false, "paginate": { "start": 0, "rows": 1 } }
              };
              const rcsbResp = await fetch(rcsbSearchUrl, { method: 'POST', headers: { 'Content-Type': 'application/json' }, body: JSON.stringify(rcsbQuery) });
              if (rcsbResp.ok) {
                const rcsbData = await rcsbResp.json();
                if (rcsbData.result_set && rcsbData.result_set.length > 0) {
                  const pdbid = rcsbData.result_set[0].identifier;
                  const biomoleculeData = {
                    nomenclature: cleanName,
                    compoundType: 'biomolecule',
                    embedUrl: buildMolViewEmbedUrl(pdbid, 'pdbid', 'biomolecule'),
                    imageUrl: buildRCSBImageUrl(pdbid),
                    pdbid: pdbid,
                    is3D: false
                  };
                  await renderBiomolecule2D(biomoleculeData, img);
                  decrementActiveLoads();
                  return;
                }
              }
            } catch (e) { console.warn('RCSB search failed:', e.message); }
          } else if (specifiedType === 'mineral') {
            // For minerals, search COD only
            try {
              const codSearchUrl = `https://www.crystallography.net/cod/result?text=${encodeURIComponent(cleanName)}&format=json`;
              const codResp = await fetch(codSearchUrl);
              if (codResp.ok) {
                const codData = await codResp.json();
                if (codData && codData.length > 0) {
                  const codid = codData[0].file || codData[0].id;
                  const mineralData = {
                    nomenclature: cleanName,
                    compoundType: 'mineral',
                    embedUrl: `https://embed.molview.org/v1/?codid=${codid}`,
                    codid: codid,
                    is3D: true
                  };
                  await show3DViewerInline(mineralData, img);
                  decrementActiveLoads();
                  return;
                }
              }
            } catch (e) { console.warn('COD search failed:', e.message); }
          }
          // For 'compound' type or if specific search failed, fall through to PubChem
        }

        // ===== PRIORITY 0: IntegratedSearch (NO SERVER NEEDED!) =====
        // Replicates search-server.js logic - properly detects biomolecules, minerals, compounds
        // Skip if compound type was explicitly specified OR if using direct SMILES
        let searchResult = null;
        const skipIntegratedSearch = moleculeData.flags?.isDirectSmiles ||
          moleculeData.flags?.compoundType === 'compound' ||
          moleculeData.flags?.compoundType === 'smiles' ||
          moleculeData.flags?.skipIntegratedSearch ||
          moleculeData.directLookup;

        if (skipIntegratedSearch) {
          console.log('%c📝 [Client] IntegratedSearch SKIPPED (explicit compound type or directLookup)', 'color: #4CAF50; font-weight: bold;');
        }

        // Wait for IntegratedSearch to load if not immediately available
        let integratedSearchAvailable = !!window.IntegratedSearch;
        if (!integratedSearchAvailable && !skipIntegratedSearch && !smiles) {
          // Wait up to 500ms for IntegratedSearch to load
          for (let i = 0; i < 5 && !window.IntegratedSearch; i++) {
            await new Promise(r => setTimeout(r, 100));
          }
          integratedSearchAvailable = !!window.IntegratedSearch;
          if (integratedSearchAvailable) {
            console.log('%c🔍 [Client] IntegratedSearch loaded after wait', 'color: #9C27B0;');
          }
        }

        if (!skipIntegratedSearch && !smiles && integratedSearchAvailable) {
          console.log('%c🔍 [Client] Using IntegratedSearch module (no server needed!)', 'background: #9C27B0; color: white; font-weight: bold; padding: 4px;');
          try {
            // Use pending search deduplication to prevent duplicate API calls
            searchResult = await getOrCreatePendingSearch(cleanName, async () => {
              const searchOptions = {
                format: 'compact',
                searchPubChem: settings.searchPubChem !== false,
                searchRCSB: settings.searchRCSB !== false,
                searchCOD: settings.searchCOD !== false
              };
              return window.IntegratedSearch.search(cleanName, searchOptions);
            });

            if (searchResult && !searchResult.error) {
              console.log('%c✅ [Client] IntegratedSearch result:', 'color: #4CAF50; font-weight: bold;', {
                type: searchResult.primary_type,
                name: searchResult.name,
                pdbid: searchResult.pdbid,
                codid: searchResult.codid,
                smiles: searchResult.canonical_smiles ? searchResult.canonical_smiles.substring(0, 30) + '...' : 'N/A'
              });

              // ============ HANDLE BIOMOLECULE (show 2D RCSB image with 3D button) ============
              if (searchResult.primary_type === 'biomolecule' && searchResult.pdbid) {
                console.log('%c🧬 [Client] BIOMOLECULE detected - showing 2D RCSB image', 'background: #E91E63; color: white; font-weight: bold; padding: 4px;');

                // ===== CACHE THE PDBID FOR FUTURE USE =====
                await setCachedSmiles(cleanName, {
                  pdbid: searchResult.pdbid,
                  compoundType: 'biomolecule'
                });
                console.log('%c📦 [Client] Cached biomolecule:', 'color: #E91E63;', cleanName, '→', searchResult.pdbid);

                // Build moleculeData for biomolecule with image URL and embed URL for 3D
                const biomoleculeData = {
                  nomenclature: searchResult.name || cleanName,
                  compoundType: 'biomolecule',
                  embedUrl: buildMolViewEmbedUrl(searchResult.pdbid, 'pdbid', 'biomolecule'),
                  imageUrl: searchResult.image_url || buildRCSBImageUrl(searchResult.pdbid),
                  pdbid: searchResult.pdbid,
                  is3D: false  // Show 2D first, 3D via button
                };

                // Render the 2D RCSB image
                await renderBiomolecule2D(biomoleculeData, img);
                decrementActiveLoads();
                return;
              }

              // ============ GET SMILES FROM RESULT ============
              // Extract BOTH canonical and isomeric SMILES for caching
              if (searchResult.canonical_smiles || searchResult.isomeric_smiles) {
                // Store both for caching
                if (!fetchedCanonicalSmiles && searchResult.canonical_smiles) {
                  fetchedCanonicalSmiles = searchResult.canonical_smiles.replace(/\s+/g, '');
                }
                if (!fetchedIsomericSmiles && searchResult.isomeric_smiles) {
                  fetchedIsomericSmiles = searchResult.isomeric_smiles.replace(/\s+/g, '');
                }

                // Select based on stereochemistry setting
                smiles = use3DSmiles
                  ? (fetchedIsomericSmiles || fetchedCanonicalSmiles)
                  : (fetchedCanonicalSmiles || fetchedIsomericSmiles);

                console.log(`%c🧬 IntegratedSearch SMILES: using ${use3DSmiles ? 'isomeric' : 'canonical'}`, 'color: #9C27B0;',
                  { canonical: fetchedCanonicalSmiles?.substring(0, 25), isomeric: fetchedIsomericSmiles?.substring(0, 25) });
              }

              // ============ HANDLE MINERAL WITHOUT SMILES (use existing 3D viewer pipeline) ============
              if (searchResult.primary_type === 'mineral' && !smiles && searchResult.codid) {
                console.log('%c💎 [Client] MINERAL detected - using existing 3D viewer pipeline', 'background: #00BCD4; color: white; font-weight: bold; padding: 4px;');

                // ===== CACHE THE CODID FOR FUTURE USE =====
                await setCachedSmiles(cleanName, {
                  codid: searchResult.codid,
                  compoundType: 'mineral'
                });
                console.log('%c📦 [Client] Cached mineral:', 'color: #00BCD4;', cleanName, '→', searchResult.codid);

                // Build moleculeData with embedUrl for the existing show3DViewerInline function
                const mineralData = {
                  nomenclature: searchResult.name || cleanName,
                  compoundType: 'mineral',
                  embedUrl: `https://embed.molview.org/v1/?codid=${searchResult.codid}`,
                  codid: searchResult.codid,
                  formula: searchResult.formula,
                  is3D: true
                };

                // Use the existing 3D viewer pipeline which handles sizing, settings, buttons, etc.
                await show3DViewerInline(mineralData, img);
                decrementActiveLoads();
                return;
              }

            } else {
              console.warn('%c⚠️ [Client] IntegratedSearch returned no result, trying PubChem fallback...', 'color: #FF9800;', searchResult?.error);
            }
          } catch (error) {
            console.warn('%c⚠️ [Client] IntegratedSearch failed, trying PubChem fallback:', 'color: #FF9800;', error.message);
          }
        } else if (skipIntegratedSearch) {
          console.log('%c📝 [Client] IntegratedSearch skipped (explicit compound type specified)', 'color: #9C27B0;');
        } else if (smiles) {
          // Already have SMILES, no need to search
        } else if (!integratedSearchAvailable) {
          console.error('%c❌ [Client] IntegratedSearch failed to load after 500ms wait! Check if integrated-search.js is loaded correctly.', 'color: #FF0000; font-weight: bold;');
        }

        // ===== PRIORITY 1: PubChem direct API (FALLBACK) =====
        // Only used when MolView Search fails, unless molviewOnlyMode is enabled
        // Uses fetchedCanonicalSmiles and fetchedIsomericSmiles declared at parent scope
        if (!smiles && !settings.molviewOnlyMode) {
          // Works for common names, trade names, drug names, and most molecules
          // console.log('%c🌐 [Client] Priority 1: Trying PubChem direct lookup (via background)...', 'color: #0088FF; font-weight: bold;');
          try {
            const pubchemUrl = `https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/${encodeURIComponent(cleanName)}/property/CanonicalSMILES,IsomericSMILES/JSON`;
            // console.log('%c🔗 [Client] PubChem URL:', 'color: #0088FF;', pubchemUrl);
            const data = await backgroundFetchJSON(pubchemUrl);

            // console.log('%c📦 [Client] PubChem raw response:', 'color: #FF00FF; font-weight: bold;', JSON.stringify(data));

            if (data && data.PropertyTable && data.PropertyTable.Properties && data.PropertyTable.Properties[0]) {
              const props = data.PropertyTable.Properties[0];
              console.log('%c📦 [Client] props:', 'color: #00FFFF;', props);

              // Extract BOTH canonical and isomeric SMILES for proper caching
              fetchedCanonicalSmiles = props.CanonicalSMILES || props.ConnectivitySMILES;
              fetchedIsomericSmiles = props.IsomericSMILES || props.SMILES;

              // Remove any spaces from SMILES (PubChem sometimes adds them for readability)
              if (fetchedCanonicalSmiles) {
                fetchedCanonicalSmiles = fetchedCanonicalSmiles.replace(/\s+/g, '');
              }
              if (fetchedIsomericSmiles) {
                fetchedIsomericSmiles = fetchedIsomericSmiles.replace(/\s+/g, '');
              }

              // Select the appropriate SMILES based on stereochemistry setting
              if (use3DSmiles) {
                smiles = fetchedIsomericSmiles || fetchedCanonicalSmiles;
                console.log('%c🔬 [Client] Using IsomericSMILES (stereochemistry ON):', 'color: #4CAF50;', smiles);
              } else {
                smiles = fetchedCanonicalSmiles || fetchedIsomericSmiles;
                console.log('%c🔬 [Client] Using CanonicalSMILES (stereochemistry OFF):', 'color: #FF9800;', smiles);
              }

              console.log('%c✅ [Client] PubChem direct SUCCESS - got both SMILES:', 'color: #00FF00; font-weight: bold;',
                { canonical: fetchedCanonicalSmiles?.substring(0, 30), isomeric: fetchedIsomericSmiles?.substring(0, 30) });
            }
          } catch (pubchemError) {
            // console.warn('⚠️ [Client] PubChem direct failed:', pubchemError.message);
          }

          // ===== PRIORITY 2: PubChem Autocomplete (for generic/class names like "sphingomyelin") =====
          // "sphingomyelin" → "Sphingomyelin 16:0", "phosphatidylcholine" → first match
          if (!smiles) {
            // console.log('%c🔎 [Client] Priority 2: Trying PubChem autocomplete for class name...', 'color: #FF8800; font-weight: bold;');
            try {
              const autocompleteUrl = `https://pubchem.ncbi.nlm.nih.gov/rest/autocomplete/compound/${encodeURIComponent(cleanName)}/json?limit=5`;
              const autocompleteData = await backgroundFetchJSON(autocompleteUrl);

              if (autocompleteData && autocompleteData.dictionary_terms && autocompleteData.dictionary_terms.compound && autocompleteData.dictionary_terms.compound.length > 0) {
                // Get the first matching compound name
                const bestMatch = autocompleteData.dictionary_terms.compound[0];
                // console.log('%c🔗 [Client] Found related compound:', 'color: #FF8800;', bestMatch);

                // Now lookup SMILES for the best match - get BOTH types
                const matchUrl = `https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/${encodeURIComponent(bestMatch)}/property/CanonicalSMILES,IsomericSMILES/JSON`;
                const matchData = await backgroundFetchJSON(matchUrl);

                if (matchData && matchData.PropertyTable && matchData.PropertyTable.Properties && matchData.PropertyTable.Properties[0]) {
                  const props = matchData.PropertyTable.Properties[0];

                  // Extract both SMILES types for caching
                  if (!fetchedCanonicalSmiles) {
                    fetchedCanonicalSmiles = (props.CanonicalSMILES || props.ConnectivitySMILES)?.replace(/\s+/g, '');
                  }
                  if (!fetchedIsomericSmiles) {
                    fetchedIsomericSmiles = (props.IsomericSMILES || props.SMILES)?.replace(/\s+/g, '');
                  }

                  // Select based on setting
                  smiles = use3DSmiles
                    ? (fetchedIsomericSmiles || fetchedCanonicalSmiles)
                    : (fetchedCanonicalSmiles || fetchedIsomericSmiles);

                  console.log('%c✅ [Client] PubChem autocomplete SUCCESS:', 'color: #00FF00; font-weight: bold;',
                    { canonical: fetchedCanonicalSmiles?.substring(0, 25), isomeric: fetchedIsomericSmiles?.substring(0, 25) });
                }
              }
            } catch (autocompleteError) {
              // console.warn('⚠️ [Client] PubChem autocomplete failed:', autocompleteError.message);
            }
          }

          // ===== PRIORITY 3: Try getting CID first, then SMILES from CID =====
          if (!smiles) {
            // console.log('%c🆔 [Client] Priority 3: Trying CID lookup...', 'color: #9c27b0; font-weight: bold;');
            try {
              const cid = await getPubChemCID(cleanName);
              if (cid) {
                // console.log('%c✅ [Client] Found CID:', 'color: #9c27b0;', cid);

                // Get SMILES from CID - extract BOTH types
                const cidUrl = `https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/${cid}/property/CanonicalSMILES,IsomericSMILES/JSON`;
                const cidData = await backgroundFetchJSON(cidUrl);

                if (cidData && cidData.PropertyTable && cidData.PropertyTable.Properties && cidData.PropertyTable.Properties[0]) {
                  const props = cidData.PropertyTable.Properties[0];

                  // Extract both SMILES types for caching
                  if (!fetchedCanonicalSmiles) {
                    fetchedCanonicalSmiles = (props.CanonicalSMILES || props.ConnectivitySMILES)?.replace(/\s+/g, '');
                  }
                  if (!fetchedIsomericSmiles) {
                    fetchedIsomericSmiles = (props.IsomericSMILES || props.SMILES)?.replace(/\s+/g, '');
                  }

                  // Select based on setting
                  smiles = use3DSmiles
                    ? (fetchedIsomericSmiles || fetchedCanonicalSmiles)
                    : (fetchedCanonicalSmiles || fetchedIsomericSmiles);

                  console.log('%c✅ [Client] CID lookup SUCCESS:', 'color: #9c27b0; font-weight: bold;',
                    { canonical: fetchedCanonicalSmiles?.substring(0, 25), isomeric: fetchedIsomericSmiles?.substring(0, 25) });
                }
              }
            } catch (cidError) {
              // console.warn('⚠️ [Client] CID lookup failed:', cidError.message);
            }
          }
        } // End of if (!smiles && !settings.molviewOnlyMode)
      }

      if (!smiles) {
        throw new Error(`Could not find SMILES for "${compoundName}" - not found in PubChem or invalid compound name`);
      }

      // ===== VALIDATE SMILES BEFORE USING =====
      // Skip validation if user explicitly specified +smiles flag
      // SMILES should not look like a compound name (no words like "acid", "ol", etc.)
      if (!moleculeData.flags?.isDirectSmiles) {
        const looksLikeName = /\b(acid|ol|ene|ane|ine|one|ide|ate|ose|yl|methyl|ethyl|propyl|butyl|amino|hydroxy|chloro|bromo|iodo|nitro|tartaric|glucose|fructose|butanol|ethanol|methanol|propanol|benzene|phenyl|insulin|hemoglobin|quartz)\b/i.test(smiles);
        const hasNomenclaturePattern = /^\(?[RSrs],?[RSrs]?\)?-?\d*-?[a-zA-Z]{4,}/i.test(smiles); // Like "(R)-2-butanol" or "(R,R)-tartaric"

        if (looksLikeName || hasNomenclaturePattern) {
          console.error('%c❌ SMILES VALIDATION FAILED - This looks like a compound name, not SMILES!', 'background: #f44336; color: white; font-weight: bold; padding: 4px;');
          console.error('Invalid SMILES:', smiles);
          console.error('Original nomenclature:', compoundName);
          throw new Error(`Invalid SMILES received for "${compoundName}" - got "${smiles}" which appears to be a compound name, not SMILES notation`);
        }
      } else {
        console.log('%c✅ Skipping SMILES validation - user specified +smiles flag', 'color: #4CAF50;');
      }

      // ===== CACHE THE SMILES FOR FUTURE USE =====
      // Store BOTH canonical and isomeric SMILES for consistent switching
      if (moleculeData.nomenclature && !moleculeData.flags?.isDirectSmiles) {
        const existingCache = await getCachedSmiles(moleculeData.nomenclature);
        // Only cache if we don't have existing data OR if we have new data from PubChem fetch
        if (!existingCache || (!existingCache.canonicalSmiles && fetchedCanonicalSmiles)) {
          // Use fetched SMILES if available (from PubChem direct), otherwise derive from current smiles
          const canonicalToCache = fetchedCanonicalSmiles || stripStereochemistry(smiles);
          const isomericToCache = fetchedIsomericSmiles || smiles;

          await setCachedSmiles(moleculeData.nomenclature, {
            canonicalSmiles: canonicalToCache,
            isomericSmiles: isomericToCache,
            compoundType: moleculeData.flags?.compoundType || 'compound'
          });
          console.log('%c💾 [Client] Cached both SMILES:', 'color: #9C27B0;',
            { canonical: canonicalToCache?.substring(0, 25), isomeric: isomericToCache?.substring(0, 25) });
        }
      }

      console.log('%c🎨 Rendering with SmilesDrawer:', 'background: #4CAF50; color: white; font-weight: bold; padding: 4px;');
      console.log(`%c📊 SMILES: ${smiles}`, 'background: #FF9800; color: white; font-weight: bold; padding: 4px;');
      console.log(`%c📝 For compound: ${compoundName}`, 'color: #9c27b0;');

      // Get flags from molecule data
      const flags = moleculeData.flags || moleculeData.options || {};

      // Check if AI Flag Control is enabled
      const enableAIFlagControl = settings.enableAIFlagControl === true;
      console.log('%c🔧 AI Flag Control in renderClientSide:', 'color: #FF6B00; font-weight: bold;', enableAIFlagControl, 'raw value:', settings.enableAIFlagControl);

      // Determine rendering options based on flagsLocked and enableAIFlagControl:
      // - If enableAIFlagControl is FALSE: Ignore all flags, only use popup settings
      // - If flagsLocked && !useDefaults: ONLY use explicit flags (start from false/off)
      // - If flagsLocked && useDefaults: Start with popup settings, then apply flags
      // - If !flagsLocked: Start with popup settings, then apply flags
      const useOnlyFlags = enableAIFlagControl && flags.flagsLocked && !flags.useDefaults;
      const shouldApplyFlags = enableAIFlagControl && (flags.flagsLocked || Object.keys(flags).length > 0);

      let showCarbons, aromaticCircles, showMethyl, showAtomNumbers;
      let flipVertical, flipHorizontal, addHydrogens, showImplicitH, compactDrawing;
      let customSize = null;
      let customRotation = null;

      if (!enableAIFlagControl) {
        // AI Flag Control is disabled - ignore all flags, use only popup settings
        console.log('%c🚫 AI Flag Control DISABLED - using popup settings only', 'color: #FF9800; font-weight: bold;');
        showCarbons = settings.m2cfShowCarbons === true;
        aromaticCircles = settings.m2cfAromaticCircles === true;
        showMethyl = settings.m2cfShowMethyls === true;
        showAtomNumbers = settings.m2cfAtomNumbers === true;
        flipVertical = settings.m2cfFlipVertical === true;
        flipHorizontal = settings.m2cfFlipHorizontal === true;
        addHydrogens = settings.m2cfAddH2 === true;
        showImplicitH = settings.m2cfShowImplicitH !== false;
        compactDrawing = settings.m2cfCompactDrawing === true;
      } else if (useOnlyFlags) {
        // Flags are locked without +d: start with all OFF, only apply explicit flags
        console.log('%c🔒 Flags locked (no defaults) - using explicit flags only', 'color: #FF5722; font-weight: bold;');
        showCarbons = flags.showCarbons === true;
        aromaticCircles = flags.aromaticCircles === true;
        showMethyl = flags.showMethyls === true;
        showAtomNumbers = flags.atomNumbers === true;
        flipVertical = flags.flipVertical === true;
        flipHorizontal = flags.flipHorizontal === true;
        addHydrogens = flags.addHydrogens === true || flags.addH2 === true;
        showImplicitH = flags.showImplicitHydrogens !== false && flags.showImplicitH !== false;
        compactDrawing = settings.m2cfCompactDrawing === true;  // Keep compact from settings (not a flag)
      } else {
        // Use popup settings as base, then apply flags on top (if enableAIFlagControl is true)
        showCarbons = settings.m2cfShowCarbons === true;
        aromaticCircles = settings.m2cfAromaticCircles === true;
        showMethyl = settings.m2cfShowMethyls === true;
        showAtomNumbers = settings.m2cfAtomNumbers === true;
        flipVertical = settings.m2cfFlipVertical === true;
        flipHorizontal = settings.m2cfFlipHorizontal === true;
        addHydrogens = settings.m2cfAddH2 === true;
        showImplicitH = settings.m2cfShowImplicitH !== false;
        compactDrawing = settings.m2cfCompactDrawing === true;

        // Apply per-molecule flag overrides if present and AI Flag Control is enabled
        if (shouldApplyFlags) {
          if (flags.showCarbons !== undefined) showCarbons = flags.showCarbons;
          if (flags.aromaticCircles !== undefined) aromaticCircles = flags.aromaticCircles;
          if (flags.showMethyls !== undefined) showMethyl = flags.showMethyls;
          if (flags.atomNumbers !== undefined) showAtomNumbers = flags.atomNumbers;
          if (flags.addHydrogens !== undefined) addHydrogens = flags.addHydrogens;
          if (flags.addH2 !== undefined) addHydrogens = flags.addH2;
          if (flags.showImplicitHydrogens !== undefined) showImplicitH = flags.showImplicitHydrogens;
          if (flags.showImplicitH !== undefined) showImplicitH = flags.showImplicitH;
          if (flags.flipHorizontal !== undefined) flipHorizontal = flags.flipHorizontal;
          if (flags.flipVertical !== undefined) flipVertical = flags.flipVertical;
        }
      }

      // Size and rotation always from flags if specified (even when AI Flag Control is disabled for visual adjustments)
      if (flags.size !== undefined) customSize = flags.size;
      if (flags.rotation !== undefined) customRotation = flags.rotation;


      console.log('%c🎌 Flag overrides:', 'color: #FF6B6B;', flags);

      // Calculate intelligent size scaling based on molecule complexity
      let intelligentSizeMultiplier = 1;

      if (settings.sdScaleByWeight && smiles) {
        // Count heavy atoms (non-hydrogen) in SMILES as a proxy for molecule size
        // Remove stereochemistry markers, charges, and counts
        const cleanSmiles = smiles.replace(/[@\/\\]/g, '').replace(/\[[^\]]+\]/g, 'X');
        const heavyAtomMatches = cleanSmiles.match(/[A-Z]/g);
        const heavyAtomCount = heavyAtomMatches ? heavyAtomMatches.length : 20;

        // Logarithmic scaling: sqrt gives gentle curve
        // Small molecules (10 atoms): sqrt(10)/sqrt(20) = 0.71 → ~1.0x (stays small)
        // Medium molecules (50 atoms): sqrt(50)/sqrt(20) = 1.58 → ~1.6x
        // Large molecules (200 atoms): sqrt(200)/sqrt(20) = 3.16 → ~3.2x
        // This makes large molecules relatively bigger without making small ones tiny
        intelligentSizeMultiplier = Math.sqrt(heavyAtomCount) / Math.sqrt(20);

        // Clamp between 0.7x and 4x to avoid extremes
        intelligentSizeMultiplier = Math.max(0.7, Math.min(4, intelligentSizeMultiplier));

        console.log(`📐 Intelligent sizing: ${heavyAtomCount} atoms → ${intelligentSizeMultiplier.toFixed(2)}x scale`);
      }

      // Calculate dimensions (default, can be overridden by +s flag or intelligent scaling)
      const baseWidth = 300;
      const baseHeight = 240;
      const sizeMultiplier = customSize || intelligentSizeMultiplier;
      const renderWidth = Math.round(baseWidth * sizeMultiplier);
      const renderHeight = Math.round(baseHeight * sizeMultiplier);

      // Get theme from settings (user selection) or auto-detect
      const theme = settings.sdTheme || (isDarkModeEnabled() ? 'dark' : 'light');

      // Get gradient colors setting (solidBondColors: false = gradient, true = solid)
      const useGradientColors = settings.sdGradientColors === true;

      // Get bond thickness from settings (default 1.0)
      const bondThicknessMultiplier = settings.sdBondThickness || 1.0;

      // Configure SmilesDrawer options
      // These match YOUR CUSTOM smiles-drawer.min.js options:
      // showCarbons, showAromaticRings, showHydrogens, atomNumbering, solidBondColors, terminalCarbons
      const smilesDrawerOptions = {
        width: renderWidth,
        height: renderHeight,
        bondThickness: 2.0 * bondThicknessMultiplier,
        bondLength: 25,
        shortBondLength: 0.85,
        bondSpacing: 6,
        atomVisualization: 'default',
        isomeric: true,
        debug: false,
        // YOUR CUSTOM OPTIONS:
        showCarbons: showCarbons,               // Display C labels
        terminalCarbons: showMethyl,            // Show CH3 labels at terminal carbons
        explicitHydrogens: addHydrogens,        // Legacy option (may not be used)
        showHydrogens: addHydrogens,            // Show H atoms as separate nodes with bonds (branching)
        showImplicitHydrogens: showImplicitH,   // Show H counts in labels (H₂O, CH₃, OH, NH₂)
        showAromaticRings: aromaticCircles,     // Draw circles inside aromatic rings
        atomNumbering: showAtomNumbers,         // Number all atoms
        solidBondColors: !useGradientColors,    // false = gradient colors, true = solid
        compactDrawing: compactDrawing,           // Compact mode (linear text strings)
        overlapSensitivity: 0.42,
        overlapResolutionIterations: 1,
        fontSizeLarge: 11,
        fontSizeSmall: 3,
        padding: 10,
        experimentalSSSR: false,
        themes: {
          dark: {
            C: '#ffffff', O: '#ff6b6b', N: '#4dabf7', F: '#51cf66',
            CL: '#20c997', BR: '#fd7e14', I: '#be4bdb', P: '#fd7e14',
            S: '#fcc419', B: '#f59f00', SI: '#f59f00', H: '#aaaaaa',
            ATOM_NUMBER: '#ffd700',  // Yellow for dark theme (good contrast)
            BACKGROUND: 'transparent'
          },
          light: {
            C: '#222222', O: '#e74c3c', N: '#3498db', F: '#27ae60',
            CL: '#16a085', BR: '#d35400', I: '#8e44ad', P: '#d35400',
            S: '#f1c40f', B: '#e67e22', SI: '#e67e22', H: '#666666',
            ATOM_NUMBER: '#2196F3',  // Blue for light theme (good contrast)
            BACKGROUND: 'transparent'
          }
        }
      };

      console.log('%c⚙️ SmilesDrawer Options:', 'color: #9C27B0; font-weight: bold;', {
        showCarbons: smilesDrawerOptions.showCarbons,
        terminalCarbons: smilesDrawerOptions.terminalCarbons,
        showHydrogens: smilesDrawerOptions.showHydrogens,
        showAromaticRings: smilesDrawerOptions.showAromaticRings,
        atomNumbering: smilesDrawerOptions.atomNumbering,
        solidBondColors: smilesDrawerOptions.solidBondColors,
        compactDrawing: smilesDrawerOptions.compactDrawing
      });

      console.log('%c🎨 Rendering SVG directly with SmilesDrawer', 'color: #9c88ff;');

      // Use queued parsing to limit concurrent operations (prevents crashes with many molecules)
      // Also reuse SvgDrawer instance to reduce memory usage
      queueSmilesDrawerParse(smiles, function (tree) {
        // Track tempDiv for cleanup - MUST always be cleaned up even on errors
        let tempDiv = null;
        let tempId = null;

        try {
          // Count atoms to determine molecule complexity
          const atomCount = smiles.replace(/[^A-Z]/g, '').length; // Rough estimate
          console.log(`Molecule has approximately ${atomCount} atoms`);

          // Adjust canvas size based on complexity
          let adjustedWidth = renderWidth;
          let adjustedHeight = renderHeight;

          if (atomCount > 100) {
            // Large molecule like insulin - use bigger canvas
            // Use square root for sub-linear scaling (parabolic curve - rate of increase decreases)
            // This prevents molecules from becoming excessively large
            const complexityFactor = Math.min(Math.sqrt(atomCount / 50) * 1.5, 3); // Max 3x with sqrt scaling
            adjustedWidth = Math.ceil(renderWidth * complexityFactor);
            adjustedHeight = Math.ceil(renderHeight * complexityFactor);
            console.log(`Large molecule detected, expanding canvas to ${adjustedWidth}x${adjustedHeight} (factor: ${complexityFactor.toFixed(2)}x)`);

            // Update options
            smilesDrawerOptions.width = adjustedWidth;
            smilesDrawerOptions.height = adjustedHeight;
          }

          // Safety check - ensure document.body exists
          if (!document.body) {
            console.error('[ChemRenderer] document.body not available');
            decrementActiveLoads();
            return;
          }

          // Draw to a temporary container first (offscreen to avoid layout thrashing)
          tempId = `temp-${Date.now()}-${Math.random().toString(36).substr(2, 9)}`;
          tempDiv = document.createElement('div');
          tempDiv.style.cssText = 'position:absolute;left:-9999px;top:-9999px;visibility:hidden;';
          tempDiv.innerHTML = `<svg id="${tempId}"></svg>`;
          document.body.appendChild(tempDiv);

          // Reuse global SvgDrawer instance (memory optimization)
          const finalDrawer = getOrCreateSvgDrawer(smilesDrawerOptions);
          finalDrawer.draw(tree, tempId, theme);

          // Get the rendered SVG
          const renderedSvg = document.getElementById(tempId);
          if (renderedSvg) {
            // Get the actual bounding box of the content
            try {
              const bbox = renderedSvg.getBBox();
              const padding = 10; // Small padding around content
              let actualWidth = Math.ceil(bbox.width + padding * 2);
              let actualHeight = Math.ceil(bbox.height + padding * 2);

              // Smart scaling for display
              const maxDisplayWidth = 600;   // Max display width
              const maxDisplayHeight = 500;  // Max display height (matches hard limits)
              const minDisplayWidth = 180;   // Min display width
              const minDisplayHeight = 120;  // Min display height

              let displayWidth = actualWidth;
              let displayHeight = actualHeight;

              // Scale down if too large (insulin case)
              if (displayWidth > maxDisplayWidth || displayHeight > maxDisplayHeight) {
                const scale = Math.min(maxDisplayWidth / displayWidth, maxDisplayHeight / displayHeight);
                displayWidth = Math.ceil(displayWidth * scale);
                displayHeight = Math.ceil(displayHeight * scale);
                console.log(`Scaling down to ${displayWidth}x${displayHeight} for display`);
              }

              // Scale up if too small (benzene case)
              if (displayWidth < minDisplayWidth && displayHeight < minDisplayHeight) {
                const scale = Math.min(minDisplayWidth / displayWidth, minDisplayHeight / displayHeight, 2.5);
                displayWidth = Math.ceil(displayWidth * scale);
                displayHeight = Math.ceil(displayHeight * scale);
                console.log(`Scaling up to ${displayWidth}x${displayHeight} for visibility`);
              }

              actualWidth = displayWidth;
              actualHeight = displayHeight;

              // Update SVG to fit content with smart scaling
              renderedSvg.setAttribute('width', actualWidth);
              renderedSvg.setAttribute('height', actualHeight);
              renderedSvg.setAttribute('viewBox', `${bbox.x - padding} ${bbox.y - padding} ${bbox.width + padding * 2} ${bbox.height + padding * 2}`);

              console.log('Final SVG size:', actualWidth, 'x', actualHeight);
            } catch (e) {
              console.warn('Could not get bbox:', e);
            }

            // Convert SVG to string and create data URL
            const svgString = new XMLSerializer().serializeToString(renderedSvg);

            // Use URL encoding instead of base64 to avoid UTF-8 issues
            const svgDataUrl = 'data:image/svg+xml;charset=utf-8,' + encodeURIComponent(svgString);

            // Note: tempDiv cleanup is handled in finally block to ensure it always happens

            // Create img element just like moleculeviewer/mol2chemfig does
            const svgImg = document.createElement('img');
            svgImg.src = svgDataUrl;
            svgImg.alt = compoundName;
            svgImg.className = 'molecule-diagram';

            // Store data attributes for instant re-rendering
            svgImg.dataset.moleculeViewer = 'true';
            svgImg.dataset.smiles = smiles;
            svgImg.dataset.nomenclature = compoundName;
            svgImg.dataset.renderWidth = String(renderWidth);
            svgImg.dataset.renderHeight = String(renderHeight);
            svgImg.dataset.compoundType = moleculeData.flags?.compoundType || 'compound';

            // Style to ensure no extra space
            svgImg.style.display = 'block';
            svgImg.style.margin = '0';
            svgImg.style.padding = '0';
            svgImg.style.border = 'none';

            // Apply transforms if needed
            if (flipVertical || flipHorizontal || customRotation) {
              let transformParts = [];
              if (flipVertical) transformParts.push('scaleY(-1)');
              if (flipHorizontal) transformParts.push('scaleX(-1)');
              if (customRotation) transformParts.push(`rotate(${customRotation}deg)`);
              svgImg.style.transform = transformParts.join(' ');
            }

            // Use the same wrapper function as mol2chemfig and moleculeviewer!
            const wrapper = wrapImageWithRotationContainer(svgImg, customRotation, compoundName, moleculeData);

            // Replace original img element
            const parent = img.parentNode;
            if (parent) {
              parent.replaceChild(wrapper, img);
              // console.log('%c✅ SVG rendered as img with standard wrapper', 'color: green; font-weight: bold;');
            }
          }
        } catch (renderError) {
          console.error('%c❌ Error during SVG rendering:', 'color: red;', renderError);
        } finally {
          // CRITICAL: Always clean up tempDiv from document.body to prevent memory leaks
          if (tempDiv && tempDiv.parentNode) {
            try {
              tempDiv.parentNode.removeChild(tempDiv);
            } catch (cleanupError) {
              // Ignore cleanup errors
            }
          }
        }
      }, function (err) {
        console.error('%c❌ SmilesDrawer parse error:', 'color: red;', err);
        console.error('SMILES that failed:', smiles);
        console.error('Error details:', err.message, err.stack);

        // Create error message
        const errorImg = document.createElement('img');
        errorImg.src = 'data:image/svg+xml;base64,PHN2ZyB3aWR0aD0iMzAwIiBoZWlnaHQ9IjgwIiB4bWxucz0iaHR0cDovL3d3dy53My5vcmcvMjAwMC9zdmciPjxyZWN0IHdpZHRoPSIzMDAiIGhlaWdodD0iODAiIGZpbGw9IiNmZmYiLz48dGV4dCB4PSIxMCIgeT0iMzAiIGZpbGw9IiNmMDAiIGZvbnQtZmFtaWx5PSJtb25vc3BhY2UiIGZvbnQtc2l6ZT0iMTIiPkVycm9yOiBGYWlsZWQgdG8gcGFyc2UgU01JTEVTPC90ZXh0Pjx0ZXh0IHg9IjEwIiB5PSI1MCIgZmlsbD0iIzY2NiIgZm9udC1mYW1pbHk9Im1vbm9zcGFjZSIgZm9udC1zaXplPSIxMCI+Q29tcG91bmQ6ICcgKyBjb21wb3VuZE5hbWUgKyAnPC90ZXh0Pjwvc3ZnPg==';
        errorImg.alt = `Error: ${err.message || 'Failed to parse SMILES'}`;
        errorImg.className = 'molecule-diagram';

        try {
          const parent = img.parentNode;
          if (parent) {
            parent.replaceChild(errorImg, img);
          }
        } catch (replaceError) {
          console.error('%c❌ Error replacing element:', 'color: red;', replaceError);
        }
      });

      decrementActiveLoads();
      return;
    } catch (error) {
      console.error('%c❌ Client-side rendering error:', 'color: red; font-weight: bold;', error);
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
      // Get viewer size from settings
      const viewerSize = settings.viewer3DSize || 'normal';

      // Define size dimensions
      const sizeDimensions = {
        'small': { width: 200, height: 150 },
        'medium': { width: 300, height: 250 },
        'normal': { width: 400, height: 350 },
        'large': { width: 600, height: 450 },
        'xlarge': { width: 800, height: 600 }
      };

      const dimensions = sizeDimensions[viewerSize] || sizeDimensions['normal'];
      console.log('%c📐 3D Viewer Size:', 'color: #ff6b6b; font-weight: bold;', viewerSize, dimensions);

      // Create container for 3D viewer with appropriate sizing (increased for better WebGL rendering)
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
          // Protein/biomolecule - check if bioassembly is enabled
          if (settings.molviewBioAssembly) {
            // Build URL based on user's selected bio assembly viewer
            const pdbId = moleculeData.pdbid.toLowerCase();
            const viewer = settings.bioAssemblyViewer || 'molstar';
            const graphics = settings.bioAssemblyGraphics || 'balanced';
            const provider = settings.bioAssemblyPdbProvider || 'rcsb';

            switch (viewer) {
              case 'molstar-me':
                // Mol* Mesoscale Explorer - for very large structures, simplified representation
                viewerUrl = `https://molstar.org/me/viewer/?pdb=${pdbId}&hide-controls=1&graphics-mode=${graphics}&pdb-provider=${provider}`;
                console.log('%c🧬 Using Mol* Mesoscale Explorer:', 'color: #FF5722; font-weight: bold;', viewerUrl);
                break;

              case 'molview':
                // MolView - falls back to asymmetric unit (no true bio assembly support)
                viewerUrl = buildMolViewEmbedUrl(moleculeData.pdbid, 'pdbid', 'biomolecule');
                console.log('%c🧬 Using MolView (asymmetric unit):', 'color: #9C27B0; font-weight: bold;', viewerUrl);
                break;

              case 'molstar':
              default:
                // Mol* Standard Viewer - best chain coloring, loads bio assembly
                viewerUrl = `https://molstar.org/viewer/?pdb=${pdbId}&assembly-id=1&hide-controls=1`;
                console.log('%c🧬 Using Mol* Standard Viewer:', 'color: #E91E63; font-weight: bold;', viewerUrl);
                break;
            }
          } else {
            // Use MolView embed for asymmetric unit (default behavior)
            // MolView is our default - Mol* is only used if MolView fails (via timeout fallback)
            viewerUrl = buildMolViewEmbedUrl(moleculeData.pdbid, 'pdbid', 'biomolecule');
            console.log('%c🧬 Built MolView URL for protein (asymmetric unit):', 'color: #9C27B0; font-weight: bold;', viewerUrl);
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
        viewer3DIframe.setAttribute('sandbox', 'allow-scripts allow-same-origin allow-popups allow-forms');
      } else {
        // For compounds, build URL based on user's selected source
        let viewerSource = settings.viewer3DSource || '3dmol';
        // console.log('%c🔍 Current 3D viewer source setting:', 'color: #ff6b6b; font-weight: bold;', viewerSource);
        // console.log('%c📦 All settings:', 'color: #4ecdc4;', settings);

        switch (viewerSource) {
          case 'molview':
            // MolView Embed - use SMILES directly to avoid 404 errors on molecules without 3D SDF data
            // Get SMILES from moleculeData, img dataset, or lookup
            let smilesMolview = moleculeData.smiles || (img && img.dataset.smiles);

            // If no SMILES in data, try to look it up
            if (!smilesMolview && compoundName) {
              console.log('%c🔍 Looking up SMILES for MolView...', 'color: #0066cc;');
              const bridgeResult = await smilesBridge(compoundName);
              if (bridgeResult && bridgeResult.smiles) {
                smilesMolview = bridgeResult.smiles;
              }
            }

            if (smilesMolview) {
              // Use buildMolViewEmbedUrl to get consistent URL with settings
              viewerUrl = buildMolViewEmbedUrl(smilesMolview, 'smiles', 'compound');
              console.log('%c📍 MolView Embed URL (SMILES):', 'color: #0066cc; font-weight: bold;', viewerUrl);
            } else {
              console.warn('%c⚠️ SMILES not found for MolView, using 3Dmol.js fallback', 'color: orange;');
              viewerSource = '3dmol';
            }
            break;

          case 'molstar': {
            // Use new MolView.com Mol* renderer (full app) via iframe
            const cidMolstar = moleculeData.cid || await getPubChemCID(compoundName).catch(() => null);
            const pdbMolstar = moleculeData.pdbid;
            const codMolstar = moleculeData.codid;

            const loadParam = pdbMolstar ? `pdb:${pdbMolstar}` : codMolstar ? `cod:${codMolstar}` : cidMolstar ? `cid:${cidMolstar}` : compoundName;
            viewerUrl = `https://app.molview.com/?load=${encodeURIComponent(loadParam)}`;
            viewer3DIframe.setAttribute('allow', 'accelerometer; autoplay; clipboard-write; encrypted-media; gyroscope; picture-in-picture; web-share; fullscreen; xr-spatial-tracking');
            viewer3DIframe.setAttribute('sandbox', 'allow-scripts allow-same-origin allow-popups allow-forms');
            console.log('%c📍 Mol* (MolView.com) URL:', 'color: #4caf50; font-weight: bold;', viewerUrl);
            break;
          }

          case 'pubchem':
            // PubChem direct embed is unreliable due to cross-origin restrictions and lack of a clean widget URL.
            // We fallback to 3Dmol.js which fetches data FROM PubChem but renders it locally and cleanly.
            console.log('%c⚠️ PubChem direct embed is deprecated. Using 3Dmol.js (local renderer) instead.', 'color: orange; font-weight: bold;');
          // Fallthrough to 3dmol

          case '3dmol':
          default:
            // 3Dmol.js - client-side viewer with custom styling
            let cid3dmol = null;
            try {
              cid3dmol = await getPubChemCID(compoundName);
            } catch (e) {
              console.warn('CID lookup failed, proceeding with name only:', e);
            }

            // Use our custom 3Dmol.js viewer with user's style preferences
            const style3D = settings.viewer3DStyle || 'stick:sphere';

            // Get per-style settings for stick and sphere radius
            const styleSettings = settings.viewer3DStyleSettings || {};
            const currentStyleSettings = styleSettings[style3D] || {};
            const stickRadius = currentStyleSettings.stickRadius || '0.15';
            const sphereRadius = currentStyleSettings.sphereRadius || '0.3';
            const autoRotate = settings.viewer3DAutoRotate !== false;
            const bgColor = encodeURIComponent(settings.viewer3DBgColor || '#1a1a2e');

            viewerUrl = chrome.runtime.getURL(
              `3dmol-viewer.html?name=${encodeURIComponent(compoundName)}&style=${style3D}&stickRadius=${stickRadius}&sphereRadius=${sphereRadius}&autoRotate=${autoRotate}&bgColor=${bgColor}`
            );

            if (cid3dmol) {
              viewerUrl += `&cid=${cid3dmol}`;
            }

            console.log('%c📍 3Dmol.js Custom Viewer URL:', 'color: #0066cc; font-weight: bold;', viewerUrl);
            console.log('%c🎨 Style settings:', 'color: #9c88ff;', { style: style3D, stickRadius, sphereRadius, autoRotate });
            break;
        }
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

        // Set a timeout - if MolView doesn't seem to work within 8 seconds, fallback to Mol*
        const molviewTimeout = setTimeout(switchToMolstar, 8000);

        // Listen for MolView error messages (postMessage)
        const errorHandler = (event) => {
          // Check if message is from MolView and indicates an error
          if (event.origin === 'https://embed.molview.org' && !fallbackTriggered) {
            const data = event.data;
            if (data && (data.error || data.type === 'error' ||
              (typeof data === 'string' && data.toLowerCase().includes('error')) ||
              (typeof data === 'string' && data.toLowerCase().includes('failed')))) {
              clearTimeout(molviewTimeout);
              console.log('%c❌ MolView error detected', 'color: red; font-weight: bold;', data);
              switchToMolstar();
              window.removeEventListener('message', errorHandler);
            }
          }
        };
        window.addEventListener('message', errorHandler);

        // NOTE: We do NOT clear the timeout on 'load' event anymore
        // because MolView fires 'load' even when it fails to load the structure
        // The timeout will run for the full 8 seconds unless an error message is detected

        // Clean up error listener after timeout period
        setTimeout(() => {
          window.removeEventListener('message', errorHandler);
        }, 10000); // Clean up after 10 seconds
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
      // Fallback to regular image - respect pubchemDirectFetch setting
      let fallbackUrl;
      if (settings.pubchemDirectFetch) {
        // Try direct fetch
        const cid = await getPubChemCID(compoundName);
        if (cid) {
          fallbackUrl = `https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/${cid}/PNG?image_size=large&record_type=2d`;
        } else {
          fallbackUrl = `${PUBCHEM_API}/pubchem/img/${encodeURIComponent(compoundName)}?size=large`;
        }
      } else {
        fallbackUrl = `${PUBCHEM_API}/pubchem/img/${encodeURIComponent(compoundName)}?size=large`;
      }

      img.src = fallbackUrl;
      img.classList.add('molecule-fadein');
      img.classList.remove('molecule-loading');
      decrementActiveLoads();
    }
  }

  // Duplicate getPubChemCID removed - using the robust version defined at the top of the file


  // Helper function to apply sharpening filter using convolution kernel
  function applySharpenFilter(ctx, width, height) {
    const imageData = ctx.getImageData(0, 0, width, height);
    const data = imageData.data;
    const outputData = ctx.createImageData(width, height);
    const output = outputData.data;

    // Sharpening kernel (3x3 convolution matrix)
    // This enhances edges and makes bond lines crisper
    const kernel = [
      0, -1, 0,
      -1, 5, -1,
      0, -1, 0
    ];

    // Apply convolution
    for (let y = 1; y < height - 1; y++) {
      for (let x = 1; x < width - 1; x++) {
        const idx = (y * width + x) * 4;

        // Skip transparent pixels (don't sharpen background)
        if (data[idx + 3] === 0) {
          output[idx] = data[idx];
          output[idx + 1] = data[idx + 1];
          output[idx + 2] = data[idx + 2];
          output[idx + 3] = data[idx + 3];
          continue;
        }

        // Apply kernel to RGB channels
        let r = 0, g = 0, b = 0;
        for (let ky = -1; ky <= 1; ky++) {
          for (let kx = -1; kx <= 1; kx++) {
            const pixelIdx = ((y + ky) * width + (x + kx)) * 4;
            const kernelIdx = (ky + 1) * 3 + (kx + 1);
            const weight = kernel[kernelIdx];

            r += data[pixelIdx] * weight;
            g += data[pixelIdx + 1] * weight;
            b += data[pixelIdx + 2] * weight;
          }
        }

        // Clamp values to 0-255 range
        output[idx] = Math.max(0, Math.min(255, r));
        output[idx + 1] = Math.max(0, Math.min(255, g));
        output[idx + 2] = Math.max(0, Math.min(255, b));
        output[idx + 3] = data[idx + 3]; // Preserve alpha
      }
    }

    // Copy edges without sharpening to avoid artifacts
    for (let x = 0; x < width; x++) {
      // Top edge
      const topIdx = x * 4;
      output[topIdx] = data[topIdx];
      output[topIdx + 1] = data[topIdx + 1];
      output[topIdx + 2] = data[topIdx + 2];
      output[topIdx + 3] = data[topIdx + 3];

      // Bottom edge
      const bottomIdx = ((height - 1) * width + x) * 4;
      output[bottomIdx] = data[bottomIdx];
      output[bottomIdx + 1] = data[bottomIdx + 1];
      output[bottomIdx + 2] = data[bottomIdx + 2];
      output[bottomIdx + 3] = data[bottomIdx + 3];
    }

    for (let y = 0; y < height; y++) {
      // Left edge
      const leftIdx = y * width * 4;
      output[leftIdx] = data[leftIdx];
      output[leftIdx + 1] = data[leftIdx + 1];
      output[leftIdx + 2] = data[leftIdx + 2];
      output[leftIdx + 3] = data[leftIdx + 3];

      // Right edge
      const rightIdx = (y * width + width - 1) * 4;
      output[rightIdx] = data[rightIdx];
      output[rightIdx + 1] = data[rightIdx + 1];
      output[rightIdx + 2] = data[rightIdx + 2];
      output[rightIdx + 3] = data[rightIdx + 3];
    }

    return outputData;
  }

  // Helper function to auto-crop canvas to molecule bounds
  function autoCropCanvas(canvas) {
    const ctx = canvas.getContext('2d');
    const imageData = ctx.getImageData(0, 0, canvas.width, canvas.height);
    const data = imageData.data;

    // Find bounding box of non-transparent pixels
    let minX = canvas.width;
    let minY = canvas.height;
    let maxX = 0;
    let maxY = 0;

    for (let y = 0; y < canvas.height; y++) {
      for (let x = 0; x < canvas.width; x++) {
        const idx = (y * canvas.width + x) * 4;
        const alpha = data[idx + 3];

        // If pixel is not transparent
        if (alpha > 0) {
          if (x < minX) minX = x;
          if (x > maxX) maxX = x;
          if (y < minY) minY = y;
          if (y > maxY) maxY = y;
        }
      }
    }

    // Add padding around the molecule (10px on each side)
    const padding = 10;
    minX = Math.max(0, minX - padding);
    minY = Math.max(0, minY - padding);
    maxX = Math.min(canvas.width - 1, maxX + padding);
    maxY = Math.min(canvas.height - 1, maxY + padding);

    // Calculate cropped dimensions
    const croppedWidth = maxX - minX + 1;
    const croppedHeight = maxY - minY + 1;

    // If the cropped area is too small or invalid, return original canvas
    if (croppedWidth <= 0 || croppedHeight <= 0 || croppedWidth > canvas.width || croppedHeight > canvas.height) {
      return canvas;
    }

    // Create new canvas with cropped dimensions
    const croppedCanvas = document.createElement('canvas');
    croppedCanvas.width = croppedWidth;
    croppedCanvas.height = croppedHeight;
    const croppedCtx = croppedCanvas.getContext('2d');

    // Copy the cropped region to new canvas
    croppedCtx.drawImage(canvas, minX, minY, croppedWidth, croppedHeight, 0, 0, croppedWidth, croppedHeight);

    return croppedCanvas;
  }


  // Helper function to load molecule image - ALWAYS uses client-side SmilesDrawer
  // Uses IntegratedSearch for compound lookup (no server needed!)
  function loadMoleculeImage(img) {
    console.log('%c🎨 [loadMoleculeImage] ALWAYS using CLIENT-SIDE SmilesDrawer + IntegratedSearch', 'background: #9C27B0; color: white; font-weight: bold; padding: 4px;');

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
        // Check if this is a chemistry image (MoleculeViewer, PubChem, or Mol2chemfig)
        if (img.classList.contains('molecule-viewer') || img.classList.contains('molecule-pubchem') || img.classList.contains('molecule-legacy')) {
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

    const textNodes = [];
    let combinedText = '';

    // Collect all text nodes and build combined text
    for (let j = 0; j < element.childNodes.length; j++) {
      const node = element.childNodes[j];
      if (node.nodeType === 3) { // Text node
        textNodes.push(node);
        combinedText += node.nodeValue;
      }
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

  const chemMatches = result.match(/\bchem:([^:]+):/g);

  console.log('[ChemRenderer] 🔍 DEBUG: chemMatches result:', chemMatches);

  if (chemMatches) {
    // SAFETY: Limit number of matches to prevent excessive processing
    if (chemMatches.length > 100) {
      console.warn('[ChemRenderer] Too many chem: patterns detected, limiting to first 100');
      chemMatches.length = 100;
    }

    log.debug(`  Found ${chemMatches.length} chem:text: patterns`);
    console.log('[ChemRenderer] ✅ Found patterns:', chemMatches);
    log.debug(`  Renderer engine: ${settings.rendererEngine}, 3D Viewer enabled: ${settings.enable3DViewer}`);

    result = result.replace(/\bchem:([^:]+):/g, (match, content) => {
      try {
        const content_trimmed = content.trim();
        if (!content_trimmed || content_trimmed.length > 500) {
          return match; // Empty or too long, skip
        }

        log.debug(`  Converting molecule: ${content_trimmed}`);
        window.chemRendererPerformance.recordStructure();

        // Extract compound name (before any flags)
        let compoundName = content_trimmed;


        // Process flags if any are detected (always enabled by default)
        let flagOverrides = { useDefaults: false }; // Default empty flags
        let currentSettings = settings; // Default to original settings

        // Check if there are flags to process OR if using new syntax (type=value)
        // New syntax: chem:smiles=CCO:, chem:mol=benzene:, chem:biomol=insulin:, chem:mineral=quartz:
        const hasNewSyntax = content_trimmed.includes('=');

        // Flag detection: Only treat +/- as flags if they appear AFTER the compound name
        // e.g., "benzene+c" has flags, but "(R)-butan-2-ol" does NOT have flags
        // Flags should be: /flags or +flag at end or -flag at end
        // A hyphen INSIDE a name like "butan-2-ol" or "(R)-" is NOT a flag
        const hasSlashFlags = content_trimmed.includes('/');
        // Check for trailing flags: must have + or - followed by single letters at the end
        // e.g., "benzene+c+n" or "histamine-c"
        const hasTrailingFlags = /[+\-][a-zA-Z](?:[+\-][a-zA-Z])*$/.test(content_trimmed) &&
          !content_trimmed.match(/^[a-zA-Z0-9()\[\]-]+$/); // Exclude pure compound names with internal hyphens
        const hasFlags = hasSlashFlags || hasTrailingFlags || hasNewSyntax;

        if (hasNewSyntax || hasFlags) {
          try {
            // Check if AI Flag Control is enabled
            const enableAIFlagControl = settings.enableAIFlagControl === true;
            console.log('%c🔧 AI Flag Control setting:', 'color: #FF6B00; font-weight: bold;', enableAIFlagControl, 'raw value:', settings.enableAIFlagControl);
            
            flagOverrides = window.parseChemFlags('chem:' + content_trimmed + ':', enableAIFlagControl);
            console.log('%c🏴 Parsed flags:', 'color: #FF6B00; font-weight: bold;', flagOverrides);

            // Use moleculeName from new syntax if available
            if (flagOverrides.moleculeName) {
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
        // 1. chem:smiles=CCO: OR +smiles → Direct SMILES, send to SmilesDrawer
        // 2. chem:biomol=insulin: OR +biomolecule/+protein → Search RCSB PDB
        // 3. chem:mineral=quartz: OR +mineral → Search COD
        // 4. chem:mol=benzene: OR +compound → Search PubChem for compounds
        // 5. No type specified → Use IntegratedSearch for auto-detection

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
            // Direct PubChem CID lookup for compounds
            console.log('%c🧪 DIRECT CID LOOKUP:', 'background: #4CAF50; color: white; font-weight: bold; padding: 4px;', idValue);
            const pubchemData = {
              cid: idValue,
              name: displayName,
              type: 'cid',
              isPubChem: true,
              directLookup: true,
              flags: {
                ...flagsForRender,
                compoundType: 'compound',
                skipIntegratedSearch: true
              }
            };
            const encoded = btoa(JSON.stringify(pubchemData));
            const converted = `<img src="" alt="chem" class="molecule-diagram molecule-pubchem" data-molecule-viewer="${encoded}" data-loaded="false" style="display: inline-block; height: auto; margin: 0 12px 8px 0; vertical-align: middle; padding: 4px; border-radius: 4px; background-color: transparent;">`;
            log.debug(`  📤 Sending CID to PubChem: ${idValue}`);
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
        // +3d flag: Use 3D viewer (PubChem-based)
        if (flagOverrides.is3D) {
          const pubchemData = {
            // Use actualSmiles if available, otherwise use compoundName for lookup
            ...(actualSmiles ? { smiles: actualSmiles } : { nomenclature: compoundName }),
            name: compoundName, // Display name for the tag
            type: actualSmiles ? 'smiles' : 'nomenclature',
            isPubChem: true,
            show3D: true,  // Signal to show 3D viewer immediately
            auto3D: true,  // Hint to viewer to auto-select best 3D mode
            flags: flagsForRender
          };

          const encoded = btoa(JSON.stringify(pubchemData));
          const converted = `<img src="" alt="chem" class="molecule-diagram molecule-pubchem" data-molecule-viewer="${encoded}" data-loaded="false" data-show-3d="true" style="display: inline-block; height: auto; margin: 0 12px 8px 0; vertical-align: middle; padding: 4px; border-radius: 4px; background-color: transparent;">`;

          log.debug(`  📤 Sending ${actualSmiles ? 'SMILES' : 'nomenclature'} to 3D viewer: ${compoundName}`);
          return converted;
        }

        // mol= syntax: Force PubChem renderer - Direct lookup, no IntegratedSearch!
        // This is triggered by chem:mol=benzene: which sets isPubchem=true and compoundType='compound'
        if (flagOverrides.isPubchem || specifiedType === 'compound') {
          console.log('%c🧪 DIRECT PUBCHEM PATH (mol= syntax)', 'background: #4CAF50; color: white; font-weight: bold; padding: 4px;');
          const pubchemData = {
            ...(actualSmiles ? { smiles: actualSmiles } : { nomenclature: compoundName }),
            name: compoundName,
            type: actualSmiles ? 'smiles' : 'nomenclature',
            isPubChem: true,
            directLookup: true,  // Signal to skip IntegratedSearch
            flags: {
              ...flagsForRender,
              compoundType: 'compound',  // Ensure compound type is set
              skipIntegratedSearch: true  // Ensure IntegratedSearch is skipped
            }
          };

          const encoded = btoa(JSON.stringify(pubchemData));
          const converted = `<img src="" alt="chem" class="molecule-diagram molecule-pubchem" data-molecule-viewer="${encoded}" data-loaded="false" style="display: inline-block; height: auto; margin: 0 12px 8px 0; vertical-align: middle; padding: 4px; border-radius: 4px; background-color: transparent;">`;

          log.debug(`  📤 Sending ${actualSmiles ? 'SMILES' : 'nomenclature'} DIRECTLY to PubChem (mol= syntax): ${compoundName}`);
          return converted;
        }

        // Check if using PubChem renderer (which supports 3D viewer)
        if (currentSettings.rendererEngine === 'pubchem') {
          // Use PubChem - supports 3D viewer when enabled
          const pubchemData = {
            ...(actualSmiles ? { smiles: actualSmiles } : { nomenclature: compoundName }),
            name: compoundName,
            type: actualSmiles ? 'smiles' : 'nomenclature',
            isPubChem: true,
            flags: flagsForRender
          };

          const encoded = btoa(JSON.stringify(pubchemData));
          const converted = `<img src="" alt="chem" class="molecule-diagram molecule-pubchem" data-molecule-viewer="${encoded}" data-loaded="false" style="display: inline-block; height: auto; margin: 0 12px 8px 0; vertical-align: middle; padding: 4px; border-radius: 4px; background-color: transparent;">`;

          log.debug(`  📤 Sending ${actualSmiles ? 'SMILES' : 'nomenclature'} to PubChem: ${compoundName}`);
          return converted;
        }

        // Default: Send to MoleculeViewer (uses IntegratedSearch for type detection)
        const rotation = currentSettings.m2cfRotate || currentSettings.mvRotate || 0;
        const scale = flagOverrides.size || 1.0;

        const moleculeViewerData = {
          // If actualSmiles available (from named syntax or direct): use smiles property
          // Otherwise: use nomenclature and let IntegratedSearch detect the type
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


  // Pattern 1: \ce{H2O} → convert to Unicode subscripts (avoids ChatGPT's LaTeX conflict)
  if (settings.renderMhchem) {
    log.debug('🧪 Applying Pattern 1: \\ce{...} for mhchem');
    const matches1 = result.match(/\\ce\{([^}]+)\}/g);
    if (matches1) {
      log.debug(`  Found ${matches1.length} \\ce{} patterns`);
      patternMatches.push(...matches1);

      result = result.replace(/\\ce\{([^}]+)\}/g, (match, content) => {
        log.debug(`  Converting mhchem: ${match}`);
        window.chemRendererPerformance.recordFormula();
        const converted = convertChemistry(content);
        log.debug(`  To Unicode: ${converted}`);
        return converted;
      });
    }

    // Pattern 2: ce{Na+} (without backslash) → Unicode conversion
    log.debug('🧪 Applying Pattern 2: ce{...} (no backslash)');
    const matches2 = text.match(/\bce\{([^}]+)\}/g);
    if (matches2) {
      log.debug(`  Found ${matches2.length} ce{} patterns`);
      patternMatches.push(...matches2);

      result = result.replace(/\bce\{([^}]+)\}/g, (match, content) => {
        log.debug(`  Converting mhchem (no backslash): ${match}`);
        window.chemRendererPerformance.recordFormula();
        const converted = convertChemistry(content);
        log.debug(`  To Unicode: ${converted}`);
        return converted;
      });
    }
  }

  // Pattern 5: Handle reaction arrows with conditions: →[Na/ether]
  // This makes reaction arrows more visible/formatted
  log.debug('🧪 Applying Pattern 5: →[...] reaction conditions');
  const matches5 = text.match(/→\[[^\]]+\]/g);
  if (matches5) {
    log.debug(`  Found ${matches5.length} reaction condition patterns`);
    patternMatches.push(...matches5);
    // Just highlight these - convert to bold/styled version
    result = result.replace(/→\[([^\]]+)\]/g, (match, condition) => {
      log.debug(`  Reaction condition found: ${match}`);
      // Keep the arrow and condition visible (users see: → [condition]
      return '→ [<strong>' + condition + '</strong>]';
    });
  }

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
        target.classList.contains('chemtex-image-loading-indicator')
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

    // For PNG images (like PubChem images), we can't easily change colors like we do with SVGs
    // But we can trigger a re-render if needed by reloading the images
    const pubchemImages = document.querySelectorAll('img.pubchem-diagram');
    pubchemImages.forEach(img => {
      // For dark mode, if background was removed, we might want to re-process the image
      // to potentially adjust how it appears in dark mode
      // For now, we'll just log that dark mode changed
      console.log(`🌓 Dark mode ${isDark ? 'enabled' : 'disabled'} - consider reprocessing PubChem image:`, img.src);
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

log.success('🎉 Content script fully loaded and ready!');
log.info('✨ Chemistry formulas will be rendered via Unicode and CodeCogs API');
log.info('ℹ️  Run window.chemRendererDebug.getLogs() to see all logs');
log.info('ℹ️  Run window.chemRendererDebug.testFormulas() to test formula detection');
log.info('ℹ️  Run window.chemRendererDebug.scanPage() to manually trigger a page scan');

// Check IntegratedSearch status
if (window.IntegratedSearch) {
  log.success('✅ IntegratedSearch module loaded');
} else {
  log.error('❌ IntegratedSearch module NOT loaded! Check integrated-search.js');
}
// ============================================
// CONTEXT MENU HANDLER
// ============================================
chrome.runtime.onMessage.addListener((request, sender, sendResponse) => {
  // Handler for "Render as Molecule" (treats text as chemical compound name)
  // Goes directly to PubChem for CID/SMILES lookup (no IntegratedSearch auto-detection)
  if (request.type === 'INSPECT_MOLECULE') {
    const text = request.text?.trim();
    if (!text) return;

    console.log('%c🧪 Rendering selection as Molecule (PubChem):', 'color: #2196F3; font-weight: bold;', text);

    // Get the current selection to find where to insert the image
    const selection = window.getSelection();
    if (!selection.rangeCount) return;

    const range = selection.getRangeAt(0);

    // Create image element (same structure as standard renderer)
    const img = document.createElement('img');
    img.className = 'molecule-diagram molecule-viewer';
    img.alt = text;
    img.title = `Searching PubChem: ${text}`;
    img.src = ''; // Empty src, will be filled by loader

    // Create molecule data - forces compound type (PubChem lookup)
    // compoundType='compound' tells renderClientSide to use PubChem directly
    const moleculeData = {
      nomenclature: text,
      type: 'nomenclature',
      options: {
        width: 400,
        height: 300,
        scale: 1,
        rotation: 0
      },
      isMoleculeViewer: true,
      flags: {
        compoundType: 'compound',  // Force PubChem lookup (not auto-detect)
        isDirectSmiles: false
      }
    };

    // Set dataset
    img.dataset.moleculeViewer = btoa(JSON.stringify(moleculeData));
    img.dataset.loaded = 'false';
    img.dataset.rotation = '0';
    img.dataset.compoundType = 'compound';


    // Style it to look good inline
    img.style.cssText = `
      display: inline-block;
      vertical-align: middle;
      margin: 0 4px;
      min-width: 100px;
      min-height: 80px;
      background: rgba(128,128,128,0.1);
      border-radius: 4px;
      padding: 4px;
    `;

    // Replace selected text with the image
    range.deleteContents();
    range.insertNode(img);

    // Clear selection
    selection.removeAllRanges();

    // Load the molecule using the same pipeline as chem:: syntax
    console.log('%c🔄 Triggering molecule render via _loadMoleculeImage...', 'color: #2196F3;');
    if (window._loadMoleculeImage) {
      window._loadMoleculeImage(img);
    } else {
      console.error('❌ _loadMoleculeImage not found on window object');
      img.alt = 'Error: Renderer not ready';
    }
  }

  // Handler for "Render as SMILES" (treats text directly as SMILES string)
  // Uses the SAME pipeline as chem:: syntax, but with isDirectSmiles flag
  // This ensures caching, persistence, and consistent UI
  if (request.type === 'RENDER_SMILES') {
    const text = request.text?.trim();
    if (!text) return;

    console.log('%c🧪 Rendering selection as SMILES (same pipeline as chem::):', 'color: #4CAF50; font-weight: bold;', text);

    // Get the current selection to find where to insert the image
    const selection = window.getSelection();
    if (!selection.rangeCount) return;

    const range = selection.getRangeAt(0);

    // Create image element (same structure as standard renderer)
    const img = document.createElement('img');
    img.className = 'molecule-diagram molecule-viewer';
    img.alt = text;  // Just the SMILES, not "SMILES: CCO"
    img.title = `SMILES: ${text}`;
    img.src = ''; // Empty src, will be filled by loader

    // Create molecule data - treats text directly as SMILES
    // isDirectSmiles flag tells renderClientSide to skip API lookup
    const moleculeData = {
      smiles: text,  // Pre-set the SMILES so it skips lookup
      nomenclature: text,  // Just the SMILES string (not "SMILES: CCO")
      type: 'smiles',

      options: {
        width: 400,
        height: 300,
        scale: 1,
        rotation: 0
      },
      isMoleculeViewer: true,
      flags: {
        compoundType: 'compound',
        isDirectSmiles: true  // This tells renderer to use text as SMILES directly
      }
    };

    // Set dataset
    img.dataset.moleculeViewer = btoa(JSON.stringify(moleculeData));
    img.dataset.smiles = text;  // Also set directly for caching
    img.dataset.loaded = 'false';
    img.dataset.rotation = '0';
    img.dataset.compoundType = 'compound';

    // Style it to look good inline
    img.style.cssText = `
      display: inline-block;
      vertical-align: middle;
      margin: 0 4px;
      min-width: 100px;
      min-height: 80px;
      background: rgba(128,128,128,0.1);
      border-radius: 4px;
      padding: 4px;
    `;

    // Replace selected text with the image
    range.deleteContents();
    range.insertNode(img);

    // Clear selection
    selection.removeAllRanges();

    // Load the molecule using the same pipeline as chem:: syntax
    console.log('%c🔄 Triggering SMILES render via _loadMoleculeImage...', 'color: #4CAF50;');
    if (window._loadMoleculeImage) {
      window._loadMoleculeImage(img);
    } else {
      console.error('❌ _loadMoleculeImage not found on window object');
      img.alt = 'Error: Renderer not ready';
    }
  }

  // Handler for "Render as Biomolecule" (goes directly to RCSB PDB)
  // Uses the SAME pipeline as chem:: syntax, but with compoundType='biomolecule' flag
  // This ensures caching, persistence, and consistent UI
  if (request.type === 'RENDER_BIOMOLECULE') {
    const text = request.text?.trim();
    if (!text) return;

    console.log('%c🧬 Rendering selection as Biomolecule (same pipeline as chem::):', 'color: #E91E63; font-weight: bold;', text);

    // Get the current selection to find where to insert the image
    const selection = window.getSelection();
    if (!selection.rangeCount) return;

    const range = selection.getRangeAt(0);

    // Create image element (same structure as standard renderer)
    const img = document.createElement('img');
    img.className = 'molecule-diagram molecule-viewer';
    img.alt = text;
    img.title = `Searching RCSB PDB: ${text}`;
    img.src = ''; // Empty src, will be filled by loader

    // Create molecule data - forces biomolecule type
    // compoundType='biomolecule' tells renderClientSide to query RCSB directly
    const moleculeData = {
      nomenclature: text,
      type: 'nomenclature',
      options: {
        width: 400,
        height: 300,
        scale: 1,
        rotation: 0
      },
      isMoleculeViewer: true,
      flags: {
        compoundType: 'biomolecule',  // Force biomolecule type - goes to RCSB
        isDirectSmiles: false
      }
    };

    // Set dataset
    img.dataset.moleculeViewer = btoa(JSON.stringify(moleculeData));
    img.dataset.loaded = 'false';
    img.dataset.rotation = '0';
    img.dataset.compoundType = 'biomolecule';

    // Style it to look good inline
    img.style.cssText = `
      display: inline-block;
      vertical-align: middle;
      margin: 0 4px;
      min-width: 100px;
      min-height: 80px;
      background: rgba(128,128,128,0.1);
      border-radius: 4px;
      padding: 4px;
    `;

    // Replace selected text with the image
    range.deleteContents();
    range.insertNode(img);

    // Clear selection
    selection.removeAllRanges();

    // Load the molecule using the same pipeline as chem:: syntax
    console.log('%c🔄 Triggering biomolecule render via _loadMoleculeImage...', 'color: #E91E63;');
    if (window._loadMoleculeImage) {
      window._loadMoleculeImage(img);
    } else {
      console.error('❌ _loadMoleculeImage not found on window object');
      img.alt = 'Error: Renderer not ready';
    }
  }

  // Handler for "Render as Mineral" (goes directly to COD crystal database)
  // Uses the SAME pipeline as chem:: syntax, but with compoundType='mineral' flag
  // This ensures caching, persistence, and consistent UI
  if (request.type === 'RENDER_MINERAL') {
    const text = request.text?.trim();
    if (!text) return;

    console.log('%c💎 Rendering selection as Mineral (same pipeline as chem::):', 'color: #00BCD4; font-weight: bold;', text);

    // Get the current selection to find where to insert the image
    const selection = window.getSelection();
    if (!selection.rangeCount) return;

    const range = selection.getRangeAt(0);

    // Create image element (same structure as standard renderer)
    const img = document.createElement('img');
    img.className = 'molecule-diagram molecule-viewer';
    img.alt = text;
    img.title = `Searching COD: ${text}`;
    img.src = ''; // Empty src, will be filled by loader

    // Create molecule data - forces mineral type
    // compoundType='mineral' tells renderClientSide to query COD directly
    const moleculeData = {
      nomenclature: text,
      type: 'nomenclature',
      options: {
        width: 400,
        height: 300,
        scale: 1,
        rotation: 0
      },
      isMoleculeViewer: true,
      flags: {
        compoundType: 'mineral',  // Force mineral type - goes to COD
        isDirectSmiles: false
      }
    };

    // Set dataset
    img.dataset.moleculeViewer = btoa(JSON.stringify(moleculeData));
    img.dataset.loaded = 'false';
    img.dataset.rotation = '0';
    img.dataset.compoundType = 'mineral';

    // Style it to look good inline
    img.style.cssText = `
      display: inline-block;
      vertical-align: middle;
      margin: 0 4px;
      min-width: 100px;
      min-height: 80px;
      background: rgba(128,128,128,0.1);
      border-radius: 4px;
      padding: 4px;
    `;

    // Replace selected text with the image
    range.deleteContents();
    range.insertNode(img);

    // Clear selection
    selection.removeAllRanges();

    // Load the molecule using the same pipeline as chem:: syntax
    console.log('%c🔄 Triggering mineral render via _loadMoleculeImage...', 'color: #00BCD4;');
    if (window._loadMoleculeImage) {
      window._loadMoleculeImage(img);
    } else {
      console.error('❌ _loadMoleculeImage not found on window object');
      img.alt = 'Error: Renderer not ready';
    }
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
        // Force PubChem lookup for compounds (no auto-detect)
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
        // Force RCSB PDB lookup
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
    const existingDialog = document.getElementById('chemtex-flag-editor');
    if (existingDialog) existingDialog.remove();

    // Create the dialog
    const dialog = document.createElement('div');
    dialog.id = 'chemtex-flag-editor';
    dialog.innerHTML = `
      <div style="position: fixed; top: 0; left: 0; right: 0; bottom: 0; background: rgba(0,0,0,0.5); z-index: 999999; display: flex; align-items: center; justify-content: center;">
        <div style="background: #fff; border-radius: 12px; padding: 24px; min-width: 350px; max-width: 450px; box-shadow: 0 8px 32px rgba(0,0,0,0.3); font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', Roboto, sans-serif;">
          <h3 style="margin: 0 0 8px 0; font-size: 18px; color: #333;">🎛️ ChemTex Flag Editor</h3>
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

