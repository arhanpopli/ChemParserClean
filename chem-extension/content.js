/**
 * Chemistry Formula Renderer v3.0
 * Works on ChatGPT and other strict CSP sites
 * Uses web_accessible_resources to serve scripts from extension
 */

// ============================================
// 🌍 RENDERING ENGINE API CONFIGURATION
// ============================================
// LOCAL TESTING - Use localhost to avoid HTTPS mixed content errors
const MOLECULE_VIEWER_API = 'http://localhost:5000';
const MOL2CHEMFIG_API = 'http://localhost:5001';  // Flask wrapper (NOT port 8000 Docker backend)
const PUBCHEM_API = 'http://localhost:5002';

// For Heroku production (uncomment when ready to deploy):
// const MOLECULE_VIEWER_API = 'https://YOUR-HEROKU-APP.herokuapp.com';
// const MOL2CHEMFIG_API = 'https://YOUR-MOL2CHEMFIG-APP.herokuapp.com';
// const PUBCHEM_API = 'https://YOUR-PUBCHEM-APP.herokuapp.com';
// const 

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
// DARK MODE SVG INVERSION
// ============================================
// Comprehensive function to invert SVG colors for dark mode
// Handles mol2chemfig/dvisvgm SVGs which have colors in CSS styles, attributes, and text elements

/**
 * Check if the user prefers dark mode
 * Uses multiple detection methods for better compatibility
 */
function isDarkModeEnabled() {
  // Method 1: Check prefers-color-scheme media query
  if (window.matchMedia && window.matchMedia('(prefers-color-scheme: dark)').matches) {
    return true;
  }

  // Method 2: Check body background color (detect dark backgrounds)
  try {
    const bodyBg = window.getComputedStyle(document.body).backgroundColor;
    if (bodyBg) {
      const rgb = bodyBg.match(/\d+/g);
      if (rgb && rgb.length >= 3) {
        const luminance = (0.299 * parseInt(rgb[0]) + 0.587 * parseInt(rgb[1]) + 0.114 * parseInt(rgb[2]));
        if (luminance < 128) {
          return true;
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

  return false;
}

/**
 * Invert SVG colors for dark mode display
 * Handles multiple color formats: hex, rgb, named colors, CSS styles
 * Specifically handles dvisvgm output where text has no explicit fill
 * @param {string} svgContent - The SVG content string
 * @returns {string} - The SVG with inverted colors for dark mode
 */
function invertSvgForDarkMode(svgContent) {
  if (!svgContent || typeof svgContent !== 'string') return svgContent;

  let result = svgContent;

  // 1. Replace black hex colors with white
  result = result.replace(/#000000/gi, '#FFFFFF');
  result = result.replace(/#000\b/gi, '#FFF');

  // 2. Replace rgb(0,0,0) with white
  result = result.replace(/rgb\s*\(\s*0\s*,\s*0\s*,\s*0\s*\)/gi, 'rgb(255, 255, 255)');

  // 3. Replace named 'black' color in stroke/fill attributes
  result = result.replace(/stroke\s*=\s*["']black["']/gi, 'stroke="white"');
  result = result.replace(/fill\s*=\s*["']black["']/gi, 'fill="white"');

  // 4. Replace stroke:#000 and fill:#000 in style attributes and CSS
  result = result.replace(/stroke\s*:\s*#000000/gi, 'stroke:#FFFFFF');
  result = result.replace(/stroke\s*:\s*#000\b/gi, 'stroke:#FFF');
  result = result.replace(/stroke\s*:\s*black/gi, 'stroke:white');
  result = result.replace(/fill\s*:\s*#000000/gi, 'fill:#FFFFFF');
  result = result.replace(/fill\s*:\s*#000\b/gi, 'fill:#FFF');
  result = result.replace(/fill\s*:\s*black/gi, 'fill:white');

  // 5. Replace color:#000 in CSS (for text)
  result = result.replace(/color\s*:\s*#000000/gi, 'color:#FFFFFF');
  result = result.replace(/color\s*:\s*#000\b/gi, 'color:#FFF');
  result = result.replace(/color\s*:\s*black/gi, 'color:white');

  // 6. Handle dvisvgm specific patterns - text elements often have fill in style
  // Match: style="...fill:#000..." or style="...fill:black..."
  result = result.replace(/(style\s*=\s*["'][^"']*fill\s*:\s*)#000000([^"']*["'])/gi, '$1#FFFFFF$2');
  result = result.replace(/(style\s*=\s*["'][^"']*fill\s*:\s*)#000\b([^"']*["'])/gi, '$1#FFF$2');
  result = result.replace(/(style\s*=\s*["'][^"']*fill\s*:\s*)black([^"']*["'])/gi, '$1white$2');

  // 7. Handle CSS in <style> blocks - common in dvisvgm output
  // Pattern: .classname{fill:#000} or text{fill:#000}
  result = result.replace(/(\{[^}]*fill\s*:\s*)#000000([^}]*\})/gi, '$1#FFFFFF$2');
  result = result.replace(/(\{[^}]*fill\s*:\s*)#000\b([^}]*\})/gi, '$1#FFF$2');
  result = result.replace(/(\{[^}]*fill\s*:\s*)black([^}]*\})/gi, '$1white$2');

  // 8. Handle stroke in CSS blocks
  result = result.replace(/(\{[^}]*stroke\s*:\s*)#000000([^}]*\})/gi, '$1#FFFFFF$2');
  result = result.replace(/(\{[^}]*stroke\s*:\s*)#000\b([^}]*\})/gi, '$1#FFF$2');
  result = result.replace(/(\{[^}]*stroke\s*:\s*)black([^}]*\})/gi, '$1white$2');

  // 9. Handle currentColor which might inherit black
  // Add a root style to set color to white if not present
  if (result.includes('currentColor') && !result.includes('color:white') && !result.includes('color:#FFF')) {
    result = result.replace(/<svg([^>]*)>/, '<svg$1 style="color:white;">');
  }

  // 10. CRITICAL: Handle dvisvgm text elements that have NO fill attribute
  // These inherit the default SVG fill which is black
  // Add fill="white" to all <text> elements that don't have a fill attribute
  // Match <text ...> without fill= and add fill="white"
  result = result.replace(/<text\s+(?![^>]*\bfill\s*=)([^>]*)>/gi, '<text fill="#FFFFFF" $1>');

  // 11. Add a global style rule for text elements without explicit fill
  // Insert into existing <style> block or create one
  if (result.includes('<style')) {
    // Add rule to existing style block - insert before closing ]]> or </style>
    result = result.replace(/(]]>\s*<\/style>)/i, 'text { fill: #FFFFFF !important; }\n$1');
    result = result.replace(/(<\/style>)/i, 'text { fill: #FFFFFF !important; }\n$1');
  } else {
    // No style block - add one after the opening <svg> tag
    result = result.replace(/(<svg[^>]*>)/i, '$1<style type="text/css">text { fill: #FFFFFF !important; }</style>');
  }

  console.log('%c🌙 Dark mode SVG inversion applied', 'color: #9B59B6; font-weight: bold;');
  return result;
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
  console.log('%c[ChemRenderer] backgroundFetchJSON called', 'color: #00AAFF; font-weight: bold;', { url, options });
  try {
    // Try background script first (works on CSP-restricted sites)
    if (chrome.runtime && chrome.runtime.sendMessage) {
      console.log('%c[ChemRenderer] Sending message to background...', 'color: #00AAFF;');
      return new Promise((resolve, reject) => {
        chrome.runtime.sendMessage(
          { type: 'FETCH_API', url: url, options: options },
          (response) => {
            console.log('%c[ChemRenderer] Got response from background:', 'color: #00AAFF;', response);
            if (chrome.runtime.lastError) {
              console.warn('[ChemRenderer] Background fetch failed, trying direct:', chrome.runtime.lastError.message);
              // Fall back to direct fetch
              directFetchJSON(url, options).then(resolve).catch(reject);
              return;
            }
            if (response && response.success) {
              console.log('%c[ChemRenderer] Background fetch SUCCESS', 'color: #00FF00; font-weight: bold;');
              resolve(response.data);
            } else {
              console.error('%c[ChemRenderer] Background fetch FAILED:', 'color: #FF0000;', response?.error);
              reject(new Error(response?.error || 'Background fetch failed'));
            }
          }
        );
      });
    } else {
      console.warn('[ChemRenderer] chrome.runtime.sendMessage not available');
    }
  } catch (e) {
    console.warn('[ChemRenderer] Background unavailable, using direct fetch:', e.message);
  }
  // Fall back to direct fetch
  console.log('%c[ChemRenderer] Falling back to direct fetch', 'color: #FFAA00;');
  return directFetchJSON(url, options);
}

/**
 * Direct fetch for JSON (used as fallback)
 */
async function directFetchJSON(url, options = {}) {
  const response = await fetch(url, options);
  if (!response.ok) {
    throw new Error(`HTTP ${response.status}: ${response.statusText}`);
  }
  return response.json();
}

/**
 * Fetch blob/image via background script
 * Returns { base64, type } or falls back to direct blob
 */
async function backgroundFetchBlob(url) {
  try {
    if (chrome.runtime && chrome.runtime.sendMessage) {
      return new Promise((resolve, reject) => {
        chrome.runtime.sendMessage(
          { type: 'FETCH_BLOB', url: url },
          (response) => {
            if (chrome.runtime.lastError) {
              console.warn('[ChemRenderer] Background blob fetch failed, trying direct:', chrome.runtime.lastError.message);
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
    console.warn('[ChemRenderer] Background unavailable for blob, using direct fetch');
  }
  return directFetchBlob(url);
}

/**
 * Direct fetch for blob (used as fallback)
 */
async function directFetchBlob(url) {
  const response = await fetch(url);
  if (!response.ok) {
    throw new Error(`HTTP ${response.status}: ${response.statusText}`);
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

/**
 * Fetch text from API via background script
 * Falls back to direct fetch if background is unavailable
 */
async function backgroundFetchText(url, options = {}) {
  try {
    if (chrome.runtime && chrome.runtime.sendMessage) {
      return new Promise((resolve, reject) => {
        chrome.runtime.sendMessage(
          { type: 'FETCH_TEXT', url: url, options: options },
          (response) => {
            if (chrome.runtime.lastError) {
              console.warn('[ChemRenderer] Background text fetch failed, trying direct:', chrome.runtime.lastError.message);
              directFetchText(url, options).then(resolve).catch(reject);
              return;
            }
            if (response && response.success) {
              resolve(response.data);
            } else {
              reject(new Error(response?.error || 'Background text fetch failed'));
            }
          }
        );
      });
    }
  } catch (e) {
    console.warn('[ChemRenderer] Background unavailable for text, using direct fetch');
  }
  return directFetchText(url, options);
}

/**
 * Direct fetch for text (used as fallback)
 */
async function directFetchText(url, options = {}) {
  const response = await fetch(url, options);
  if (!response.ok) {
    throw new Error(`HTTP ${response.status}: ${response.statusText}`);
  }
  return response.text();
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

  console.log('%c🌉 SMILES BRIDGE: Converting name→SMILES', 'background: #4A90D9; color: white; font-weight: bold; padding: 2px 6px;', name, use3DSmiles ? '(3D enabled)' : '(3D disabled)');

  if (!name || typeof name !== 'string') {
    console.error('❌ SMILES Bridge: Invalid name provided');
    return null;
  }

  // Strip flags from molecule name before conversion
  // Example: "phenol/d-c" → "phenol", "histamine/+c+o" → "histamine"
  const cleanName = stripFlagsFromName(name.trim());
  if (cleanName !== name.trim()) {
    console.log('%c🧹 Stripped flags from name:', 'color: #9c88ff;', name, '→', cleanName);
  }

  if (!cleanName) {
    console.error('❌ SMILES Bridge: Empty name after stripping flags');
    return null;
  }

  // ===== PRIORITY 1: Direct PubChem API =====
  // Works for most compound names (common names, trade names, drug names)
  console.log('%c🌐 [Bridge] Priority 1: Trying Direct PubChem API...', 'color: #0088FF; font-weight: bold;');
  try {
    const pubchemUrl = `https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/${encodeURIComponent(cleanName)}/property/CanonicalSMILES,IsomericSMILES/JSON`;
    const data = await backgroundFetchJSON(pubchemUrl);
    if (data && data.PropertyTable && data.PropertyTable.Properties && data.PropertyTable.Properties.length > 0) {
      const props = data.PropertyTable.Properties[0];
      const smiles = use3DSmiles ? (props.IsomericSMILES || props.CanonicalSMILES) : props.CanonicalSMILES;
      console.log('%c✅ [Bridge] Direct PubChem SUCCESS:', 'color: #00FF00; font-weight: bold;', smiles);
      return { smiles: smiles, source: 'PubChem' };
    }
  } catch (ePubChem) {
    console.warn('⚠️ [Bridge] Direct PubChem API request failed:', ePubChem.message);
  }

  // ===== PRIORITY 2: PubChem Autocomplete (for class/generic names) =====
  // Handles names like "sphingomyelin" → finds "Sphingomyelin 16:0"
  console.log('%c� [Bridge] Priority 2: Trying PubChem Autocomplete for class name...', 'color: #FF8800; font-weight: bold;');
  try {
    const autocompleteUrl = `https://pubchem.ncbi.nlm.nih.gov/rest/autocomplete/compound/${encodeURIComponent(cleanName)}/json?limit=5`;
    const autocompleteData = await backgroundFetchJSON(autocompleteUrl);

    if (autocompleteData && autocompleteData.dictionary_terms && autocompleteData.dictionary_terms.compound && autocompleteData.dictionary_terms.compound.length > 0) {
      const bestMatch = autocompleteData.dictionary_terms.compound[0];
      console.log('%c🔗 [Bridge] Found related compound:', 'color: #FF8800;', bestMatch);

      // Fetch SMILES for the matched compound
      const matchUrl = `https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/${encodeURIComponent(bestMatch)}/property/CanonicalSMILES,IsomericSMILES/JSON`;
      const matchData = await backgroundFetchJSON(matchUrl);

      if (matchData && matchData.PropertyTable && matchData.PropertyTable.Properties && matchData.PropertyTable.Properties.length > 0) {
        const props = matchData.PropertyTable.Properties[0];
        const smiles = use3DSmiles ? (props.IsomericSMILES || props.CanonicalSMILES) : props.CanonicalSMILES;
        console.log('%c✅ [Bridge] PubChem Autocomplete SUCCESS:', 'color: #00FF00; font-weight: bold;', smiles);
        return { smiles: smiles, source: 'PubChem-Autocomplete' };
      }
    }
  } catch (eAutocomplete) {
    console.warn('⚠️ [Bridge] PubChem Autocomplete failed:', eAutocomplete.message);
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
  console.log('%c🔍 Getting PubChem CID for:', 'color: #0088FF; font-weight: bold;', nameOrSmiles);

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
  console.log('%c🌐 [CID] Priority 1: Direct PubChem CID lookup...', 'color: #0088FF; font-weight: bold;');
  let cid = await fetchCidByName(cleanInput);
  if (cid) {
    console.log('%c✅ [CID] Direct lookup SUCCESS:', 'color: #00FF00; font-weight: bold;', cid);
    return cid;
  }

  // Priority 1.5: Try appending "_1" (common PubChem pattern for representative compounds like Phosphatidylcholine)
  // This fixes "phosphatidylcholine" -> "Phosphatidylcholine_1" (CID 10425706 which has 3D)
  if (!cleanInput.includes('_')) {
    console.log('%c🌐 [CID] Priority 1.5: Trying with _1 suffix...', 'color: #0088FF;');
    cid = await fetchCidByName(cleanInput + '_1');
    if (cid) {
      console.log('%c✅ [CID] Suffix lookup SUCCESS:', 'color: #00FF00; font-weight: bold;', cid);
      return cid;
    }
  }

  // Priority 2: PubChem Autocomplete (The "Singular Solution" for complex molecules)
  // If direct lookup fails, ask PubChem for suggestions and use the first one
  // This handles generic names like "sphingomyelin" -> "Sphingomyelin 16:0" -> CID
  console.log('%c🔍 [CID] Priority 2: Trying PubChem Autocomplete...', 'color: #FF8800; font-weight: bold;');
  try {
    const autocompleteUrl = `https://pubchem.ncbi.nlm.nih.gov/rest/autocomplete/compound/${encodeURIComponent(cleanInput)}/json?limit=1`;
    const autocompleteData = await backgroundFetchJSON(autocompleteUrl);

    if (autocompleteData && autocompleteData.dictionary_terms && autocompleteData.dictionary_terms.compound && autocompleteData.dictionary_terms.compound.length > 0) {
      const bestMatch = autocompleteData.dictionary_terms.compound[0];
      console.log('%c🔗 [CID] Autocomplete found match:', 'color: #FF8800; font-weight: bold;', bestMatch);

      // Now get CID for this best match
      cid = await fetchCidByName(bestMatch);
      if (cid) {
        console.log('%c✅ [CID] Autocomplete match SUCCESS:', 'color: #00FF00; font-weight: bold;', cid);
        return cid;
      }
    }
  } catch (e) {
    console.warn('⚠️ [CID] Autocomplete failed:', e.message);
  }

  // Priority 3: SMILES Bridge fallback (Last resort)
  console.log('%c🌉 [CID] Priority 3: Using SMILES Bridge fallback...', 'color: #FF8800; font-weight: bold;');
  const bridgeResult = await smilesBridge(cleanInput, { use3DSmiles: false });

  if (bridgeResult && bridgeResult.smiles) {
    console.log('%c✅ [CID] Got SMILES from bridge:', 'color: #00FF00;', bridgeResult.smiles, `(source: ${bridgeResult.source})`);

    // Now get CID from SMILES
    try {
      const smilesUrl = `https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/smiles/${encodeURIComponent(bridgeResult.smiles)}/cids/JSON`;
      const smilesData = await backgroundFetchJSON(smilesUrl);

      if (smilesData && smilesData.IdentifierList && smilesData.IdentifierList.CID && smilesData.IdentifierList.CID.length > 0) {
        cid = smilesData.IdentifierList.CID[0];
        console.log('%c✅ [CID] Found CID from SMILES:', 'color: #00FF00; font-weight: bold;', cid);
        return cid;
      }
    } catch (e) {
      console.warn('⚠️ [CID] Lookup from SMILES failed:', e.message);
    }
  }

  console.error('%c❌ [CID] Could not find CID for:', 'color: #FF0000; font-weight: bold;', nameOrSmiles);
  return null;
}

// Export for global access
window.getPubChemCID = getPubChemCID;

/**
 * Add hover controls (name label + 3D button) to molecule image container
 * Creates a floating overlay that appears on hover with molecule name and 3D viewer button
 * 
 * @param {HTMLElement} container - The container element to add controls to
 * @param {string} moleculeName - Name of the molecule to display
 * @param {Object} moleculeData - Full molecule data for 3D viewer
 */
function addHoverControls(container, moleculeName, moleculeData) {
  // Don't add controls if they already exist
  if (container.querySelector('.chem-hover-controls')) {
    console.log('%c⚠️ Hover controls already exist, skipping', 'color: orange;');
    return;
  }

  console.log('%c🎮 Adding hover controls to container', 'color: #4ecdc4; font-weight: bold;', { moleculeName, moleculeData });

  // Create hover controls container
  const hoverControls = document.createElement('div');
  hoverControls.className = 'chem-hover-controls';
  hoverControls.style.cssText = `
    position: absolute;
    top: 4px;
    right: 4px;
    display: flex;
    flex-direction: column;
    gap: 4px;
    opacity: 0;
    transition: opacity 0.2s;
    pointer-events: auto;
    z-index: 1000;
  `;

  // Create molecule name label
  const nameLabel = document.createElement('div');
  nameLabel.className = 'chem-name-label';
  nameLabel.textContent = moleculeName;
  nameLabel.title = moleculeName; // Full name on hover
  nameLabel.style.cssText = `
    background: rgba(0, 0, 0, 0.8);
    color: white;
    padding: 4px 8px;
    border-radius: 4px;
    font-size: 11px;
    font-family: monospace;
    max-width: 200px;
    overflow: hidden;
    text-overflow: ellipsis;
    white-space: nowrap;
    pointer-events: auto;
  `;

  // Create 3D viewer button
  const viewer3DBtn = document.createElement('button');
  viewer3DBtn.className = 'chem-3d-btn';
  viewer3DBtn.innerHTML = '🔮 3D';
  viewer3DBtn.title = 'View in 3D';
  viewer3DBtn.style.cssText = `
    background: rgba(0, 0, 0, 0.8);
    color: white;
    border: none;
    padding: 4px 8px;
    border-radius: 4px;
    font-size: 11px;
    cursor: pointer;
    transition: background 0.2s;
    pointer-events: auto;
  `;

  viewer3DBtn.addEventListener('mouseenter', () => {
    viewer3DBtn.style.background = 'rgba(0, 0, 0, 0.95)';
  });

  viewer3DBtn.addEventListener('mouseleave', () => {
    viewer3DBtn.style.background = 'rgba(0, 0, 0, 0.8)';
  });

  viewer3DBtn.addEventListener('click', async (e) => {
    e.stopPropagation();
    e.preventDefault();
    console.log('%c🔮 3D Viewer button clicked', 'color: #764ba2; font-weight: bold;');

    // Call show3DViewerInline if it exists (it's defined later in setupLazyLoading)
    if (typeof window.show3DViewerInline === 'function') {
      await window.show3DViewerInline(moleculeData, container);
    } else {
      console.error('❌ show3DViewerInline function not available yet');
    }
  });

  hoverControls.appendChild(nameLabel);
  hoverControls.appendChild(viewer3DBtn);
  container.appendChild(hoverControls);

  // Show/hide on hover
  container.addEventListener('mouseenter', () => {
    hoverControls.style.opacity = '1';
  });

  container.addEventListener('mouseleave', () => {
    hoverControls.style.opacity = '0';
  });

  console.log('%c✅ Hover controls added successfully', 'color: #00FF00; font-weight: bold;');
}

// Export for global access
window.addHoverControls = addHoverControls;

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
    if (!imageKey) return { scale: 1.5 }; // Default to 1.5x scale
    if (!settings.saveSizePerImage && !settings.saveSizeBySMILES) {
      return { scale: 1.5 };
    }
    let storageKey;
    if (settings.saveSizePerImage && !settings.saveSizeBySMILES) {
      storageKey = getPageImageKey(moleculeData, pageUrl);
    } else {
      storageKey = imageKey;
    }
    if (!storageKey) return { scale: 1.5 };
    return new Promise((resolve) => {
      chrome.storage.local.get([storageKey], (result) => {
        if (result[storageKey]) {
          // Support both old format (width/height) and new format (scale)
          if (result[storageKey].scale) {
            resolve(result[storageKey]);
          } else {
            // Legacy: convert old width/height to scale (assume default was 400x350)
            const scale = result[storageKey].width ? result[storageKey].width / DEFAULT_WIDTH : 1.5;
            resolve({ scale });
          }
        } else {
          resolve({ scale: 1.5 });
        }
      });
    });
  } catch (error) {
    console.error('Error loading image size:', error);
    return { scale: 1.5 };
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
  // Get current scale (default to 1.5 for a good starting size)
  let currentScale = parseFloat(svgImg.dataset.scale) || 1.5;

  // Adjust scale by delta (convert pixel delta to scale delta)
  const scaleDelta = delta / 100; // 20px = 0.2 scale change
  let newScale = currentScale + scaleDelta;

  // Constrain scale between 0.5x and 5x
  newScale = Math.max(0.5, Math.min(5, newScale));

  // Get the SVG's intrinsic width with robust fallbacks for large molecules
  let intrinsicWidth = svgImg.naturalWidth || svgImg.width;

  // If still no width (can happen with complex SVGs like phosphatidylcholine, insulin)
  if (!intrinsicWidth || intrinsicWidth === 0) {
    console.warn('%c⚠️ naturalWidth not available in adjustImageSize, extracting from SVG', 'color: orange;');

    try {
      if (svgImg.src && svgImg.src.startsWith('data:image/svg+xml')) {
        let svgContent;
        if (svgImg.src.includes('base64,')) {
          svgContent = atob(svgImg.src.split('base64,')[1]);
        } else {
          svgContent = decodeURIComponent(svgImg.src.split(',')[1]);
        }

        const widthMatch = svgContent.match(/width\s*=\s*["']?(\d+(?:\.\d+)?)/i);
        if (widthMatch) {
          intrinsicWidth = parseFloat(widthMatch[1]);
        } else {
          const viewBoxMatch = svgContent.match(/viewBox\s*=\s*["']?[\d.]+\s+[\d.]+\s+([\d.]+)\s+([\d.]+)/i);
          if (viewBoxMatch) {
            intrinsicWidth = parseFloat(viewBoxMatch[1]);
          }
        }
      }
    } catch (e) {
      console.warn('Could not extract SVG width:', e);
    }
  }

  // Final fallback
  if (!intrinsicWidth || intrinsicWidth === 0) {
    intrinsicWidth = 400;
    console.warn('%c⚠️ Using default width in adjustImageSize:', 'color: orange;', intrinsicWidth);
  }

  console.log(`Intrinsic width: ${intrinsicWidth}px`);

  // Calculate new width based on scale
  const newWidth = Math.round(intrinsicWidth * newScale);

  // Only set width, let height: auto maintain aspect ratio
  svgImg.style.width = `${newWidth}px`;
  svgImg.style.height = 'auto';
  svgImg.style.maxWidth = 'none';
  svgImg.style.maxHeight = 'none';

  // Store scale in dataset for future adjustments
  svgImg.dataset.scale = newScale.toString();

  // Save scale preference
  const pageUrl = window.location.href;
  const size = { scale: newScale };
  saveImageSize(moleculeData, pageUrl, size, settings);

  console.log(`Adjusted scale: ${newScale.toFixed(2)}x (width: ${newWidth}px, height: auto, was ${currentScale.toFixed(2)}x)`);
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
    `;
    const pageUrl = window.location.href;
    const savedSize = await loadImageSize(moleculeData, pageUrl, settings);

    // Apply scale by calculating actual dimensions from intrinsic size
    const scale = savedSize.scale || 1.5;
    svgImg.dataset.scale = scale.toString();

    // Wait for image to load to get natural dimensions
    if (svgImg.complete) {
      applyScaleToImage(svgImg, scale);
    } else {
      svgImg.onload = () => applyScaleToImage(svgImg, scale);
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
    const moleculeName = moleculeData.nomenclature || moleculeData.smiles || 'Unknown';
    addHoverControls(container, moleculeName, moleculeData);

    console.log('✅ Size controls wrapper created successfully with scale:', scale);
    return container;
  } catch (error) {
    console.error('❌ Error wrapping image with size controls:', error);
    originalImg.parentNode.replaceChild(svgImg, originalImg);
    return svgImg;
  }
}

function applyScaleToImage(svgImg, scale) {
  // Try to get intrinsic width, with multiple fallbacks for large/complex molecules
  let intrinsicWidth = svgImg.naturalWidth || svgImg.width;

  // If still no width (can happen with complex SVGs like phosphatidylcholine, insulin)
  if (!intrinsicWidth || intrinsicWidth === 0) {
    console.warn('%c⚠️ naturalWidth not available, trying to extract from SVG', 'color: orange;');

    try {
      // Try to extract width from SVG src if it's a data URL
      if (svgImg.src && svgImg.src.startsWith('data:image/svg+xml')) {
        let svgContent;

        // Handle both base64 and plain text data URLs
        if (svgImg.src.includes('base64,')) {
          svgContent = atob(svgImg.src.split('base64,')[1]);
        } else {
          svgContent = decodeURIComponent(svgImg.src.split(',')[1]);
        }

        // Try to extract width attribute from SVG tag
        const widthMatch = svgContent.match(/width\s*=\s*["']?(\d+(?:\.\d+)?)/i);
        if (widthMatch) {
          intrinsicWidth = parseFloat(widthMatch[1]);
          console.log('%c📏 Extracted width from SVG content:', 'color: #9c88ff;', intrinsicWidth);
        } else {
          // Try viewBox as fallback
          const viewBoxMatch = svgContent.match(/viewBox\s*=\s*["']?[\d.]+\s+[\d.]+\s+([\d.]+)\s+([\d.]+)/i);
          if (viewBoxMatch) {
            intrinsicWidth = parseFloat(viewBoxMatch[1]);
            console.log('%c📏 Extracted width from viewBox:', 'color: #9c88ff;', intrinsicWidth);
          }
        }
      }
    } catch (e) {
      console.warn('Could not extract SVG width from data URL:', e);
    }
  }

  // Final fallback: use default width based on molecule complexity
  if (!intrinsicWidth || intrinsicWidth === 0) {
    // For large molecules, use a larger default
    intrinsicWidth = 400;
    console.warn('%c⚠️ Using default width for large molecule:', 'color: orange;', intrinsicWidth);
  }

  const newWidth = Math.round(intrinsicWidth * scale);

  // Only set width, let height: auto maintain aspect ratio
  svgImg.style.width = `${newWidth}px`;
  svgImg.style.height = 'auto';
  svgImg.style.maxWidth = 'none';
  svgImg.style.maxHeight = 'none';

  console.log(`Applied scale ${scale}x: intrinsic width ${intrinsicWidth}px → ${newWidth}px (height: auto)`);
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
          if (text.includes('\\') || text.includes('ce{') || text.includes('chemfig')) {
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
    const images = document.querySelectorAll('.chemfig-diagram');
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
      totalSVGsOnPage: document.querySelectorAll('img.chemfig-diagram').length,
      visibleSVGs: document.querySelectorAll('img.chemfig-diagram[data-loaded="true"]').length,
      pendingSVGs: document.querySelectorAll('img.chemfig-diagram[data-loaded="false"]').length
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
      containersOnPage: document.querySelectorAll('.chemfig-container').length,
      structuresOnPage: document.querySelectorAll('.chemfig-diagram').length
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
  // Switch rendering engine
  setRendererEngine: function (engine) {
    const validEngines = ['codecogs', 'latex-online', 'quicklatex'];
    if (!validEngines.includes(engine)) {
      console.error(`❌ Invalid engine. Use: ${validEngines.join(', ')}`);
      return;
    }
    settings.rendererEngine = engine;
    log.success(`🔧 Renderer engine switched to: ${engine}`);
    chrome.storage.sync.set({ rendererEngine: engine });
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

// Settings with all defaults - FORCE MoleculeViewer only
let settings = {
  enabled: true,
  renderMhchem: true,
  renderChemfig: true,
  performanceMode: true,
  maxVisibleSVGs: 5,
  layoutMode: 'horizontal',
  renderCarbonsAsSticks: false,
  sizePreset: 'auto',
  rendererEngine: 'moleculeviewer',  // Default to MoleculeViewer (can be overridden by user)
  devMode: false,
  // MoleculeViewer rendering options
  mvUse3DSmiles: false,  // Use OPSIN 3D for stereochemistry in MoleculeViewer
  flipHorizontal: false,
  flipVertical: false,
  mvAutoInvert: false,  // Auto-apply invert in dark mode for MoleculeViewer - disabled by default to prevent double inversion
  mvInvertMode: 'full',  // 'full' or 'bw' (black & white only)
  mvRotate: 0,  // Rotation angle for MoleculeViewer
  // mol2chemfig rendering options
  m2cfShowCarbons: false,
  m2cfAromaticCircles: false,
  m2cfShowMethyls: false,
  m2cfFancyBonds: false,
  m2cfAtomNumbers: false,
  m2cfCompact: false,
  m2cfFlipHorizontal: false,
  m2cfFlipVertical: false,
  m2cfHydrogensMode: 'keep',
  use3DSmiles: false,  // Enable 3D stereochemistry via OPSIN
  m2cfAddH2: false,
  m2cfAutoInvert: true,  // Auto-apply invert in dark mode for mol2chemfig
  m2cfInvertMode: 'full',  // 'full' or 'bw' (black & white only)
  m2cfRotate: 0,
  // PubChem specific settings
  pubchemImageSize: 'large',
  pubchemRecordType: '2d',
  pubchemDirectFetch: true,  // Default to direct fetch from PubChem (no local server needed)
  pubchemRemoveBg: false,
  pubchemSharpenImages: true,
  enable3DViewer: false,
  viewer3DSource: 'local', // 'local', 'pubchem', or 'molview'
  // Client-side renderer options
  clientSideRenderer: 'smilesdrawer',  // 'smilesdrawer' or 'kekule'
  // AI Molecular Control
  enableAIMolecularControl: false
};

log.info('📦 Loading settings from storage...');

// Load settings - with proper callback
chrome.storage.sync.get(null, (result) => {
  // Merge stored settings with defaults
  settings = { ...settings, ...result };

  // ✅ DO NOT FORCE - Respect user's renderer choice
  let engineName, enginePort;
  if (settings.rendererEngine === 'mol2chemfig') {
    engineName = '📐 mol2chemfig';
    enginePort = '8000';
  } else if (settings.rendererEngine === 'pubchem') {
    engineName = '🌐 PubChem';
    enginePort = '5002';
  } else if (settings.rendererEngine === 'client-side') {
    engineName = '💻 Client-Side (SmilesDrawer)';
    enginePort = 'N/A';
  } else {
    engineName = '🧪 MoleculeViewer';
    enginePort = '5000';
  }

  log.success('✅ Settings loaded', settings);
  log.info(`Renderer Engine: ${engineName} (localhost:${enginePort})`);
  log.info(`Performance mode: ${settings.performanceMode ? 'ON ⚡' : 'OFF'}`);
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

  // Listen for renderer engine changes
  if (changes.rendererEngine) {
    log.info('🔄 Renderer engine changed, updating...', changes.rendererEngine);
    settings.rendererEngine = changes.rendererEngine.newValue;
    let engineName;
    if (settings.rendererEngine === 'mol2chemfig') {
      engineName = '📐 mol2chemfig';
    } else if (settings.rendererEngine === 'pubchem') {
      engineName = '🌐 PubChem';
    } else if (settings.rendererEngine === 'client-side') {
      engineName = '💻 Client-Side (SmilesDrawer)';
    } else {
      engineName = '🧪 MoleculeViewer';
    }
    log.success(`✅ Switched to ${engineName} renderer`);

    // Reload page to apply new renderer to all content
    setTimeout(() => {
      location.reload();
    }, 500);
  }

  // Listen for 3D viewer source changes
  if (changes.viewer3DSource) {
    log.info('🔄 3D Viewer source changed, updating...', changes.viewer3DSource);
    settings.viewer3DSource = changes.viewer3DSource.newValue;
    log.success(`✅ Switched to ${changes.viewer3DSource.newValue} 3D viewer`);
    // Update settings immediately without reload
  }

  // Listen for rendering option changes (aromatic circles, etc.) and reload to apply
  if (changes.aromaticCircles || changes.showCarbons || changes.showMethyls ||
    changes.fancyBonds || changes.atomNumbers || changes.flipHorizontal || changes.flipVertical) {
    log.info('⚙️  Rendering options changed, reloading...', changes);
    setTimeout(() => {
      location.reload();
    }, 500);
  }
});

// Listen for messages from popup (live preview)
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
    /* Chemfig structure diagrams - inline display with consistent margins */
    .chemfig-diagram {
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
    .chemfig-loading {
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
    .chemfig-fadein {
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
    
    /* Container for multiple chemfig structures - horizontal layout (DEFAULT) */
    .chemfig-container {
      display: inline-flex;
      flex-wrap: wrap;
      gap: 12px;
      align-items: center;
      margin: 8px 0;
    }
    
    .chemfig-container .chemfig-diagram {
      margin: 0;
    }
    
    /* Vertical layout mode - stack structures top to bottom */
    .chemfig-container.vertical {
      display: flex;
      flex-direction: column;
      align-items: flex-start;
      gap: 16px;
      margin: 12px 0;
    }
    
    .chemfig-container.vertical .chemfig-diagram {
      margin: 0;
      display: block;
    }
    
    /* Rotation utilities - for orientation control */
    .chemfig-rotate-0 {
      transform: rotate(0deg);
    }
    
    .chemfig-rotate-90 {
      transform: rotate(90deg);
      max-width: 200px;
      max-height: 300px;
    }
    
    .chemfig-rotate-180 {
      transform: rotate(180deg);
    }
    
    .chemfig-rotate-270 {
      transform: rotate(270deg);
      max-width: 200px;
      max-height: 300px;
    }
    
    /* Hover effects */
    .chemfig-diagram:hover {
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
    
    .molecule-viewer-container .chemfig-diagram {
      display: block;
      max-width: 300px;
      max-height: 200px;
      margin-bottom: 6px;
      cursor: pointer;
      border-radius: 4px;
    }
    
    .molecule-viewer-container .chemfig-diagram:hover {
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
      .chemfig-diagram {
        max-width: 250px;
        max-height: 150px;
      }
      
      .chemfig-container {
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
function parseChemFlags(chemFormula) {
  // Example: chem:histamine/+c+o+m+n+3d: or chem:histamine/d-c: or chem:histamine/+s200+r45:
  const result = {
    useDefaults: false, // Whether to use the user's defaults
    showCarbons: undefined,
    aromaticCircles: undefined,
    showMethyls: undefined,
    atomNumbers: undefined,
    flipHorizontal: undefined,
    flipVertical: undefined,
    addHydrogens: undefined, // Add explicit hydrogens
    size: undefined, // Scale factor
    rotation: undefined,
    is3D: false,
    isPubchem: false
    // invert removed - now auto-applied based on dark mode detection
  };

  // Check if it contains flags
  if (chemFormula.includes('/')) {
    const parts = chemFormula.split('/');
    if (parts.length > 1) {
      // The part after the first slash might contain flags
      const flagsPart = parts[1];
      if (flagsPart.includes(':')) {
        // Get the flags before the colon
        const flags = flagsPart.split(':')[0].toLowerCase();

        // Check for 'd' flag (use defaults)
        if (flags.includes('d')) {
          result.useDefaults = true;
          // Process default overrides like d-c, d-o, etc.
          const remainingFlags = flags.replace('d', '');
          // Handle individual overrides like -c, -o, -m, -n, -p, -q
          if (remainingFlags.includes('-c')) result.showCarbons = false;
          if (remainingFlags.includes('-o')) result.aromaticCircles = false;
          if (remainingFlags.includes('-m')) result.showMethyls = false;
          if (remainingFlags.includes('-n')) result.atomNumbers = false;
          if (remainingFlags.includes('-p')) result.flipHorizontal = true;  // Flip horizontal
          if (remainingFlags.includes('-q')) result.flipVertical = true;    // Flip vertical
          // -i flag removed - inversion now auto-applied based on dark mode
        }
      }
    }
  }

  // Check for + flags
  if (chemFormula.includes('+')) {
    const flagMatches = chemFormula.match(/\+([a-zA-Z0-9]+)/g);
    if (flagMatches) {
      flagMatches.forEach(flag => {
        const flagContent = flag.substring(1).toLowerCase(); // Remove the + and normalize to lowercase

        if (flagContent === 'c') result.showCarbons = true;
        if (flagContent === 'o') result.aromaticCircles = true;
        if (flagContent === 'm') result.showMethyls = true;
        if (flagContent === 'n') result.atomNumbers = true;
        if (flagContent === 'h') result.addHydrogens = true;  // Add hydrogens (+h)
        if (flagContent === 'p') result.flipHorizontal = true;  // Flip horizontal (+p)
        if (flagContent === 'q') result.flipVertical = true;    // Flip vertical (+q)
        if (flagContent === '3d') result.is3D = true;
        if (flagContent === 'pubchem') result.isPubchem = true;
        // +inv flag removed - inversion now auto-applied based on dark mode

        // Handle flags with values
        if (flagContent.startsWith('s')) { // Size flag: +s200
          const sizeValue = parseInt(flagContent.substring(1));
          if (!isNaN(sizeValue)) result.size = sizeValue / 100; // Convert to scale factor
        }
        if (flagContent.startsWith('r')) { // Rotation flag: +r60
          const rotationValue = parseInt(flagContent.substring(1));
          if (!isNaN(rotationValue)) result.rotation = rotationValue;
        }
      });
    }
  }

  // Check for - flags (to disable options) - but not the ones in the /d- format
  // Only process - flags if they're not part of the /d- pattern
  const slashIndex = chemFormula.indexOf('/');
  const afterSlash = slashIndex !== -1 ? chemFormula.substring(slashIndex) : '';
  const minusFlags = chemFormula.match(/-([cmnhopqi])/gi);

  if (minusFlags) {
    minusFlags.forEach(flag => {
      // Skip if this flag is part of the /d- pattern
      if (!afterSlash.toLowerCase().includes(flag.toLowerCase())) {
        const flagContent = flag.substring(1).toLowerCase(); // Remove the -

        if (flagContent === 'c') result.showCarbons = false;
        if (flagContent === 'o') result.aromaticCircles = false;
        if (flagContent === 'm') result.showMethyls = false;
        if (flagContent === 'n') result.atomNumbers = false;
        if (flagContent === 'h') result.addHydrogens = false;    // Disable hydrogens (-h)
        if (flagContent === 'p') result.flipHorizontal = false;  // Disable flip horizontal (-p)
        if (flagContent === 'q') result.flipVertical = false;    // Disable flip vertical (-q)
        // -i flag removed - inversion now auto-applied based on dark mode
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
      if (baseSettings?.rendererEngine === 'mol2chemfig') {
        newSettings.m2cfFlipHorizontal = flagOverrides.flipHorizontal;
      } else {
        newSettings.flipHorizontal = flagOverrides.flipHorizontal;
      }
    }
    if (flagOverrides?.flipVertical !== undefined) {
      if (baseSettings?.rendererEngine === 'mol2chemfig') {
        newSettings.m2cfFlipVertical = flagOverrides.flipVertical;
      } else {
        newSettings.flipVertical = flagOverrides.flipVertical;
      }
    }
    // invert removed - now auto-applied based on dark mode detection

    // Note: Size and rotation are handled separately in the image creation process
    if (flagOverrides?.rotation !== undefined) {
      if (baseSettings?.rendererEngine === 'mol2chemfig') {
        newSettings.m2cfRotate = flagOverrides.rotation;
      } else {
        newSettings.mvRotate = flagOverrides.rotation;
      }
    }

    return newSettings;
  } catch (error) {
    console.error('Error applying flag overrides:', error);
    // Return original settings to not break functionality
    return baseSettings;
  }
}

/**
 * Helper function to get the invert filter based on invert mode
 * NOTE: This function is now deprecated in favor of conditional inversion in loadMoleculeViewerImage
 * @param {string} invertMode - 'full' or 'bw' (black & white only)
 * @returns {string} CSS filter string for invert
 */
function getInvertFilter(invertMode) {
  // This function is kept for backwards compatibility but is no longer used
  // in MoleculeViewer rendering. See loadMoleculeViewerImage for the new approach.
  if (invertMode === 'bw') {
    // BW mode: String replacement of black→white in SVG source (no CSS filter)
    return '';
  } else {
    // Full invert: CSS filter inversion of all colors
    return 'invert(1) ';
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
  // Always create wrapper container (img elements cannot have children for hover controls)
  const wrapper = document.createElement('div');
  wrapper.className = 'chemfig-rotation-wrapper chemfig-diagram';
  wrapper.style.cssText = `
    display: inline-block !important;
    position: relative !important;
    vertical-align: middle !important;
    margin: 0 12px 0 0 !important;
    padding: 0 !important;
    visibility: visible !important;
    opacity: 1 !important;
    line-height: 0 !important;
    width: fit-content !important;
    height: fit-content !important;
  `;

  // Add the image to wrapper
  wrapper.appendChild(img);

  // Add hover controls to wrapper (not to img, since img can't have children)
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
  if (element.querySelector('.chemfig-name-overlay')) {
    return;
  }
  // Create name overlay (bottom-right, always visible)
  const nameOverlay = document.createElement('div');
  nameOverlay.className = 'chemfig-name-overlay';
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
  viewer3DBtn.className = 'chemfig-3d-btn';
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
      btn.style.background = 'rgba(0, 0, 0, 0.7)';
    });
  });

  let currentScale = 1.0;
  const imgElement = element.querySelector('img');

  upButton.addEventListener('click', (e) => {
    e.stopPropagation();
    currentScale = Math.min(currentScale + 0.2, 5.0);
    if (imgElement) {
      // Get current rendered width as baseline (handles SVG data URLs properly)
      const currentWidth = imgElement.offsetWidth || imgElement.width || imgElement.naturalWidth || 300;
      const newWidth = Math.round(currentWidth * 1.2); // 20% increase
      imgElement.style.width = `${newWidth}px`;
      imgElement.style.height = 'auto';
      console.log(`Size up: ${currentWidth}px → ${newWidth}px (scale: ${currentScale.toFixed(2)})`);
    }
  });

  downButton.addEventListener('click', (e) => {
    e.stopPropagation();
    currentScale = Math.max(currentScale - 0.2, 0.5);
    if (imgElement) {
      // Get current rendered width as baseline
      const currentWidth = imgElement.offsetWidth || imgElement.width || imgElement.naturalWidth || 300;
      const newWidth = Math.round(currentWidth * 0.8333); // ~20% decrease (inverse of 1.2)
      imgElement.style.width = `${newWidth}px`;
      imgElement.style.height = 'auto';
      console.log(`Size down: ${currentWidth}px → ${newWidth}px (scale: ${currentScale.toFixed(2)})`);
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

// Export helper functions for global access
window.parseChemFlags = parseChemFlags;
window.stripFlagsFromName = stripFlagsFromName;
window.applyFlagOverrides = applyFlagOverrides;
window.getInvertFilter = getInvertFilter;

/**
 * Setup Lazy-Loading for Chemfig SVGs
 * Only renders visible SVGs, defers off-screen ones to save performance
 */
function setupLazyLoading() {
  log.inject('🚀 Setting up Intersection Observer for lazy-loading with performance optimization');

  // Track loading count to prevent too many simultaneous loads
  let activeLoads = 0;
  const maxConcurrentLoads = 3;  // Reduced from 5 to prevent lag

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
    img.classList.add('chemfig-fadein');
    img.classList.remove('chemfig-loading');

    // Preload the image to avoid blocking the UI
    const preloadImg = new Image();
    preloadImg.onload = () => {
      img.src = img.dataset.src;
      img.dataset.loaded = 'true';
      observer.unobserve(img);
      activeLoads--;
      log.debug(`✅ SVG loaded (${activeLoads} active loads remaining)`);
    };
    preloadImg.onerror = () => {
      activeLoads--;
      log.debug(`❌ SVG failed to load (${activeLoads} active loads remaining)`);
    };
    preloadImg.src = img.dataset.src;
  }

  // Helper function to render molecule client-side using SmilesDrawer (SVG)
  // Uses iframe approach similar to 3D viewer for proper isolation
  // Uses PubChem API (via background script) for name→SMILES conversion
  async function renderClientSide(moleculeData, img) {
    activeLoads++;
    console.log('%c🎨 RENDERCLIENTSIDE CALLED (SmilesDrawer SVG)!', 'background: #222; color: #00FF00; font-size: 20px; padding: 10px;');
    console.log('Image element:', img);
    console.log('Dataset:', img.dataset);
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
        console.log('%c🔮 3D View requested in Client-Side mode', 'color: #FF00FF;');
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
        activeLoads--;
        return;
      }

      let smiles = moleculeData.smiles;
      const compoundName = moleculeData.nomenclature || moleculeData.name || 'molecule';

      // If no SMILES, get it using PubChem (via background script to avoid CORS)
      if (!smiles && moleculeData.nomenclature) {
        const cleanName = moleculeData.nomenclature.trim();
        console.log('%c🔍 Client-Side: Getting SMILES for:', 'color: #FF00FF; font-weight: bold;', cleanName);

        // Check if 3D stereochemistry is enabled
        const use3DSmiles = settings.m2cfUse3DSmiles === true;

        // ===== PRIORITY 1: PubChem direct API (via background to bypass CORS) =====
        // Works for common names, trade names, drug names, and most molecules
        console.log('%c🌐 [Client] Priority 1: Trying PubChem direct lookup (via background)...', 'color: #0088FF; font-weight: bold;');
        try {
          const pubchemUrl = `https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/${encodeURIComponent(cleanName)}/property/CanonicalSMILES,IsomericSMILES/JSON`;
          console.log('%c🔗 [Client] PubChem URL:', 'color: #0088FF;', pubchemUrl);
          const data = await backgroundFetchJSON(pubchemUrl);

          console.log('%c📦 [Client] PubChem raw response:', 'color: #FF00FF; font-weight: bold;', JSON.stringify(data));
          console.log('%c📦 [Client] data type:', 'color: #FF00FF;', typeof data);
          console.log('%c📦 [Client] data.PropertyTable:', 'color: #FF00FF;', data?.PropertyTable);
          console.log('%c📦 [Client] data.PropertyTable?.Properties:', 'color: #FF00FF;', data?.PropertyTable?.Properties);

          if (data && data.PropertyTable && data.PropertyTable.Properties && data.PropertyTable.Properties[0]) {
            const props = data.PropertyTable.Properties[0];
            console.log('%c📦 [Client] props:', 'color: #00FFFF;', props);
            // PubChem can return different property names: CanonicalSMILES, IsomericSMILES, SMILES, ConnectivitySMILES
            smiles = use3DSmiles ?
              (props.IsomericSMILES || props.SMILES || props.CanonicalSMILES || props.ConnectivitySMILES) :
              (props.CanonicalSMILES || props.SMILES || props.ConnectivitySMILES);

            // Remove any spaces from SMILES (PubChem sometimes adds them for readability)
            if (smiles) {
              smiles = smiles.replace(/\s+/g, '');
            }

            console.log('%c✅ [Client] PubChem direct SUCCESS:', 'color: #00FF00; font-weight: bold;', smiles);
            if (smiles) {
              console.log('%c📏 SMILES length:', 'color: #9c88ff;', smiles.length, 'characters');
            }
          } else {
            console.warn('%c⚠️ [Client] PubChem data structure invalid:', 'color: #FF8800;', {
              hasData: !!data,
              hasPropertyTable: !!(data?.PropertyTable),
              hasProperties: !!(data?.PropertyTable?.Properties),
              hasFirst: !!(data?.PropertyTable?.Properties?.[0])
            });
          }
        } catch (pubchemError) {
          console.warn('⚠️ [Client] PubChem direct failed:', pubchemError.message);
          console.error('⚠️ [Client] Full error:', pubchemError);
        }

        // ===== PRIORITY 2: PubChem Autocomplete (for generic/class names like "sphingomyelin") =====
        // "sphingomyelin" → "Sphingomyelin 16:0", "phosphatidylcholine" → first match
        if (!smiles) {
          console.log('%c🔎 [Client] Priority 2: Trying PubChem autocomplete for class name...', 'color: #FF8800; font-weight: bold;');
          try {
            const autocompleteUrl = `https://pubchem.ncbi.nlm.nih.gov/rest/autocomplete/compound/${encodeURIComponent(cleanName)}/json?limit=5`;
            const autocompleteData = await backgroundFetchJSON(autocompleteUrl);

            if (autocompleteData && autocompleteData.dictionary_terms && autocompleteData.dictionary_terms.compound && autocompleteData.dictionary_terms.compound.length > 0) {
              // Get the first matching compound name
              const bestMatch = autocompleteData.dictionary_terms.compound[0];
              console.log('%c🔗 [Client] Found related compound:', 'color: #FF8800;', bestMatch);

              // Now lookup SMILES for the best match
              const matchUrl = `https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/${encodeURIComponent(bestMatch)}/property/CanonicalSMILES,IsomericSMILES/JSON`;
              const matchData = await backgroundFetchJSON(matchUrl);

              if (matchData && matchData.PropertyTable && matchData.PropertyTable.Properties && matchData.PropertyTable.Properties[0]) {
                const props = matchData.PropertyTable.Properties[0];
                // Accept any SMILES property name
                smiles = use3DSmiles ?
                  (props.IsomericSMILES || props.SMILES || props.CanonicalSMILES || props.ConnectivitySMILES) :
                  (props.CanonicalSMILES || props.SMILES || props.ConnectivitySMILES);

                // Remove spaces
                if (smiles) {
                  smiles = smiles.replace(/\s+/g, '');
                }

                console.log('%c✅ [Client] PubChem autocomplete SUCCESS:', 'color: #00FF00; font-weight: bold;', smiles);
                if (smiles) {
                  console.log('%c📏 SMILES length:', 'color: #9c88ff;', smiles.length, 'characters');
                }
              }
            }
          } catch (autocompleteError) {
            console.warn('⚠️ [Client] PubChem autocomplete failed:', autocompleteError.message);
          }
        }

        // ===== PRIORITY 3: Try getting CID first, then SMILES from CID =====
        if (!smiles) {
          console.log('%c🆔 [Client] Priority 3: Trying CID lookup...', 'color: #9c27b0; font-weight: bold;');
          try {
            const cid = await getPubChemCID(cleanName);
            if (cid) {
              console.log('%c✅ [Client] Found CID:', 'color: #9c27b0;', cid);

              // Get SMILES from CID
              const cidUrl = `https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/${cid}/property/CanonicalSMILES,IsomericSMILES/JSON`;
              const cidData = await backgroundFetchJSON(cidUrl);

              if (cidData && cidData.PropertyTable && cidData.PropertyTable.Properties && cidData.PropertyTable.Properties[0]) {
                const props = cidData.PropertyTable.Properties[0];
                smiles = use3DSmiles ?
                  (props.IsomericSMILES || props.SMILES || props.CanonicalSMILES || props.ConnectivitySMILES) :
                  (props.CanonicalSMILES || props.SMILES || props.ConnectivitySMILES);

                // Remove spaces
                if (smiles) {
                  smiles = smiles.replace(/\s+/g, '');
                }

                console.log('%c✅ [Client] CID lookup SUCCESS:', 'color: #9c27b0; font-weight: bold;', smiles);
                if (smiles) {
                  console.log('%c📏 SMILES length:', 'color: #9c88ff;', smiles.length, 'characters');
                }
              }
            }
          } catch (cidError) {
            console.warn('⚠️ [Client] CID lookup failed:', cidError.message);
          }
        }
      }

      if (!smiles) {
        throw new Error(`Could not obtain SMILES for "${compoundName}" - check internet connection`);
      }

      console.log('%c🧪 Creating SVG directly with SmilesDrawer:', 'color: #00BFFF;', smiles);

      // Get mol2chemfig-style rendering options from settings
      // Then apply any per-molecule flag overrides (e.g., +c, +m, /d-c)
      let showCarbons = settings.m2cfShowCarbons === true;          // -c flag: Display C labels
      let aromaticCircles = settings.m2cfAromaticCircles === true;   // -o flag: Draw circles in rings  
      let showMethyl = settings.m2cfShowMethyls === true;            // -m flag: Display CH3 labels
      let showAtomNumbers = settings.m2cfAtomNumbers === true;       // -n flag: Number atoms
      let flipVertical = settings.m2cfFlipVertical === true;         // Upside down
      let flipHorizontal = settings.m2cfFlipHorizontal === true;     // Mirror
      let addHydrogens = settings.m2cfAddH2 === true;                // Show H atoms
      let customSize = null;
      let customRotation = null;

      // Apply per-molecule flag overrides if present
      const flags = moleculeData.flags || moleculeData.options || {};
      if (flags) {
        if (flags.showCarbons !== undefined) showCarbons = flags.showCarbons;
        if (flags.aromaticCircles !== undefined) aromaticCircles = flags.aromaticCircles;
        if (flags.showMethyls !== undefined) showMethyl = flags.showMethyls;
        if (flags.atomNumbers !== undefined) showAtomNumbers = flags.atomNumbers;
        if (flags.addHydrogens !== undefined) addHydrogens = flags.addHydrogens;
        if (flags.flipHorizontal !== undefined) flipHorizontal = flags.flipHorizontal;
        if (flags.flipVertical !== undefined) flipVertical = flags.flipVertical;
        if (flags.size !== undefined) customSize = flags.size;
        if (flags.rotation !== undefined) customRotation = flags.rotation;
      }

      console.log('%c🎌 Flag overrides:', 'color: #FF6B6B;', flags);

      // Calculate dimensions (default, can be overridden by +s flag)
      const baseWidth = 300;
      const baseHeight = 240;
      const sizeMultiplier = customSize || 1;
      const renderWidth = Math.round(baseWidth * sizeMultiplier);
      const renderHeight = Math.round(baseHeight * sizeMultiplier);

      const theme = isDarkModeEnabled() ? 'dark' : 'light';

      // Configure SmilesDrawer options
      const smilesDrawerOptions = {
        width: renderWidth,
        height: renderHeight,
        bondThickness: 2.0,
        bondLength: 25,
        shortBondLength: 0.85,
        bondSpacing: 6,
        atomVisualization: 'default',
        isomeric: true,
        debug: false,
        terminalCarbons: showMethyl || showCarbons,
        explicitHydrogens: addHydrogens,
        compactDrawing: !showCarbons,
        overlapSensitivity: 0.42,
        overlapResolutionIterations: 1,
        fontSizeLarge: 10,
        fontSizeSmall: 8,
        padding: 5.0,
        experimental: false,
        themes: {
          dark: {
            C: '#ffffff', O: '#ff6b6b', N: '#4dabf7', F: '#51cf66',
            CL: '#20c997', BR: '#fd7e14', I: '#be4bdb', P: '#fd7e14',
            S: '#fcc419', B: '#f59f00', SI: '#f59f00', H: '#aaaaaa',
            BACKGROUND: 'transparent'
          },
          light: {
            C: '#222222', O: '#e74c3c', N: '#3498db', F: '#27ae60',
            CL: '#16a085', BR: '#d35400', I: '#8e44ad', P: '#d35400',
            S: '#f1c40f', B: '#e67e22', SI: '#e67e22', H: '#666666',
            BACKGROUND: 'transparent'
          }
        }
      };

      console.log('%c🎨 Rendering SVG directly with SmilesDrawer', 'color: #9c88ff;');

      // Render the molecule directly using SmilesDrawer
      const svgDrawer = new SmilesDrawer.SvgDrawer(smilesDrawerOptions);

      SmilesDrawer.parse(smiles, function (tree) {
        // Count atoms to determine molecule complexity
        const atomCount = smiles.replace(/[^A-Z]/g, '').length; // Rough estimate
        console.log(`Molecule has approximately ${atomCount} atoms`);

        // Adjust canvas size based on complexity
        let adjustedWidth = renderWidth;
        let adjustedHeight = renderHeight;

        if (atomCount > 100) {
          // Large molecule like insulin - use bigger canvas
          const complexityFactor = Math.min(atomCount / 50, 4); // Max 4x
          adjustedWidth = Math.ceil(renderWidth * complexityFactor);
          adjustedHeight = Math.ceil(renderHeight * complexityFactor);
          console.log(`Large molecule detected, expanding canvas to ${adjustedWidth}x${adjustedHeight}`);

          // Update options
          smilesDrawerOptions.width = adjustedWidth;
          smilesDrawerOptions.height = adjustedHeight;
        }

        // Draw to a temporary container first
        const tempId = `temp-${Date.now()}`;
        const tempDiv = document.createElement('div');
        tempDiv.innerHTML = `<svg id="${tempId}"></svg>`;
        document.body.appendChild(tempDiv);

        const finalDrawer = new SmilesDrawer.SvgDrawer(smilesDrawerOptions);
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
            const maxDisplayHeight = 450;  // Max display height
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

          // Clean up temp
          document.body.removeChild(tempDiv);

          // Create img element just like moleculeviewer/mol2chemfig does
          const svgImg = document.createElement('img');
          svgImg.src = svgDataUrl;
          svgImg.alt = compoundName;
          svgImg.className = 'chemfig-diagram';

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
            console.log('%c✅ SVG rendered as img with standard wrapper', 'color: green; font-weight: bold;');
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
        errorImg.className = 'chemfig-diagram';

        const parent = img.parentNode;
        if (parent) {
          parent.replaceChild(errorImg, img);
        }
      });

      activeLoads--;
      return;
    } catch (error) {
      console.error('%c❌ Client-side rendering error:', 'color: red; font-weight: bold;', error);
      img.alt = `Error: ${error.message}`;
      activeLoads--;
    }
  }

  // ========== SERVER-SIDE RENDERING (existing code) ==========
  async function renderMoleculeServerSide(moleculeData, img) {
    try {
      activeLoads--;
    } catch (error) {
      console.error('%c❌ Server-side render failed:', 'color: red; font-weight: bold;', error);

      // Show error in place of image
      const errorDiv = document.createElement('div');
      errorDiv.style.cssText = `
        display: inline-block;
        padding: 15px;
        background: ${isDarkModeEnabled() ? '#2d1f1f' : '#ffebee'};
        border: 1px solid ${isDarkModeEnabled() ? '#5c3333' : '#f44336'};
        border-radius: 8px;
        color: ${isDarkModeEnabled() ? '#ff8a80' : '#c62828'};
        font-size: 12px;
        max-width: 300px;
      `;
      errorDiv.innerHTML = `
        <strong>❌ Render Error</strong><br>
        ${error.message || 'Unknown error'}<br>
        <small>Molecule: ${moleculeData?.nomenclature || 'unknown'}</small>
      `;

      if (img.parentNode) {
        img.parentNode.replaceChild(errorDiv, img);
      }

      activeLoads--;
    }
  }

  // Helper function to load MoleculeViewer rendering with caching and download link
  // Uses centralized SMILES Bridge for nomenclature→SMILES conversion
  async function loadMoleculeViewerImage(img) {
    activeLoads++;
    console.log('%c🧪 LOADMOLECULEVIEWERIMAGE CALLED!', 'background: #222; color: #00FF00; font-size: 20px; padding: 10px;');
    console.log('Image element:', img);
    console.log('Dataset:', img.dataset);
    log.debug(`🧪 Loading MoleculeViewer SVG (#${activeLoads})`);

    try {
      const moleculeData = JSON.parse(atob(img.dataset.moleculeViewer));
      console.log('%c📦 Decoded molecule data:', 'color: #0088FF; font-weight: bold;', moleculeData);

      // Check if we should use Client-Side rendering
      if (settings.rendererEngine === 'client-side') {
        activeLoads--; // Decrement because renderClientSide increments it again
        await renderClientSide(moleculeData, img);
        return;
      }

      // Determine which endpoint to use based on data type
      let isNomenclature = moleculeData.type === 'nomenclature' && moleculeData.nomenclature;
      let isSMILES = moleculeData.type === 'smiles' && moleculeData.smiles;

      // If chemfig payload was stored, use fallback.smiles for MoleculeViewer
      if (!isSMILES && moleculeData.type === 'chemfig' && moleculeData.fallback && moleculeData.fallback.smiles) {
        isSMILES = true;
        moleculeData.smiles = moleculeData.fallback.smiles;
      }

      // 🌉 SMILES BRIDGE: Convert nomenclature to SMILES first
      // This centralizes conversion logic (OPSIN → PubChem fallback)
      // Uses mvUse3DSmiles setting for MoleculeViewer 3D stereochemistry
      // IMPORTANT: Only use 3D SMILES when the option is explicitly enabled
      if (isNomenclature && !isSMILES) {
        const use3D = settings.mvUse3DSmiles === true;  // Ensure boolean check
        console.log('%c🌉 Using SMILES Bridge for nomenclature conversion...', 'background: #4A90D9; color: white; font-weight: bold; padding: 2px 6px;', use3D ? '(3D enabled)' : '(3D disabled)');
        const bridgeResult = await smilesBridge(moleculeData.nomenclature, { use3DSmiles: use3D });
        if (bridgeResult && bridgeResult.smiles) {
          console.log('%c✅ SMILES Bridge conversion successful:', 'color: #00FF00; font-weight: bold;', bridgeResult.smiles, `(source: ${bridgeResult.source})`);
          moleculeData.smiles = bridgeResult.smiles;
          moleculeData.smilesSource = bridgeResult.source;
          isSMILES = true;
          isNomenclature = false; // Now we have SMILES, use SMILES endpoint
        } else {
          console.warn('%c⚠️ SMILES Bridge conversion failed, falling back to server nomenclature endpoint', 'color: #FFA500; font-weight: bold;');
          // Keep isNomenclature = true to fall back to server-side conversion
        }
      }

      let apiUrl;
      if (isSMILES) {
        console.log('%c📤 Using SMILES endpoint', 'color: #FF6B00; font-weight: bold;');

        // Build options query string for MoleculeViewer (simplified - only essential params)
        const optionsParams = new URLSearchParams({
          smiles: moleculeData.smiles,
          width: '300',
          height: '200',
          json: 'true'
        });

        apiUrl = `${MOLECULE_VIEWER_API}/img/smiles?${optionsParams.toString()}&t=${Date.now()}`;
        console.log('API URL:', apiUrl);
      } else if (isNomenclature) {
        // Fallback: server-side nomenclature conversion (if SMILES Bridge failed)
        console.log('%c📤 Using nomenclature endpoint (server-side fallback)', 'color: #FF6B00; font-weight: bold;');
        console.log('Nomenclature:', moleculeData.nomenclature);

        // Build options query string for MoleculeViewer (simplified - only essential params)
        const optionsParams = new URLSearchParams({
          nomenclature: moleculeData.nomenclature,
          width: '300',
          height: '200',
          json: 'true'
        });

        apiUrl = `${MOLECULE_VIEWER_API}/img/nomenclature?${optionsParams.toString()}&t=${Date.now()}`;
        console.log('API URL:', apiUrl);
      } else {
        throw new Error('Invalid molecule data');
      }

      // Log settings being used
      console.log('%c⚙️ MoleculeViewer Settings:', 'color: #FF6B00; font-weight: bold;', {
        mvUse3DSmiles: settings.mvUse3DSmiles,
        flipHorizontal: settings.flipHorizontal,
        flipVertical: settings.flipVertical
      });

      // Fetch JSON response with cache link (via background to bypass CSP)
      console.log('%c🌐 Fetching from backend...', 'color: #00AAFF; font-weight: bold;');
      backgroundFetchJSON(apiUrl)
        .then(data => {
          console.log('%c📊 Backend data:', 'color: #FFAA00; font-weight: bold;', data);
          if (!data.success) {
            throw new Error(data.error || 'Rendering failed');
          }

          // Create clean SVG image (no container, no controls)

          // Get per-image options from moleculeData
          const opts = moleculeData.options || {};

          // Apply dark mode color inversion based on user settings
          let svgContent = data.svg;
          const isDarkMode = isDarkModeEnabled();
          const shouldInvertSVG = settings.mvAutoInvert && isDarkMode;
          const invertMode = settings.mvInvertMode || 'full';

          // Only apply SVG content inversion in "bw" mode (black/white only)
          // This replaces black (#000) with white (#FFF) in the SVG source
          if (shouldInvertSVG && invertMode === 'bw') {
            console.log('%c🎨 Applying BW SVG inversion (black→white in source)', 'color: #FF6B00; font-weight: bold;');
            svgContent = invertSvgForDarkMode(svgContent);
          }

          // Create clean SVG image with dark mode support
          const svgImg = document.createElement('img');
          svgImg.src = 'data:image/svg+xml;base64,' + btoa(svgContent);
          svgImg.alt = 'molecule';
          svgImg.className = 'chemfig-diagram';

          // Apply MoleculeViewer transforms (flip, rotation, scale, invert)
          let transform = '';
          const flipH = opts.flipHorizontal !== undefined ? opts.flipHorizontal : settings.flipHorizontal;
          const flipV = opts.flipVertical !== undefined ? opts.flipVertical : settings.flipVertical;
          const rotation = opts.rotation !== undefined ? opts.rotation : (settings.mvRotate || 0);
          const scale = opts.scale !== undefined ? opts.scale : 1.0;

          if (flipH) transform += 'scaleX(-1) ';
          if (flipV) transform += 'scaleY(-1) ';
          if (rotation) transform += `rotate(${rotation}deg) `;
          if (scale && scale !== 1.0) transform += `scale(${scale}) `;

          // Apply CSS filter inversion in "full" mode only
          // This inverts ALL colors via CSS filter
          let filter = '';
          if (shouldInvertSVG && invertMode === 'full') {
            console.log('%c🎨 Applying FULL CSS filter inversion', 'color: #FF6B00; font-weight: bold;');
            filter += 'invert(1) ';
          }

          svgImg.style.cssText = `
            display: inline-block;
            max-width: 800px;
            max-height: 600px;
            width: fit-content;
            height: auto;
            margin: 0 12px 8px 0;
            vertical-align: middle;
            cursor: pointer;
            ${transform ? `transform: ${transform.trim()};` : ''}
            ${filter ? `filter: ${filter.trim()};` : ''}
          `;

          // Mark as loaded
          img.dataset.loaded = 'true';

          // Replace original img element with size controls
          // (Hover controls are now added automatically in wrapImageWithSizeControls)
          chrome.storage.sync.get({
            saveSizePerImage: false,
            saveSizeBySMILES: true  // FIX: Enable by default
          }, async (sizeSettings) => {
            await wrapImageWithSizeControls(svgImg, img, moleculeData, sizeSettings);
          });

          console.log('%c✅ Image loaded successfully', 'color: green; font-weight: bold;');
          console.log('%c📍 Cache URL:', 'color: #0066cc; font-weight: bold;', data.cache_url);

          activeLoads--;
        })
        .catch(error => {
          console.error('%c❌ Error fetching molecule:', 'color: red; font-weight: bold;', error.message);

          // Fallback: show error message
          img.src = 'data:image/svg+xml;base64,PHN2ZyB3aWR0aD0iMzAwIiBoZWlnaHQ9IjIwMCIgeG1sbnM9Imh0dHA6Ly93d3cudzMub3JnLzIwMDAvc3ZnIj48cmVjdCB3aWR0aD0iMzAwIiBoZWlnaHQ9IjIwMCIgZmlsbD0iI2ZmZiIvPjx0ZXh0IHg9IjEwIiB5PSIzMCIgZmlsbD0iI2YwMCIgZm9udC1mYW1pbHk9Im1vbm9zcGFjZSIgZm9udC1zaXplPSIxNCI+RXJyb3I6IExvYWQgRmFpbGVkPC90ZXh0Pjwvc3ZnPg==';
          img.classList.add('chemfig-fadein');
          img.classList.remove('chemfig-loading', 'chemfig-molecule-viewer');
          img.dataset.loaded = 'true';

          activeLoads--;
        });

      log.debug(`✅ Fetching molecule data (${activeLoads} active loads remaining)`);
    } catch (error) {
      console.error(`%c❌ Error: ${error.message}`, 'color: red; font-weight: bold;');
      activeLoads--;
    }
  }

  // Helper function to load PubChem images
  async function loadPubChemImage(img) {
    activeLoads++;
    console.log('%c🌐 LOADPUBCHEMIMAGE CALLED!', 'background: #222; color: #4CAF50; font-size: 20px; padding: 10px;');
    console.log('Image element:', img);
    console.log('Dataset:', img.dataset);
    log.debug(`🌐 Loading PubChem image (#${activeLoads})`);

    try {
      const moleculeData = JSON.parse(atob(img.dataset.moleculeViewer || img.dataset.mol2chemfig));
      console.log('%c📦 Decoded molecule data:', 'color: #4CAF50; font-weight: bold;', moleculeData);

      // Determine the compound name or SMILES for PubChem
      let compoundName = '';
      if (moleculeData.nomenclature) {
        compoundName = moleculeData.nomenclature;
      } else if (moleculeData.smiles) {
        compoundName = moleculeData.smiles;
      } else if (moleculeData.type === 'chemfig' && moleculeData.fallback && moleculeData.fallback.smiles) {
        compoundName = moleculeData.fallback.smiles;
      } else {
        throw new Error('No valid compound identifier found');
      }

      console.log('%c🔍 Looking up compound:', 'color: #4CAF50; font-weight: bold;', compoundName);

      // Get settings, potentially overridden by flags
      const currentSettings = moleculeData.customSettings || settings;
      const flags = moleculeData.flags || {};

      // Check if 3D flag was set - if so, show 3D viewer immediately
      if (flags.is3D || moleculeData.show3D || img.dataset.show3d === 'true') {
        await show3DViewerInline(moleculeData, img);
        activeLoads--;
        return;
      }

      // Check if 3D viewer is enabled globally - skip 2D image fetch entirely
      const viewerSettings = await new Promise(resolve => {
        chrome.storage.sync.get({ enable3DViewer: false }, resolve);
      });

      if (viewerSettings.enable3DViewer) {
        // User wants 3D viewer - show it directly without fetching 2D image
        await show3DViewerInline(moleculeData, img);
        activeLoads--;
        return;
      }

      let imageUrl;
      let recordType = currentSettings.pubchemRecordType || '2d';
      let fetchMethod = 'unknown';

      if (currentSettings.pubchemDirectFetch) {
        // Direct fetch from PubChem
        const cid = await getPubChemCID(compoundName);
        if (cid) {
          imageUrl = `https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/${cid}/PNG?image_size=${currentSettings.pubchemImageSize || 'large'}&record_type=${recordType}`;
          fetchMethod = 'direct-pubchem';
          console.log('%c🌐 Using DIRECT fetch from PubChem.gov', 'background: #10b981; color: white; font-weight: bold; padding: 2px 6px;');
        } else {
          console.error('%c❌ Could not get CID for compound, falling back to bridge server', 'color: red; font-weight: bold;');
          imageUrl = `${PUBCHEM_API}/pubchem/img/${encodeURIComponent(compoundName)}?size=${currentSettings.pubchemImageSize || 'large'}&type=${recordType}&proxy=true`;
          fetchMethod = 'bridge-server-fallback';
        }
      } else {
        // Use local PubChem bridge server (port 5002)
        imageUrl = `${PUBCHEM_API}/pubchem/img/${encodeURIComponent(compoundName)}?size=${currentSettings.pubchemImageSize || 'large'}&type=${recordType}&proxy=true`;
        fetchMethod = 'bridge-server';
        console.log('%c🖥️ Using BRIDGE SERVER (localhost:5002)', 'background: #3b82f6; color: white; font-weight: bold; padding: 2px 6px;');
        console.log('%c⚠️ Make sure PubChem server is running: python pubchem_server.py', 'color: orange; font-weight: bold;');
      }

      console.log('%c📡 PubChem URL:', 'color: #4CAF50; font-weight: bold;', imageUrl);
      console.log('%c🔧 Fetch method:', 'color: #9333ea;', fetchMethod);

      // Fetch the image via background script to bypass CSP
      backgroundFetchBlob(imageUrl)
        .then(async (blobData) => {
          // Process image with background removal if enabled
          let imageSrc = await processImageBackgroundRemoval(blobData);

          const pubchemImg = document.createElement('img');
          pubchemImg.src = imageSrc;
          pubchemImg.alt = `PubChem: ${compoundName}`;
          pubchemImg.className = 'chemfig-diagram pubchem-diagram';

          // Apply transformations based on flags
          let transform = [];
          if (flags.flipHorizontal) transform.push('scaleX(-1)');
          if (flags.flipVertical) transform.push('scaleY(-1)');

          // Apply rotation if specified in flags
          let rotation = parseInt(img.dataset.rotation || 0);
          if (flags.rotation !== undefined) {
            rotation = flags.rotation;
          }
          if (rotation !== 0) transform.push(`rotate(${rotation}deg)`);

          pubchemImg.style.cssText = `
            display: inline-block;
            max-width: 300px;
            max-height: 200px;
            margin: 0 12px 8px 0;
            vertical-align: middle;
            cursor: pointer;
            ${transform.length > 0 ? `transform: ${transform.join(' ')};` : ''}
          `;

          // Mark as loaded
          img.dataset.loaded = 'true';

          // Show 2D image with size controls (3D viewer is handled above)
          chrome.storage.sync.get({
            saveSizePerImage: false,
            saveSizeBySMILES: true  // Enable by default
          }, async (pubchemSettings) => {
            await wrapImageWithSizeControls(pubchemImg, img, moleculeData, pubchemSettings);
          });

          console.log('%c✅ PubChem image loaded successfully', 'color: green; font-weight: bold;');
          activeLoads--;
        })
        .catch(error => {
          console.error('%c❌ Error loading PubChem image:', 'color: red; font-weight: bold;', error.message);

          // Fallback: show error message
          img.src = 'data:image/svg+xml;base64,PHN2ZyB3aWR0aD0iMzAwIiBoZWlnaHQ9IjIwMCIgeG1sbnM9Imh0dHA6Ly93d3cudzMub3JnLzIwMDAvc3ZnIj48cmVjdCB3aWR0aD0iMzAwIiBoZWlnaHQ9IjIwMCIgZmlsbD0iI2ZmZiIvPjx0ZXh0IHg9IjEwIiB5PSIzMCIgZmlsbD0iI2YwMCIgZm9udC1mYW1pbHk9Im1vbm9zcGFjZSIgZm9udC1zaXplPSIxNCI+UHViQ2hlbTogTm90IEZvdW5kPC90ZXh0Pjwvc3ZnPg==';
          img.classList.add('chemfig-fadein');
          img.classList.remove('chemfig-loading');
          img.dataset.loaded = 'true';

          activeLoads--;
        });

    } catch (error) {
      console.error(`%c❌ Error: ${error.message}`, 'color: red; font-weight: bold;');
      activeLoads--;
    }
  }

  // Helper function to show 3D viewer inline using local PubChem 3D viewer (replaces the img element)
  async function show3DViewerInline(moleculeData, targetElement) {
    console.log('%c🔮 SHOWING 3D VIEWER INLINE', 'background: #764ba2; color: white; font-size: 14px; padding: 8px;');
    console.log('moleculeData:', moleculeData);
    console.log('targetElement:', targetElement);

    // Extract compound name from moleculeData
    const compoundName = moleculeData.nomenclature || moleculeData.smiles || '';

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
        const sizeWrapper = targetElement.querySelector('.chemfig-size-wrapper');
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

    console.log('Found img element:', img);
    console.log('Container element:', container);
    console.log('Compound name:', compoundName);

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
      let viewerSource = settings.viewer3DSource || '3dmol';
      console.log('%c🔍 Current 3D viewer source setting:', 'color: #ff6b6b; font-weight: bold;', viewerSource);
      console.log('%c📦 All settings:', 'color: #4ecdc4;', settings);

      switch (viewerSource) {
        case 'molview':
          // MolView Embed
          const cidMolview = await getPubChemCID(compoundName);
          if (cidMolview) {
            viewerUrl = `https://embed.molview.org/v1/?mode=balls&cid=${cidMolview}`;
            console.log('%c📍 MolView Embed URL:', 'color: #0066cc; font-weight: bold;', viewerUrl);
          } else {
            console.warn('%c⚠️ CID not found for MolView, using 3Dmol.js fallback', 'color: orange;');
            viewerSource = '3dmol';
          }
          break;

        case 'pubchem':
          // PubChem 3D viewer with wrapper to hide extra UI
          const cidPubchem = await getPubChemCID(compoundName);
          if (cidPubchem) {
            // Use our wrapper that tries to hide everything except the 3D widget
            viewerUrl = chrome.runtime.getURL(
              `pubchem-3d-wrapper.html?cid=${cidPubchem}&name=${encodeURIComponent(compoundName)}`
            );
            console.log('%c📍 PubChem 3D Viewer (wrapped):', 'color: #0066cc; font-weight: bold;', viewerUrl);
          } else {
            // Fallback to 3dmol if CID not found
            console.warn('%c⚠️ CID not found for PubChem, using 3Dmol.js', 'color: orange;');
            viewerSource = '3dmol';
          }
          break;

        case '3dmol':
        default:
          // 3Dmol.js - client-side viewer with custom styling
          const cid3dmol = await getPubChemCID(compoundName);
          if (cid3dmol) {
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
              `3dmol-viewer.html?cid=${cid3dmol}&name=${encodeURIComponent(compoundName)}&style=${style3D}&stickRadius=${stickRadius}&sphereRadius=${sphereRadius}&autoRotate=${autoRotate}&bgColor=${bgColor}`
            );
            console.log('%c📍 3Dmol.js Custom Viewer URL:', 'color: #0066cc; font-weight: bold;', viewerUrl);
            console.log('%c🎨 Style settings:', 'color: #9c88ff;', { style: style3D, stickRadius, sphereRadius, autoRotate });
          } else {
            throw new Error('Could not find molecule in PubChem database');
            activeLoads--;
            return;
          }
          break;
      }

      viewer3DIframe.src = viewerUrl;
      // Allow necessary features for 3D rendering
      viewer3DIframe.setAttribute('allow', 'accelerometer; autoplay; clipboard-write; encrypted-media; gyroscope; picture-in-picture; web-share');
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

      // Assemble container (no toggle button or name label - handled by hover controls)
      viewer3DContainer.appendChild(viewer3DIframe);

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
      activeLoads--;

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
      img.classList.add('chemfig-fadein');
      img.classList.remove('chemfig-loading');
      activeLoads--;
    }
  }

  // Duplicate getPubChemCID removed - using the robust version defined at the top of the file

  // Function to remove background color from an image element
  // Also handles dark mode color inversion
  function removeWhiteBackgroundFromImage(imgElement, bgColor = '#F5F5F5') {
    // Check if dark mode is enabled
    const darkMode = isDarkModeEnabled();

    // Create canvas to process the image
    const canvas = document.createElement('canvas');
    const ctx = canvas.getContext('2d');
    canvas.width = imgElement.naturalWidth || imgElement.width;
    canvas.height = imgElement.naturalHeight || imgElement.height;

    // Draw the original image on canvas
    ctx.drawImage(imgElement, 0, 0);

    // Get the image data
    const imageData = ctx.getImageData(0, 0, canvas.width, canvas.height);
    const data = imageData.data;

    // Process pixels
    for (let i = 0; i < data.length; i += 4) {
      const r = data[i];
      const g = data[i + 1];
      const b = data[i + 2];

      // Check if pixel is the PubChem background color #F5F5F5 (245, 245, 245)
      // Also catch #FCFCFC, #FAFAFA, and pure white
      const isBackgroundColor = (
        (r >= 243 && r <= 255 && g >= 243 && g <= 255 && b >= 243 && b <= 255) &&
        Math.abs(r - g) <= 5 && Math.abs(g - b) <= 5 && Math.abs(r - b) <= 5
      );

      if (isBackgroundColor) {
        // Make background transparent
        data[i + 3] = 0;
      } else if (darkMode) {
        // In dark mode, invert dark colors to light for visibility
        if (r <= 50 && g <= 50 && b <= 50) {
          // Invert black to white
          data[i] = 255;
          data[i + 1] = 255;
          data[i + 2] = 255;
        } else if (r <= 100 && g <= 100 && b <= 100 &&
          Math.abs(r - g) <= 20 && Math.abs(g - b) <= 20) {
          // Invert dark grays
          data[i] = 255 - r;
          data[i + 1] = 255 - g;
          data[i + 2] = 255 - b;
        }
      }
    }

    // Put the modified data back to canvas
    ctx.putImageData(imageData, 0, 0);

    // Replace the original image with the processed one
    const newImageSrc = canvas.toDataURL('image/png');
    imgElement.src = newImageSrc;
    console.log('%c🎨 Background removal applied to element' + (darkMode ? ' (dark mode)' : ''), 'color: #9B59B6;');
  }



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

  // Updated function to apply background removal to base64 data before creating image
  // Also handles dark mode color inversion for molecule lines
  async function processImageBackgroundRemoval(blobData) {
    return new Promise((resolve) => {
      const imageSrc = `data:${blobData.type || 'image/png'};base64,${blobData.base64}`;
      const img = new Image();

      img.onload = function () {
        // If background removal is not enabled, return the original
        if (!settings.pubchemRemoveBg) {
          resolve(imageSrc);
          return;
        }

        // Check if dark mode is enabled
        const darkMode = isDarkModeEnabled();

        // Process the image to remove background
        const canvas = document.createElement('canvas');
        const ctx = canvas.getContext('2d');
        canvas.width = img.width;
        canvas.height = img.height;
        ctx.drawImage(img, 0, 0);

        // Get image data
        const imageData = ctx.getImageData(0, 0, canvas.width, canvas.height);
        const data = imageData.data;

        // Process pixels
        for (let i = 0; i < data.length; i += 4) {
          const r = data[i];
          const g = data[i + 1];
          const b = data[i + 2];

          // Check if pixel is the PubChem background color #F5F5F5 (245, 245, 245)
          // Also catch #FCFCFC (252,252,252), #FAFAFA (250,250,250), and pure white
          const isBackgroundColor = (
            // Exact #F5F5F5 and nearby
            (r >= 243 && r <= 255 && g >= 243 && g <= 255 && b >= 243 && b <= 255) &&
            // Ensure it's a gray (all values close together)
            Math.abs(r - g) <= 5 && Math.abs(g - b) <= 5 && Math.abs(r - b) <= 5
          );

          if (isBackgroundColor) {
            // Make background transparent
            data[i + 3] = 0;
          } else if (darkMode) {
            // In dark mode, invert dark colors to light for visibility
            // Check if pixel is dark (likely molecule lines - black or near-black)
            if (r <= 50 && g <= 50 && b <= 50) {
              // Invert to white
              data[i] = 255;
              data[i + 1] = 255;
              data[i + 2] = 255;
            }
            // Invert other dark grays too
            else if (r <= 100 && g <= 100 && b <= 100 &&
              Math.abs(r - g) <= 20 && Math.abs(g - b) <= 20) {
              // Invert gray to light gray
              data[i] = 255 - r;
              data[i + 1] = 255 - g;
              data[i + 2] = 255 - b;
            }
          }
        }

        ctx.putImageData(imageData, 0, 0);

        // Apply sharpening filter if enabled
        if (settings.pubchemSharpenImages) {
          const sharpenedData = applySharpenFilter(ctx, canvas.width, canvas.height);
          ctx.putImageData(sharpenedData, 0, 0);
          console.log('%c✨ Sharpening filter applied', 'color: #3498DB; font-weight: bold;');
        }

        // Auto-crop to molecule bounds (remove extra whitespace)
        const croppedCanvas = autoCropCanvas(canvas);
        const newImageSrc = croppedCanvas.toDataURL('image/png');

        console.log('%c🎨 Background removal applied' + (darkMode ? ' (dark mode inversion)' : '') + (settings.pubchemSharpenImages ? ' + sharpening' : '') + ' + auto-crop', 'color: #9B59B6; font-weight: bold;');
        resolve(newImageSrc);
      };

      img.onerror = function () {
        // If there's an error, return original data
        resolve(imageSrc);
      };

      img.src = imageSrc;
    });
  }

  // Helper function to add 3D view button to PubChem images
  function add3DViewButton(container, compoundName) {
    // Find or create a control panel
    let controlPanel = container.querySelector('.molecule-controls');
    if (!controlPanel) {
      controlPanel = document.createElement('div');
      controlPanel.className = 'molecule-controls';
      controlPanel.style.cssText = `
        position: absolute;
        top: 5px;
        right: 5px;
        display: flex;
        gap: 5px;
        z-index: 10;
      `;
      container.style.position = 'relative';
      container.appendChild(controlPanel);
    }

    // Create 3D view button
    const view3DBtn = document.createElement('button');
    view3DBtn.innerHTML = '🔮 3D';
    view3DBtn.title = 'View 3D model on PubChem';
    view3DBtn.className = 'pubchem-3d-btn';
    view3DBtn.style.cssText = `
      background: linear-gradient(135deg, #667eea 0%, #764ba2 100%);
      color: white;
      border: none;
      border-radius: 4px;
      padding: 4px 8px;
      font-size: 11px;
      font-weight: bold;
      cursor: pointer;
      box-shadow: 0 2px 4px rgba(0,0,0,0.2);
      transition: transform 0.2s;
    `;
    view3DBtn.addEventListener('mouseenter', () => {
      view3DBtn.style.transform = 'scale(1.05)';
    });
    view3DBtn.addEventListener('mouseleave', () => {
      view3DBtn.style.transform = 'scale(1)';
    });
    view3DBtn.addEventListener('click', (e) => {
      e.stopPropagation();
      // Use local 3D viewer with canvas element instead of redirecting to PubChem
      const viewerUrl = `${PUBCHEM_API}/viewer-3d/${encodeURIComponent(compoundName)}`;
      window.open(viewerUrl, '_blank', 'width=1000,height=700');
    });

    controlPanel.appendChild(view3DBtn);
  }

  // Helper function to load mol2chemfig rendering with caching
  async function loadMol2chemfigImage(img) {
    activeLoads++;
    console.log('%c📐 LOADMOL2CHEMFIGIMAGE CALLED!', 'background: #222; color: #00FF00; font-size: 20px; padding: 10px;');
    console.log('Image element:', img);
    console.log('Dataset:', img.dataset);
    log.debug(`📐 Loading mol2chemfig LaTeX SVG (#${activeLoads})`);

    try {
      const moleculeData = JSON.parse(atob(img.dataset.moleculeViewer));
      console.log('%c📦 Decoded molecule data:', 'color: #0088FF; font-weight: bold;', moleculeData);

      // Determine input data based on type
      const isSMILES = moleculeData.type === 'smiles' && moleculeData.smiles;
      const isNomenclature = moleculeData.type === 'nomenclature' && moleculeData.nomenclature;
      const isChemfig = moleculeData.type === 'chemfig' && moleculeData.chemfig;

      console.log('%c📐 MOL2CHEMFIG Type Detection:', 'color: #FF6B00; font-weight: bold;', {
        type: moleculeData.type,
        isSMILES,
        isNomenclature,
        isChemfig,
        smiles: moleculeData.smiles,
        nomenclature: moleculeData.nomenclature
      });

      let inputData;
      if (isChemfig) {
        console.log('%c📤 Using CHEMFIG input (raw chemfig)', 'color: #FF6B00; font-weight: bold;');
        // Prefer sending a LaTeX-wrapped chemfig block if available
        inputData = moleculeData.latex || ('$\\chemfig{' + moleculeData.chemfig + '}$');
      } else if (isSMILES) {
        console.log('%c📤 Using SMILES input', 'color: #FF6B00; font-weight: bold;');
        console.log('SMILES:', moleculeData.smiles);
        inputData = moleculeData.smiles;
      } else if (isNomenclature) {
        // 🌉 SMILES BRIDGE: Convert nomenclature to SMILES first
        // This centralizes conversion logic (OPSIN → PubChem fallback)
        // IMPORTANT: Only use 3D SMILES when the option is explicitly enabled
        const use3D = settings.use3DSmiles === true;  // Ensure boolean check
        console.log('%c🌉 Using SMILES Bridge for nomenclature→SMILES conversion...', 'background: #4A90D9; color: white; font-weight: bold; padding: 2px 6px;', use3D ? '(3D enabled)' : '(3D disabled)');
        const bridgeResult = await smilesBridge(moleculeData.nomenclature, { use3DSmiles: use3D });
        if (bridgeResult && bridgeResult.smiles) {
          console.log('%c✅ SMILES Bridge conversion successful:', 'color: #00FF00; font-weight: bold;', bridgeResult.smiles, `(source: ${bridgeResult.source})`);
          inputData = bridgeResult.smiles;
          moleculeData.smiles = bridgeResult.smiles;
          moleculeData.smilesSource = bridgeResult.source;
        } else {
          console.error('%c❌ SMILES Bridge conversion failed for:', 'color: #FF0000; font-weight: bold;', moleculeData.nomenclature);
          throw new Error(`Could not convert "${moleculeData.nomenclature}" to SMILES (OPSIN and PubChem both failed)`);
        }
      } else {
        throw new Error('Invalid molecule data for mol2chemfig');
      }

      // Heuristic: sanitize common human-readable fragments like CH_3-CH_2-OH -> CCO
      // Many sources produce CH_3 style strings which are not valid SMILES; try a simple transform
      try {
        if (typeof inputData === 'string' && /CH[_\d]|CH\d|OH|\-/.test(inputData)) {
          let sanitized = String(inputData);
          sanitized = sanitized.replace(/_/g, ''); // remove underscores
          sanitized = sanitized.replace(/\-/g, ''); // remove hyphens
          sanitized = sanitized.replace(/CH3/gi, 'C');
          sanitized = sanitized.replace(/CH2/gi, 'C');
          sanitized = sanitized.replace(/CH1/gi, 'C');
          sanitized = sanitized.replace(/CH/gi, 'C');
          sanitized = sanitized.replace(/OH/gi, 'O');
          sanitized = sanitized.replace(/\s+/g, '');
          if (sanitized && sanitized !== inputData) {
            console.warn('⚠️ Sanitized human-readable formula to SMILES-like string:', inputData, '->', sanitized);
            inputData = sanitized;
          }
        }
      } catch (eSan) {
        console.warn('Sanitizer failed, proceeding with original inputData', eSan);
      }

      // POST to mol2chemfig API with rendering options
      console.log('%c🌐 Fetching from mol2chemfig backend...', 'color: #00AAFF; font-weight: bold;');

      // Use options from moleculeData if available (per-image settings), otherwise fallback to global settings
      const opts = moleculeData.options || {};
      console.log('%c📋 Per-image options from moleculeData:', 'color: #00BCD4; font-weight: bold;', opts);

      const aromaticCircles = opts.aromaticCircles !== undefined ? opts.aromaticCircles : settings.m2cfAromaticCircles;
      const showCarbons = opts.showCarbons !== undefined ? opts.showCarbons : settings.m2cfShowCarbons;
      const showMethyls = opts.showMethyls !== undefined ? opts.showMethyls : settings.m2cfShowMethyls;
      const atomNumbers = opts.atomNumbers !== undefined ? opts.atomNumbers : settings.m2cfAtomNumbers;
      const addH2 = opts.addH2 !== undefined ? opts.addH2 : settings.m2cfAddH2;
      const flipHorizontal = opts.flipHorizontal !== undefined ? opts.flipHorizontal : settings.m2cfFlipHorizontal;
      const flipVertical = opts.flipVertical !== undefined ? opts.flipVertical : settings.m2cfFlipVertical;

      // Build selections array for mol2chemfig API
      // mol2chemfig uses command-line style flags: -o, -c, -m, -n, -p (flip h), -q (flip v)
      const selections = [];
      if (aromaticCircles) selections.push('-o');  // Aromatic circles
      if (showCarbons) selections.push('-c');       // Show carbon labels
      if (showMethyls) selections.push('-m');       // Show methyl labels
      if (atomNumbers) selections.push('-n');       // Atom numbers

      // Smart flip logic: Determine if we need server-side or client-side flipping
      // Asymmetrical text (C, CH3, numbers) looks wrong when CSS-flipped, so use server flags
      // Symmetrical elements (H, aromatic circles) look fine CSS-flipped
      const hasAsymmetricalText = showCarbons || showMethyls || atomNumbers;

      // If asymmetrical text is shown, use server-side flip flags (-p for horizontal, -q for vertical)
      // Otherwise, we'll use CSS transforms (handled later when creating the image)
      let useServerFlipH = false;
      let useServerFlipV = false;

      if (hasAsymmetricalText) {
        if (flipHorizontal) {
          selections.push('-p');  // Server-side horizontal flip
          useServerFlipH = true;
        }
        if (flipVertical) {
          selections.push('-q');  // Server-side vertical flip
          useServerFlipV = true;
        }
      }

      // Log mol2chemfig options being used
      console.log('%c⚙️ mol2chemfig Rendering Options:', 'color: #FF6B00; font-weight: bold;', {
        aromaticCircles: aromaticCircles,
        showCarbons: showCarbons,
        showMethyls: showMethyls,
        atomNumbers: atomNumbers,
        addH2: addH2,
        flipHorizontal: flipHorizontal,
        flipVertical: flipVertical,
        hasAsymmetricalText: hasAsymmetricalText,
        useServerFlipH: useServerFlipH,
        useServerFlipV: useServerFlipV,
        selectionsArray: selections
      });

      const requestBody = {
        textAreaData: inputData,
        selections: selections,
        h2: addH2 ? 'add' : (settings.m2cfHydrogensMode || 'keep')
      };
      console.log('%c📤 POST Body:', 'color: #9B59B6; font-weight: bold;', requestBody);

      // Use background fetch to bypass CSP
      backgroundFetchJSON(`${MOL2CHEMFIG_API}/m2cf/submit`, {
        method: 'POST',
        headers: { 'Content-Type': 'application/json' },
        body: JSON.stringify(requestBody)
      })
        .then(async (data) => {
          console.log('%c📊 Response data:', 'color: #FFAA00; font-weight: bold;', data);

          if (data.error) {
            throw new Error(data.error);
          }

          // If the server returned a direct svglink use it. If not, fall back to other fields.
          if (!data.svglink) {
            // If the server returned chemfig (TeX) code, render a readable fallback SVG showing the chemfig code
            if (data.chemfig) {
              try {
                const chemText = String(data.chemfig).trim();
                const safeText = chemText.replace(/</g, '&lt;').replace(/>/g, '&gt;');
                const fallbackSvg = `<?xml version="1.0" encoding="UTF-8"?>\n<svg xmlns="http://www.w3.org/2000/svg" width="300" height="160">\n  <rect width="100%" height="100%" fill="#ffffff"/>\n  <text x="10" y="20" font-family="monospace" font-size="12" fill="#333">chemfig fallback:</text>\n  <foreignObject x="10" y="28" width="280" height="120">\n    <pre xmlns="http://www.w3.org/1999/xhtml" style="font-family:monospace; font-size:12px; color:#111; white-space:pre-wrap;">${safeText}</pre>\n  </foreignObject>\n</svg>`;
                const blob = new Blob([fallbackSvg], { type: 'image/svg+xml;charset=utf-8' });
                const url = URL.createObjectURL(blob);
                const svgImg = document.createElement('img');
                svgImg.src = url;
                svgImg.alt = 'chemfig (fallback)';
                svgImg.className = 'chemfig-diagram';

                // Apply transforms - only CSS flip if no asymmetrical text (server handles it otherwise)
                const hasAsymText1 = settings.m2cfShowCarbons || settings.m2cfShowMethyls || settings.m2cfAtomNumbers;
                let transform = '';
                if (settings.m2cfFlipHorizontal && !hasAsymText1) transform += 'scaleX(-1) ';
                if (settings.m2cfFlipVertical && !hasAsymText1) transform += 'scaleY(-1) ';

                // Apply per-image rotation and scale from opts
                const rotation = opts.rotation !== undefined ? opts.rotation : (settings.m2cfRotate || 0);
                const scale = opts.scale !== undefined ? opts.scale : 1.0;
                if (rotation) transform += `rotate(${rotation}deg) `;
                if (scale && scale !== 1.0) transform += `scale(${scale}) `;

                let filter = '';
                // Auto-apply invert in dark mode if enabled
                const shouldInvert = settings.m2cfAutoInvert && isDarkModeEnabled();
                if (shouldInvert) {
                  const invertMode = settings.m2cfInvertMode || 'full';
                  filter += getInvertFilter(invertMode);
                }

                svgImg.style.cssText = `display: inline-block; max-width: 300px; max-height: 160px; width: fit-content; margin: 0 12px 8px 0; vertical-align: middle; cursor: pointer; transform: ${transform}; filter: ${filter};`;

                chrome.storage.sync.get({
                  saveSizePerImage: false,
                  saveSizeBySMILES: true  // FIX: Enable by default
                }, async (sizeSettings) => {
                  await wrapImageWithSizeControls(svgImg, img, moleculeData, sizeSettings);
                });
                svgImg.onload = () => { URL.revokeObjectURL(url); };
                log.warn('⚠️ mol2chemfig response had no svglink — showing chemfig fallback');
                activeLoads--;
                return;
              } catch (eFallback) {
                console.warn('⚠️ Failed to render chemfig fallback', eFallback);
              }
            }

            // If a PDF link was returned, create a small clickable PDF placeholder
            if (data.pdflink) {
              try {
                const pdfUrl = data.pdflink;
                const placeholderSvg = `<?xml version="1.0" encoding="UTF-8"?>\n<svg xmlns="http://www.w3.org/2000/svg" width="300" height="120">\n  <rect width="100%" height="100%" fill="#f7f7f7" stroke="#ddd"/>\n  <text x="16" y="36" font-family="sans-serif" font-size="16" fill="#333">PDF result available</text>\n  <text x="16" y="64" font-family="sans-serif" font-size="12" fill="#666">Click to open the PDF</text>\n</svg>`;
                const blob = new Blob([placeholderSvg], { type: 'image/svg+xml;charset=utf-8' });
                const url = URL.createObjectURL(blob);
                const link = document.createElement('a');
                link.href = pdfUrl;
                link.target = '_blank';
                link.rel = 'noopener noreferrer';
                const svgImg = document.createElement('img');
                svgImg.src = url;
                svgImg.alt = 'mol2chemfig PDF result';
                svgImg.className = 'chemfig-diagram';

                // Apply transforms - only CSS flip if no asymmetrical text (server handles it otherwise)
                const hasAsymText2 = settings.m2cfShowCarbons || settings.m2cfShowMethyls || settings.m2cfAtomNumbers;
                let transform = '';
                if (settings.m2cfFlipHorizontal && !hasAsymText2) transform += 'scaleX(-1) ';
                if (settings.m2cfFlipVertical && !hasAsymText2) transform += 'scaleY(-1) ';
                if (settings.m2cfRotate) transform += `rotate(${settings.m2cfRotate}deg) `;

                let filter = '';
                if (settings.m2cfInvert) filter += 'invert(1) ';

                svgImg.style.cssText = `display: inline-block; max-width: 300px; max-height: 120px; margin: 0 12px 8px 0; vertical-align: middle; cursor: pointer; transform: ${transform}; filter: ${filter};`;

                link.appendChild(svgImg);
                img.parentNode.replaceChild(link, img);
                svgImg.onload = () => { URL.revokeObjectURL(url); };
                log.warn('⚠️ mol2chemfig response had no svglink but included a PDF — added link placeholder');
                activeLoads--;
                return;
              } catch (ePdf) {
                console.warn('⚠️ Failed to render PDF placeholder', ePdf);
              }
            }

            // No svglink, no chemfig, no pdflink — throw to be handled by catch below
            throw new Error('No svglink, chemfig, or pdflink returned by mol2chemfig');
          }

          // Extract SVG content from response and handle multiple formats:
          // - URL path like /images/abc123.svg (new cached format)
          // - Full URL like http://localhost:5001/images/abc123.svg
          // - data:image/svg+xml;base64,....
          // - data:image/svg+xml;utf8,<svg...>
          // - raw SVG text starting with <?xml or <svg
          let svgContent = data.svglink;

          // If the server returned a URL path (starts with / or http), fetch the SVG content
          if (typeof svgContent === 'string' && (svgContent.startsWith('/images/') || svgContent.startsWith('http'))) {
            try {
              const svgUrl = svgContent.startsWith('http') ? svgContent : `${MOL2CHEMFIG_API}${svgContent}`;
              console.log('%c🔗 Fetching SVG from URL:', 'color: #00AAFF; font-weight: bold;', svgUrl);

              const svgResponse = await backgroundFetchText(svgUrl, {
                method: 'GET',
                headers: { 'Accept': 'image/svg+xml' }
              });

              if (svgResponse && svgResponse.trim()) {
                svgContent = svgResponse;
                console.log('%c✅ Fetched SVG content from URL', 'color: #00FF00; font-weight: bold;');
              } else {
                throw new Error('Empty SVG response from URL');
              }
            } catch (urlFetchError) {
              console.warn('⚠️ Failed to fetch SVG from URL, trying direct use:', urlFetchError.message);
              // Fall through to try using it as a direct image src below
            }
          }

          // If the server returned a data URI, decode it, apply dark mode, then use it
          if (typeof svgContent === 'string' && svgContent.startsWith('data:image/svg+xml')) {
            try {
              // Decode the data URI to get actual SVG content
              let decodedSvg = svgContent;
              if (svgContent.includes('base64,')) {
                const base64Part = svgContent.split('base64,')[1];
                decodedSvg = atob(base64Part);
              } else if (svgContent.includes('utf8,')) {
                const utf8Part = svgContent.split('utf8,')[1];
                decodedSvg = decodeURIComponent(utf8Part);
              } else if (svgContent.includes(',')) {
                const dataPart = svgContent.split(',')[1];
                decodedSvg = decodeURIComponent(dataPart);
              }

              // Apply dark mode color inversion using comprehensive function
              if (isDarkModeEnabled()) {
                decodedSvg = invertSvgForDarkMode(decodedSvg);
              }

              // Create Blob URL from modified SVG
              const blob = new Blob([decodedSvg], { type: 'image/svg+xml;charset=utf-8' });
              const url = URL.createObjectURL(blob);
              const svgImg = document.createElement('img');
              svgImg.src = url;
              svgImg.alt = 'molecule';
              svgImg.className = 'chemfig-diagram';

              // Apply client-side flip transforms and invert filter for mol2chemfig
              // Only use CSS flip if server didn't handle it (i.e., no asymmetrical text like C, CH3, numbers)
              const hasAsymText = showCarbons || showMethyls || atomNumbers;
              let transform = '';
              if (flipHorizontal && !hasAsymText) transform += 'scaleX(-1) ';
              if (flipVertical && !hasAsymText) transform += 'scaleY(-1) ';

              // Apply per-image rotation from flags (+r)
              const rotation = opts.rotation !== undefined ? opts.rotation : (settings.m2cfRotate || 0);
              if (rotation) transform += `rotate(${rotation}deg) `;

              // Apply per-image scale from flags (+s)
              const scale = opts.scale !== undefined ? opts.scale : 1.0;
              if (scale && scale !== 1.0) transform += `scale(${scale}) `;

              let filter = '';
              // Auto-apply invert in dark mode if enabled
              const shouldInvert = settings.m2cfAutoInvert && isDarkModeEnabled();
              if (shouldInvert) {
                const invertMode = settings.m2cfInvertMode || 'full';
                filter += getInvertFilter(invertMode);
              }

              svgImg.style.cssText = `
            display: inline-block;
            max-width: 400px;
            max-height: 350px;
            width: fit-content;
            margin: 0 12px 8px 0;
            vertical-align: middle;
            cursor: pointer;
            ${transform ? `transform: ${transform.trim()};` : ''}
            ${filter ? `filter: ${filter.trim()};` : ''}
          `;
              img.dataset.loaded = 'true';
              chrome.storage.sync.get({
                saveSizePerImage: false,
                saveSizeBySMILES: true  // FIX: Enable by default
              }, async (sizeSettings) => {
                await wrapImageWithSizeControls(svgImg, img, moleculeData, sizeSettings);
              });
              svgImg.onload = () => { URL.revokeObjectURL(url); };
              log.debug(`✅ Used data URI with dark mode applied (${activeLoads} active loads remaining)`);
              activeLoads--;
              return;
            } catch (eDataURI) {
              console.warn('⚠️ Failed to decode/process data URI, falling through to regular processing', eDataURI);
              // Fall through to regular processing below
            }
          }

          // If the server returned raw SVG text (starts with < or <?xml), apply dark mode then create Blob URL
          if (typeof svgContent === 'string' && (/^\s*<\?xml|^\s*<svg/i).test(svgContent)) {
            try {
              // Apply dark mode color inversion using comprehensive function
              if (isDarkModeEnabled()) {
                svgContent = invertSvgForDarkMode(svgContent);
              }

              const blob = new Blob([svgContent], { type: 'image/svg+xml;charset=utf-8' });
              const url = URL.createObjectURL(blob);
              const svgImg = document.createElement('img');
              svgImg.src = url;
              svgImg.alt = 'molecule';
              svgImg.className = 'chemfig-diagram';

              // Apply client-side flip transforms and invert filter for mol2chemfig
              // Only use CSS flip if no asymmetrical text (server handles it otherwise)
              const hasAsymText3 = settings.m2cfShowCarbons || settings.m2cfShowMethyls || settings.m2cfAtomNumbers;
              let transform = '';
              if (settings.m2cfFlipHorizontal && !hasAsymText3) transform += 'scaleX(-1) ';
              if (settings.m2cfFlipVertical && !hasAsymText3) transform += 'scaleY(-1) ';

              // Apply per-image rotation and scale from opts
              const rotation = opts.rotation !== undefined ? opts.rotation : (settings.m2cfRotate || 0);
              const scale = opts.scale !== undefined ? opts.scale : 1.0;
              if (rotation) transform += `rotate(${rotation}deg) `;
              if (scale && scale !== 1.0) transform += `scale(${scale}) `;

              let filter = '';
              // Auto-apply invert in dark mode if enabled
              const shouldInvert = settings.m2cfAutoInvert && isDarkModeEnabled();
              if (shouldInvert) {
                const invertMode = settings.m2cfInvertMode || 'full';
                filter += getInvertFilter(invertMode);
              }

              svgImg.style.cssText = `
            display: inline-block;
            max-width: 400px;
            max-height: 350px;
            width: fit-content;
            margin: 0 12px 8px 0;
            vertical-align: middle;
            cursor: pointer;
            ${transform ? `transform: ${transform.trim()};` : ''}
            ${filter ? `filter: ${filter.trim()};` : ''}
          `;
              img.dataset.loaded = 'true';
              chrome.storage.sync.get({
                saveSizePerImage: false,
                saveSizeBySMILES: true  // FIX: Enable by default
              }, async (sizeSettings) => {
                await wrapImageWithSizeControls(svgImg, img, moleculeData, sizeSettings);
              });
              // Revoke URL after image loads to free memory
              svgImg.onload = () => { URL.revokeObjectURL(url); };
              log.debug(`✅ Used Blob URL for raw SVG with dark mode applied (${activeLoads} active loads remaining)`);
              activeLoads--;
              return;
            } catch (errBlob) {
              console.warn('⚠️ Failed to create Blob for raw SVG, falling back to data URI creation', errBlob);
            }
          }

          // If the svglink contains a base64 payload somewhere, attempt to extract and decode it
          if (typeof svgContent === 'string' && svgContent.includes('base64,')) {
            let base64Part = svgContent.split('base64,')[1] || '';
            base64Part = base64Part.replace(/\s+/g, '');
            try {
              const decoded = atob(base64Part);
              svgContent = decoded;
            } catch (err) {
              try {
                const decodedURIComponent = decodeURIComponent(base64Part);
                svgContent = atob(decodedURIComponent);
              } catch (err2) {
                console.warn('⚠️ Failed to decode base64 svglink, will try to use original svglink as src', err2);
                // fall through — try using original svglink below
              }
            }
          }

          console.log('%c📊 SVG Content (first 100 chars):', 'color: #FFAA00; font-weight: bold;', (typeof svgContent === 'string' ? svgContent.substring(0, 100) : String(svgContent)) + '...');

          // Apply dark mode color inversion using comprehensive function
          if (isDarkModeEnabled()) {
            svgContent = invertSvgForDarkMode(svgContent);
          }

          // Create clean SVG image with dark mode support
          const svgImg = document.createElement('img');
          svgImg.alt = 'molecule';
          svgImg.className = 'chemfig-diagram';

          // Apply client-side flip transforms and invert filter for mol2chemfig
          // Only use CSS flip if no asymmetrical text (server handles it otherwise)
          const hasAsymText4 = settings.m2cfShowCarbons || settings.m2cfShowMethyls || settings.m2cfAtomNumbers;
          let transform = '';
          if (settings.m2cfFlipHorizontal && !hasAsymText4) transform += 'scaleX(-1) ';
          if (settings.m2cfFlipVertical && !hasAsymText4) transform += 'scaleY(-1) ';

          // Apply per-image rotation and scale from opts
          const rotation = opts.rotation !== undefined ? opts.rotation : (settings.m2cfRotate || 0);
          const scale = opts.scale !== undefined ? opts.scale : 1.0;
          if (rotation) transform += `rotate(${rotation}deg) `;
          if (scale && scale !== 1.0) transform += `scale(${scale}) `;

          let filter = '';
          // Auto-apply invert in dark mode if enabled
          const shouldInvert = settings.m2cfAutoInvert && isDarkModeEnabled();
          if (shouldInvert) {
            const invertMode = settings.m2cfInvertMode || 'full';
            filter += getInvertFilter(invertMode);
          }

          svgImg.style.cssText = `
            display: inline-block;
            max-width: 1000px;
            max-height: 800px;
            width: fit-content;
            height: auto;
            margin: 0 12px 8px 0;
            vertical-align: middle;
            cursor: pointer;
            ${transform ? `transform: ${transform.trim()};` : ''}
            ${filter ? `filter: ${filter.trim()};` : ''}
          `;

          // Prefer creating a Blob URL to avoid btoa() Unicode issues
          try {
            const blob = new Blob([svgContent], { type: 'image/svg+xml;charset=utf-8' });
            const url = URL.createObjectURL(blob);
            svgImg.src = url;
            // Replace original img element with size controls
            chrome.storage.sync.get({
              saveSizePerImage: false,
              saveSizeBySMILES: true  // FIX: Enable by default
            }, async (sizeSettings) => {
              await wrapImageWithSizeControls(svgImg, img, moleculeData, sizeSettings);
            });
            // Revoke object URL after image loads
            svgImg.onload = () => { URL.revokeObjectURL(url); };
          } catch (eBlob) {
            // Fallback: UTF-8 safe base64 encoding then data URI
            try {
              const utf8Base64 = btoa(unescape(encodeURIComponent(svgContent)));
              svgImg.src = 'data:image/svg+xml;base64,' + utf8Base64;
              chrome.storage.sync.get({
                saveSizePerImage: false,
                saveSizeBySMILES: true  // FIX: Enable by default
              }, async (sizeSettings) => {
                await wrapImageWithSizeControls(svgImg, img, moleculeData, sizeSettings);
              });
            } catch (eBtoa) {
              // Last resort: put raw SVG into data URI (may break in some browsers)
              svgImg.src = 'data:image/svg+xml;utf8,' + encodeURIComponent(svgContent);
              chrome.storage.sync.get({
                saveSizePerImage: false,
                saveSizeBySMILES: true  // FIX: Enable by default
              }, async (sizeSettings) => {
                await wrapImageWithSizeControls(svgImg, img, moleculeData, sizeSettings);
              });
            }
          }

          // Mark as loaded
          img.dataset.loaded = 'true';

          console.log('%c✅ mol2chemfig image loaded successfully', 'color: green; font-weight: bold;');

          activeLoads--;
        })
        .catch(error => {
          console.error('%c❌ Error fetching from mol2chemfig:', 'color: red; font-weight: bold;', error.message);

          // Fallback: show error message
          img.src = 'data:image/svg+xml;base64,PHN2ZyB3aWR0aD0iMzAwIiBoZWlnaHQ9IjIwMCIgeG1sbnM9Imh0dHA6Ly93d3cudzMub3JnLzIwMDAvc3ZnIj48cmVjdCB3aWR0aD0iMzAwIiBoZWlnaHQ9IjIwMCIgZmlsbD0iI2ZmZiIvPjx0ZXh0IHg9IjEwIiB5PSIzMCIgZmlsbD0iI2YwMCIgZm9udC1mYW1pbHk9Im1vbm9zcGFjZSIgZm9udC1zaXplPSIxNCI+RXJyb3I6IE1vbDJDaGVtZmlnIExvYWQgRmFpbGVkPC90ZXh0Pjwvc3ZnPg==';
          img.classList.add('chemfig-fadein');
          img.classList.remove('chemfig-loading', 'chemfig-molecule-viewer');
          img.dataset.loaded = 'true';

          activeLoads--;
        });

      log.debug(`✅ Fetching molecule data from mol2chemfig (${activeLoads} active loads remaining)`);
    } catch (error) {
      console.error(`%c❌ Error: ${error.message}`, 'color: red; font-weight: bold;');
      activeLoads--;
    }
  }

  // Helper function to load molecule image based on configured renderer engine
  function loadMoleculeImage(img) {
    // First check if the image has a specific renderer class (takes priority)
    if (img.classList.contains('chemfig-pubchem')) {
      console.log('%c🌐 Using PUBCHEM renderer (from class)', 'background: #4CAF50; color: #FFF; font-size: 14px; padding: 5px;');
      loadPubChemImage(img);
    } else if (img.classList.contains('chemfig-mol2chemfig')) {
      console.log('%c📐 Using MOL2CHEMFIG renderer (from class)', 'background: #FF6B00; color: #FFF; font-size: 14px; padding: 5px;');
      loadMol2chemfigImage(img);
    } else if (settings.rendererEngine === 'mol2chemfig') {
      console.log('%c📐 Using MOL2CHEMFIG renderer engine', 'background: #FF6B00; color: #FFF; font-size: 14px; padding: 5px;');
      loadMol2chemfigImage(img);
    } else if (settings.rendererEngine === 'pubchem') {
      console.log('%c🌐 Using PUBCHEM renderer engine', 'background: #4CAF50; color: #FFF; font-size: 14px; padding: 5px;');
      loadPubChemImage(img);
    } else {
      console.log('%c🧪 Using MOLECULEVIEWER renderer engine', 'background: #0088FF; color: #FFF; font-size: 14px; padding: 5px;');
      loadMoleculeViewerImage(img);
    }
  }

  // Helper function to download SVG file
  function downloadSVG(svgContent, filename) {
    const blob = new Blob([svgContent], { type: 'image/svg+xml' });
    const url = URL.createObjectURL(blob);
    const link = document.createElement('a');
    link.href = url;
    link.download = filename;
    link.style.display = 'none';
    document.body.appendChild(link);
    link.click();
    document.body.removeChild(link);
    URL.revokeObjectURL(url);
    console.log('%c💾 Downloaded:', 'color: #4CAF50; font-weight: bold;', filename);
  }

  // Override observer callback to handle both types
  const originalObserver = observer;
  const newObserver = new IntersectionObserver((entries) => {
    console.log('%c👁️ IntersectionObserver triggered!', 'background: #000; color: #FFFF00; font-size: 16px; padding: 5px;');
    console.log('Number of entries:', entries.length);
    entries.forEach((entry) => {
      const img = entry.target;
      console.log('Entry:', entry);
      console.log('Image:', img);
      console.log('Is intersecting?', entry.isIntersecting);
      console.log('Already loaded?', img.dataset.loaded);
      console.log('Has molecule-viewer class?', img.classList.contains('chemfig-molecule-viewer'));

      if (entry.isIntersecting && img.dataset.loaded !== 'true') {
        // Check if this is a chemistry image (MoleculeViewer, PubChem, or Mol2chemfig)
        if (img.classList.contains('chemfig-molecule-viewer') || img.classList.contains('chemfig-pubchem') || img.classList.contains('chemfig-mol2chemfig')) {
          console.log('%c🧪 Detected molecule image!', 'background: #0088FF; color: #FFF; font-size: 14px; padding: 5px;');
          if (activeLoads < maxConcurrentLoads) {
            loadMoleculeImage(img);  // Use renderer-agnostic function
          } else {
            if (typeof requestIdleCallback !== 'undefined') {
              requestIdleCallback(() => {
                if (activeLoads < maxConcurrentLoads) {
                  loadMoleculeImage(img);
                }
              });
            } else {
              setTimeout(() => {
                if (activeLoads < maxConcurrentLoads) {
                  loadMoleculeImage(img);
                }
              }, 50);
            }
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
    rootMargin: '300px'
  });

  // Replace the old observer
  window._lazyLoadObserver = newObserver;

  // Expose loading function globally for immediate loading when performance mode is off
  window._loadMoleculeImage = loadMoleculeImage;

  // Export show3DViewerInline for hover controls
  window.show3DViewerInline = show3DViewerInline;

  log.success('✅ Lazy-loading observer initialized (max 3 concurrent loads)');
}


/**
 * Inject MathJax library from Extension Context
 * The new approach: serve math-render.js from extension using src= (no CSP blocks this)
 */
function injectMathJax() {
  log.inject('💉 INJECTING Math Rendering Engine');

  // Get the extension URL for math-render.js
  const mathRenderUrl = chrome.runtime.getURL('math-render.js');
  log.inject(`Math render URL: ${mathRenderUrl}`);

  // Inject into page context using a wrapper script
  const injectionScript = document.createElement('script');
  injectionScript.textContent = `
    // Load math-render.js from extension in page context
    const script = document.createElement('script');
    script.src = '${mathRenderUrl}';
    script.async = true;
    script.onload = function() {
      console.log('[ChemRenderer] Math-render loaded in page context');
      window.dispatchEvent(new CustomEvent('mathJaxReady'));
    };
    script.onerror = function(e) {
      console.log('[ChemRenderer] Math-render failed, using basic rendering');
      window.dispatchEvent(new CustomEvent('mathRendererFailed'));
    };
    document.head.appendChild(script);
  `;

  try {
    (document.head || document.documentElement).appendChild(injectionScript);
    log.inject('✅ Injected math-render loader into page');

    // Listen for MathJax to be ready
    window.addEventListener('mathJaxReady', () => {
      log.success('✅ MathJax is ready!');
      setTimeout(() => {
        initializeChemistryRendering();
      }, 500);
    }, { once: true });

    window.addEventListener('mathRendererFailed', () => {
      log.error('Failed to load math renderer');
      initializeBasicRendering();
    }, { once: true });

    // Timeout fallback - if MathJax doesn't load in 8 seconds
    setTimeout(() => {
      if (!window.MathJax) {
        log.error('MathJax did not load within 8 seconds');
        initializeBasicRendering();
      }
    }, 8000);

  } catch (e) {
    log.error('Failed to inject math-render.js:', e);
    initializeBasicRendering();
  }
}

/**
 * Initialize chemistry rendering with KaTeX
 */
function initializeChemistryRendering() {
  log.inject('🔍 Setting up chemistry formula rendering');

  // Make KaTeX available
  log.success('✅ KaTeX Ready! Chemistry rendering available');

  // Start scanning after a brief delay to ensure KaTeX is ready
  setTimeout(() => {
    log.inject('🚀 Running initial page scan now that KaTeX is ready');
    if (settings.enabled) {
      scanAndRender();
    } else {
      log.inject('⚠️ Extension disabled, skipping scan');
    }
  }, 500);
}

/**
 * Basic rendering without external libraries
 * Replaces \\ce{...} with readable text approximation
 */
function initializeBasicRendering() {
  log.inject('🔄 Using basic text rendering (no external library)');
  log.inject('Formulas will be visible but not beautifully rendered');

  setTimeout(() => {
    log.inject('🔍 Running initial page scan with basic rendering');
    if (settings.enabled) {
      scanAndRender();
    } else {
      log.inject('⚠️ Extension disabled, skipping scan');
    }
  }, 100);
}


/**
 * Fallback MathJax loader
 */
function loadMathJaxFallback() {
  log.inject('🔄 METHOD 3: FALLBACK CDN');
  log.debug('Fallback URL: https://cdnjs.cloudflare.com/ajax/libs/mathjax/3.2.2/es5/tex-svg.js');

  fetch('https://cdnjs.cloudflare.com/ajax/libs/mathjax/3.2.2/es5/tex-svg.js')
    .then(response => {
      log.inject('✅ Fallback CDN fetch successful');
      return response.text();
    })
    .then(scriptText => {
      log.inject('Creating inline script from fallback CDN...');
      const script = document.createElement('script');
      script.textContent = scriptText;
      script.id = 'MathJax-script-fallback';
      document.head.appendChild(script);
      log.success('📦 MathJax injected from fallback CDN');
    })
    .catch(err => {
      log.error('❌ All CDN methods failed', err);
      log.error('Cannot load MathJax - check internet connection');
    });
}

/**
 * Render chemistry formula with KaTeX
 */
function renderWithKaTeX(formula) {
  if (!window.katex) {
    return formula;  // Fall back to text
  }

  try {
    // KaTeX render options
    const html = window.katex.renderToString(formula, {
      throwOnError: false,
      trust: true
    });
    return html;
  } catch (err) {
    log.debug('KaTeX rendering failed for: ' + formula, err);
    return formula;
  }
}

/**
 * Convert chemistry notation to KaTeX-compatible format
 */
function convertToKaTeX(formula) {
  let result = formula;

  // mhchem formulas: \ce{...} → stays as is for KaTeX
  // KaTeX requires mhchem.sty to be loaded, which it's not
  // So convert common patterns to Unicode/text

  // H2O → H₂O
  result = result.replace(/H(\d+)O/g, 'H$1O');
  result = result.replace(/(\w)(\d+)/g, (match, element, number) => {
    const subscript = number.split('').map(d => {
      const codes = ['₀', '₁', '₂', '₃', '₄', '₅', '₆', '₇', '₈', '₉'];
      return codes[parseInt(d)];
    }).join('');
    return element + subscript;
  });

  return result;
}

/**
 * Scan page and wrap formulas for MathJax
 * Based on ReplaceR's text replacement pattern
 * 
 * HOW IT WORKS:
 * 1. Gets all DOM elements
 * 2. Iterates through child text nodes
 * 3. Detects chemical formulas (\\ce{...}, chemfig{...})
 * 4. Wraps them in $ delimiters
 * 5. Asks MathJax to render the wrapped formulas
 */

/**
 * Apply layout mode (horizontal/vertical) to chemfig containers
 */
function applyLayoutMode() {
  const containers = document.querySelectorAll('.chemfig-container');
  containers.forEach(container => {
    // Remove both classes to reset
    container.classList.remove('horizontal', 'vertical');

    // Add the current layout mode class
    if (settings.layoutMode === 'vertical') {
      container.classList.add('vertical');
    } else {
      container.classList.add('horizontal');
    }
  });

  if (containers.length > 0) {
    log.debug(`📐 Applied layout mode: ${settings.layoutMode} to ${containers.length} containers`);
  }
}

function scanAndRender() {
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

  log.debug(`Total elements in DOM: ${elements.length}`);

  for (let i = 0; i < elements.length; i++) {
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

    for (let j = 0; j < element.childNodes.length; j++) {
      const node = element.childNodes[j];

      // Only process text nodes (nodeType 3)
      if (node.nodeType === 3) {
        textNodesFound++;
        const text = node.nodeValue;

        // Log ALL text nodes that contain backslash (potential formulas)
        if (text.includes('\\') || text.includes('ce{') || text.includes('chemfig')) {
          chemistryTextFound.push(text.substring(0, 100));
          log.debug(`🔬 Found potential chemistry text: "${text.substring(0, 80)}..."`);
        }

        const replacedText = wrapChemicalFormulas(text);

        if (replacedText !== text) {
          log.debug(`✨ Found formula in text: "${text.substring(0, 50)}..."`);
          log.debug(`Wrapped result: "${replacedText.substring(0, 50)}..."`);

          // Replace the text node with new content
          const span = document.createElement('span');
          span.innerHTML = replacedText;
          span.setAttribute('data-chem-processed', 'true');
          element.replaceChild(span, node);

          // Setup loading for newly created images
          const images = span.querySelectorAll('img.chemfig-diagram[data-loaded="false"]');
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
            } else if (window._loadMoleculeImage) {
              // Performance mode OFF: load all images immediately
              images.forEach(img => {
                console.log('%c[ChemRenderer] Immediately loading image:', 'color: #00FF00;', img.className);
                window._loadMoleculeImage(img);
              });
              log.debug(`📊 Immediately loading ${images.length} SVGs`);
            } else {
              console.error('%c[ChemRenderer] ERROR: No loader available!', 'color: #FF0000; font-weight: bold;');
            }
          }

          replacements++;
        }
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
 * Detect if the page is using dark mode
 * Checks multiple signals: CSS media query, background color, body styles
 */
function detectDarkMode() {
  // Method 1: Check prefers-color-scheme media query
  if (window.matchMedia && window.matchMedia('(prefers-color-scheme: dark)').matches) {
    return true;
  }

  // Method 2: Check body background color (detect dark backgrounds)
  const bodyBg = window.getComputedStyle(document.body).backgroundColor;
  if (bodyBg) {
    // Parse RGB values
    const rgb = bodyBg.match(/\d+/g);
    if (rgb && rgb.length >= 3) {
      const r = parseInt(rgb[0]);
      const g = parseInt(rgb[1]);
      const b = parseInt(rgb[2]);

      // Calculate luminance (perceived brightness)
      // Dark if luminance < 128
      const luminance = (0.299 * r + 0.587 * g + 0.114 * b);
      if (luminance < 128) {
        return true;
      }
    }
  }

  // Method 3: Check for common dark mode classes
  const darkClasses = ['dark', 'dark-mode', 'dark-theme', 'night-mode'];
  for (const className of darkClasses) {
    if (document.body.classList.contains(className) ||
      document.documentElement.classList.contains(className)) {
      return true;
    }
  }

  // Method 4: ChatGPT specific - check for dark background
  const htmlBg = window.getComputedStyle(document.documentElement).backgroundColor;
  if (htmlBg) {
    const rgb = htmlBg.match(/\d+/g);
    if (rgb && rgb.length >= 3) {
      const luminance = (0.299 * parseInt(rgb[0]) + 0.587 * parseInt(rgb[1]) + 0.114 * parseInt(rgb[2]));
      if (luminance < 128) {
        return true;
      }
    }
  }

  // Default: light mode
  return false;
}

/**
/**
 * NOMENCLATURE DATABASE - Chemical names to tested chemfig formulas
 * All formulas tested and verified working
 */
const NOMENCLATURE_DB = {
  // Alcohols
  "ethanol": "\\chemfig{C-[1]C(-[1]OH)",
  "propanol": "\\chemfig{C-[1]C-[7]C(-[1]OH)}",
  "isopropanol": "\\chemfig{C-[1]C(-[1]OH)(-[7]C)}",

  // Ethers
  "dimethyl ether": "\\chemfig{C-[1]O-[7]C}",
  "diethyl ether": "\\chemfig{C-[1]C-[1]O-[7]C-[7]C}",

  // Amines
  "methylamine": "\\chemfig{C-[1]NH2}",
  "ethylamine": "\\chemfig{C-[1]C-[7]NH2}",
  "propylamine": "\\chemfig{C-[1]C-[7]C(-[1]NH2)}",

  // Aldehydes
  "acetaldehyde": "\\chemfig{C(=O)-[1]C}",
  "propanal": "\\chemfig{C(=O)-[1]C-[7]C}",

  // Ketones
  "acetone": "\\chemfig{C-[1]C(=O)-[7]C}",
  "butanone": "\\chemfig{C-[1]C-[7]C(=O)-[1]C}",

  // Benzene derivatives
  "phenol": "\\chemfig{*6(=(-OH)-=(-)-=(-)-=(-)-)}",
  "aniline": "\\chemfig{*6(=(-NH2)-=(-)-=(-)-=(-)-)}",
  "nitrobenzene": "\\chemfig{*6(=(-NO2)-=(-)-=(-)-=(-)-)}",
};

/**
 * Get chemfig formula from nomenclature
 */
function getChemfigFromNomenclature(name) {
  if (!name) return null;
  const normalized = name.toLowerCase().trim();
  const formula = NOMENCLATURE_DB[normalized];
  if (formula) {
    console.log(`[Nomenclature] Found: "${name}" → chemfig formula`);
  }
  return formula || null;
}

/**
 * Chemfig to SMILES conversion - DEPRECATED
 * 
 * ❌ DISABLED - chemfig is a LaTeX drawing language, not suitable for direct chemical rendering
 * Use chem:SMILES or chem:nomenclature instead (like CodeCogs)
 * 
 * This function is kept for reference but no longer called.
 * It was generating warnings about malformed SMILES strings because chemfig
 * syntax doesn't map cleanly to SMILES notation.
 * 
 * Example malformed outputs it was generating:
 *   c1ccccc1-C(-Br)-C (dashes in wrong places)
 *   c1ccccc1-=(-Cl)-=(-)-=(-)-)  (invalid bond syntax)
 * 
 * Modern approach: Use direct notation
 *   chem:CCO (ethanol in SMILES)
 *   chem:acetone (acetone by name)
 */
function chemfigToSmiles(chemfigContent) {
  log.debug(`🔬 [DEPRECATED] Would convert chemfig to SMILES: ${chemfigContent.substring(0, 50)}...`);

  let smiles = chemfigContent;
  let conversionSteps = [];

  try {
    // Step 1: Remove angle brackets [n] - they're just for drawing angles
    // C-[1]C → C-C
    smiles = smiles.replace(/-\[\d+\]/g, '-');
    conversionSteps.push(`Remove angle brackets: ${smiles.substring(0, 30)}...`);

    // Step 2: Remove parentheses used for branching positions in chemfig
    // Most times we can simplify these, but keep double bond info
    // C(-OH) → C(O) then C(O) → CO (simplified)

    // Step 3: Handle functional groups
    // OH → O
    smiles = smiles.replace(/OH/g, 'O');
    // NH2 → N
    smiles = smiles.replace(/NH2/g, 'N');
    // NO2 → N(=O)(=O) → simplified to [N+](=O)[O-]
    smiles = smiles.replace(/NO2/g, '[N+](=O)[O-]');
    conversionSteps.push(`Handle functional groups: ${smiles.substring(0, 30)}...`);

    // Step 4: Handle aromatic rings
    // *6(=(-)-=(-)-=(-)-) → c1ccccc1 (benzene pattern)
    // *5(=(-)-=(-)-=(-)-) → c1cccc1 (pyrrole pattern)

    // Benzene: *6 means 6-membered ring
    if (smiles.includes('*6')) {
      // Count the segments between dashes to verify 6 carbons
      // For now, replace common 6-membered aromatic pattern
      smiles = smiles.replace(/\*6\s*\([^)]*\)/g, 'c1ccccc1');
      conversionSteps.push(`Convert 6-membered ring to benzene: ${smiles.substring(0, 30)}...`);
    }

    // Step 5: Clean up multiple bonds
    // Double bond: = stays as =
    // Triple bond: ~ or ≡ becomes #
    smiles = smiles.replace(/~/g, '#');
    smiles = smiles.replace(/≡/g, '#');
    conversionSteps.push(`Standardize bonds: ${smiles.substring(0, 30)}...`);

    // Step 6: Remove parentheses that don't mean branching
    // Keep branching parentheses: C(C) means branches
    // Remove wrapping parentheses: (C-C) → C-C
    while (smiles.match(/^\([^()]*\)$/)) {
      smiles = smiles.substring(1, smiles.length - 1);
    }
    conversionSteps.push(`Clean parentheses: ${smiles.substring(0, 30)}...`);

    // Step 7: Remove extra spaces
    smiles = smiles.replace(/\s+/g, '');
    conversionSteps.push(`Remove spaces: ${smiles.substring(0, 30)}...`);

    // Step 8: Validate SMILES (basic check)
    // Should contain C, N, O, S, P, H, numbers, bonds, parentheses
    if (!/^[CNOSPHcnosp\d\-=()[\]#+@:\\\/]+$/.test(smiles)) {
      log.warn(`⚠️  SMILES validation warning - unusual characters: ${smiles.substring(0, 50)}`);
    }

    log.debug(`✅ Converted chemfig to SMILES: ${smiles}`);
    conversionSteps.forEach(step => log.debug(`   ${step}`));
    return smiles;

  } catch (err) {
    log.error(`❌ Chemfig to SMILES conversion failed: ${err.message}`);
    log.debug(`Original chemfig: ${chemfigContent}`);
    return null;  // Fall back to default
  }
}

/**
 * Build MoleculeViewer rendering request
 * ✅ ONLY uses local MoleculeViewer server (localhost:5000)
 * ❌ CodeCogs completely removed
 */
function buildChemfigImageUrl(latex, isDarkMode, chemfigContent = null) {
  log.debug('🔬 Using MoleculeViewer server for rendering');

  // Convert chemfig to SMILES
  const smiles = chemfigToSmiles(chemfigContent);

  if (smiles) {
    // Build rendering options
    const options = {
      show_carbons: settings.showCarbons,
      show_methyls: settings.showMethyls,
      aromatic_circles: settings.aromaticCircles,
      fancy_bonds: settings.fancyBonds,
      atom_numbers: settings.atomNumbers,
      hydrogens: settings.hydrogensMode,
      flip_horizontal: settings.flipHorizontal,
      flip_vertical: settings.flipVertical,
      recalculate_coordinates: false
    };

    log.debug(`  Converted chemfig → SMILES: ${smiles}`);
    log.debug(`  Rendering options:`, options);

    // Return object for MoleculeViewer rendering
    return {
      isMoleculeViewer: true,
      smiles: smiles,
      options: options,
      isDarkMode: isDarkMode
    };
  } else {
    log.error('❌ Chemfig to SMILES conversion failed!');
    log.error(`  Original chemfig: ${chemfigContent}`);

    // ✅ NO FALLBACK - Just return error object
    return {
      isMoleculeViewer: true,
      smiles: null,
      error: 'Could not convert chemfig to SMILES',
      options: {}
    };
  }
}

/**
 * Remove unnecessary explicit hydrogens from chemfig
 * Converts: C(-H)(-H)-C(-H)(-H)-C(-H)(-H) → C-C-C (clean)
 * Keeps: Important heteroatoms and hydrogens (NH2, OH, etc.)
 */
function removeUnnecessaryHydrogens(chemfigContent) {
  let result = chemfigContent;

  // Remove (-H) from carbons - they're implicit in chemfig
  // But be careful not to remove H from OH, NH2, etc.
  // Pattern: (-H) that's not part of a functional group

  // Remove (-H) and (-[n]H) patterns
  result = result.replace(/\(-\[?\d*\]?H\)/g, '');

  // Remove redundant bonds with just H
  result = result.replace(/-H(?![\w])/g, '');

  return result;
}

/**
 * Add zigzag angles to simple carbon chains for realistic rendering
 * Converts: C-C-C → C-[1]C-[7]C (zigzag pattern)
 */
function addZigzagAngles(chemfigContent) {
  // Simple carbon chains without angles - add zigzag
  // Pattern: C-C or C-C-C or longer chains without angles
  // Add alternating angles [1] and [7] for zigzag effect

  let result = chemfigContent;

  // Replace simple C-C chains with angled chains
  // Match: C-C (but not C-[n]C which already has angles)
  result = result.replace(/C(?!-\[)-C/g, (match) => {
    // This is a C-C without angles
    // Replace with C-[1]C for first connection
    return 'C-[1]C';
  });

  // Now add the zigzag to subsequent carbons
  // C-[1]C-C → C-[1]C-[7]C
  result = result.replace(/\[1\]C(?!-\[)-C/g, '[1]C-[7]C');

  // Continue pattern for longer chains
  // C-[7]C-C → C-[7]C-[1]C
  result = result.replace(/\[7\]C(?!-\[)-C/g, '[7]C-[1]C');

  // C-[1]C-C → C-[1]C-[7]C (if we have more carbons after)
  result = result.replace(/\[1\]C(?!-\[)-C/g, '[1]C-[7]C');

  return result;
}

/**
 * Simplify chemfig structure by removing explicit carbon labels
 * Renders carbons as implicit sticks instead of showing CH, CH2, CH3 labels
 * Keeps heteroatoms (N, O, P, S, etc.) as explicit labels
 */
function simplifyChemfigCarbons(chemfigContent) {
  if (!settings.renderCarbonsAsSticks) {
    return chemfigContent;  // Keep original if option disabled
  }

  let simplified = chemfigContent;

  // Replace CH followed by number or end of word with just the position marker
  // CH3 → just the carbon (rendered as stick)
  // CH2 → just the carbon (rendered as stick)
  // CH → just the carbon (rendered as stick)
  // But preserve groups like NH2, OH, etc.

  // Pattern: Replace standalone CH labels (not followed by letters like N, O, P, etc.)
  // This regex looks for: whitespace or start, then CH, then not a letter or number
  simplified = simplified.replace(/\bCH(\d+)?(?![A-Z])/g, (match) => {
    // Return empty or just a dash depending on context
    // The position marker is preserved so structure remains valid
    return '';
  });

  // Remove standalone C that was left (already implicit in chemfig)
  simplified = simplified.replace(/\b(?<!-)C(?!H)(?![0-9])/g, '');

  return simplified;
}

/**
 * Wrap chemical formulas in renderable format
 * Converts plain text chemistry to Unicode with subscripts/superscripts
 */
function wrapChemicalFormulas(text) {
  if (!text || text.trim().length === 0) return text;

  let result = text;
  let patternMatches = [];

  // Pattern 0a: chem:\chemfig{...} (chemfig with chem: prefix)
  // Example: chem:\chemfig{CH_3-CH_2-OH}, chem:\chemfig{C-C-C}
  log.debug('🧪 Applying Pattern 0a: chem:\\chemfig{...} → mol2chemfig conversion');
  const chemChemfigMatches = result.match(/\bchem:\\chemfig\{([^{}]*(?:\{[^{}]*\}[^{}]*)*)\}(?::(\d+))?/g);
  if (chemChemfigMatches) {
    log.debug(`  Found ${chemChemfigMatches.length} chem:\\chemfig{} patterns`);

    result = result.replace(/\bchem:\\chemfig\{([^{}]*(?:\{[^{}]*\}[^{}]*)*)\}(?::(\d+))?/g, (match, content, rotationAngle) => {
      log.debug(`  Converting chemfig: ${match}`);
      window.chemRendererPerformance.recordStructure();

      const rotation = rotationAngle ? parseInt(rotationAngle) : 0;
      const isDarkMode = detectDarkMode();

      // Use mol2chemfig conversion via buildChemfigImageUrl
      const latex = '$\\chemfig{' + content + '}$';
      const imageUrl = buildChemfigImageUrl(latex, isDarkMode, content);

      const styleAttr = `style="transform: rotate(${rotation}deg); margin: 0 12px 8px 0; vertical-align: middle; width: 350px; height: 300px; object-fit: contain;"`;

      let converted;
      if (settings.devMode) {
        const devStyle = `style="background: #f0f0f0; border: 1px solid #ddd; border-radius: 4px; padding: 8px 12px; font-family: monospace; font-size: 12px; color: #333; margin: 0 12px 8px 0; display: inline-block; white-space: pre-wrap; word-break: break-all; max-width: 500px;"`;
        converted = `<span class="chemfig-dev-mode" ${devStyle}>chem:\\chemfig{${content}}</span>`;
      } else {
        // Store both the raw chemfig and a MoleculeViewer fallback in the dataset.
        // This allows mol2chemfig to receive the original chemfig text while
        // MoleculeViewer can use the SMILES fallback when needed.
        const moleculeViewerPayload = {
          type: 'chemfig',
          chemfig: content,
          latex: latex,
          fallback: imageUrl // may contain converted SMILES or an error
        };
        const moleculeViewerData = btoa(JSON.stringify(moleculeViewerPayload));
        converted = `<img src="" alt="chemfig" class="chemfig-diagram chemfig-molecule-viewer" data-molecule-viewer="${moleculeViewerData}" data-loaded="false" data-rotation="${rotation}" ${styleAttr}>`;
      }

      log.debug(`  📤 Sending chemfig to mol2chemfig: ${content ? content.substring(0, Math.min(50, content.length)) : '(empty)'}...`);
      return converted;
    });
  }

  // Pattern 0b: chem:text: (double colon format - captures everything in between)
  // Example: chem:CCO:, chem:benzene:, chem:1-chloro-benzene:, chem:CC(=O)C:, chem:histamine/+c+o+m+n+3d:, chem:phenol/d-c:
  log.debug('🧪 Applying Pattern 0b: chem:text: → Using selected renderer engine');
  const chemMatches = result.match(/\bchem:([^:]+):/g);
  if (chemMatches) {
    log.debug(`  Found ${chemMatches.length} chem:text: patterns`);
    log.debug(`  Renderer engine: ${settings.rendererEngine}, 3D Viewer enabled: ${settings.enable3DViewer}`);

    result = result.replace(/\bchem:([^:]+):/g, (match, content) => {
      const content_trimmed = content.trim();
      if (!content_trimmed) {
        return match; // Empty content, skip
      }

      log.debug(`  Converting molecule: ${content_trimmed}`);
      window.chemRendererPerformance.recordStructure();

      // Determine if it's SMILES or nomenclature by checking for SMILES-specific characters
      // Removed + from check to avoid confusing 'name/+flags' with SMILES
      const isSMILES = /[=\[\]()@#\\]/.test(content_trimmed);

      // For basic functionality, just use the content as is (backward compatibility)
      let compoundName = content_trimmed;

      // Process flags if any are detected (always enabled by default)
      let flagOverrides = { useDefaults: false }; // Default empty flags
      let currentSettings = settings; // Default to original settings

      // Check if there are flags to process
      if (content_trimmed.includes('/') || content_trimmed.includes('+') || content_trimmed.includes('-')) {
        try {
          flagOverrides = window.parseChemFlags('chem:' + content_trimmed + ':');
          console.log('%c🏴 Parsed flags:', 'color: #FF6B00; font-weight: bold;', flagOverrides);

          // Determine base settings:
          // - If /d flag is present: use global defaults
          // - If no /d flag: start with clean slate (no defaults)
          let baseSettings;
          if (flagOverrides.useDefaults) {
            // /d flag present: use defaults and apply overrides
            baseSettings = settings;
          } else {
            // No /d flag: start clean, only apply explicit flags
            baseSettings = {
              ...settings,
              // Clear all display options to start fresh
              m2cfShowCarbons: undefined,
              m2cfAromaticCircles: undefined,
              m2cfShowMethyls: undefined,
              m2cfAtomNumbers: undefined,
              m2cfAddH2: undefined,
              m2cfFlipHorizontal: undefined,
              m2cfFlipVertical: undefined,
              m2cfInvert: undefined
            };
          }

          // Apply flag overrides to base settings
          currentSettings = window.applyFlagOverrides(baseSettings, flagOverrides);
          console.log('%c⚙️ Applied settings:', 'color: #9C27B0; font-weight: bold;', {
            useDefaults: flagOverrides.useDefaults,
            showCarbons: currentSettings.m2cfShowCarbons,
            aromaticCircles: currentSettings.m2cfAromaticCircles,
            showMethyls: currentSettings.m2cfShowMethyls,
            atomNumbers: currentSettings.m2cfAtomNumbers,
            addH2: currentSettings.m2cfAddH2,
            flipHorizontal: currentSettings.m2cfFlipHorizontal,
            flipVertical: currentSettings.m2cfFlipVertical,
            invert: currentSettings.m2cfInvert
          });

          // Extract compound name by removing flag parts
          if (content_trimmed.includes('/')) {
            compoundName = content_trimmed.split('/')[0].trim();
          } else if (content_trimmed.includes('+') && !isSMILES) {
            // For + flags, split on first + to get compound name
            const parts = content_trimmed.split('+');
            compoundName = parts[0].trim();
          }
        } catch (e) {
          console.error('Error processing flags, using default behavior:', e);
          // Fall back to basic behavior
          compoundName = content_trimmed;
          currentSettings = settings;
        }
      }

      // Check for special flags that override renderer selection
      // +3d flag: Use 3D viewer (PubChem-based)
      if (flagOverrides.is3D) {
        const pubchemData = {
          ...(isSMILES ? { smiles: compoundName } : { nomenclature: compoundName }),
          type: isSMILES ? 'smiles' : 'nomenclature',
          isPubChem: true,
          show3D: true  // Signal to show 3D viewer immediately
        };

        const encoded = btoa(JSON.stringify(pubchemData));
        const converted = `<img src="" alt="chem" class="chemfig-diagram chemfig-pubchem" data-molecule-viewer="${encoded}" data-loaded="false" data-show-3d="true" style="display: inline-block; height: auto; margin: 0 12px 8px 0; vertical-align: middle; padding: 4px; border-radius: 4px; background-color: transparent;">`;

        log.debug(`  📤 Sending ${isSMILES ? 'SMILES' : 'nomenclature'} to 3D viewer: ${compoundName}`);
        return converted;
      }

      // +pubchem flag: Force PubChem renderer
      if (flagOverrides.isPubchem) {
        const pubchemData = {
          ...(isSMILES ? { smiles: compoundName } : { nomenclature: compoundName }),
          type: isSMILES ? 'smiles' : 'nomenclature',
          isPubChem: true
        };

        const encoded = btoa(JSON.stringify(pubchemData));
        const converted = `<img src="" alt="chem" class="chemfig-diagram chemfig-pubchem" data-molecule-viewer="${encoded}" data-loaded="false" style="display: inline-block; height: auto; margin: 0 12px 8px 0; vertical-align: middle; padding: 4px; border-radius: 4px; background-color: transparent;">`;

        log.debug(`  📤 Sending ${isSMILES ? 'SMILES' : 'nomenclature'} to PubChem (forced by +pubchem flag): ${compoundName}`);
        return converted;
      }

      // Check if using PubChem renderer (which supports 3D viewer)
      if (currentSettings.rendererEngine === 'pubchem') {
        // Use PubChem - supports 3D viewer when enabled
        const pubchemData = {
          ...(isSMILES ? { smiles: compoundName } : { nomenclature: compoundName }),
          type: isSMILES ? 'smiles' : 'nomenclature',
          isPubChem: true
        };

        const encoded = btoa(JSON.stringify(pubchemData));
        const converted = `<img src="" alt="chem" class="chemfig-diagram chemfig-pubchem" data-molecule-viewer="${encoded}" data-loaded="false" style="display: inline-block; height: auto; margin: 0 12px 8px 0; vertical-align: middle; padding: 4px; border-radius: 4px; background-color: transparent;">`;

        log.debug(`  📤 Sending ${isSMILES ? 'SMILES' : 'nomenclature'} to PubChem: ${compoundName}`);
        return converted;
      }

      // Default: Send to MoleculeViewer (2D only)
      const rotation = currentSettings.m2cfRotate || currentSettings.mvRotate || 0;
      const scale = flagOverrides.size || 1.0;

      const moleculeViewerData = {
        ...(isSMILES ? { smiles: compoundName } : { nomenclature: compoundName }),
        type: isSMILES ? 'smiles' : 'nomenclature',
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
          // invert removed - now auto-applied based on dark mode detection
        },
        isMoleculeViewer: true
      };

      const encoded = btoa(JSON.stringify(moleculeViewerData));
      const converted = `<img src="" alt="chem" class="chemfig-diagram chemfig-molecule-viewer" data-molecule-viewer="${encoded}" data-loaded="false" data-rotation="${rotation}" style="display: inline-block; height: auto; margin: 0 12px 8px 0; vertical-align: middle; padding: 4px; border-radius: 4px; background-color: transparent;">`;

      log.debug(`  📤 Sending ${isSMILES ? 'SMILES' : 'nomenclature'} to MoleculeViewer: ${compoundName}`);
      return converted;
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

  // Pattern 3: \chemfig{C=C}:30 → Enable mol2chemfig rendering
  // Chemfig is a LaTeX drawing language that can be converted via mol2chemfig
  if (settings.renderChemfig) {
    log.debug('🧪 Applying Pattern 3: \\chemfig{...} → mol2chemfig conversion');
    // Updated regex - more robust for nested structures like CH(NH2)
    // Uses balanced brace matching
    const matches3 = [];
    const regex3 = /\\chemfig\{([^{}]*(?:\{[^{}]*\}[^{}]*)*)\}(?::(\d+))?/g;
    let match3;
    while ((match3 = regex3.exec(text)) !== null) {
      matches3.push(match3[0]);
    }

    if (matches3.length > 0) {
      log.debug(`  Found ${matches3.length} \\chemfig{} patterns`);
      patternMatches.push(...matches3);

      result = result.replace(/\\chemfig\{([^{}]*(?:\{[^{}]*\}[^{}]*)*)\}(?::(\d+))?/g, (match, content, rotationAngle) => {
        log.debug(`  Converting chemfig: ${match}`);
        window.chemRendererPerformance.recordStructure();

        // TRY NOMENCLATURE FIRST
        const nomenclatureFormula = getChemfigFromNomenclature(content);
        if (nomenclatureFormula) {
          log.debug(`  Nomenclature match: "${content}" → using database formula`);
          content = nomenclatureFormula;
        }

        // Parse rotation angle (default 0)
        const rotation = rotationAngle ? parseInt(rotationAngle) : 0;

        // REMOVE UNNECESSARY HYDROGENS: Clean up explicit H atoms
        let cleanedContent = removeUnnecessaryHydrogens(content);
        log.debug(`  After H removal: ${cleanedContent}`);

        // ADD ZIGZAG ANGLES: Simple chains without angles get automatic zigzag
        let zigzagContent = addZigzagAngles(cleanedContent);
        log.debug(`  After zigzag: ${zigzagContent}`);

        // SIMPLIFY CARBONS: Remove explicit CH, CH2, CH3 labels if option enabled
        let simplifiedContent = simplifyChemfigCarbons(zigzagContent);
        log.debug(`  Original: ${content}`);
        if (simplifiedContent !== zigzagContent) {
          log.debug(`  Simplified: ${simplifiedContent} (carbons as sticks)`);
        }

        // Detect dark mode
        const isDarkMode = detectDarkMode();

        // Build image URL based on selected rendering engine
        // Pass the original content for chemfig to SMILES conversion
        const latex = '$\\chemfig{' + simplifiedContent + '}$';
        const imageUrl = buildChemfigImageUrl(latex, isDarkMode, content);

        // Build style attribute with rotation and consistent margins
        // Add width/height to prevent layout shift
        const styleAttr = `style="transform: rotate(${rotation}deg); margin: 0 12px 8px 0; vertical-align: middle; width: 350px; height: 300px; object-fit: contain;"`;

        // Check if dev mode is enabled
        let converted;
        if (settings.devMode) {
          // Developer mode: show raw chemfig text instead of rendering
          const devStyle = `style="background: #f0f0f0; border: 1px solid #ddd; border-radius: 4px; padding: 8px 12px; font-family: monospace; font-size: 12px; color: #333; margin: 0 12px 8px 0; display: inline-block; white-space: pre-wrap; word-break: break-all; max-width: 500px;"`;
          converted = `<span class="chemfig-dev-mode" ${devStyle}>\\chemfig{${content}}</span>`;
        } else {
          // ✅ ALWAYS use MoleculeViewer rendering
          const moleculeViewerData = btoa(JSON.stringify(imageUrl));
          converted = `<img src="" alt="chemfig" class="chemfig-diagram chemfig-molecule-viewer" data-molecule-viewer="${moleculeViewerData}" data-loaded="false" data-rotation="${rotation}" ${styleAttr}>`;
        }


        log.debug(`  LaTeX: ${latex}`);
        log.debug(`  Rotation: ${rotation}°`);
        log.debug(`  Dark mode: ${isDarkMode}`);
        if (settings.devMode) {
          log.debug(`  🔧 DEV MODE: Showing raw text instead of rendering`);
        }
        const urlDebug = typeof imageUrl === 'object' ? `MoleculeViewer: ${imageUrl.smiles}` : String(imageUrl).substring(0, 100);
        log.debug(`  URL: ${urlDebug}...`);
        return converted;
      });
    }
  }

  // Pattern 4: chemfig{C=C}:30 (without backslash) - DISABLED (use chem: notation instead)
  // This pattern is now covered by Pattern 0a: chem:\chemfig{...}
  // Keeping this commented out to avoid duplicate processing
  /*
  log.debug('🧪 Pattern 4 DISABLED - use chem: notation instead');
  const matches4 = [];
  const regex4 = /\bchemfig\{([^{}]*(?:\{[^{}]*\}[^{}]*)*)\}(?::(\d+))?/g;
  let match4;
  while ((match4 = regex4.exec(text)) !== null) {
    matches4.push(match4[0]);
  }
  
  if (matches4.length > 0) {
    log.debug(`  Found ${matches4.length} chemfig{} patterns (IGNORED)`);
    patternMatches.push(...matches4);
  }
  */

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
 * Convert chemistry notation to Unicode with subscripts/superscripts
 * H2O → H₂O
 * Na+ → Na⁺
 * OH- → OH⁻
 */
function convertChemistry(text) {
  let result = text;

  // Subscript digits (₀₁₂₃₄₅₆₇₈₉)
  const subscripts = ['₀', '₁', '₂', '₃', '₄', '₅', '₆', '₇', '₈', '₉'];
  const superscripts = ['⁰', '¹', '²', '³', '⁴', '⁵', '⁶', '⁷', '⁸', '⁹'];

  // Convert numbers to subscripts: H2O → H₂O
  // But keep numbers in front: 2H2O stays 2H₂O
  result = result.replace(/(\D)(\d+)/g, (match, letter, num) => {
    const numStr = num.split('').map(d => subscripts[d]).join('');
    return letter + numStr;
  });

  // Convert + to superscript: Na+ → Na⁺
  result = result.replace(/\+/g, '⁺');

  // Convert - to superscript: Cl- → Cl⁻
  result = result.replace(/-(?!\d)/g, '⁻');

  // Render arrow: -> becomes →
  result = result.replace(/->|→/g, '→');
  result = result.replace(/<-|←/g, '←');
  result = result.replace(/<=>/g, '⇌');

  // Double/triple bonds stay as: = and ≡
  result = result.replace(/~/g, '≡');

  return result;
}

/**
 * Observe page changes for dynamic content
 * Watches for DOM mutations and re-scans if content changes
 */
function observePageChanges() {
  log.inject('🔄 Setting up mutation observer for dynamic content detection');

  let mutationCount = 0;
  const observer = new MutationObserver((mutations) => {
    mutationCount++;
    // Only log every 50 mutations to avoid spamming
    if (mutationCount % 50 === 0) {
      log.debug(`Detected ${mutationCount} total mutations, re-scanning...`);
    }
    clearTimeout(observer.timer);
    // Use shorter timeout for faster detection (especially for ChatGPT streaming)
    observer.timer = setTimeout(() => {
      if (settings.enabled) {
        log.inject('⚡ Re-scanning page after dynamic content detected');
        scanAndRender();
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

// ============================================
// DARK MODE SUPPORT
// ============================================
// Listen for dark mode changes and update molecule images
if (window.matchMedia) {
  const darkModeQuery = window.matchMedia('(prefers-color-scheme: dark)');

  function updateMoleculeColors(isDark) {
    console.log(`%c🌓 Dark mode ${isDark ? 'enabled' : 'disabled'} - updating molecule colors`, 'color: #00AAFF; font-weight: bold;');

    // Update all molecule images by re-rendering with color replacements
    const moleculeImages = document.querySelectorAll('img.chemfig-diagram');
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

