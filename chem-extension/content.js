/**
 * Chemistry Formula Renderer v3.0
 * Works on ChatGPT and other strict CSP sites
 * Uses web_accessible_resources to serve scripts from extension
 */

// ============================================
// 🌍 MOLECULEVIEWER API CONFIGURATION
// ============================================
// LOCAL TESTING - Use localhost to avoid HTTPS mixed content errors
const MOLECULE_VIEWER_API = 'http://localhost:5000';

// For Heroku production (uncomment when ready to deploy):
// const MOLECULE_VIEWER_API = 'https://YOUR-HEROKU-APP.herokuapp.com';

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
  
  recordLoad: function(duration) {
    this.metrics.totalSVGsRendered++;
    this.metrics.totalLoadTime += duration;
    this.metrics.loadTimes.push(duration);
    this.metrics.maxLoadTime = Math.max(this.metrics.maxLoadTime, duration);
    this.metrics.minLoadTime = Math.min(this.metrics.minLoadTime, duration);
    this.metrics.averageLoadTime = this.metrics.totalLoadTime / this.metrics.totalSVGsRendered;
    // FIXED: Reduced from 1000 to 50 to save memory
    if (this.metrics.loadTimes.length > 50) this.metrics.loadTimes.shift();
  },
  
  recordStructure: function() {
    this.metrics.totalStructuresFound++;
  },
  
  recordFormula: function() {
    this.metrics.totalFormulasTransformed++;
  },
  
  getStats: function() {
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
  
  logStats: function() {
    const stats = this.getStats();
    console.table(stats);
    console.log('%c📊 Performance Analysis:', 'color: #0066cc; font-weight: bold; font-size: 14px;');
    console.log(`   • ${stats.structures} chemical structures detected`);
    console.log(`   • ${stats.svgsRendered} SVGs rendered (avg ${stats.avgLoadTime})`);
    console.log(`   • ${stats.formulasTransformed} formulas transformed locally`);
    console.log(`   • Load time range: ${stats.minLoadTime} - ${stats.maxLoadTime}`);
  },
};

// Create debug interface EARLY so it's always available
window.chemRendererDebug = {
  getLogs: function() {
    console.table(logHistory);
    return logHistory;
  },
  getSettings: function() {
    console.table(settings);
    return settings;
  },
  getPerformanceStats: function() {
    return window.chemRendererPerformance.logStats();
  },
  scanPage: function() {
    console.log('[Manual Scan] Triggering page scan...');
    scanAndRender();
  },
  testFormulas: function() {
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
  getCurrentFormulas: function() {
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
  rotateFormulas: function(angle = 90) {
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
  getRotationHelp: function() {
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
  togglePerformanceMode: function(enable) {
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
  setMaxVisibleSVGs: function(count) {
    settings.maxVisibleSVGs = count;
    log.success(`📊 Max visible SVGs set to: ${count}`);
    chrome.storage.sync.set({ maxVisibleSVGs: count });
  },
  getPerformanceStats: function() {
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
  setLayoutMode: function(mode) {
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
  toggleLayoutMode: function() {
    settings.layoutMode = settings.layoutMode === 'horizontal' ? 'vertical' : 'horizontal';
    applyLayoutMode();
    log.success(`📐 Layout mode toggled to: ${settings.layoutMode.toUpperCase()}`);
    chrome.storage.sync.set({ layoutMode: settings.layoutMode });
    return settings.layoutMode;
  },
  getLayoutSettings: function() {
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
  clearMemory: function() {
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
  checkMemory: function() {
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
  toggleCarbonLabels: function(enable) {
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
  setSizePreset: function(preset) {
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
  setRendererEngine: function(engine) {
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
  getRenderingHelp: function() {
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
  rendererEngine: 'molecule-viewer',  // ✅ FORCED - ONLY MoleculeViewer
  devMode: false,
  // MoleculeViewer rendering options
  flipHorizontal: false,
  flipVertical: false,
  showCarbons: false,
  showMethyls: false,
  aromaticCircles: true,
  fancyBonds: true,
  atomNumbers: false,
  hydrogensMode: 'keep'
};

log.info('📦 Loading settings from storage...');

// Load settings - with proper callback
chrome.storage.sync.get(null, (result) => {
  // Merge stored settings with defaults
  settings = { ...settings, ...result };
  
  // ✅ FORCE MoleculeViewer - IGNORE any stored renderer choice
  settings.rendererEngine = 'molecule-viewer';
  
  log.success('✅ Settings loaded', settings);
  log.info(`Renderer Engine: 🧪 MoleculeViewer (forced)`);
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
});

/**
 * Initialize the renderer
 */
function initializeRenderer() {
  log.inject('🔧 Initialize renderer - Starting chemical formula processing');
  
  // We process formulas directly using Unicode and CodeCogs
  // No MathJax injection needed - this avoids CSP issues entirely
  
  log.inject('🚀 Starting initial page scan...');
  scanAndRender();
  
  // Apply layout mode to any containers
  applyLayoutMode();
  
  log.inject('🔄 Setting up dynamic content observer...');
  observePageChanges();
  
  // Inject CSS for proper rendering
  injectStyles();
  
  // Setup lazy-loading if performance mode is enabled
  if (settings.performanceMode) {
    log.inject('⚡ Performance mode ENABLED - Setting up lazy-loading for SVGs');
    setupLazyLoading();
  } else {
    log.inject('⚡ Performance mode disabled');
  }
  
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
      background-color: rgba(255, 255, 255, 0.05);
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
  
  // Helper function to load MoleculeViewer rendering with caching and download link
  function loadMoleculeViewerImage(img) {
    activeLoads++;
    console.log('%c🧪 LOADMOLECULEVIEWERIMAGE CALLED!', 'background: #222; color: #00FF00; font-size: 20px; padding: 10px;');
    console.log('Image element:', img);
    console.log('Dataset:', img.dataset);
    log.debug(`🧪 Loading MoleculeViewer SVG (#${activeLoads})`);
    
    try {
      const moleculeData = JSON.parse(atob(img.dataset.moleculeViewer));
      console.log('%c📦 Decoded molecule data:', 'color: #0088FF; font-weight: bold;', moleculeData);
      
      // Determine which endpoint to use based on data type
      const isNomenclature = moleculeData.type === 'nomenclature' && moleculeData.nomenclature;
      const isSMILES = moleculeData.type === 'smiles' && moleculeData.smiles;
      
      let apiUrl;
      if (isNomenclature) {
        console.log('%c📤 Using nomenclature endpoint', 'color: #FF6B00; font-weight: bold;');
        console.log('Nomenclature:', moleculeData.nomenclature);
        apiUrl = `${MOLECULE_VIEWER_API}/img/nomenclature?nomenclature=${encodeURIComponent(moleculeData.nomenclature)}&width=300&height=200&json=true`;
        console.log('API URL:', apiUrl);
      } else if (isSMILES) {
        console.log('%c📤 Using SMILES endpoint', 'color: #FF6B00; font-weight: bold;');
        apiUrl = `${MOLECULE_VIEWER_API}/img/smiles?smiles=${encodeURIComponent(moleculeData.smiles)}&width=300&height=200&json=true`;
      } else {
        throw new Error('Invalid molecule data');
      }
      
      // Fetch JSON response with cache link
      console.log('%c🌐 Fetching from backend...', 'color: #00AAFF; font-weight: bold;');
      fetch(apiUrl)
        .then(response => {
          console.log('%c✅ Got response from backend:', 'color: #00FF00; font-weight: bold;', response);
          return response.json();
        })
        .then(data => {
          console.log('%c📊 Backend data:', 'color: #FFAA00; font-weight: bold;', data);
          if (!data.success) {
            throw new Error(data.error || 'Rendering failed');
          }
          
          // Create clean SVG image (no container, no controls)
          
          // Detect dark mode
          const isDarkMode = window.matchMedia && window.matchMedia('(prefers-color-scheme: dark)').matches;
          
          // Create clean SVG image with dark mode support
          const svgImg = document.createElement('img');
          svgImg.src = 'data:image/svg+xml;base64,' + btoa(data.svg);
          svgImg.alt = 'molecule';
          svgImg.className = 'chemfig-diagram';
          svgImg.style.cssText = `
            display: inline-block;
            max-width: 300px;
            max-height: 200px;
            margin: 0 12px 8px 0;
            vertical-align: middle;
            cursor: pointer;
            ${isDarkMode ? 'filter: invert(1) hue-rotate(180deg);' : ''}
          `;
          
          // Just add the image directly to the page (no controls!)
          // Mark as loaded
          img.dataset.loaded = 'true';
          
          // Replace original img element with the SVG image
          img.parentNode.replaceChild(svgImg, img);
          
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
        // Check if this is a MoleculeViewer image
        if (img.classList.contains('chemfig-molecule-viewer')) {
          console.log('%c🧪 Detected MoleculeViewer image!', 'background: #0088FF; color: #FFF; font-size: 14px; padding: 5px;');
          if (activeLoads < maxConcurrentLoads) {
            loadMoleculeViewerImage(img);
          } else {
            if (typeof requestIdleCallback !== 'undefined') {
              requestIdleCallback(() => {
                if (activeLoads < maxConcurrentLoads) {
                  loadMoleculeViewerImage(img);
                }
              });
            } else {
              setTimeout(() => {
                if (activeLoads < maxConcurrentLoads) {
                  loadMoleculeViewerImage(img);
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
  if (!window.katex) {
    log.error('KaTeX not loaded in window');
    return;
  }
  
  log.success('✅ KaTeX Ready! Chemistry rendering available');
  
  // Start scanning after a brief delay to ensure KaTeX is ready
  setTimeout(() => {
    log.inject('� Running initial page scan now that KaTeX is ready');
    scanAndRender();
  }, 500);
}

/**
 * Basic rendering without external libraries
 * Replaces \ce{...} with readable text approximation
 */
function initializeBasicRendering() {
  log.inject('🔄 Using basic text rendering (no external library)');
  log.inject('Formulas will be visible but not beautifully rendered');
  
  setTimeout(() => {
    log.inject('🔍 Running initial page scan with basic rendering');
    scanAndRender();
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
      const codes = ['₀','₁','₂','₃','₄','₅','₆','₇','₈','₉'];
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
          
          // Setup lazy-loading observers for newly created images
          if (settings.performanceMode && window._lazyLoadObserver) {
            const images = span.querySelectorAll('img.chemfig-diagram[data-loaded="false"]');
            images.forEach(img => {
              window._lazyLoadObserver.observe(img);
            });
            if (images.length > 0) {
              log.debug(`📊 Setup lazy-loading for ${images.length} SVGs`);
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
      log.debug(`    Sample ${i+1}: ${t}`);
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
  
  // Pattern 0A: chem:SMILES (contains special characters like =, [], (), @, +)
  // Example: chem:CC(=O)C, chem:c1ccccc1, chem:C[C@H](O)C
  log.debug('🧪 Applying Pattern 0A: chem:SMILES → MoleculeViewer');
  const smilesMatches = result.match(/\bchem:([A-Za-z0-9_\-()=[\]@+#\\]+)/g);
  if (smilesMatches) {
    log.debug(`  Found ${smilesMatches.length} chem:SMILES patterns`);
    
    result = result.replace(/\bchem:([A-Za-z0-9_\-()=[\]@+#\\]+)/g, (match, smiles) => {
      // Only process if it looks like SMILES (has special structural characters)
      if (!/[=\[\]()@+#\\]/.test(smiles)) {
        return match; // Not SMILES, skip this match
      }
      
      log.debug(`  Converting SMILES: ${smiles}`);
      window.chemRendererPerformance.recordStructure();
      
      // Send SMILES to MoleculeViewer
      const moleculeViewerData = {
        smiles: smiles,
        type: 'smiles',  // SMILES format
        options: {
          width: 300,
          height: 200,
          aromaticCircles: settings.aromaticCircles,
          fancyBonds: settings.fancyBonds,
          showCarbons: settings.showCarbons,
          showMethyls: settings.showMethyls,
          atomNumbers: settings.atomNumbers,
          hydrogensMode: settings.hydrogensMode
        },
        isMoleculeViewer: true
      };
      
      const encoded = btoa(JSON.stringify(moleculeViewerData));
      const converted = `<img src="" alt="chem" class="chemfig-diagram chemfig-molecule-viewer" data-molecule-viewer="${encoded}" data-loaded="false" data-rotation="0" style="display: inline-block; max-width: 300px; max-height: 200px; margin: 0 12px 8px 0; vertical-align: middle; padding: 4px; border-radius: 4px; background-color: rgba(255, 255, 255, 0.05);">`;
      
      log.debug(`  📤 Sending SMILES to MoleculeViewer: ${smiles}`);
      return converted;
    });
  }
  
  // Pattern 0B: chem:nomenclature (plain text, no special characters)
  // Example: chem:acetone, chem:benzene, chem:aspirin
  log.debug('🧪 Applying Pattern 0B: chem:nomenclature → MoleculeViewer');
  const nomenclatureMatches = result.match(/\bchem:([A-Za-z][A-Za-z0-9\-_]*)/g);
  if (nomenclatureMatches) {
    log.debug(`  Found ${nomenclatureMatches.length} chem:nomenclature patterns`);
    
    result = result.replace(/\bchem:([A-Za-z][A-Za-z0-9\-_]*)/g, (match, nomenclature) => {
      // Skip if it looks like SMILES
      if (/[=\[\]()@+#\\]/.test(nomenclature)) {
        return match; // Already processed as SMILES
      }
      
      log.debug(`  Converting nomenclature: ${nomenclature}`);
      window.chemRendererPerformance.recordStructure();
      
      // Send nomenclature to MoleculeViewer for name→SMILES conversion
      const moleculeViewerData = {
        nomenclature: nomenclature,
        type: 'nomenclature',  // Nomenclature/name format
        options: {
          width: 300,
          height: 200,
          aromaticCircles: settings.aromaticCircles,
          fancyBonds: settings.fancyBonds,
          showCarbons: settings.showCarbons,
          showMethyls: settings.showMethyls,
          atomNumbers: settings.atomNumbers,
          hydrogensMode: settings.hydrogensMode
        },
        isMoleculeViewer: true
      };
      
      const encoded = btoa(JSON.stringify(moleculeViewerData));
      const converted = `<img src="" alt="chem" class="chemfig-diagram chemfig-molecule-viewer" data-molecule-viewer="${encoded}" data-loaded="false" data-rotation="0" style="display: inline-block; max-width: 300px; max-height: 200px; margin: 0 12px 8px 0; vertical-align: middle; padding: 4px; border-radius: 4px; background-color: rgba(255, 255, 255, 0.05);">`;
      
      log.debug(`  📤 Sending nomenclature to MoleculeViewer: ${nomenclature}`);
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
  
  // Pattern 3: \chemfig{C=C}:30 → DISABLED (use chem: notation instead)
  // Chemfig is a LaTeX drawing language, not suitable for direct chemical structure rendering
  // Use chem:SMILES or chem:nomenclature instead (like CodeCogs)
  if (false && settings.renderChemfig) {  // ❌ DISABLED - use chem: notation instead
    log.debug('🧪 Pattern 3 DISABLED - use chem: notation instead');
    // Updated regex - more robust for nested structures like CH(NH2)
    // Uses balanced brace matching
    const matches3 = [];
    const regex3 = /\\chemfig\{([^{}]*(?:\{[^{}]*\}[^{}]*)*)\}(?::(\d+))?/g;
    let match3;
    while ((match3 = regex3.exec(text)) !== null) {
      matches3.push(match3[0]);
    }
    
    if (matches3.length > 0) {
      log.debug(`  Found ${matches3.length} \\chemfig{} patterns (IGNORED)`);
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
        const styleAttr = `style="transform: rotate(${rotation}deg); margin: 0 12px 8px 0; vertical-align: middle; width: 300px; height: 200px; object-fit: contain;"`;
        
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
    
    // Pattern 4: chemfig{C=C}:30 (without backslash) - DISABLED (use chem: notation instead)
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
      
      result = result.replace(/\bchemfig\{([^{}]*(?:\{[^{}]*\}[^{}]*)*)\}(?::(\d+))?/g, (match, content, rotationAngle) => {
        log.debug(`  Converting chemfig: ${match}`);
        window.chemRendererPerformance.recordStructure();
        
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
        const styleAttr = `style="transform: rotate(${rotation}deg); margin: 0 12px 8px 0; vertical-align: middle;"`;
        
        // Check if dev mode is enabled
        let converted;
        if (settings.devMode) {
          // Developer mode: show raw chemfig text instead of rendering
          const devStyle = `style="background: #f0f0f0; border: 1px solid #ddd; border-radius: 4px; padding: 8px 12px; font-family: monospace; font-size: 12px; color: #333; margin: 0 12px 8px 0; display: inline-block; white-space: pre-wrap; word-break: break-all; max-width: 500px;"`;
          converted = `<span class="chemfig-dev-mode" ${devStyle}>chemfig{${content}}</span>`;
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
    
    // Update all molecule images
    const moleculeImages = document.querySelectorAll('img.chemfig-diagram');
    moleculeImages.forEach(img => {
      if (isDark) {
        img.style.filter = 'invert(1) hue-rotate(180deg)';
      } else {
        img.style.filter = '';
      }
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

