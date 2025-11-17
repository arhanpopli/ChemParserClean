/**
 * Image Size Controls for Chemistry Renderer
 * Adds up/down arrow controls to adjust molecule image sizes
 */

// ============================================
// SIZE CONTROL CONFIGURATION
// ============================================
const SIZE_STEP = 20; // Pixels to increase/decrease per click
const MIN_SIZE = 100; // Minimum size in pixels
const MAX_SIZE = 800; // Maximum size in pixels
const DEFAULT_WIDTH = 300;
const DEFAULT_HEIGHT = 200;

// ============================================
// SIZE STORAGE MANAGEMENT
// ============================================

/**
 * Generate a unique key for an image based on its SMILES or content
 */
function getImageKey(moleculeData) {
  if (!moleculeData) return null;

  // Use SMILES as the primary key if available
  if (moleculeData.smiles) {
    return `smiles:${moleculeData.smiles}`;
  }

  // Use nomenclature as fallback
  if (moleculeData.nomenclature) {
    return `nomenclature:${moleculeData.nomenclature}`;
  }

  return null;
}

/**
 * Generate a unique key for a specific image on a specific page
 */
function getPageImageKey(moleculeData, pageUrl) {
  const imageKey = getImageKey(moleculeData);
  if (!imageKey) return null;

  // Include page URL for per-page storage
  return `page:${pageUrl}:${imageKey}`;
}

/**
 * Load saved size for an image
 */
async function loadImageSize(moleculeData, pageUrl, settings) {
  try {
    const imageKey = getImageKey(moleculeData);
    if (!imageKey) return { width: DEFAULT_WIDTH, height: DEFAULT_HEIGHT };

    // Check if size saving is enabled
    if (!settings.saveSizePerImage && !settings.saveSizeBySMILES) {
      return { width: DEFAULT_WIDTH, height: DEFAULT_HEIGHT };
    }

    // Determine which storage key to use
    let storageKey;
    if (settings.saveSizePerImage && !settings.saveSizeBySMILES) {
      // Save per page
      storageKey = getPageImageKey(moleculeData, pageUrl);
    } else {
      // Save globally by SMILES (saveSizeBySMILES is enabled)
      storageKey = imageKey;
    }

    if (!storageKey) return { width: DEFAULT_WIDTH, height: DEFAULT_HEIGHT };

    // Load from chrome.storage
    return new Promise((resolve) => {
      chrome.storage.local.get([storageKey], (result) => {
        if (result[storageKey]) {
          resolve(result[storageKey]);
        } else {
          resolve({ width: DEFAULT_WIDTH, height: DEFAULT_HEIGHT });
        }
      });
    });
  } catch (error) {
    console.error('Error loading image size:', error);
    return { width: DEFAULT_WIDTH, height: DEFAULT_HEIGHT };
  }
}

/**
 * Save size for an image
 */
async function saveImageSize(moleculeData, pageUrl, size, settings) {
  try {
    const imageKey = getImageKey(moleculeData);
    if (!imageKey) return;

    // Check if size saving is enabled
    if (!settings.saveSizePerImage && !settings.saveSizeBySMILES) {
      return;
    }

    // Determine which storage key to use
    let storageKey;
    if (settings.saveSizePerImage && !settings.saveSizeBySMILES) {
      // Save per page
      storageKey = getPageImageKey(moleculeData, pageUrl);
    } else {
      // Save globally by SMILES
      storageKey = imageKey;
    }

    if (!storageKey) return;

    // Save to chrome.storage
    const data = {};
    data[storageKey] = size;
    chrome.storage.local.set(data, () => {
      console.log(`Saved size for ${storageKey}:`, size);
    });
  } catch (error) {
    console.error('Error saving image size:', error);
  }
}

// ============================================
// SIZE CONTROL UI
// ============================================

/**
 * Create size control buttons
 */
function createSizeControls(container, svgImg, moleculeData, settings) {
  // Create controls container
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
    z-index: 10;
  `;

  // Create up arrow button
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

  // Create down arrow button
  const downButton = document.createElement('button');
  downButton.className = 'chem-size-btn chem-size-down';
  downButton.innerHTML = '▼';
  downButton.title = 'Decrease size';
  downButton.style.cssText = upButton.style.cssText;

  // Add hover effects
  [upButton, downButton].forEach(btn => {
    btn.addEventListener('mouseenter', () => {
      btn.style.background = 'rgba(0, 0, 0, 0.9)';
    });
    btn.addEventListener('mouseleave', () => {
      btn.style.background = 'rgba(0, 0, 0, 0.7)';
    });
  });

  // Add click handlers
  upButton.addEventListener('click', (e) => {
    e.stopPropagation();
    adjustImageSize(container, svgImg, moleculeData, SIZE_STEP, settings);
  });

  downButton.addEventListener('click', (e) => {
    e.stopPropagation();
    adjustImageSize(container, svgImg, moleculeData, -SIZE_STEP, settings);
  });

  // Append buttons to controls
  controlsDiv.appendChild(upButton);
  controlsDiv.appendChild(downButton);

  // Show controls on hover
  container.addEventListener('mouseenter', () => {
    controlsDiv.style.opacity = '1';
  });
  container.addEventListener('mouseleave', () => {
    controlsDiv.style.opacity = '0';
  });

  return controlsDiv;
}

/**
 * Adjust image size
 */
function adjustImageSize(container, svgImg, moleculeData, delta, settings) {
  // Get current size
  const currentWidth = parseInt(svgImg.style.maxWidth) || DEFAULT_WIDTH;
  const currentHeight = parseInt(svgImg.style.maxHeight) || DEFAULT_HEIGHT;

  // Calculate new size
  let newWidth = currentWidth + delta;
  let newHeight = currentHeight + Math.round(delta * (currentHeight / currentWidth));

  // Clamp to min/max
  newWidth = Math.max(MIN_SIZE, Math.min(MAX_SIZE, newWidth));
  newHeight = Math.max(MIN_SIZE, Math.min(MAX_SIZE, newHeight));

  // Apply new size
  svgImg.style.maxWidth = `${newWidth}px`;
  svgImg.style.maxHeight = `${newHeight}px`;

  // Save size if enabled
  const pageUrl = window.location.href;
  const size = { width: newWidth, height: newHeight };
  saveImageSize(moleculeData, pageUrl, size, settings);

  console.log(`Adjusted size: ${newWidth}x${newHeight}`);
}

/**
 * Wrap an image with a container and add size controls
 */
async function wrapImageWithSizeControls(svgImg, originalImg, moleculeData, settings) {
  try {
    // Create container
    const container = document.createElement('div');
    container.className = 'chem-image-container';
    container.style.cssText = `
      position: relative;
      display: inline-block;
      margin: 0 12px 8px 0;
      vertical-align: middle;
    `;

    // Load saved size
    const pageUrl = window.location.href;
    const savedSize = await loadImageSize(moleculeData, pageUrl, settings);

    // Apply saved size to image
    svgImg.style.maxWidth = `${savedSize.width}px`;
    svgImg.style.maxHeight = `${savedSize.height}px`;

    // Store molecule data on the image for later reference
    if (moleculeData) {
      svgImg.dataset.moleculeData = JSON.stringify(moleculeData);
    }

    // Create size controls
    const controls = createSizeControls(container, svgImg, moleculeData, settings);

    // Insert container before original image
    originalImg.parentNode.insertBefore(container, originalImg);

    // Move SVG image into container
    container.appendChild(svgImg);

    // Add controls to container
    container.appendChild(controls);

    // Remove original image
    originalImg.remove();

    return container;
  } catch (error) {
    console.error('Error wrapping image with size controls:', error);
    // Fallback: just replace the image directly
    originalImg.parentNode.replaceChild(svgImg, originalImg);
    return svgImg;
  }
}
