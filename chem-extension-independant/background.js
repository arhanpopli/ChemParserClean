/**
 * Background Service Worker for Chemistry Formula Renderer
 * Handles fetch requests from content scripts to bypass CSP restrictions
 */

console.log('[ChemistryLaTeX] Background service worker loaded');

// Listen for messages from content scripts
chrome.runtime.onMessage.addListener((request, sender, sendResponse) => {
  if (request.type === 'FETCH_API') {
    // Handle API fetch requests (JSON)
    console.log('[Background] FETCH_API request:', request.url);
    handleFetchRequest(request.url, request.options)
      .then(response => {
        console.log('[Background] FETCH_API success');
        sendResponse({ success: true, data: response });
      })
      .catch(error => {
        console.error('[Background] FETCH_API error:', error.message);
        sendResponse({ success: false, error: error.message });
      });
    return true; // Keep the message channel open for async response
  }

  if (request.type === 'FETCH_BLOB') {
    // Handle blob fetch requests (for images)
    console.log('[Background] FETCH_BLOB request:', request.url);
    handleBlobRequest(request.url)
      .then(response => {
        console.log('[Background] FETCH_BLOB success');
        sendResponse({ success: true, data: response });
      })
      .catch(error => {
        console.error('[Background] FETCH_BLOB error:', error.message);
        sendResponse({ success: false, error: error.message });
      });
    return true;
  }

  if (request.type === 'FETCH_TEXT') {
    // Handle text fetch requests
    console.log('[Background] FETCH_TEXT request:', request.url);
    handleTextRequest(request.url, request.options)
      .then(response => {
        console.log('[Background] FETCH_TEXT success');
        sendResponse({ success: true, data: response });
      })
      .catch(error => {
        console.error('[Background] FETCH_TEXT error:', error.message);
        sendResponse({ success: false, error: error.message });
      });
    return true;
  }
});

/**
 * Fetch JSON data from API
 */
async function handleFetchRequest(url, options = {}) {
  try {
    console.log('[Background] handleFetchRequest called with:', { url, options });

    // Build fetch options - handle body for POST requests
    const fetchOptions = {
      method: options.method || 'GET',
      headers: {
        'Accept': 'application/json',
        ...options.headers
      }
    };

    // Only add body for non-GET requests
    if (options.body && options.method && options.method.toUpperCase() !== 'GET') {
      fetchOptions.body = options.body;
    }

    console.log('[Background] Final fetch options:', fetchOptions);

    const response = await fetch(url, fetchOptions);
    console.log('[Background] Response status:', response.status, response.statusText);

    if (!response.ok) {
      const errorText = await response.text();
      console.error('[Background] Error response body:', errorText);
      throw new Error(`HTTP ${response.status}: ${response.statusText}`);
    }

    const data = await response.json();
    console.log('[Background] Response data received, keys:', Object.keys(data));
    return data;
  } catch (error) {
    console.error('[Background] Fetch error:', error);
    throw error;
  }
}

/**
 * Fetch blob data (images) and convert to base64
 */
async function handleBlobRequest(url) {
  try {
    const response = await fetch(url);

    if (!response.ok) {
      throw new Error(`HTTP ${response.status}: ${response.statusText}`);
    }

    const blob = await response.blob();
    const base64 = await blobToBase64(blob);

    return {
      base64: base64,
      type: blob.type
    };
  } catch (error) {
    console.error('[Background] Blob fetch error:', error);
    throw error;
  }
}

/**
 * Fetch text data from API
 * Service workers have full network access via host_permissions: "<all_urls>"
 * They are NOT restricted by the extension's CSP connect-src directive
 */
async function handleTextRequest(url, options = {}) {
  console.log('[Background] handleTextRequest called with:', { url, options });

  try {
    // Validate URL
    let parsedUrl;
    try {
      parsedUrl = new URL(url);
      console.log('[Background] Parsed URL:', parsedUrl.href, 'Protocol:', parsedUrl.protocol);
    } catch (urlError) {
      console.error('[Background] Invalid URL:', url, urlError);
      throw new Error(`Invalid URL: ${url}`);
    }

    // Build clean fetch options - don't blindly spread options to avoid issues
    const fetchOptions = {
      method: options.method || 'GET',
      headers: {
        'Accept': 'text/plain, */*',
        ...(options.headers || {})
      },
      // Important: Use cors mode for cross-origin requests
      mode: 'cors',
      // Don't send credentials for cross-origin requests (avoids CORS preflight issues)
      credentials: 'omit'
    };

    // Only add body for non-GET requests
    if (options.body && fetchOptions.method.toUpperCase() !== 'GET') {
      fetchOptions.body = options.body;
    }

    console.log('[Background] Fetching with options:', JSON.stringify(fetchOptions));

    const response = await fetch(url, fetchOptions);
    console.log('[Background] Response received:', {
      status: response.status,
      statusText: response.statusText,
      ok: response.ok,
      type: response.type,
      url: response.url,
      headers: Object.fromEntries([...response.headers.entries()])
    });

    if (!response.ok) {
      const errorBody = await response.text().catch(() => '(could not read error body)');
      console.error('[Background] HTTP Error response body:', errorBody);
      throw new Error(`HTTP ${response.status}: ${response.statusText} - ${errorBody.substring(0, 200)}`);
    }

    const text = await response.text();
    console.log('[Background] Text response received, length:', text.length, 'preview:', text.substring(0, 150));

    // Validate that we got SVG content
    if (text.includes('<?xml') || text.includes('<svg')) {
      console.log('[Background] ✅ Valid SVG detected');
    } else {
      console.warn('[Background] ⚠️ Response may not be valid SVG');
    }

    return text;
  } catch (error) {
    console.error('[Background] Text fetch error:', {
      message: error.message,
      name: error.name,
      stack: error.stack,
      url: url
    });
    throw error;
  }
}

/**
 * Convert blob to base64 string
 */
function blobToBase64(blob) {
  return new Promise((resolve, reject) => {
    const reader = new FileReader();
    reader.onloadend = () => {
      // Remove the data URL prefix to get just the base64
      const base64 = reader.result.split(',')[1];
      resolve(base64);
    };
    reader.onerror = reject;
    reader.readAsDataURL(blob);
  });
}

// Listen for tab updates
chrome.tabs.onUpdated.addListener((tabId, changeInfo, tab) => {
  if (changeInfo.status === 'complete') {
    console.log('[ChemistryLaTeX] Tab updated:', tab.url);
  }
});

// Create context menu items
chrome.runtime.onInstalled.addListener(() => {
  // ===== TEXT SELECTION CONTEXT MENUS =====
  // Option 1: Render as Molecule (treats text as chemical name/nomenclature)
  chrome.contextMenus.create({
    id: "inspect-molecule",
    title: "Render as Molecule: '%s'",
    contexts: ["selection"]
  });

  // Option 2: Render as SMILES (treats text directly as SMILES string)
  chrome.contextMenus.create({
    id: "render-smiles",
    title: "Render as SMILES: '%s'",
    contexts: ["selection"]
  });

  // Option 3: Render as Biomolecule (goes directly to RCSB PDB)
  chrome.contextMenus.create({
    id: "render-biomolecule",
    title: "Render as Biomolecule: '%s'",
    contexts: ["selection"]
  });

  // Option 4: Render as Mineral (goes directly to COD crystal database)
  chrome.contextMenus.create({
    id: "render-mineral",
    title: "Render as Mineral: '%s'",
    contexts: ["selection"]
  });

  // Option 5: Render as IUPAC (goes to OPSIN for precise nomenclature)
  chrome.contextMenus.create({
    id: "render-iupac",
    title: "Render as IUPAC: '%s'",
    contexts: ["selection"]
  });

  // ===== IMAGE RE-RENDER CONTEXT MENUS =====
  // Parent menu for image re-rendering options
  chrome.contextMenus.create({
    id: "rerender-parent",
    title: "ChemistryLaTeX: Re-render as...",
    contexts: ["image"]
  });

  // Re-render as Molecule (PubChem lookup)
  chrome.contextMenus.create({
    id: "rerender-molecule",
    parentId: "rerender-parent",
    title: "Molecule (PubChem)",
    contexts: ["image"]
  });


  // Re-render as SMILES
  chrome.contextMenus.create({
    id: "rerender-smiles",
    parentId: "rerender-parent",
    title: "SMILES (direct render)",
    contexts: ["image"]
  });

  // Re-render as Biomolecule
  chrome.contextMenus.create({
    id: "rerender-biomolecule",
    parentId: "rerender-parent",
    title: "Biomolecule (RCSB PDB)",
    contexts: ["image"]
  });

  // Re-render as Mineral
  chrome.contextMenus.create({
    id: "rerender-mineral",
    parentId: "rerender-parent",
    title: "Mineral (COD)",
    contexts: ["image"]
  });

  // Separator
  chrome.contextMenus.create({
    id: "separator-flags",
    parentId: "rerender-parent",
    type: "separator",
    contexts: ["image"]
  });

  // Edit Flags (opens dialog to toggle individual flags)
  chrome.contextMenus.create({
    id: "edit-flags",
    parentId: "rerender-parent",
    title: "Edit Flags...",
    contexts: ["image"]
  });
});



// Handle SETTINGS_CHANGED messages from popup - broadcast to all tabs
chrome.runtime.onMessage.addListener((request, sender, sendResponse) => {
  if (request.type === 'SETTINGS_CHANGED') {
    console.log('[Background] Broadcasting settings change to all tabs:', request.settings, 'Changed keys:', request.changedKeys);

    // Send to all tabs
    chrome.tabs.query({}, (tabs) => {
      tabs.forEach(tab => {
        // Skip chrome:// and edge:// URLs
        if (tab.url && !tab.url.startsWith('chrome://') && !tab.url.startsWith('edge://')) {
          chrome.tabs.sendMessage(tab.id, {
            type: 'APPLY_SETTINGS',
            settings: request.settings,
            changedSettings: request.changedKeys || []  // Forward the list of changed setting keys
          }).catch(() => {
            // Tab might not have content script loaded, ignore
          });
        }
      });
    });

    sendResponse({ success: true });
    return true;
  }

  // Handle RELOAD_ALL_IMAGES from popup - broadcast to all tabs
  if (request.type === 'RELOAD_ALL_IMAGES') {
    console.log('[Background] Broadcasting reload all images to all tabs');

    chrome.tabs.query({}, (tabs) => {
      tabs.forEach(tab => {
        if (tab.url && !tab.url.startsWith('chrome://') && !tab.url.startsWith('edge://')) {
          chrome.tabs.sendMessage(tab.id, {
            type: 'RELOAD_ALL_IMAGES'
          }).catch(() => {
            // Tab might not have content script loaded, ignore
          });
        }
      });
    });

    sendResponse({ success: true });
    return true;
  }

  // Handle CLEAR_CACHE from popup - broadcast to all tabs
  if (request.type === 'CLEAR_CACHE') {
    console.log('[Background] Broadcasting clear cache to all tabs');

    chrome.tabs.query({}, (tabs) => {
      tabs.forEach(tab => {
        if (tab.url && !tab.url.startsWith('chrome://') && !tab.url.startsWith('edge://')) {
          chrome.tabs.sendMessage(tab.id, {
            type: 'CLEAR_CACHE'
          }).catch(() => {
            // Tab might not have content script loaded, ignore
          });
        }
      });
    });

    sendResponse({ success: true });
    return true;
  }
});

// Handle context menu clicks
chrome.contextMenus.onClicked.addListener((info, tab) => {
  if (info.menuItemId === "inspect-molecule") {
    const selectedText = info.selectionText;
    console.log('[Background] Inspecting molecule (nomenclature):', selectedText);

    // Send message to content script to handle the inspection
    chrome.tabs.sendMessage(tab.id, {
      type: 'INSPECT_MOLECULE',
      text: selectedText
    });
  }

  if (info.menuItemId === "render-smiles") {
    const selectedText = info.selectionText;
    console.log('[Background] Rendering as SMILES:', selectedText);

    // Send message to content script to handle SMILES rendering
    chrome.tabs.sendMessage(tab.id, {
      type: 'RENDER_SMILES',
      text: selectedText
    });
  }

  if (info.menuItemId === "render-biomolecule") {
    const selectedText = info.selectionText;
    console.log('[Background] Rendering as Biomolecule (RCSB PDB):', selectedText);

    // Send message to content script to handle biomolecule rendering
    chrome.tabs.sendMessage(tab.id, {
      type: 'RENDER_BIOMOLECULE',
      text: selectedText
    });
  }

  if (info.menuItemId === "render-mineral") {
    const selectedText = info.selectionText;
    console.log('[Background] Rendering as Mineral (COD):', selectedText);

    // Send message to content script to handle mineral rendering
    chrome.tabs.sendMessage(tab.id, {
      type: 'RENDER_MINERAL',
      text: selectedText
    });
  }

  if (info.menuItemId === "render-iupac") {
    const selectedText = info.selectionText;
    console.log('[Background] Rendering as IUPAC (OPSIN):', selectedText);

    // Send message to content script to handle IUPAC rendering
    chrome.tabs.sendMessage(tab.id, {
      type: 'RENDER_IUPAC',
      text: selectedText
    });
  }

  // ===== IMAGE RE-RENDER HANDLERS =====
  // These handle right-click on images to change rendering type
  if (info.menuItemId === "rerender-molecule") {
    console.log('[Background] Re-rendering image as Molecule (auto-detect)');
    chrome.tabs.sendMessage(tab.id, {
      type: 'RERENDER_IMAGE',
      srcUrl: info.srcUrl,
      renderAs: 'molecule'
    });
  }

  if (info.menuItemId === "rerender-smiles") {
    console.log('[Background] Re-rendering image as SMILES');
    chrome.tabs.sendMessage(tab.id, {
      type: 'RERENDER_IMAGE',
      srcUrl: info.srcUrl,
      renderAs: 'smiles'
    });
  }

  if (info.menuItemId === "rerender-biomolecule") {
    console.log('[Background] Re-rendering image as Biomolecule');
    chrome.tabs.sendMessage(tab.id, {
      type: 'RERENDER_IMAGE',
      srcUrl: info.srcUrl,
      renderAs: 'biomolecule'
    });
  }

  if (info.menuItemId === "rerender-mineral") {
    console.log('[Background] Re-rendering image as Mineral');
    chrome.tabs.sendMessage(tab.id, {
      type: 'RERENDER_IMAGE',
      srcUrl: info.srcUrl,
      renderAs: 'mineral'
    });
  }

  // Edit Flags handler - opens dialog to toggle flags
  if (info.menuItemId === "edit-flags") {
    console.log('[Background] Opening flag editor for image');
    chrome.tabs.sendMessage(tab.id, {
      type: 'EDIT_FLAGS',
      srcUrl: info.srcUrl
    });
  }
});

