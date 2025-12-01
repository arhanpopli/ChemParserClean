/**
 * Background Service Worker for Chemistry Formula Renderer
 * Handles fetch requests from content scripts to bypass CSP restrictions
 */

console.log('[ChemRenderer] Background service worker loaded');

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
      url: response.url
    });

    if (!response.ok) {
      const errorBody = await response.text().catch(() => '(could not read error body)');
      console.error('[Background] HTTP Error response body:', errorBody);
      throw new Error(`HTTP ${response.status}: ${response.statusText}`);
    }

    const text = await response.text();
    console.log('[Background] Text response received, length:', text.length, 'preview:', text.substring(0, 100));
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
    console.log('[ChemRenderer] Tab updated:', tab.url);
  }
});

// Create context menu item
chrome.runtime.onInstalled.addListener(() => {
  chrome.contextMenus.create({
    id: "inspect-molecule",
    title: "Render as Molecule: '%s'",
    contexts: ["selection"]
  });
});

// Handle context menu clicks
chrome.contextMenus.onClicked.addListener((info, tab) => {
  if (info.menuItemId === "inspect-molecule") {
    const selectedText = info.selectionText;
    console.log('[Background] Inspecting molecule:', selectedText);

    // Send message to content script to handle the inspection
    chrome.tabs.sendMessage(tab.id, {
      type: 'INSPECT_MOLECULE',
      text: selectedText
    });
  }
});
