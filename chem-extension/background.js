/**
 * Background Worker for Chemistry Renderer
 * Handles service worker operations and proxies fetch requests to localhost servers
 */

console.log('[ChemRenderer] Background worker loaded');

// Listen for tab updates
chrome.tabs.onUpdated.addListener((tabId, changeInfo, tab) => {
  if (changeInfo.status === 'complete') {
    console.log('[ChemRenderer] Tab updated:', tab.url);
  }
});

// Handle messages from content script to proxy fetch requests
chrome.runtime.onMessage.addListener((request, sender, sendResponse) => {
  if (request.action === 'fetchUrl') {
    console.log('[ChemRenderer] Proxying fetch request:', request.url);

    // Make the fetch request from the background script (not subject to content script CSP)
    fetch(request.url, {
      method: request.method || 'GET',
      headers: request.headers || {},
      body: request.body || undefined
    })
      .then(response => {
        console.log('[ChemRenderer] Fetch response status:', response.status);

        // Determine content type
        const contentType = response.headers.get('content-type');

        if (contentType && contentType.includes('application/json')) {
          // Return JSON
          return response.json().then(data => {
            sendResponse({
              success: true,
              status: response.status,
              ok: response.ok,
              data: data,
              type: 'json'
            });
          });
        } else if (contentType && (contentType.includes('image/') || contentType.includes('application/octet-stream'))) {
          // Return blob as base64
          return response.blob().then(blob => {
            return new Promise((resolve) => {
              const reader = new FileReader();
              reader.onloadend = () => {
                sendResponse({
                  success: true,
                  status: response.status,
                  ok: response.ok,
                  data: reader.result, // base64 data URL
                  type: 'blob',
                  contentType: contentType
                });
                resolve();
              };
              reader.readAsDataURL(blob);
            });
          });
        } else {
          // Return text
          return response.text().then(text => {
            sendResponse({
              success: true,
              status: response.status,
              ok: response.ok,
              data: text,
              type: 'text'
            });
          });
        }
      })
      .catch(error => {
        console.error('[ChemRenderer] Fetch error:', error);
        sendResponse({
          success: false,
          error: error.message
        });
      });

    // Return true to indicate we'll respond asynchronously
    return true;
  }
});
