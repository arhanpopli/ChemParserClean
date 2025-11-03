/**
 * Background Worker for Chemistry Renderer
 * Handles service worker operations
 */

console.log('[ChemRenderer] Background worker loaded');

// Listen for tab updates
chrome.tabs.onUpdated.addListener((tabId, changeInfo, tab) => {
  if (changeInfo.status === 'complete') {
    console.log('[ChemRenderer] Tab updated:', tab.url);
  }
});

// Handle header modifications (Manifest V3 doesn't support this directly in background.js)
// Instead, we'll use declarativeNetRequest or modify at content script level
