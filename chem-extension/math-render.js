/**
 * MathJax Renderer - Loaded from Extension Context
 * This script runs in the extension context, not the page context
 * So it can fetch MathJax from CDN without CSP restrictions
 */

console.log('[ChemRenderer] math-render.js loaded');

// Load MathJax configuration
window.MathJax = {
  tex: {
    inlineMath: [['$', '$']],
    displayMath: [['$$', '$$']],
    processEscapes: true
  },
  svg: {
    fontCache: 'global'
  },
  startup: {
    pageReady() {
      console.log('[ChemRenderer] MathJax page ready');
      return MathJax.typesetPromise();
    }
  }
};

// Load MathJax library
const script = document.createElement('script');
script.src = 'https://cdn.jsdelivr.net/npm/mathjax@3/es5/tex-svg.js';
script.async = true;

script.onload = function() {
  console.log('[ChemRenderer] ✅ MathJax loaded successfully from CDN');
  window.__mathJaxReady = true;
  // Signal to content script that MathJax is ready
  window.dispatchEvent(new CustomEvent('mathJaxReady'));
};

script.onerror = function(e) {
  console.error('[ChemRenderer] Failed to load MathJax:', e);
  // Try fallback CDN
  console.log('[ChemRenderer] Trying fallback CDN...');
  script.src = 'https://cdnjs.cloudflare.com/ajax/libs/mathjax/3.2.2/es5/tex-svg.min.js';
  document.head.appendChild(script);
};

document.head.appendChild(script);
