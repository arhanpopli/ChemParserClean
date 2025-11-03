/**
 * CSP Bypass Helper
 * Executes scripts in page context to avoid CSP restrictions
 * Injects MathJax and other libraries directly into the page
 */

// Run MathJax loading in page context
function injectMathJaxWithBypass() {
  console.log('[ChemRenderer] [CSP-BYPASS] Attempting to inject MathJax...');
  
  const script = document.createElement('script');
  script.textContent = `
    (function() {
      console.log('[ChemRenderer] [PAGE-CONTEXT] Loading MathJax...');
      
      // Load MathJax config
      window.MathJax = {
        tex: {
          inlineMath: [['$', '$']],
          displayMath: [['$$', '$$']]
        },
        svg: {
          fontCache: 'global'
        }
      };
      
      // Create and inject script tag for MathJax
      var script = document.createElement('script');
      script.async = true;
      script.src = 'https://cdn.jsdelivr.net/npm/mathjax@3/es5/tex-svg.js';
      
      script.onload = function() {
        console.log('[ChemRenderer] [PAGE-CONTEXT] MathJax loaded successfully');
        window.__mathJaxReady = true;
        if (window.__chemRendererCallback) {
          window.__chemRendererCallback();
        }
      };
      
      script.onerror = function() {
        console.error('[ChemRenderer] [PAGE-CONTEXT] Failed to load MathJax from CDN');
        // Try alternate CDN
        script.src = 'https://cdnjs.cloudflare.com/ajax/libs/mathjax/3.2.2/es5/tex-svg.min.js';
        document.head.appendChild(script);
      };
      
      document.head.appendChild(script);
    })();
  `;
  
  try {
    (document.head || document.documentElement).appendChild(script);
    script.remove();
    console.log('[ChemRenderer] [CSP-BYPASS] MathJax injection script sent to page context');
    return true;
  } catch (e) {
    console.error('[ChemRenderer] [CSP-BYPASS] Failed to inject:', e);
    return false;
  }
}

// Alternative: Use fetch + blob URLs
function injectMathJaxViaBlob() {
  console.log('[ChemRenderer] [CSP-BYPASS] Attempting blob URL injection...');
  
  fetch('https://cdn.jsdelivr.net/npm/mathjax@3/es5/tex-svg.js')
    .then(response => response.text())
    .then(code => {
      const blob = new Blob([code], { type: 'application/javascript' });
      const blobUrl = URL.createObjectURL(blob);
      
      const script = document.createElement('script');
      script.src = blobUrl;
      script.async = true;
      
      script.onload = () => {
        console.log('[ChemRenderer] [CSP-BYPASS] MathJax loaded via blob URL');
        window.__mathJaxReady = true;
      };
      
      script.onerror = () => {
        console.error('[ChemRenderer] [CSP-BYPASS] Failed to load from blob URL');
      };
      
      document.head.appendChild(script);
    })
    .catch(err => {
      console.error('[ChemRenderer] [CSP-BYPASS] Blob fetch failed:', err);
    });
}

// Method 3: Chrome extension context injection
function injectFromExtensionContext() {
  console.log('[ChemRenderer] [CSP-BYPASS] Using extension context injection...');
  
  // Create iframe with minimal CSP
  const iframe = document.createElement('iframe');
  iframe.style.display = 'none';
  iframe.sandbox.add('allow-scripts');
  iframe.sandbox.add('allow-same-origin');
  iframe.src = 'about:blank';
  
  iframe.onload = function() {
    const iframeDoc = iframe.contentDocument || iframe.contentWindow.document;
    const script = iframeDoc.createElement('script');
    script.textContent = `
      // This runs in iframe with less restrictive CSP
      const mjScript = document.createElement('script');
      mjScript.src = 'https://cdn.jsdelivr.net/npm/mathjax@3/es5/tex-svg.js';
      mjScript.async = true;
      mjScript.onload = () => {
        console.log('[ChemRenderer] MathJax loaded in iframe');
        parent.window.__mathJaxReady = true;
      };
      document.head.appendChild(mjScript);
    `;
    iframeDoc.body.appendChild(script);
  };
  
  document.body.appendChild(iframe);
}

// Export functions for content script to use
window.__chemRendererBypass = {
  injectMathJaxWithBypass,
  injectMathJaxViaBlob,
  injectFromExtensionContext
};

console.log('[ChemRenderer] [CSP-BYPASS] Helper module loaded');
