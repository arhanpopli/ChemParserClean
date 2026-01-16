// Expose svg-path-properties to the popup page without a bundler.
// The vendored module is CommonJS; we load it via an IIFE-style require shim.
// This file runs in an extension page (self origin), so it's safe to attach to window.

(function () {
  // Minimal CommonJS loader for the single vendored file.
  function __requireSvgPathProperties() {
    // eslint-disable-next-line no-undef
    const module = { exports: {} };
    const exports = module.exports;

    // BEGIN vendored module
    /* eslint-disable */
    // The content below is injected at build-time by copy; keep as a single require so we can update easily.
    // In Chromium extension pages, `fetch` is allowed for self resources, but synchronous XHR is blocked in MV3.
    // So we embed the module source directly in this wrapper via a build step.
    /* eslint-enable */

    return module.exports;
  }

  // We don't inline automatically here; instead we expect window.SvgPathProperties to be set by the inlined build.
  // If it's not, fall back to a lightweight path sampler using native SVG APIs.

  if (window.SvgPathProperties) return;

  window.SvgPathProperties = {
    // API compatibility with svg-path-properties
    // Returns object with getTotalLength() and getPointAtLength().
    // This fallback uses SVGPathElement directly.
    svgPathProperties: function (d) {
      const svgNS = 'http://www.w3.org/2000/svg';
      const path = document.createElementNS(svgNS, 'path');
      path.setAttribute('d', d);
      return {
        getTotalLength: function () {
          return path.getTotalLength();
        },
        getPointAtLength: function (len) {
          const p = path.getPointAtLength(len);
          return { x: p.x, y: p.y };
        }
      };
    }
  };
})();
