/**
 * SVG collider utilities for the popup molecules.
 *
 * Goals:
 * - Build compound bodies from SVG paths (tighter than circles/rects)
 * - Stay dependency-light (uses `window.SvgPathProperties` fallback backed by native SVGPathElement)
 *
 * Coordinate system:
 * - Input vertices in SVG viewBox space.
 * - Caller provides scaling to screen pixels.
 */

(function () {
  const EPS = 1e-6;

  function parseSvgString(svgString) {
    const parser = new DOMParser();
    const doc = parser.parseFromString(svgString, 'image/svg+xml');
    const svg = doc.querySelector('svg');
    return svg;
  }

  function parseViewBox(svgEl) {
    const vb = (svgEl.getAttribute('viewBox') || '').trim();
    if (!vb) {
      const w = parseFloat(svgEl.getAttribute('width') || '100');
      const h = parseFloat(svgEl.getAttribute('height') || '100');
      return { minX: 0, minY: 0, width: w, height: h };
    }
    const parts = vb.split(/[ ,]+/).map(Number);
    return { minX: parts[0], minY: parts[1], width: parts[2], height: parts[3] };
  }

  function samplePathD(d, targetPointCount) {
    // Prefer svg-path-properties API, but fallback provided by wrapper.
    const props = window.SvgPathProperties?.svgPathProperties
      ? window.SvgPathProperties.svgPathProperties(d)
      : null;

    // If path is invalid, return empty.
    if (!props) return [];

    const total = props.getTotalLength();
    if (!isFinite(total) || total <= EPS) return [];

    const count = Math.max(8, targetPointCount || 32);
    const pts = [];
    for (let i = 0; i < count; i++) {
      const t = (i / count) * total;
      const p = props.getPointAtLength(t);
      pts.push({ x: p.x, y: p.y });
    }

    return pts;
  }

  // Andrew monotonic chain convex hull. Returns hull in CCW order.
  function convexHull(points) {
    if (!points || points.length < 3) return points || [];

    // Dedup a bit
    const pts = points
      .filter(p => p && isFinite(p.x) && isFinite(p.y))
      .map(p => ({ x: p.x, y: p.y }))
      .sort((a, b) => (a.x - b.x) || (a.y - b.y));

    function cross(o, a, b) {
      return (a.x - o.x) * (b.y - o.y) - (a.y - o.y) * (b.x - o.x);
    }

    const lower = [];
    for (const p of pts) {
      while (lower.length >= 2 && cross(lower[lower.length - 2], lower[lower.length - 1], p) <= 0) {
        lower.pop();
      }
      lower.push(p);
    }

    const upper = [];
    for (let i = pts.length - 1; i >= 0; i--) {
      const p = pts[i];
      while (upper.length >= 2 && cross(upper[upper.length - 2], upper[upper.length - 1], p) <= 0) {
        upper.pop();
      }
      upper.push(p);
    }

    upper.pop();
    lower.pop();
    const hull = lower.concat(upper);

    return hull;
  }

  function centroid(points) {
    if (!points.length) return { x: 0, y: 0 };
    let cx = 0, cy = 0;
    for (const p of points) {
      cx += p.x;
      cy += p.y;
    }
    return { x: cx / points.length, y: cy / points.length };
  }

  function translate(points, dx, dy) {
    return points.map(p => ({ x: p.x + dx, y: p.y + dy }));
  }

  // Simple Ramer–Douglas–Peucker polyline simplification.
  function simplifyRDP(points, epsilon) {
    if (!points || points.length < 3) return points || [];
    const eps = epsilon || 0.75; // in SVG units

    function distToSegment(p, a, b) {
      const vx = b.x - a.x;
      const vy = b.y - a.y;
      const wx = p.x - a.x;
      const wy = p.y - a.y;

      const c1 = vx * wx + vy * wy;
      if (c1 <= 0) return Math.hypot(p.x - a.x, p.y - a.y);
      const c2 = vx * vx + vy * vy;
      if (c2 <= c1) return Math.hypot(p.x - b.x, p.y - b.y);
      const t = c1 / c2;
      const px = a.x + t * vx;
      const py = a.y + t * vy;
      return Math.hypot(p.x - px, p.y - py);
    }

    function rdp(start, end, out) {
      let maxDist = 0;
      let index = -1;
      for (let i = start + 1; i < end; i++) {
        const d = distToSegment(points[i], points[start], points[end]);
        if (d > maxDist) {
          maxDist = d;
          index = i;
        }
      }
      if (maxDist > eps && index !== -1) {
        rdp(start, index, out);
        out.pop();
        rdp(index, end, out);
      } else {
        out.push(points[start]);
        out.push(points[end]);
      }
    }

    const out = [];
    rdp(0, points.length - 1, out);

    // Dedup consecutive duplicates
    const dedup = [out[0]];
    for (let i = 1; i < out.length; i++) {
      const prev = dedup[dedup.length - 1];
      const cur = out[i];
      if (Math.hypot(cur.x - prev.x, cur.y - prev.y) > EPS) dedup.push(cur);
    }
    return dedup;
  }

  /**
   * Build compound polygons from SVG by grouping paths.
   *
   * Strategy:
   * - For each significant path, sample points -> convex hull -> simplify.
   * - Create a Matter body from each hull (part).
   * - Combine into compound body.
   *
   * Options:
   * - maxParts: cap number of path parts used
   * - samplesPerPath: int
   * - simplifyEps: float in SVG units
   */
  function buildCompoundFromSvgString(svgString, options) {
    const opts = options || {};

    const svg = parseSvgString(svgString);
    if (!svg) return null;

    const vb = parseViewBox(svg);

    // collect all paths
    const paths = Array.from(svg.querySelectorAll('path'));
    if (!paths.length) {
      return {
        viewBox: vb,
        parts: [
          // fallback: rectangle box
          [
            { x: vb.minX, y: vb.minY },
            { x: vb.minX + vb.width, y: vb.minY },
            { x: vb.minX + vb.width, y: vb.minY + vb.height },
            { x: vb.minX, y: vb.minY + vb.height }
          ]
        ]
      };
    }

    const samplesPerPath = opts.samplesPerPath || 36;
    const simplifyEps = opts.simplifyEps ?? 0.8;
    const maxParts = opts.maxParts || 12;

    // Rank paths by bbox area so we keep big/important paths first.
    const ranked = paths
      .map(p => {
        const d = p.getAttribute('d') || '';
        const pts = samplePathD(d, samplesPerPath);
        if (!pts.length) return null;
        let minX = Infinity, minY = Infinity, maxX = -Infinity, maxY = -Infinity;
        for (const q of pts) {
          minX = Math.min(minX, q.x);
          minY = Math.min(minY, q.y);
          maxX = Math.max(maxX, q.x);
          maxY = Math.max(maxY, q.y);
        }
        const area = (maxX - minX) * (maxY - minY);
        return { d, pts, area };
      })
      .filter(Boolean)
      .sort((a, b) => b.area - a.area)
      .slice(0, maxParts);

    if (!ranked.length) return null;

    // For each path, compute hull and simplify.
    const rawParts = [];
    for (const item of ranked) {
      const hull = convexHull(item.pts);
      if (hull.length < 3) continue;
      // close the loop for simplification stability
      const closed = hull.concat([hull[0]]);
      const simplified = simplifyRDP(closed, simplifyEps);
      // remove last if equals first
      const simp = simplified.length > 1 ? simplified.slice(0, -1) : simplified;
      if (simp.length >= 3) rawParts.push(simp);
    }

    if (!rawParts.length) return null;

    // Normalize each part around its own centroid so Matter.js can build them well.
    const parts = rawParts.map(poly => {
      const c = centroid(poly);
      return translate(poly, -c.x, -c.y);
    });

    // Provide per-part offsets so caller can place bodies correctly.
    const partCenters = rawParts.map(poly => centroid(poly));

    return {
      viewBox: vb,
      parts,
      partCenters
    };
  }

  window.MoleculeCollider = {
    buildCompoundFromSvgString
  };
})();
