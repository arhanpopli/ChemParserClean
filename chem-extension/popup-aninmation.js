/**
 * Advanced Molecule Animation for ChemTex
 * Uses SmilesDrawer to generate real molecule structures.
 * Implements an Atom-Level Physics Engine for realistic bonding/collision simulation.
 */

(function () {
    // --- Configuration ---
    var CONFIG = {
        SMILES_LIST: [
            'CN1C=NC2=C1C(=O)N(C(=O)N2C)C', // Caffeine
            'CC(=O)QC1=CC=CC=C1C(=O)O',     // Aspirin (Fixed Q->O typo just in case, or similar)
            'CC(=O)OC1=CC=CC=C1C(=O)O',     // Aspirin Correct
            'C1=CC=C(C=C1)C(=O)NO',         // Benzamide
            'CCO',                          // Ethanol
            'C(=O)(O)C1=CC=CC=C1',          // Benzoic acid
            'C1CCCCC1',                     // Cyclohexane
            'c1ccccc1',                     // Benzene
            'C1=CN=CN1',                    // Pyrimidine
            'C',                            // Methane (Test)
            'O'                             // Water (Test)
        ],
        ATOM_RADIUS: 10,    // Physics radius
        DRAW_SCALE: 0.7,    // Visual scale
        PHYSICS_SUBSTEPS: 3,
        RESTITUTION: 0.6,
        FRICTION: 0.99,
        ROTATION_FRICTION: 0.98
    };

    var container = document.getElementById('moleculeContainer');
    if (!container) return;

    var WIDTH = 360;
    var HEIGHT = 600;

    function updateDims() {
        WIDTH = window.innerWidth;
        HEIGHT = window.innerHeight;
    }
    window.addEventListener('resize', updateDims);
    updateDims();

    // --- Physics Engine Classes ---

    class MolecularBody {
        constructor(type, spriteUrl, atomPoints, width, height, x, y) {
            this.type = type;
            this.width = width;
            this.height = height;

            // Visual element
            this.element = document.createElement('div');
            this.element.className = 'floating-molecule';
            this.element.style.position = 'absolute';
            this.element.style.left = '0';
            this.element.style.top = '0';
            this.element.style.width = width + 'px';
            this.element.style.height = height + 'px';
            if (spriteUrl) {
                this.element.style.backgroundImage = 'url(' + spriteUrl + ')';
            } else {
                this.element.style.backgroundColor = 'rgba(100,100,100,0.5)';
                this.element.style.borderRadius = '50%';
            }
            this.element.style.backgroundSize = 'contain';
            this.element.style.backgroundRepeat = 'no-repeat';
            this.element.style.pointerEvents = 'none';
            this.element.style.opacity = '0.6';
            this.element.style.zIndex = '-1';
            container.appendChild(this.element);

            // Physics State
            this.pos = { x: x, y: y };
            this.vel = { x: (Math.random() - 0.5) * 2, y: (Math.random() - 0.5) * 2 };
            this.angle = Math.random() * Math.PI * 2;
            this.angVel = (Math.random() - 0.5) * 0.05;
            // Mass proportional to atom count
            this.mass = (atomPoints && atomPoints.length > 0) ? atomPoints.length : 1;

            this.localAtoms = atomPoints || [{ x: 0, y: 0 }];
            this.worldAtoms = [];
            this.updateWorldPoints();
        }

        updateWorldPoints() {
            var cos = Math.cos(this.angle);
            var sin = Math.sin(this.angle);
            this.worldAtoms = this.localAtoms.map(function (p) {
                return {
                    x: this.pos.x + (p.x * cos - p.y * sin),
                    y: this.pos.y + (p.x * sin + p.y * cos)
                };
            }.bind(this));
        }

        render() {
            var tx = this.pos.x - this.width / 2;
            var ty = this.pos.y - this.height / 2;
            this.element.style.transform =
                'translate(' + tx + 'px, ' + ty + 'px) rotate(' + this.angle + 'rad)';
        }

        integrate(dt) {
            this.pos.x += this.vel.x * dt;
            this.pos.y += this.vel.y * dt;
            this.angle += this.angVel * dt;

            // Wall collisions
            var margin = 20;
            if (this.pos.x < margin) { this.pos.x = margin; this.vel.x *= -0.8; }
            if (this.pos.x > WIDTH - margin) { this.pos.x = WIDTH - margin; this.vel.x *= -0.8; }
            if (this.pos.y < margin) { this.pos.y = margin; this.vel.y *= -0.8; }
            if (this.pos.y > HEIGHT - margin) { this.pos.y = HEIGHT - margin; this.vel.y *= -0.8; }

            this.vel.x *= CONFIG.FRICTION;
            this.vel.y *= CONFIG.FRICTION;
            this.angVel *= CONFIG.ROTATION_FRICTION;

            // Constrain velocity
            var speed = Math.sqrt(this.vel.x * this.vel.x + this.vel.y * this.vel.y);
            if (speed > 4) {
                this.vel.x *= (4 / speed);
                this.vel.y *= (4 / speed);
            }

            this.updateWorldPoints();
        }
    }

    // --- Logic to Generate Molecule Data from SmilesDrawer ---

    var moleculeCache = {};

    // Robust Mock Context
    class GeometryRecorder {
        constructor() {
            this.points = [];
            this.font = '10px sans-serif';
            // Default properties
            this.fillStyle = '#000';
            this.strokeStyle = '#000';
            this.lineWidth = 1;
        }
        moveTo(x, y) { this.points.push({ x: x, y: y }); }
        lineTo(x, y) { this.points.push({ x: x, y: y }); }
        fillText(txt, x, y) { this.points.push({ x: x, y: y }); }
        arc(x, y, r, sa, ea) { this.points.push({ x: x, y: y }); }
        rect(x, y, w, h) { this.points.push({ x: x, y: y }); this.points.push({ x: x + w, y: y + h }); }

        measureText(t) { return { width: 5 }; }
        beginPath() { } stroke() { } fill() { } save() { } restore() { }
        translate() { } rotate() { } scale() { } resetTransform() { } setTransform() { }
        closePath() { } arcTo() { } bezierCurveTo() { } quadraticCurveTo() { }
        clearRect() { }
        setLineDash() { }
    }

    function processMolecule(smiles, callback) {
        if (moleculeCache[smiles]) {
            callback(moleculeCache[smiles]);
            return;
        }

        try {
            var SD = window.SmilesDrawer || window.SmiDrawer;
            if (!SD) {
                console.warn("SmilesDrawer not available. Using fallback.");
                // Provide a fallback "atom"
                setTimeout(function () {
                    callback({ spriteUrl: null, width: 30, height: 30, atomPoints: [{ x: 0, y: 0 }] });
                }, 10);
                return;
            }

            var recorder = new GeometryRecorder();
            // Try to create drawer. 
            // Note: Options needs to be robust
            var drawerOpts = { width: 200, height: 200 };
            var drawer = new SD.Drawer(drawerOpts);

            // Mock Canvas
            var mockCanvas = {
                width: 200, height: 200,
                getContext: function (t) { return recorder; },
                setAttribute: function (k, v) { this[k] = v; },
                getAttribute: function (k) { return this[k]; },
                style: {}
            };

            SD.parse(smiles, function (tree) {
                try {
                    drawer.draw(tree, mockCanvas, 'light', false);
                } catch (e) { console.warn("Mock Draw Error", e); }

                var points = recorder.points;
                if (points.length === 0) points.push({ x: 100, y: 100 });

                // Bounds
                var minX = Infinity, maxX = -Infinity, minY = Infinity, maxY = -Infinity;
                points.forEach(function (p) {
                    if (p.x < minX) minX = p.x;
                    if (p.x > maxX) maxX = p.x;
                    if (p.y < minY) minY = p.y;
                    if (p.y > maxY) maxY = p.y;
                });

                if (Math.abs(maxX - minX) < 1) { maxX += 10; minX -= 10; }
                if (Math.abs(maxY - minY) < 1) { maxY += 10; minY -= 10; }

                var width = maxX - minX + 30;
                var height = maxY - minY + 30;
                var cx = (minX + maxX) / 2;
                var cy = (minY + maxY) / 2;

                // Consolidate atoms
                var atoms = [];
                points.forEach(function (p) {
                    var exists = false;
                    for (var k = 0; k < atoms.length; k++) {
                        var dx = atoms[k].x - p.x;
                        var dy = atoms[k].y - p.y;
                        if (dx * dx + dy * dy < 100) { exists = true; break; }
                    }
                    if (!exists) atoms.push(p);
                });

                var localAtoms = atoms.map(function (p) {
                    return {
                        x: (p.x - cx) * CONFIG.DRAW_SCALE,
                        y: (p.y - cy) * CONFIG.DRAW_SCALE
                    };
                });

                // Render Real Sprite
                var finalCanvas = document.createElement('canvas');
                finalCanvas.width = width * CONFIG.DRAW_SCALE;
                finalCanvas.height = height * CONFIG.DRAW_SCALE;
                var ctx = finalCanvas.getContext('2d');

                ctx.scale(CONFIG.DRAW_SCALE, CONFIG.DRAW_SCALE);
                ctx.translate(-minX + 15, -minY + 15);

                try {
                    var realDrawer = new SD.Drawer({ width: width, height: height });
                    realDrawer.draw(tree, finalCanvas, 'light', false);
                } catch (e) {
                    console.warn("Real Draw Error", e);
                    // Fallback visual
                    ctx.beginPath();
                    ctx.arc(cx, cy, 15, 0, Math.PI * 2);
                    ctx.fillStyle = "#888";
                    ctx.fill();
                }

                var data = {
                    spriteUrl: finalCanvas.toDataURL(),
                    width: finalCanvas.width,
                    height: finalCanvas.height,
                    atomPoints: localAtoms
                };
                moleculeCache[smiles] = data;
                callback(data);

            }, function (err) {
                console.warn("Parse Error", err);
                // Fallback
                callback({ spriteUrl: null, width: 30, height: 30, atomPoints: [{ x: 0, y: 0 }] });
            });

        } catch (e) {
            console.error("ProcessMolecule Error", e);
        }
    }

    // --- Main Loop ---

    var bodies = [];

    function resolveCollision(b1, b2) {
        var collisionNormal = { x: 0, y: 0 };
        var maxDepth = 0;
        var collided = false;

        var radSum = CONFIG.ATOM_RADIUS * 2;
        var radSumSq = radSum * radSum;

        // Check bounding box first
        var dx = b1.pos.x - b2.pos.x;
        var dy = b1.pos.y - b2.pos.y;
        var distEst = Math.sqrt(dx * dx + dy * dy);
        if (distEst > (b1.width + b2.width)) return; // Too far

        // Atom check
        for (var i = 0; i < b1.worldAtoms.length; i++) {
            var p1 = b1.worldAtoms[i];
            for (var j = 0; j < b2.worldAtoms.length; j++) {
                var p2 = b2.worldAtoms[j];

                var adx = p1.x - p2.x;
                var ady = p1.y - p2.y;
                var distSq = adx * adx + ady * ady;

                if (distSq < radSumSq) {
                    var dist = Math.sqrt(distSq);
                    if (dist < 0.1) dist = 0.1;

                    var depth = radSum - dist;
                    var nx = adx / dist;
                    var ny = ady / dist;

                    collisionNormal.x += nx;
                    collisionNormal.y += ny;
                    if (depth > maxDepth) maxDepth = depth;
                    collided = true;
                }
            }
        }

        if (collided) {
            var len = Math.sqrt(collisionNormal.x * collisionNormal.x + collisionNormal.y * collisionNormal.y);
            if (len > 0) {
                collisionNormal.x /= len;
                collisionNormal.y /= len;
            }

            // Separate
            var push = maxDepth * 0.4;
            b1.pos.x += collisionNormal.x * push;
            b1.pos.y += collisionNormal.y * push;
            b2.pos.x -= collisionNormal.x * push;
            b2.pos.y -= collisionNormal.y * push;

            // Bounce
            var dvx = b1.vel.x - b2.vel.x;
            var dvy = b1.vel.y - b2.vel.y;
            var velAlongNormal = dvx * collisionNormal.x + dvy * collisionNormal.y;

            if (velAlongNormal > 0) return;

            var j = -(1 + CONFIG.RESTITUTION) * velAlongNormal;
            j /= (1 / b1.mass + 1 / b2.mass);

            var ix = j * collisionNormal.x;
            var iy = j * collisionNormal.y;

            b1.vel.x += ix / b1.mass;
            b1.vel.y += iy / b1.mass;
            b2.vel.x -= ix / b2.mass;
            b2.vel.y -= iy / b2.mass;

            // Spin
            var spin = (Math.random() - 0.5) * 0.2;
            b1.angVel += spin;
            b2.angVel -= spin;
        }
    }

    function animate() {
        var dt = 1.0;

        // Physics SUBSTEPS
        for (var n = 0; n < CONFIG.PHYSICS_SUBSTEPS; n++) {
            var subDt = dt / CONFIG.PHYSICS_SUBSTEPS;
            for (var i = 0; i < bodies.length; i++) {
                bodies[i].integrate(subDt);
            }

            for (var i = 0; i < bodies.length; i++) {
                for (var j = i + 1; j < bodies.length; j++) {
                    resolveCollision(bodies[i], bodies[j]);
                }
            }
        }

        // Render
        for (var i = 0; i < bodies.length; i++) {
            bodies[i].render();
        }

        requestAnimationFrame(animate);
    }

    function init() {
        container.innerHTML = '';
        bodies = [];

        var loadedCount = 0;
        var targetCount = 12;

        for (var i = 0; i < targetCount; i++) {
            var s = CONFIG.SMILES_LIST[i % CONFIG.SMILES_LIST.length];
            // Delay spawn slightly to avoid lag spike
            (function (smiles, index) {
                setTimeout(function () {
                    processMolecule(smiles, function (data) {
                        var x = Math.random() * (WIDTH - 100) + 50;
                        var y = Math.random() * (HEIGHT - 100) + 50;
                        var b = new MolecularBody('mol', data.spriteUrl, data.atomPoints, data.width, data.height, x, y);
                        bodies.push(b);
                    });
                }, index * 100);
            })(s, i);
        }

        animate();
    }

    // Wait for SmilesDrawer
    var attempts = 0;
    var checkInterval = setInterval(function () {
        if ((window.SmilesDrawer || window.SmiDrawer) && document.body) {
            clearInterval(checkInterval);
            init();
        }
        attempts++;
        if (attempts > 50) { // 5 seconds timeout
            clearInterval(checkInterval);
            console.warn("SmilesDrawer timed out, forcing init with fallbacks");
            init(); // Will trigger fallback paths
        }
    }, 100);

})();
