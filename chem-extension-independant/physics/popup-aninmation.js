/**
 * Bouncing Molecule Animation for ChemistryLaTeX Popup
 * Uses mol2chemfig SVGs for molecule rendering
 * Features improved physics engine with AABB collision detection
 * Sub-stepping to prevent tunneling and realistic bounce physics
 */

(function () {
    var container = document.getElementById('moleculeContainer');
    if (!container) return;

    var WIDTH = 400;
    var HEIGHT = 650;

    // Physics constants
    var SUB_STEPS = 4; // Number of physics sub-steps per frame to prevent tunneling
    var RESTITUTION = 0.85; // Bounciness factor
    var WALL_RESTITUTION = 0.9;
    var MAX_SPEED = 3.0; // Maximum velocity to prevent tunneling
    var SEPARATION_FACTOR = 1.0; // How aggressively to separate overlapping molecules

    // Molecule Data - mol2chemfig SVGs only
    var moleculeData = {
        "ethane": "<?xml version='1.0' encoding='iso-8859-1'?><svg version='1.1' baseProfile='full'              xmlns='http://www.w3.org/2000/svg'                      xmlns:rdkit='http://www.rdkit.org/xml'                      xmlns:xlink='http://www.w3.org/1999/xlink'                  xml:space='preserve'width='105px' height='100px' viewBox='0 0 105 100'><!-- END OF HEADER --><path class='bond-0 atom-0 atom-1' d='M 18.5,50.0 L 86.0,50.0' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' /></svg>",
        "propane": "<?xml version='1.0' encoding='iso-8859-1'?><svg version='1.1' baseProfile='full'              xmlns='http://www.w3.org/2000/svg'                      xmlns:rdkit='http://www.rdkit.org/xml'                      xmlns:xlink='http://www.w3.org/1999/xlink'                  xml:space='preserve'width='137px' height='100px' viewBox='0 0 137 100'><!-- END OF HEADER --><path class='bond-0 atom-0 atom-1' d='M 6.9,67.8 L 68.3,32.2' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' /><path class='bond-1 atom-1 atom-2' d='M 68.3,32.2 L 129.8,67.8' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' /><path d='M 65.3,34.0 L 68.3,32.2 L 71.4,34.0' style='fill:none;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-miterlimit:10;stroke-opacity:1;' /></svg>",
        "methanol": "<?xml version='1.0' encoding='iso-8859-1'?><svg version='1.1' baseProfile='full'              xmlns='http://www.w3.org/2000/svg'                      xmlns:rdkit='http://www.rdkit.org/xml'                      xmlns:xlink='http://www.w3.org/1999/xlink'                  xml:space='preserve'width='105px' height='100px' viewBox='0 0 105 100'><!-- END OF HEADER --><path class='bond-0 atom-0 atom-1' d='M 5.2,49.9 L 61.6,49.9' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' /><path class='atom-1' d='M 63.8 49.9Q 63.8 45.3, 66.1 42.8Q 68.4 40.2, 72.6 40.2Q 76.8 40.2, 79.1 42.8Q 81.4 45.3, 81.4 49.9Q 81.4 54.5, 79.1 57.2Q 76.8 59.8, 72.6 59.8Q 68.4 59.8, 66.1 57.2Q 63.8 54.6, 63.8 49.9M 72.6 57.6Q 75.5 57.6, 77.1 55.7Q 78.7 53.7, 78.7 49.9Q 78.7 46.2, 77.1 44.3Q 75.5 42.4, 72.6 42.4Q 69.7 42.4, 68.1 44.2Q 66.5 46.1, 66.5 49.9Q 66.5 53.8, 68.1 55.7Q 69.7 57.6, 72.6 57.6' fill='#000000'/><path class='atom-1' d='M 84.3 40.4L 86.9 40.4L 86.9 48.5L 96.7 48.5L 96.7 40.4L 99.2 40.4L 99.2 59.5L 96.7 59.5L 96.7 50.7L 86.9 50.7L 86.9 59.5L 84.3 59.5L 84.3 40.4' fill='#000000'/></svg>",
        "dimethylEther": "<?xml version='1.0' encoding='iso-8859-1'?><svg version='1.1' baseProfile='full'              xmlns='http://www.w3.org/2000/svg'                      xmlns:rdkit='http://www.rdkit.org/xml'                      xmlns:xlink='http://www.w3.org/1999/xlink'                  xml:space='preserve'width='137px' height='100px' viewBox='0 0 137 100'><!-- END OF HEADER --><path class='bond-0 atom-0 atom-1' d='M 6.9,72.8 L 56.8,44.0' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' /><path class='bond-1 atom-1 atom-2' d='M 79.9,44.0 L 129.8,72.8' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' /><path class='atom-1' d='M 59.1 37.4Q 59.1 32.6, 61.5 29.9Q 63.9 27.2, 68.3 27.2Q 72.8 27.2, 75.2 29.9Q 77.6 32.6, 77.6 37.4Q 77.6 42.3, 75.2 45.1Q 72.8 47.8, 68.3 47.8Q 63.9 47.8, 61.5 45.1Q 59.1 42.3, 59.1 37.4M 68.3 45.5Q 71.4 45.5, 73.1 43.5Q 74.7 41.4, 74.7 37.4Q 74.7 33.4, 73.1 31.5Q 71.4 29.4, 68.3 29.4Q 65.3 29.4, 63.6 31.4Q 62.0 33.4, 62.0 37.4Q 62.0 41.4, 63.6 43.5Q 65.3 45.5, 68.3 45.5' fill='#000000'/></svg>",
        "acetone": "<?xml version='1.0' encoding='iso-8859-1'?><svg version='1.1' baseProfile='full'              xmlns='http://www.w3.org/2000/svg'                      xmlns:rdkit='http://www.rdkit.org/xml'                      xmlns:xlink='http://www.w3.org/1999/xlink'                  xml:space='preserve'width='137px' height='127px' viewBox='0 0 137 127'><!-- END OF HEADER --><path class='bond-0 atom-0 atom-1' d='M 8.3,120.3 L 68.3,85.7' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' /><path class='bond-1 atom-1 atom-2' d='M 73.6,88.7 L 73.6,28.8' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' /><path class='bond-1 atom-1 atom-2' d='M 63.1,88.7 L 63.1,28.8' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' /><path class='bond-2 atom-1 atom-3' d='M 68.3,85.7 L 128.4,120.3' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' /><path d='M 65.3,87.4 L 68.3,85.7 L 71.4,87.4' style='fill:none;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-miterlimit:10;stroke-opacity:1;' /><path class='atom-2' d='M 59.3 16.3Q 59.3 11.6, 61.7 9.0Q 64.0 6.3, 68.3 6.3Q 72.7 6.3, 75.0 9.0Q 77.4 11.6, 77.4 16.3Q 77.4 21.1, 75.0 23.8Q 72.7 26.5, 68.3 26.5Q 64.0 26.5, 61.7 23.8Q 59.3 21.1, 59.3 16.3M 68.3 24.3Q 71.3 24.3, 73.0 22.3Q 74.6 20.3, 74.6 16.3Q 74.6 12.5, 73.0 10.5Q 71.3 8.6, 68.3 8.6Q 65.4 8.6, 63.7 10.5Q 62.1 12.5, 62.1 16.3Q 62.1 20.3, 63.7 22.3Q 65.4 24.3, 68.3 24.3' fill='#000000'/></svg>",
        "acetaldehyde": "<?xml version='1.0' encoding='iso-8859-1'?><svg version='1.1' baseProfile='full'              xmlns='http://www.w3.org/2000/svg'                      xmlns:rdkit='http://www.rdkit.org/xml'                      xmlns:xlink='http://www.w3.org/1999/xlink'                  xml:space='preserve'width='137px' height='100px' viewBox='0 0 137 100'><!-- END OF HEADER --><path class='bond-0 atom-0 atom-1' d='M 6.9,59.7 L 65.2,26.0' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' /><path class='bond-1 atom-1 atom-2' d='M 65.2,26.0 L 112.6,53.3' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' /><path class='bond-1 atom-1 atom-2' d='M 65.2,37.7 L 107.6,62.1' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' /><path d='M 62.3,27.7 L 65.2,26.0 L 67.6,27.4' style='fill:none;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-miterlimit:10;stroke-opacity:1;' /><path class='atom-2' d='M 112.3 64.1Q 112.3 59.5, 114.6 57.0Q 116.9 54.4, 121.1 54.4Q 125.3 54.4, 127.6 57.0Q 129.9 59.5, 129.9 64.1Q 129.9 68.8, 127.6 71.4Q 125.3 74.0, 121.1 74.0Q 116.9 74.0, 114.6 71.4Q 112.3 68.8, 112.3 64.1M 121.1 71.9Q 124.0 71.9, 125.6 69.9Q 127.2 68.0, 127.2 64.1Q 127.2 60.4, 125.6 58.5Q 124.0 56.6, 121.1 56.6Q 118.2 56.6, 116.6 58.5Q 115.0 60.3, 115.0 64.1Q 115.0 68.0, 116.6 69.9Q 118.2 71.9, 121.1 71.9' fill='#000000'/></svg>",
        "isobutane": "<?xml version='1.0' encoding='iso-8859-1'?><svg version='1.1' baseProfile='full'              xmlns='http://www.w3.org/2000/svg'                      xmlns:rdkit='http://www.rdkit.org/xml'                      xmlns:xlink='http://www.w3.org/1999/xlink'                  xml:space='preserve'width='137px' height='127px' viewBox='0 0 137 127'><!-- END OF HEADER --><path class='bond-0 atom-0 atom-1' d='M 6.9,116.6 L 68.3,81.1' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' /><path class='bond-1 atom-1 atom-2' d='M 68.3,81.1 L 129.8,116.6' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' /><path class='bond-2 atom-1 atom-3' d='M 68.3,81.1 L 68.3,10.1' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' /></svg>",
        "benzene": "<?xml version='1.0' encoding='iso-8859-1'?><svg version='1.1' baseProfile='full'              xmlns='http://www.w3.org/2000/svg'                      xmlns:rdkit='http://www.rdkit.org/xml'                      xmlns:xlink='http://www.w3.org/1999/xlink'                  xml:space='preserve'width='150px' height='137px' viewBox='0 0 150 137'><!-- END OF HEADER --><path class='bond-0 atom-0 atom-1' d='M 142.5,68.3 L 108.7,126.8' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' /><path class='bond-0 atom-0 atom-1' d='M 130.8,68.3 L 102.9,116.7' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' /><path class='bond-1 atom-1 atom-2' d='M 108.7,126.8 L 41.2,126.8' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' /><path class='bond-2 atom-2 atom-3' d='M 41.2,126.8 L 7.5,68.3' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' /><path class='bond-2 atom-2 atom-3' d='M 47.1,116.7 L 19.2,68.3' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' /><path class='bond-3 atom-3 atom-4' d='M 7.5,68.3 L 41.3,9.9' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' /><path class='bond-4 atom-4 atom-5' d='M 41.3,9.9 L 108.8,9.9' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' /><path class='bond-4 atom-4 atom-5' d='M 47.1,20.0 L 102.9,20.0' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' /><path class='bond-5 atom-5 atom-0' d='M 108.8,9.9 L 142.5,68.3' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' /><path d='M 140.8,71.3 L 142.5,68.3 L 140.8,65.4' style='fill:none;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-miterlimit:10;stroke-opacity:1;' /><path d='M 110.4,123.9 L 108.7,126.8 L 105.4,126.8' style='fill:none;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-miterlimit:10;stroke-opacity:1;' /><path d='M 44.6,126.8 L 41.2,126.8 L 39.6,123.9' style='fill:none;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-miterlimit:10;stroke-opacity:1;' /><path d='M 9.2,71.3 L 7.5,68.3 L 9.2,65.4' style='fill:none;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-miterlimit:10;stroke-opacity:1;' /><path d='M 39.6,12.8 L 41.3,9.9 L 44.6,9.9' style='fill:none;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-miterlimit:10;stroke-opacity:1;' /><path d='M 105.4,9.9 L 108.8,9.9 L 110.4,12.8' style='fill:none;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-miterlimit:10;stroke-opacity:1;' /></svg>",
        "cyclohexane": "<?xml version='1.0' encoding='iso-8859-1'?><svg version='1.1' baseProfile='full'              xmlns='http://www.w3.org/2000/svg'                      xmlns:rdkit='http://www.rdkit.org/xml'                      xmlns:xlink='http://www.w3.org/1999/xlink'                  xml:space='preserve'width='150px' height='137px' viewBox='0 0 150 137'><!-- END OF HEADER --><path class='bond-0 atom-0 atom-1' d='M 142.5,68.3 L 108.7,126.8' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' /><path class='bond-1 atom-1 atom-2' d='M 108.7,126.8 L 41.2,126.8' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' /><path class='bond-2 atom-2 atom-3' d='M 41.2,126.8 L 7.5,68.3' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' /><path class='bond-3 atom-3 atom-4' d='M 7.5,68.3 L 41.3,9.9' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' /><path class='bond-4 atom-4 atom-5' d='M 41.3,9.9 L 108.8,9.9' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' /><path class='bond-5 atom-5 atom-0' d='M 108.8,9.9 L 142.5,68.3' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' /><path d='M 140.8,71.3 L 142.5,68.3 L 140.8,65.4' style='fill:none;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-miterlimit:10;stroke-opacity:1;' /><path d='M 110.4,123.9 L 108.7,126.8 L 105.4,126.8' style='fill:none;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-miterlimit:10;stroke-opacity:1;' /><path d='M 44.6,126.8 L 41.2,126.8 L 39.6,123.9' style='fill:none;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-miterlimit:10;stroke-opacity:1;' /><path d='M 9.2,71.3 L 7.5,68.3 L 9.2,65.4' style='fill:none;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-miterlimit:10;stroke-opacity:1;' /><path d='M 39.6,12.8 L 41.3,9.9 L 44.6,9.9' style='fill:none;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-miterlimit:10;stroke-opacity:1;' /><path d='M 105.4,9.9 L 108.8,9.9 L 110.4,12.8' style='fill:none;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-miterlimit:10;stroke-opacity:1;' /></svg>",
        "pyridine": "<?xml version='1.0' encoding='iso-8859-1'?><svg version='1.1' baseProfile='full'              xmlns='http://www.w3.org/2000/svg'                      xmlns:rdkit='http://www.rdkit.org/xml'                      xmlns:xlink='http://www.w3.org/1999/xlink'                  xml:space='preserve'width='150px' height='137px' viewBox='0 0 150 137'><!-- END OF HEADER --><path class='bond-0 atom-0 atom-1' d='M 142.5,68.3 L 110.5,123.8' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' /><path class='bond-0 atom-0 atom-1' d='M 131.4,68.3 L 104.9,114.2' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' /><path class='bond-1 atom-1 atom-2' d='M 110.5,123.8 L 46.4,123.8' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' /><path class='bond-2 atom-2 atom-3' d='M 46.4,123.8 L 20.8,79.5' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' /><path class='bond-2 atom-2 atom-3' d='M 52.0,114.2 L 29.2,74.7' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' /><path class='bond-3 atom-3 atom-4' d='M 20.8,57.2 L 46.4,12.9' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' /><path class='bond-4 atom-4 atom-5' d='M 46.4,12.9 L 110.5,12.9' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' /><path class='bond-4 atom-4 atom-5' d='M 52.0,22.5 L 104.9,22.5' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' /><path class='bond-5 atom-5 atom-0' d='M 110.5,12.9 L 142.5,68.3' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' /><path d='M 140.9,71.1 L 142.5,68.3 L 140.9,65.6' style='fill:none;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-miterlimit:10;stroke-opacity:1;' /><path d='M 112.1,121.1 L 110.5,123.8 L 107.3,123.8' style='fill:none;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-miterlimit:10;stroke-opacity:1;' /><path d='M 49.6,123.8 L 46.4,123.8 L 45.1,121.6' style='fill:none;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-miterlimit:10;stroke-opacity:1;' /><path d='M 45.1,15.1 L 46.4,12.9 L 49.6,12.9' style='fill:none;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-miterlimit:10;stroke-opacity:1;' /><path d='M 107.3,12.9 L 110.5,12.9 L 112.1,15.6' style='fill:none;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-miterlimit:10;stroke-opacity:1;' /><path class='atom-3' d='M 10.4 59.3L 16.3 68.9Q 16.9 69.8, 17.9 71.6Q 18.8 73.3, 18.9 73.4L 18.9 59.3L 21.3 59.3L 21.3 77.4L 18.8 77.4L 12.4 66.9Q 11.7 65.7, 10.9 64.3Q 10.1 62.9, 9.9 62.4L 9.9 77.4L 7.5 77.4L 7.5 59.3L 10.4 59.3' fill='#000000'/></svg>",
        "THF": "<?xml version='1.0' encoding='iso-8859-1'?><svg version='1.1' baseProfile='full'              xmlns='http://www.w3.org/2000/svg'                      xmlns:rdkit='http://www.rdkit.org/xml'                      xmlns:xlink='http://www.w3.org/1999/xlink'                  xml:space='preserve'width='129px' height='132px' viewBox='0 0 129 132'><!-- END OF HEADER --><path class='bond-0 atom-0 atom-1' d='M 122.5,65.6 L 81.6,9.4' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' /><path class='bond-1 atom-1 atom-2' d='M 81.6,9.4 L 15.5,30.8' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' /><path class='bond-2 atom-2 atom-3' d='M 15.5,30.8 L 15.5,88.1' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' /><path class='bond-3 atom-3 atom-4' d='M 26.8,104.0 L 81.6,121.8' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' /><path class='bond-4 atom-4 atom-0' d='M 81.6,121.8 L 122.5,65.6' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' /><path d='M 120.4,62.8 L 122.5,65.6 L 120.4,68.4' style='fill:none;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-miterlimit:10;stroke-opacity:1;' /><path d='M 83.6,12.2 L 81.6,9.4 L 78.3,10.4' style='fill:none;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-miterlimit:10;stroke-opacity:1;' /><path d='M 18.8,29.8 L 15.5,30.8 L 15.5,33.7' style='fill:none;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-miterlimit:10;stroke-opacity:1;' /><path d='M 78.9,120.9 L 81.6,121.8 L 83.6,119.0' style='fill:none;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-miterlimit:10;stroke-opacity:1;' /><path class='atom-3' d='M 6.5 100.4Q 6.5 95.7, 8.8 93.0Q 11.1 90.4, 15.5 90.4Q 19.9 90.4, 22.2 93.0Q 24.5 95.7, 24.5 100.4Q 24.5 105.2, 22.2 107.9Q 19.8 110.6, 15.5 110.6Q 11.1 110.6, 8.8 107.9Q 6.5 105.2, 6.5 100.4M 15.5 108.4Q 18.5 108.4, 20.1 106.4Q 21.7 104.4, 21.7 100.4Q 21.7 96.5, 20.1 94.6Q 18.5 92.6, 15.5 92.6Q 12.5 92.6, 10.8 94.6Q 9.2 96.5, 9.2 100.4Q 9.2 104.4, 10.8 106.4Q 12.5 108.4, 15.5 108.4' fill='#000000'/></svg>",
        "phenol": "<?xml version='1.0' encoding='iso-8859-1'?><svg version='1.1' baseProfile='full'              xmlns='http://www.w3.org/2000/svg'                      xmlns:rdkit='http://www.rdkit.org/xml'                      xmlns:xlink='http://www.w3.org/1999/xlink'                  xml:space='preserve'width='195px' height='137px' viewBox='0 0 195 137'><!-- END OF HEADER --><path class='bond-0 atom-0 atom-1' d='M 156.0,68.3 L 112.8,68.3' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' /><path class='bond-1 atom-1 atom-2' d='M 112.8,68.3 L 87.1,113.0' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' /><path class='bond-1 atom-1 atom-2' d='M 103.9,68.3 L 82.6,105.3' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' /><path class='bond-2 atom-2 atom-3' d='M 87.1,113.0 L 35.5,113.0' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' /><path class='bond-3 atom-3 atom-4' d='M 35.5,113.0 L 9.8,68.3' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' /><path class='bond-3 atom-3 atom-4' d='M 40.0,105.3 L 18.7,68.3' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' /><path class='bond-4 atom-4 atom-5' d='M 9.8,68.3 L 35.5,23.7' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' /><path class='bond-5 atom-5 atom-6' d='M 35.5,23.7 L 87.1,23.7' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' /><path class='bond-5 atom-5 atom-6' d='M 40.0,31.4 L 82.6,31.4' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' /><path class='bond-6 atom-6 atom-1' d='M 87.1,23.7 L 112.8,68.3' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' /><path d='M 88.3,110.7 L 87.1,113.0 L 84.5,113.0' style='fill:none;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-miterlimit:10;stroke-opacity:1;' /><path d='M 38.1,113.0 L 35.5,113.0 L 34.2,110.7' style='fill:none;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-miterlimit:10;stroke-opacity:1;' /><path d='M 11.0,70.6 L 9.8,68.3 L 11.0,66.1' style='fill:none;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-miterlimit:10;stroke-opacity:1;' /><path d='M 34.2,26.0 L 35.5,23.7 L 38.1,23.7' style='fill:none;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-miterlimit:10;stroke-opacity:1;' /><path d='M 84.5,23.7 L 87.1,23.7 L 88.3,26.0' style='fill:none;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-miterlimit:10;stroke-opacity:1;' /><path class='atom-0' d='M 157.7 68.4Q 157.7 64.9, 159.4 62.9Q 161.1 61.0, 164.4 61.0Q 167.6 61.0, 169.3 62.9Q 171.1 64.9, 171.1 68.4Q 171.1 71.9, 169.3 74.0Q 167.6 76.0, 164.4 76.0Q 161.1 76.0, 159.4 74.0Q 157.7 72.0, 157.7 68.4M 164.4 74.3Q 166.6 74.3, 167.8 72.8Q 169.0 71.3, 169.0 68.4Q 169.0 65.5, 167.8 64.1Q 166.6 62.6, 164.4 62.6Q 162.1 62.6, 160.9 64.1Q 159.7 65.5, 159.7 68.4Q 159.7 71.3, 160.9 72.8Q 162.1 74.3, 164.4 74.3' fill='#000000'/><path class='atom-0' d='M 173.3 61.1L 175.3 61.1L 175.3 67.3L 182.8 67.3L 182.8 61.1L 184.8 61.1L 184.8 75.7L 182.8 75.7L 182.8 69.0L 175.3 69.0L 175.3 75.7L 173.3 75.7L 173.3 61.1' fill='#000000'/></svg>",
        "toluene": "<?xml version='1.0' encoding='iso-8859-1'?><svg version='1.1' baseProfile='full'              xmlns='http://www.w3.org/2000/svg'                      xmlns:rdkit='http://www.rdkit.org/xml'                      xmlns:xlink='http://www.w3.org/1999/xlink'                  xml:space='preserve'width='195px' height='137px' viewBox='0 0 195 137'><!-- END OF HEADER --><path class='bond-0 atom-0 atom-1' d='M 184.8,68.3 L 126.4,68.3' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' /><path class='bond-1 atom-1 atom-2' d='M 126.4,68.3 L 97.2,118.9' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' /><path class='bond-1 atom-1 atom-2' d='M 116.3,68.3 L 92.2,110.1' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' /><path class='bond-2 atom-2 atom-3' d='M 97.2,118.9 L 38.9,118.9' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' /><path class='bond-3 atom-3 atom-4' d='M 38.9,118.9 L 9.8,68.3' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' /><path class='bond-3 atom-3 atom-4' d='M 44.0,110.1 L 19.9,68.3' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' /><path class='bond-4 atom-4 atom-5' d='M 9.8,68.3 L 38.9,17.8' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' /><path class='bond-5 atom-5 atom-6' d='M 38.9,17.8 L 97.3,17.8' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' /><path class='bond-5 atom-5 atom-6' d='M 44.0,26.6 L 92.2,26.6' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' /><path class='bond-6 atom-6 atom-1' d='M 97.3,17.8 L 126.4,68.3' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' /><path d='M 98.7,116.3 L 97.2,118.9 L 94.3,118.9' style='fill:none;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-miterlimit:10;stroke-opacity:1;' /><path d='M 41.8,118.9 L 38.9,118.9 L 37.5,116.3' style='fill:none;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-miterlimit:10;stroke-opacity:1;' /><path d='M 11.2,70.9 L 9.8,68.3 L 11.2,65.8' style='fill:none;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-miterlimit:10;stroke-opacity:1;' /><path d='M 37.5,20.4 L 38.9,17.8 L 41.8,17.8' style='fill:none;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-miterlimit:10;stroke-opacity:1;' /><path d='M 94.3,17.8 L 97.3,17.8 L 98.7,20.4' style='fill:none;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-miterlimit:10;stroke-opacity:1;' /></svg>",
        "catechol": "<?xml version='1.0' encoding='iso-8859-1'?><svg version='1.1' baseProfile='full'              xmlns='http://www.w3.org/2000/svg'                      xmlns:rdkit='http://www.rdkit.org/xml'                      xmlns:xlink='http://www.w3.org/1999/xlink'                  xml:space='preserve'width='176px' height='150px' viewBox='0 0 176 150'><!-- END OF HEADER --><path class='bond-0 atom-0 atom-1' d='M 38.9,122.5 L 75.8,101.1' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' /><path class='bond-1 atom-1 atom-2' d='M 75.8,101.1 L 75.8,48.6' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' /><path class='bond-1 atom-1 atom-2' d='M 83.7,96.6 L 83.7,53.2' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' /><path class='bond-2 atom-2 atom-3' d='M 75.8,48.6 L 38.9,27.3' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' /><path class='bond-3 atom-2 atom-4' d='M 75.8,48.6 L 121.3,22.4' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' /><path class='bond-4 atom-4 atom-5' d='M 121.3,22.4 L 166.8,48.6' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' /><path class='bond-4 atom-4 atom-5' d='M 121.3,31.5 L 158.9,53.2' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' /><path class='bond-5 atom-5 atom-6' d='M 166.8,48.6 L 166.8,101.1' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' /><path class='bond-6 atom-6 atom-7' d='M 166.8,101.1 L 121.3,127.4' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' /><path class='bond-6 atom-6 atom-7' d='M 158.9,96.6 L 121.3,118.3' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' /><path class='bond-7 atom-7 atom-1' d='M 121.3,127.4 L 75.8,101.1' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' /><path d='M 119.0,23.7 L 121.3,22.4 L 123.6,23.7' style='fill:none;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-miterlimit:10;stroke-opacity:1;' /><path d='M 164.5,47.3 L 166.8,48.6 L 166.8,51.2' style='fill:none;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-miterlimit:10;stroke-opacity:1;' /><path d='M 166.8,98.5 L 166.8,101.1 L 164.5,102.5' style='fill:none;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-miterlimit:10;stroke-opacity:1;' /><path d='M 123.6,126.1 L 121.3,127.4 L 119.0,126.1' style='fill:none;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-miterlimit:10;stroke-opacity:1;' /><path class='atom-0' d='M 8.8 120.1L 10.8 120.1L 10.8 126.4L 18.4 126.4L 18.4 120.1L 20.4 120.1L 20.4 134.9L 18.4 134.9L 18.4 128.1L 10.8 128.1L 10.8 134.9L 8.8 134.9L 8.8 120.1' fill='#000000'/><path class='atom-0' d='M 23.5 127.5Q 23.5 123.9, 25.3 121.9Q 27.0 119.9, 30.3 119.9Q 33.6 119.9, 35.4 121.9Q 37.1 123.9, 37.1 127.5Q 37.1 131.1, 35.4 133.1Q 33.6 135.2, 30.3 135.2Q 27.0 135.2, 25.3 133.1Q 23.5 131.1, 23.5 127.5M 30.3 133.5Q 32.6 133.5, 33.8 132.0Q 35.0 130.4, 35.0 127.5Q 35.0 124.5, 33.8 123.1Q 32.6 121.6, 30.3 121.6Q 28.1 121.6, 26.8 123.0Q 25.6 124.5, 25.6 127.5Q 25.6 130.5, 26.8 132.0Q 28.1 133.5, 30.3 133.5' fill='#000000'/><path class='atom-3' d='M 8.8 15.0L 10.8 15.0L 10.8 21.3L 18.4 21.3L 18.4 15.0L 20.4 15.0L 20.4 29.9L 18.4 29.9L 18.4 23.0L 10.8 23.0L 10.8 29.9L 8.8 29.9L 8.8 15.0' fill='#000000'/><path class='atom-3' d='M 23.5 22.4Q 23.5 18.8, 25.3 16.8Q 27.0 14.8, 30.3 14.8Q 33.6 14.8, 35.4 16.8Q 37.1 18.8, 37.1 22.4Q 37.1 26.0, 35.4 28.1Q 33.6 30.1, 30.3 30.1Q 27.0 30.1, 25.3 28.1Q 23.5 26.0, 23.5 22.4M 30.3 28.4Q 32.6 28.4, 33.8 26.9Q 35.0 25.4, 35.0 22.4Q 35.0 19.5, 33.8 18.0Q 32.6 16.5, 30.3 16.5Q 28.1 16.5, 26.8 18.0Q 25.6 19.5, 25.6 22.4Q 25.6 25.4, 26.8 26.9Q 28.1 28.4, 30.3 28.4' fill='#000000'/></svg>",
        "naphthalene": "<?xml version='1.0' encoding='iso-8859-1'?><svg version='1.1' baseProfile='full'              xmlns='http://www.w3.org/2000/svg'                      xmlns:rdkit='http://www.rdkit.org/xml'                      xmlns:xlink='http://www.w3.org/1999/xlink'                  xml:space='preserve'width='215px' height='150px' viewBox='0 0 215 150'><!-- END OF HEADER --><path class='bond-0 atom-0 atom-1' d='M 10.8,102.9 L 10.8,47.1' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' /><path class='bond-0 atom-0 atom-1' d='M 19.1,98.0 L 19.1,52.0' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' /><path class='bond-1 atom-1 atom-2' d='M 10.8,47.1 L 59.0,19.3' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' /><path class='bond-2 atom-2 atom-3' d='M 59.0,19.3 L 107.2,47.1' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' /><path class='bond-2 atom-2 atom-3' d='M 59.0,28.9 L 98.9,52.0' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' /><path class='bond-3 atom-3 atom-4' d='M 107.2,47.1 L 155.5,19.3' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' /><path class='bond-4 atom-4 atom-5' d='M 155.5,19.3 L 203.8,47.1' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' /><path class='bond-4 atom-4 atom-5' d='M 155.5,28.9 L 195.4,52.0' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' /><path class='bond-5 atom-5 atom-6' d='M 203.8,47.1 L 203.8,102.9' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' /><path class='bond-6 atom-6 atom-7' d='M 203.8,102.9 L 155.5,130.7' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' /><path class='bond-6 atom-6 atom-7' d='M 195.4,98.0 L 155.5,121.1' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' /><path class='bond-7 atom-7 atom-8' d='M 155.5,130.7 L 107.2,102.9' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' /><path class='bond-8 atom-8 atom-9' d='M 107.2,102.9 L 59.0,130.7' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' /><path class='bond-8 atom-8 atom-9' d='M 98.9,98.0 L 59.0,121.1' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' /><path class='bond-9 atom-9 atom-0' d='M 59.0,130.7 L 10.8,102.9' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' /><path class='bond-10 atom-8 atom-3' d='M 107.2,102.9 L 107.2,47.1' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' /><path d='M 10.8,100.1 L 10.8,102.9 L 13.2,104.3' style='fill:none;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-miterlimit:10;stroke-opacity:1;' /><path d='M 10.8,49.9 L 10.8,47.1 L 13.2,45.7' style='fill:none;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-miterlimit:10;stroke-opacity:1;' /><path d='M 56.6,20.7 L 59.0,19.3 L 61.4,20.7' style='fill:none;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-miterlimit:10;stroke-opacity:1;' /><path d='M 153.1,20.7 L 155.5,19.3 L 157.9,20.7' style='fill:none;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-miterlimit:10;stroke-opacity:1;' /><path d='M 201.3,45.7 L 203.8,47.1 L 203.8,49.9' style='fill:none;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-miterlimit:10;stroke-opacity:1;' /><path d='M 203.8,100.1 L 203.8,102.9 L 201.3,104.3' style='fill:none;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-miterlimit:10;stroke-opacity:1;' /><path d='M 157.9,129.3 L 155.5,130.7 L 153.1,129.3' style='fill:none;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-miterlimit:10;stroke-opacity:1;' /><path d='M 61.4,129.3 L 59.0,130.7 L 56.6,129.3' style='fill:none;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-miterlimit:10;stroke-opacity:1;' /></svg>",
        "anisole": "<?xml version='1.0' encoding='iso-8859-1'?><svg version='1.1' baseProfile='full'              xmlns='http://www.w3.org/2000/svg'                      xmlns:rdkit='http://www.rdkit.org/xml'                      xmlns:xlink='http://www.w3.org/1999/xlink'                  xml:space='preserve'width='222px' height='145px' viewBox='0 0 222 145'><!-- END OF HEADER --><path class='bond-0 atom-0 atom-1' d='M 11.1,56.0 L 39.6,86.8' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' /><path class='bond-1 atom-1 atom-2' d='M 57.5,94.5 L 102.4,84.4' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' /><path class='bond-2 atom-2 atom-3' d='M 102.4,84.4 L 118.8,31.7' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' /><path class='bond-2 atom-2 atom-3' d='M 111.7,82.3 L 125.3,38.7' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' /><path class='bond-3 atom-3 atom-4' d='M 118.8,31.7 L 172.7,19.5' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' /><path class='bond-4 atom-4 atom-5' d='M 172.7,19.5 L 210.1,60.1' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' /><path class='bond-4 atom-4 atom-5' d='M 169.8,28.7 L 200.8,62.2' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' /><path class='bond-5 atom-5 atom-6' d='M 210.1,60.1 L 193.7,112.8' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' /><path class='bond-6 atom-6 atom-7' d='M 193.7,112.8 L 139.8,125.0' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' /><path class='bond-6 atom-6 atom-7' d='M 187.2,105.8 L 142.7,115.8' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' /><path class='bond-7 atom-7 atom-2' d='M 139.8,125.0 L 102.4,84.4' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' /><path d='M 118.0,34.3 L 118.8,31.7 L 121.5,31.1' style='fill:none;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-miterlimit:10;stroke-opacity:1;' /><path d='M 170.0,20.1 L 172.7,19.5 L 174.5,21.6' style='fill:none;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-miterlimit:10;stroke-opacity:1;' /><path d='M 208.2,58.1 L 210.1,60.1 L 209.3,62.7' style='fill:none;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-miterlimit:10;stroke-opacity:1;' /><path d='M 194.5,110.2 L 193.7,112.8 L 191.0,113.4' style='fill:none;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-miterlimit:10;stroke-opacity:1;' /><path d='M 142.5,124.4 L 139.8,125.0 L 138.0,122.9' style='fill:none;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-miterlimit:10;stroke-opacity:1;' /><path class='atom-1' d='M 41.4 96.6Q 41.4 92.8, 43.2 90.7Q 45.1 88.6, 48.5 88.6Q 52.0 88.6, 53.9 90.7Q 55.7 92.8, 55.7 96.6Q 55.7 100.4, 53.8 102.5Q 52.0 104.7, 48.5 104.7Q 45.1 104.7, 43.2 102.5Q 41.4 100.4, 41.4 96.6M 48.5 102.9Q 50.9 102.9, 52.2 101.3Q 53.5 99.7, 53.5 96.6Q 53.5 93.5, 52.2 92.0Q 50.9 90.4, 48.5 90.4Q 46.2 90.4, 44.9 91.9Q 43.6 93.5, 43.6 96.6Q 43.6 99.7, 44.9 101.3Q 46.2 102.9, 48.5 102.9' fill='#000000'/></svg>",
        "dioxane": "<?xml version='1.0' encoding='iso-8859-1'?><svg version='1.1' baseProfile='full'              xmlns='http://www.w3.org/2000/svg'                      xmlns:rdkit='http://www.rdkit.org/xml'                      xmlns:xlink='http://www.w3.org/1999/xlink'                  xml:space='preserve'width='150px' height='137px' viewBox='0 0 150 137'><!-- END OF HEADER --><path class='bond-0 atom-0 atom-1' d='M 135.8,68.2 L 105.4,120.9' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' /><path class='bond-1 atom-1 atom-2' d='M 105.4,120.9 L 54.5,120.9' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' /><path class='bond-2 atom-2 atom-3' d='M 38.4,110.2 L 14.2,68.2' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' /><path class='bond-3 atom-3 atom-4' d='M 14.2,68.2 L 44.6,15.6' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' /><path class='bond-4 atom-4 atom-5' d='M 44.6,15.6 L 95.5,15.6' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' /><path class='bond-5 atom-5 atom-0' d='M 111.7,26.5 L 135.8,68.2' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' /><path d='M 134.3,70.8 L 135.8,68.2 L 134.6,66.1' style='fill:none;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-miterlimit:10;stroke-opacity:1;' /><path d='M 106.9,118.2 L 105.4,120.9 L 102.9,120.9' style='fill:none;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-miterlimit:10;stroke-opacity:1;' /><path d='M 15.4,70.3 L 14.2,68.2 L 15.7,65.6' style='fill:none;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-miterlimit:10;stroke-opacity:1;' /><path d='M 43.1,18.2 L 44.6,15.6 L 47.1,15.6' style='fill:none;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-miterlimit:10;stroke-opacity:1;' /><path class='atom-2' d='M 36.7 120.9Q 36.7 116.8, 38.7 114.5Q 40.8 112.2, 44.6 112.2Q 48.4 112.2, 50.5 114.5Q 52.5 116.8, 52.5 120.9Q 52.5 125.1, 50.4 127.5Q 48.4 129.9, 44.6 129.9Q 40.8 129.9, 38.7 127.5Q 36.7 125.1, 36.7 120.9M 44.6 127.9Q 47.2 127.9, 48.6 126.2Q 50.1 124.4, 50.1 120.9Q 50.1 117.5, 48.6 115.8Q 47.2 114.1, 44.6 114.1Q 42.0 114.1, 40.5 115.8Q 39.1 117.5, 39.1 120.9Q 39.1 124.4, 40.5 126.2Q 42.0 127.9, 44.6 127.9' fill='#000000'/><path class='atom-5' d='M 97.5 15.6Q 97.5 11.5, 99.5 9.2Q 101.6 6.9, 105.4 6.9Q 109.2 6.9, 111.3 9.2Q 113.3 11.5, 113.3 15.6Q 113.3 19.8, 111.2 22.2Q 109.2 24.5, 105.4 24.5Q 101.6 24.5, 99.5 22.2Q 97.5 19.8, 97.5 15.6M 105.4 22.6Q 108.0 22.6, 109.4 20.8Q 110.9 19.1, 110.9 15.6Q 110.9 12.2, 109.4 10.5Q 108.0 8.8, 105.4 8.8Q 102.8 8.8, 101.3 10.5Q 99.9 12.2, 99.9 15.6Q 99.9 19.1, 101.3 20.8Q 102.8 22.6, 105.4 22.6' fill='#000000'/></svg>",
        "aceticAcid": "<?xml version='1.0' encoding='iso-8859-1'?><svg version='1.1' baseProfile='full'              xmlns='http://www.w3.org/2000/svg'                      xmlns:rdkit='http://www.rdkit.org/xml'                      xmlns:xlink='http://www.w3.org/1999/xlink'                  xml:space='preserve'width='137px' height='127px' viewBox='0 0 137 127'><!-- END OF HEADER --><path class='bond-0 atom-0 atom-1' d='M 6.9,106.6 L 56.9,77.7' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' /><path class='bond-1 atom-1 atom-2' d='M 61.2,80.2 L 61.2,30.3' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' /><path class='bond-1 atom-1 atom-2' d='M 52.6,80.2 L 52.6,30.3' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' /><path class='bond-2 atom-1 atom-3' d='M 56.9,77.7 L 97.6,101.1' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' /><path d='M 54.4,79.1 L 56.9,77.7 L 58.9,78.8' style='fill:none;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-miterlimit:10;stroke-opacity:1;' /><path class='atom-2' d='M 49.4 19.9Q 49.4 16.0, 51.3 13.8Q 53.3 11.6, 56.9 11.6Q 60.5 11.6, 62.5 13.8Q 64.4 16.0, 64.4 19.9Q 64.4 23.9, 62.5 26.2Q 60.5 28.4, 56.9 28.4Q 53.3 28.4, 51.3 26.2Q 49.4 23.9, 49.4 19.9M 56.9 26.5Q 59.4 26.5, 60.8 24.9Q 62.1 23.2, 62.1 19.9Q 62.1 16.7, 60.8 15.1Q 59.4 13.4, 56.9 13.4Q 54.4 13.4, 53.1 15.1Q 51.7 16.7, 51.7 19.9Q 51.7 23.2, 53.1 24.9Q 54.4 26.5, 56.9 26.5' fill='#000000'/><path class='atom-3' d='M 99.5 106.6Q 99.5 102.7, 101.4 100.5Q 103.3 98.3, 107.0 98.3Q 110.6 98.3, 112.5 100.5Q 114.5 102.7, 114.5 106.6Q 114.5 110.6, 112.5 112.9Q 110.6 115.1, 107.0 115.1Q 103.4 115.1, 101.4 112.9Q 99.5 110.6, 99.5 106.6M 107.0 113.3Q 109.5 113.3, 110.8 111.6Q 112.2 109.9, 112.2 106.6Q 112.2 103.4, 110.8 101.8Q 109.5 100.2, 107.0 100.2Q 104.5 100.2, 103.1 101.8Q 101.8 103.4, 101.8 106.6Q 101.8 109.9, 103.1 111.6Q 104.5 113.3, 107.0 113.3' fill='#000000'/><path class='atom-3' d='M 117.0 98.5L 119.3 98.5L 119.3 105.4L 127.6 105.4L 127.6 98.5L 129.8 98.5L 129.8 114.9L 127.6 114.9L 127.6 107.3L 119.3 107.3L 119.3 114.9L 117.0 114.9L 117.0 98.5' fill='#000000'/></svg>",
        "histamine": "<?xml version='1.0' encoding='iso-8859-1'?><svg version='1.1' baseProfile='full'              xmlns='http://www.w3.org/2000/svg'                      xmlns:rdkit='http://www.rdkit.org/xml'                      xmlns:xlink='http://www.w3.org/1999/xlink'                  xml:space='preserve'width='251px' height='130px' viewBox='0 0 251 130'><!-- END OF HEADER --><path class='bond-0 atom-0 atom-1' d='M 42.4,63.6 L 78.3,54.4' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' /><path class='bond-1 atom-1 atom-2' d='M 78.3,54.4 L 108.7,85.3' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' /><path class='bond-2 atom-2 atom-3' d='M 108.7,85.3 L 150.8,74.5' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' /><path class='bond-3 atom-3 atom-4' d='M 150.8,74.5 L 184.3,102.1' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' /><path class='bond-3 atom-3 atom-4' d='M 158.6,72.4 L 184.8,94.0' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' /><path class='bond-4 atom-4 atom-5' d='M 184.3,102.1 L 214.8,82.6' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' /><path class='bond-5 atom-5 atom-6' d='M 218.9,71.2 L 210.0,36.7' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' /><path class='bond-6 atom-6 atom-7' d='M 210.0,36.7 L 172.8,34.4' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' /><path class='bond-6 atom-6 atom-7' d='M 204.9,42.9 L 172.4,40.9' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' /><path class='bond-7 atom-7 atom-3' d='M 163.7,41.6 L 150.8,74.5' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' /><path d='M 76.5,54.8 L 78.3,54.4 L 79.8,55.9' style='fill:none;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-miterlimit:10;stroke-opacity:1;' /><path d='M 107.2,83.8 L 108.7,85.3 L 110.8,84.8' style='fill:none;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-miterlimit:10;stroke-opacity:1;' /><path d='M 182.6,100.7 L 184.3,102.1 L 185.8,101.1' style='fill:none;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-miterlimit:10;stroke-opacity:1;' /><path d='M 210.5,38.4 L 210.0,36.7 L 208.2,36.6' style='fill:none;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-miterlimit:10;stroke-opacity:1;' /><path class='atom-0' d='M 12.6 59.1L 14.2 59.1L 14.2 64.3L 20.5 64.3L 20.5 59.1L 22.2 59.1L 22.2 71.4L 20.5 71.4L 20.5 65.7L 14.2 65.7L 14.2 71.4L 12.6 71.4L 12.6 59.1' fill='#000000'/><path class='atom-0' d='M 24.6 70.9Q 24.9 70.2, 25.6 69.7Q 26.3 69.3, 27.3 69.3Q 28.5 69.3, 29.2 70.0Q 29.9 70.6, 29.9 71.8Q 29.9 73.0, 29.0 74.1Q 28.1 75.3, 26.3 76.6L 30.0 76.6L 30.0 77.5L 24.6 77.5L 24.6 76.7Q 26.1 75.7, 27.0 74.9Q 27.9 74.1, 28.3 73.3Q 28.7 72.6, 28.7 71.9Q 28.7 71.1, 28.3 70.7Q 28.0 70.2, 27.3 70.2Q 26.6 70.2, 26.2 70.5Q 25.8 70.8, 25.4 71.3L 24.6 70.9' fill='#000000'/><path class='atom-0' d='M 33.6 59.1L 37.6 65.6Q 38.0 66.2, 38.6 67.4Q 39.3 68.6, 39.3 68.6L 39.3 59.1L 40.9 59.1L 40.9 71.4L 39.3 71.4L 34.9 64.2Q 34.4 63.4, 33.9 62.5Q 33.4 61.5, 33.2 61.2L 33.2 71.4L 31.6 71.4L 31.6 59.1L 33.6 59.1' fill='#000000'/><path class='atom-5' d='M 218.2 72.6L 222.2 79.1Q 222.6 79.8, 223.2 80.9Q 223.9 82.1, 223.9 82.2L 223.9 72.6L 225.5 72.6L 225.5 84.9L 223.9 84.9L 219.5 77.8Q 219.0 77.0, 218.5 76.0Q 218.0 75.0, 217.8 74.7L 217.8 84.9L 216.2 84.9L 216.2 72.6L 218.2 72.6' fill='#000000'/><path class='atom-5' d='M 227.9 72.6L 229.6 72.6L 229.6 77.8L 235.9 77.8L 235.9 72.6L 237.6 72.6L 237.6 84.9L 235.9 84.9L 235.9 79.2L 229.6 79.2L 229.6 84.9L 227.9 84.9L 227.9 72.6' fill='#000000'/><path class='atom-7' d='M 164.0 27.9L 168.0 34.4Q 168.4 35.1, 169.0 36.2Q 169.7 37.4, 169.7 37.5L 169.7 27.9L 171.3 27.9L 171.3 40.2L 169.7 40.2L 165.3 33.1Q 164.8 32.3, 164.3 31.3Q 163.8 30.4, 163.6 30.1L 163.6 40.2L 162.0 40.2L 162.0 27.9L 164.0 27.9' fill='#000000'/></svg>",
        "amphetamine": "<?xml version='1.0' encoding='iso-8859-1'?><svg version='1.1' baseProfile='full'              xmlns='http://www.w3.org/2000/svg'                      xmlns:rdkit='http://www.rdkit.org/xml'                      xmlns:xlink='http://www.w3.org/1999/xlink'                  xml:space='preserve'width='264px' height='148px' viewBox='0 0 264 148'><!-- END OF HEADER --><path class='bond-0 atom-0 atom-1' d='M 13.2,94.7 L 62.6,77.6' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' /><path class='bond-1 atom-1 atom-2' d='M 62.6,77.6 L 70.7,35.4' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' /><path class='bond-2 atom-1 atom-3' d='M 62.6,77.6 L 102.1,111.8' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' /><path class='bond-3 atom-3 atom-4' d='M 102.1,111.8 L 151.4,94.7' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' /><path class='bond-4 atom-4 atom-5' d='M 151.4,94.7 L 161.3,43.4' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' /><path class='bond-4 atom-4 atom-5' d='M 160.0,91.7 L 168.2,49.3' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' /><path class='bond-5 atom-5 atom-6' d='M 161.3,43.4 L 210.7,26.3' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' /><path class='bond-6 atom-6 atom-7' d='M 210.7,26.3 L 250.2,60.5' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' /><path class='bond-6 atom-6 atom-7' d='M 209.0,35.2 L 241.6,63.5' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' /><path class='bond-7 atom-7 atom-8' d='M 250.2,60.5 L 240.3,111.8' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' /><path class='bond-8 atom-8 atom-9' d='M 240.3,111.8 L 190.9,128.9' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' /><path class='bond-8 atom-8 atom-9' d='M 233.5,105.9 L 192.7,120.0' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' /><path class='bond-9 atom-9 atom-4' d='M 190.9,128.9 L 151.4,94.7' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' /><path d='M 100.1,110.1 L 102.1,111.8 L 104.5,111.0' style='fill:none;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-miterlimit:10;stroke-opacity:1;' /><path d='M 160.8,46.0 L 161.3,43.4 L 163.8,42.5' style='fill:none;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-miterlimit:10;stroke-opacity:1;' /><path d='M 208.2,27.1 L 210.7,26.3 L 212.7,28.0' style='fill:none;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-miterlimit:10;stroke-opacity:1;' /><path d='M 248.2,58.8 L 250.2,60.5 L 249.7,63.1' style='fill:none;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-miterlimit:10;stroke-opacity:1;' /><path d='M 240.8,109.2 L 240.3,111.8 L 237.9,112.7' style='fill:none;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-miterlimit:10;stroke-opacity:1;' /><path d='M 193.4,128.1 L 190.9,128.9 L 189.0,127.2' style='fill:none;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-miterlimit:10;stroke-opacity:1;' /><path class='atom-2' d='M 69.2 18.9L 74.0 26.7Q 74.5 27.5, 75.3 28.9Q 76.1 30.3, 76.1 30.4L 76.1 18.9L 78.1 18.9L 78.1 33.7L 76.0 33.7L 70.8 25.1Q 70.2 24.1, 69.6 23.0Q 68.9 21.8, 68.8 21.5L 68.8 33.7L 66.8 33.7L 66.8 18.9L 69.2 18.9' fill='#000000'/><path class='atom-2' d='M 80.9 18.9L 82.9 18.9L 82.9 25.2L 90.5 25.2L 90.5 18.9L 92.5 18.9L 92.5 33.7L 90.5 33.7L 90.5 26.9L 82.9 26.9L 82.9 33.7L 80.9 33.7L 80.9 18.9' fill='#000000'/><path class='atom-2' d='M 95.4 33.2Q 95.7 32.2, 96.6 31.7Q 97.4 31.2, 98.6 31.2Q 100.1 31.2, 100.9 32.0Q 101.8 32.8, 101.8 34.2Q 101.8 35.7, 100.7 37.0Q 99.6 38.4, 97.4 40.0L 101.9 40.0L 101.9 41.1L 95.3 41.1L 95.3 40.2Q 97.2 38.9, 98.2 37.9Q 99.3 36.9, 99.8 36.1Q 100.4 35.2, 100.4 34.3Q 100.4 33.4, 99.9 32.8Q 99.4 32.3, 98.6 32.3Q 97.8 32.3, 97.3 32.6Q 96.8 32.9, 96.4 33.7L 95.4 33.2' fill='#000000'/></svg>",
        "diethylEther": "<?xml version='1.0' encoding='iso-8859-1'?><svg version='1.1' baseProfile='full'              xmlns='http://www.w3.org/2000/svg'                      xmlns:rdkit='http://www.rdkit.org/xml'                      xmlns:xlink='http://www.w3.org/1999/xlink'                  xml:space='preserve'width='215px' height='100px' viewBox='0 0 215 100'><!-- END OF HEADER --><path class='bond-0 atom-0 atom-1' d='M 10.8,59.8 L 59.0,32.0' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' /><path class='bond-1 atom-1 atom-2' d='M 59.0,32.0 L 98.2,54.6' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' /><path class='bond-2 atom-2 atom-3' d='M 116.3,54.6 L 155.5,32.0' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' /><path class='bond-3 atom-3 atom-4' d='M 155.5,32.0 L 203.8,59.8' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' /><path d='M 56.6,33.4 L 59.0,32.0 L 61.0,33.1' style='fill:none;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-miterlimit:10;stroke-opacity:1;' /><path d='M 153.5,33.1 L 155.5,32.0 L 157.9,33.4' style='fill:none;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-miterlimit:10;stroke-opacity:1;' /><path class='atom-2' d='M 100.0 59.9Q 100.0 56.1, 101.9 54.0Q 103.8 51.8, 107.2 51.8Q 110.7 51.8, 112.6 54.0Q 114.5 56.1, 114.5 59.9Q 114.5 63.7, 112.6 65.9Q 110.7 68.0, 107.2 68.0Q 103.8 68.0, 101.9 65.9Q 100.0 63.7, 100.0 59.9M 107.2 66.3Q 109.7 66.3, 110.9 64.7Q 112.3 63.0, 112.3 59.9Q 112.3 56.8, 110.9 55.2Q 109.7 53.6, 107.2 53.6Q 104.8 53.6, 103.5 55.2Q 102.2 56.7, 102.2 59.9Q 102.2 63.0, 103.5 64.7Q 104.8 66.3, 107.2 66.3' fill='#000000'/></svg>",
        "diisopropylEther": "<?xml version='1.0' encoding='iso-8859-1'?><svg version='1.1' baseProfile='full'              xmlns='http://www.w3.org/2000/svg'                      xmlns:rdkit='http://www.rdkit.org/xml'                      xmlns:xlink='http://www.w3.org/1999/xlink'                  xml:space='preserve'width='215px' height='127px' viewBox='0 0 215 127'><!-- END OF HEADER --><path class='bond-0 atom-0 atom-1' d='M 10.8,25.6 L 59.0,53.4' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' /><path class='bond-1 atom-1 atom-2' d='M 59.0,53.4 L 59.0,109.1' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' /><path class='bond-2 atom-1 atom-3' d='M 59.0,53.4 L 98.2,30.8' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' /><path class='bond-3 atom-3 atom-4' d='M 116.3,30.8 L 155.5,53.4' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' /><path class='bond-4 atom-4 atom-5' d='M 155.5,53.4 L 203.8,25.6' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' /><path class='bond-5 atom-4 atom-6' d='M 155.5,53.4 L 155.5,109.1' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' /><path class='atom-3' d='M 100.0 25.6Q 100.0 21.8, 101.9 19.7Q 103.8 17.6, 107.2 17.6Q 110.7 17.6, 112.6 19.7Q 114.5 21.8, 114.5 25.6Q 114.5 29.4, 112.6 31.6Q 110.7 33.8, 107.2 33.8Q 103.8 33.8, 101.9 31.6Q 100.0 29.5, 100.0 25.6M 107.2 32.0Q 109.7 32.0, 110.9 30.4Q 112.3 28.8, 112.3 25.6Q 112.3 22.5, 110.9 20.9Q 109.7 19.4, 107.2 19.4Q 104.8 19.4, 103.5 20.9Q 102.2 22.5, 102.2 25.6Q 102.2 28.8, 103.5 30.4Q 104.8 32.0, 107.2 32.0' fill='#000000'/></svg>",
        "MTBE": "<?xml version='1.0' encoding='iso-8859-1'?><svg version='1.1' baseProfile='full'              xmlns='http://www.w3.org/2000/svg'                      xmlns:rdkit='http://www.rdkit.org/xml'                      xmlns:xlink='http://www.w3.org/1999/xlink'                  xml:space='preserve'width='176px' height='137px' viewBox='0 0 176 137'><!-- END OF HEADER --><path class='bond-0 atom-0 atom-1' d='M 8.8,68.3 L 51.6,43.7' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' /><path class='bond-1 atom-1 atom-2' d='M 71.4,43.7 L 114.1,68.3' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' /><path class='bond-2 atom-2 atom-3' d='M 114.1,68.3 L 144.5,15.7' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' /><path class='bond-3 atom-2 atom-4' d='M 114.1,68.3 L 83.7,121.0' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' /><path class='bond-4 atom-2 atom-5' d='M 114.1,68.3 L 166.8,98.8' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' /><path class='atom-1' d='M 53.6 38.0Q 53.6 33.9, 55.6 31.5Q 57.6 29.2, 61.5 29.2Q 65.3 29.2, 67.3 31.5Q 69.4 33.9, 69.4 38.0Q 69.4 42.2, 67.3 44.6Q 65.2 46.9, 61.5 46.9Q 57.7 46.9, 55.6 44.6Q 53.6 42.2, 53.6 38.0M 61.5 45.0Q 64.1 45.0, 65.5 43.2Q 66.9 41.4, 66.9 38.0Q 66.9 34.6, 65.5 32.9Q 64.1 31.2, 61.5 31.2Q 58.8 31.2, 57.4 32.9Q 56.0 34.6, 56.0 38.0Q 56.0 41.5, 57.4 43.2Q 58.8 45.0, 61.5 45.0' fill='#000000'/></svg>",
        "aspirin": "<?xml version='1.0' encoding='iso-8859-1'?><svg version='1.1' baseProfile='full'              xmlns='http://www.w3.org/2000/svg'                      xmlns:rdkit='http://www.rdkit.org/xml'                      xmlns:xlink='http://www.w3.org/1999/xlink'                  xml:space='preserve'width='264px' height='222px' viewBox='0 0 264 222'><!-- END OF HEADER --><path class='bond-0 atom-0 atom-1' d='M 13.2,151.3 L 50.8,117.4' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' /><path class='bond-1 atom-1 atom-2' d='M 55.0,118.8 L 45.9,76.2' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' /><path class='bond-1 atom-1 atom-2' d='M 47.6,120.3 L 38.5,77.8' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' /><path class='bond-2 atom-1 atom-3' d='M 50.8,117.4 L 90.8,130.3' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' /><path class='bond-3 atom-3 atom-4' d='M 107.3,125.6 L 136.6,99.1' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' /><path class='bond-4 atom-4 atom-5' d='M 136.6,99.1 L 126.0,49.5' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' /><path class='bond-4 atom-4 atom-5' d='M 143.1,93.2 L 134.4,52.2' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' /><path class='bond-5 atom-5 atom-6' d='M 126.0,49.5 L 163.6,15.6' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' /><path class='bond-6 atom-6 atom-7' d='M 163.6,15.6 L 211.8,31.2' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' /><path class='bond-6 atom-6 atom-7' d='M 165.5,24.2 L 205.3,37.1' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' /><path class='bond-7 atom-7 atom-8' d='M 211.8,31.2 L 222.4,80.7' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' /><path class='bond-8 atom-8 atom-9' d='M 222.4,80.7 L 184.8,114.7' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' /><path class='bond-8 atom-8 atom-9' d='M 214.1,78.0 L 183.0,106.1' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' /><path class='bond-9 atom-9 atom-10' d='M 184.8,114.7 L 195.4,164.2' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' /><path class='bond-10 atom-10 atom-11' d='M 192.2,167.1 L 234.2,180.7' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' /><path class='bond-10 atom-10 atom-11' d='M 194.5,159.9 L 236.5,173.5' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' /><path class='bond-11 atom-10 atom-12' d='M 195.4,164.2 L 166.1,190.7' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' /><path class='bond-12 atom-9 atom-4' d='M 184.8,114.7 L 136.6,99.1' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' /><path d='M 48.9,119.1 L 50.8,117.4 L 52.8,118.0' style='fill:none;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-miterlimit:10;stroke-opacity:1;' /><path d='M 126.6,52.0 L 126.0,49.5 L 127.9,47.8' style='fill:none;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-miterlimit:10;stroke-opacity:1;' /><path d='M 161.8,17.3 L 163.6,15.6 L 166.1,16.4' style='fill:none;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-miterlimit:10;stroke-opacity:1;' /><path d='M 209.4,30.4 L 211.8,31.2 L 212.4,33.7' style='fill:none;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-miterlimit:10;stroke-opacity:1;' /><path d='M 221.9,78.2 L 222.4,80.7 L 220.6,82.4' style='fill:none;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-miterlimit:10;stroke-opacity:1;' /><path d='M 194.9,161.7 L 195.4,164.2 L 193.9,165.5' style='fill:none;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-miterlimit:10;stroke-opacity:1;' /><path class='atom-2' d='M 33.6 67.9Q 33.6 64.5, 35.3 62.5Q 37.0 60.6, 40.2 60.6Q 43.4 60.6, 45.1 62.5Q 46.8 64.5, 46.8 67.9Q 46.8 71.4, 45.1 73.4Q 43.4 75.3, 40.2 75.3Q 37.1 75.3, 35.3 73.4Q 33.6 71.4, 33.6 67.9M 40.2 73.7Q 42.4 73.7, 43.6 72.3Q 44.8 70.8, 44.8 67.9Q 44.8 65.1, 43.6 63.7Q 42.4 62.2, 40.2 62.2Q 38.0 62.2, 36.8 63.6Q 35.7 65.1, 35.7 67.9Q 35.7 70.8, 36.8 72.3Q 38.0 73.7, 40.2 73.7' fill='#000000'/><path class='atom-3' d='M 92.4 133.0Q 92.4 129.6, 94.1 127.7Q 95.8 125.7, 99.0 125.7Q 102.2 125.7, 103.9 127.7Q 105.6 129.6, 105.6 133.0Q 105.6 136.5, 103.9 138.5Q 102.2 140.5, 99.0 140.5Q 95.9 140.5, 94.1 138.5Q 92.4 136.5, 92.4 133.0M 99.0 138.9Q 101.2 138.9, 102.4 137.4Q 103.6 135.9, 103.6 133.0Q 103.6 130.2, 102.4 128.8Q 101.2 127.4, 99.0 127.4Q 96.8 127.4, 95.6 128.8Q 94.5 130.2, 94.5 133.0Q 94.5 135.9, 95.6 137.4Q 96.8 138.9, 99.0 138.9' fill='#000000'/><path class='atom-11' d='M 237.0 179.8Q 237.0 176.4, 238.7 174.5Q 240.4 172.6, 243.6 172.6Q 246.8 172.6, 248.5 174.5Q 250.2 176.4, 250.2 179.8Q 250.2 183.3, 248.5 185.3Q 246.8 187.3, 243.6 187.3Q 240.5 187.3, 238.7 185.3Q 237.0 183.4, 237.0 179.8M 243.6 185.7Q 245.8 185.7, 247.0 184.2Q 248.2 182.7, 248.2 179.8Q 248.2 177.0, 247.0 175.6Q 245.8 174.2, 243.6 174.2Q 241.4 174.2, 240.2 175.6Q 239.1 177.0, 239.1 179.8Q 239.1 182.7, 240.2 184.2Q 241.4 185.7, 243.6 185.7' fill='#000000'/><path class='atom-12' d='M 137.0 191.1L 139.0 191.1L 139.0 197.2L 146.3 197.2L 146.3 191.1L 148.3 191.1L 148.3 205.4L 146.3 205.4L 146.3 198.8L 139.0 198.8L 139.0 205.4L 137.0 205.4L 137.0 191.1' fill='#000000'/><path class='atom-12' d='M 151.2 198.2Q 151.2 194.7, 152.9 192.8Q 154.6 190.9, 157.8 190.9Q 161.0 190.9, 162.7 192.8Q 164.4 194.7, 164.4 198.2Q 164.4 201.7, 162.7 203.7Q 160.9 205.6, 157.8 205.6Q 154.6 205.6, 152.9 203.7Q 151.2 201.7, 151.2 198.2M 157.8 204.0Q 160.0 204.0, 161.2 202.5Q 162.4 201.1, 162.4 198.2Q 162.4 195.4, 161.2 194.0Q 160.0 192.5, 157.8 192.5Q 155.6 192.5, 154.4 193.9Q 153.2 195.4, 153.2 198.2Q 153.2 201.1, 154.4 202.5Q 155.6 204.0, 157.8 204.0' fill='#000000'/></svg>",
        "paracetamol": "<?xml version='1.0' encoding='iso-8859-1'?><svg version='1.1' baseProfile='full'              xmlns='http://www.w3.org/2000/svg'                      xmlns:rdkit='http://www.rdkit.org/xml'                      xmlns:xlink='http://www.w3.org/1999/xlink'                  xml:space='preserve'width='309px' height='155px' viewBox='0 0 309 155'><!-- END OF HEADER --><path class='bond-0 atom-0 atom-1' d='M 15.5,88.2 L 60.7,76.1' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' /><path class='bond-1 atom-1 atom-2' d='M 63.6,78.9 L 74.0,40.2' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' /><path class='bond-1 atom-1 atom-2' d='M 56.8,77.1 L 67.2,38.3' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' /><path class='bond-2 atom-1 atom-3' d='M 60.7,76.1 L 87.3,102.6' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' /><path class='bond-3 atom-3 atom-4' d='M 100.4,107.4 L 139.1,97.1' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' /><path class='bond-4 atom-4 atom-5' d='M 139.1,97.1 L 151.2,51.8' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' /><path class='bond-4 atom-4 atom-5' d='M 147.0,95.0 L 157.0,57.5' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' /><path class='bond-5 atom-5 atom-6' d='M 151.2,51.8 L 196.5,39.7' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' /><path class='bond-6 atom-6 atom-7' d='M 196.5,39.7 L 229.6,72.8' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' /><path class='bond-6 atom-6 atom-7' d='M 194.4,47.5 L 221.8,74.9' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' /><path class='bond-7 atom-7 atom-8' d='M 229.6,72.8 L 267.3,62.7' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' /><path class='bond-8 atom-7 atom-9' d='M 229.6,72.8 L 217.5,118.1' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' /><path class='bond-9 atom-9 atom-10' d='M 217.5,118.1 L 172.3,130.2' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' /><path class='bond-9 atom-9 atom-10' d='M 211.8,112.3 L 174.4,122.4' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' /><path class='bond-10 atom-10 atom-4' d='M 172.3,130.2 L 139.1,97.1' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' /><path d='M 58.5,76.7 L 60.7,76.1 L 62.0,77.4' style='fill:none;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-miterlimit:10;stroke-opacity:1;' /><path d='M 150.6,54.1 L 151.2,51.8 L 153.5,51.2' style='fill:none;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-miterlimit:10;stroke-opacity:1;' /><path d='M 194.2,40.3 L 196.5,39.7 L 198.2,41.3' style='fill:none;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-miterlimit:10;stroke-opacity:1;' /><path d='M 218.1,115.8 L 217.5,118.1 L 215.3,118.7' style='fill:none;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-miterlimit:10;stroke-opacity:1;' /><path d='M 174.5,129.6 L 172.3,130.2 L 170.6,128.5' style='fill:none;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-miterlimit:10;stroke-opacity:1;' /><path class='atom-2' d='M 66.8 30.8Q 66.8 27.6, 68.3 25.9Q 69.9 24.1, 72.8 24.1Q 75.8 24.1, 77.4 25.9Q 78.9 27.6, 78.9 30.8Q 78.9 34.0, 77.3 35.9Q 75.7 37.7, 72.8 37.7Q 69.9 37.7, 68.3 35.9Q 66.8 34.1, 66.8 30.8M 72.8 36.2Q 74.9 36.2, 76.0 34.9Q 77.1 33.5, 77.1 30.8Q 77.1 28.2, 76.0 26.9Q 74.9 25.6, 72.8 25.6Q 70.8 25.6, 69.7 26.9Q 68.6 28.2, 68.6 30.8Q 68.6 33.5, 69.7 34.9Q 70.8 36.2, 72.8 36.2' fill='#000000'/><path class='atom-3' d='M 90.9 102.6L 95.3 109.6Q 95.7 110.3, 96.4 111.5Q 97.1 112.8, 97.1 112.9L 97.1 102.6L 98.9 102.6L 98.9 115.8L 97.1 115.8L 92.4 108.1Q 91.9 107.2, 91.3 106.2Q 90.7 105.2, 90.5 104.9L 90.5 115.8L 88.8 115.8L 88.8 102.6L 90.9 102.6' fill='#000000'/><path class='atom-3' d='M 88.7 117.2L 90.5 117.2L 90.5 122.8L 97.2 122.8L 97.2 117.2L 99.0 117.2L 99.0 130.4L 97.2 130.4L 97.2 124.3L 90.5 124.3L 90.5 130.4L 88.7 130.4L 88.7 117.2' fill='#000000'/><path class='atom-8' d='M 268.8 60.7Q 268.8 57.5, 270.4 55.7Q 272.0 54.0, 274.9 54.0Q 277.9 54.0, 279.4 55.7Q 281.0 57.5, 281.0 60.7Q 281.0 63.9, 279.4 65.8Q 277.8 67.6, 274.9 67.6Q 272.0 67.6, 270.4 65.8Q 268.8 64.0, 268.8 60.7M 274.9 66.1Q 276.9 66.1, 278.0 64.7Q 279.1 63.4, 279.1 60.7Q 279.1 58.1, 278.0 56.8Q 276.9 55.5, 274.9 55.5Q 272.9 55.5, 271.8 56.8Q 270.7 58.1, 270.7 60.7Q 270.7 63.4, 271.8 64.7Q 272.9 66.1, 274.9 66.1' fill='#000000'/><path class='atom-8' d='M 283.1 54.1L 284.9 54.1L 284.9 59.8L 291.7 59.8L 291.7 54.1L 293.4 54.1L 293.4 67.4L 291.7 67.4L 291.7 61.3L 284.9 61.3L 284.9 67.4L 283.1 67.4L 283.1 54.1' fill='#000000'/></svg>",
        "ibuprofen": "<?xml version='1.0' encoding='iso-8859-1'?><svg version='1.1' baseProfile='full'              xmlns='http://www.w3.org/2000/svg'                      xmlns:rdkit='http://www.rdkit.org/xml'                      xmlns:xlink='http://www.w3.org/1999/xlink'                  xml:space='preserve'width='374px' height='183px' viewBox='0 0 374 183'><!-- END OF HEADER --><path class='bond-0 atom-0 atom-1' d='M 18.7,86.6 L 64.1,83.8' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' /><path class='bond-1 atom-1 atom-2' d='M 64.1,83.8 L 84.4,43.1' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' /><path class='bond-2 atom-1 atom-3' d='M 64.1,83.8 L 89.3,121.7' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' /><path class='bond-3 atom-3 atom-4' d='M 89.3,121.7 L 134.7,119.0' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' /><path class='bond-4 atom-4 atom-5' d='M 134.7,119.0 L 155.0,78.2' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' /><path class='bond-4 atom-4 atom-5' d='M 142.6,118.5 L 159.3,84.8' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' /><path class='bond-5 atom-5 atom-6' d='M 155.0,78.2 L 200.4,75.4' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' /><path class='bond-6 atom-6 atom-7' d='M 200.4,75.4 L 225.5,113.4' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' /><path class='bond-6 atom-6 atom-7' d='M 196.9,82.5 L 217.7,113.9' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' /><path class='bond-7 atom-7 atom-8' d='M 225.5,113.4 L 271.0,110.6' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' /><path class='bond-8 atom-8 atom-9' d='M 271.0,110.6 L 296.1,148.5' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' /><path class='bond-9 atom-8 atom-10' d='M 271.0,110.6 L 291.3,69.9' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' /><path class='bond-10 atom-10 atom-11' d='M 295.2,69.6 L 274.4,38.2' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' /><path class='bond-10 atom-10 atom-11' d='M 289.5,73.4 L 268.7,42.0' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' /><path class='bond-11 atom-10 atom-12' d='M 291.3,69.9 L 329.3,67.5' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' /><path class='bond-12 atom-7 atom-13' d='M 225.5,113.4 L 205.2,154.1' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' /><path class='bond-13 atom-13 atom-14' d='M 205.2,154.1 L 159.8,156.9' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' /><path class='bond-13 atom-13 atom-14' d='M 200.9,147.5 L 163.3,149.9' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' /><path class='bond-14 atom-14 atom-4' d='M 159.8,156.9 L 134.7,119.0' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' /><path d='M 88.0,119.9 L 89.3,121.7 L 91.5,121.6' style='fill:none;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-miterlimit:10;stroke-opacity:1;' /><path d='M 154.0,80.3 L 155.0,78.2 L 157.3,78.1' style='fill:none;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-miterlimit:10;stroke-opacity:1;' /><path d='M 198.1,75.6 L 200.4,75.4 L 201.7,77.3' style='fill:none;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-miterlimit:10;stroke-opacity:1;' /><path d='M 290.2,71.9 L 291.3,69.9 L 293.2,69.7' style='fill:none;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-miterlimit:10;stroke-opacity:1;' /><path d='M 206.3,152.1 L 205.2,154.1 L 203.0,154.3' style='fill:none;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-miterlimit:10;stroke-opacity:1;' /><path d='M 162.1,156.8 L 159.8,156.9 L 158.6,155.0' style='fill:none;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-miterlimit:10;stroke-opacity:1;' /><path class='atom-11' d='M 260.2 31.9Q 260.2 28.9, 261.7 27.1Q 263.3 25.4, 266.1 25.4Q 269.0 25.4, 270.5 27.1Q 272.1 28.9, 272.1 31.9Q 272.1 35.1, 270.5 36.9Q 269.0 38.6, 266.1 38.6Q 263.3 38.6, 261.7 36.9Q 260.2 35.1, 260.2 31.9M 266.1 37.2Q 268.1 37.2, 269.2 35.9Q 270.2 34.5, 270.2 31.9Q 270.2 29.4, 269.2 28.1Q 268.1 26.9, 266.1 26.9Q 264.2 26.9, 263.1 28.1Q 262.0 29.4, 262.0 31.9Q 262.0 34.6, 263.1 35.9Q 264.2 37.2, 266.1 37.2' fill='#000000'/><path class='atom-12' d='M 330.8 67.1Q 330.8 64.0, 332.3 62.3Q 333.8 60.6, 336.7 60.6Q 339.5 60.6, 341.1 62.3Q 342.6 64.0, 342.6 67.1Q 342.6 70.2, 341.1 72.0Q 339.5 73.8, 336.7 73.8Q 333.9 73.8, 332.3 72.0Q 330.8 70.3, 330.8 67.1M 336.7 72.3Q 338.7 72.3, 339.7 71.0Q 340.8 69.7, 340.8 67.1Q 340.8 64.6, 339.7 63.3Q 338.7 62.0, 336.7 62.0Q 334.7 62.0, 333.7 63.3Q 332.6 64.6, 332.6 67.1Q 332.6 69.7, 333.7 71.0Q 334.7 72.3, 336.7 72.3' fill='#000000'/><path class='atom-12' d='M 344.6 60.7L 346.4 60.7L 346.4 66.2L 353.0 66.2L 353.0 60.7L 354.7 60.7L 354.7 73.6L 353.0 73.6L 353.0 67.6L 346.4 67.6L 346.4 73.6L 344.6 73.6L 344.6 60.7' fill='#000000'/></svg>",
        "dopamine": "<?xml version='1.0' encoding='iso-8859-1'?><svg version='1.1' baseProfile='full'              xmlns='http://www.w3.org/2000/svg'                      xmlns:rdkit='http://www.rdkit.org/xml'                      xmlns:xlink='http://www.w3.org/1999/xlink'                  xml:space='preserve'width='308px' height='165px' viewBox='0 0 308 165'><!-- END OF HEADER --><path class='bond-0 atom-0 atom-1' d='M 44.9,76.9 L 79.9,66.2' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' /><path class='bond-1 atom-1 atom-2' d='M 79.9,66.2 L 111.2,95.6' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' /><path class='bond-2 atom-2 atom-3' d='M 111.2,95.6 L 152.3,83.1' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' /><path class='bond-3 atom-3 atom-4' d='M 152.3,83.1 L 162.0,41.3' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' /><path class='bond-3 atom-3 atom-4' d='M 159.4,81.0 L 167.4,46.4' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' /><path class='bond-4 atom-4 atom-5' d='M 162.0,41.3 L 203.1,28.9' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' /><path class='bond-5 atom-5 atom-6' d='M 203.1,28.9 L 234.4,58.2' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' /><path class='bond-5 atom-5 atom-6' d='M 201.4,36.1 L 227.3,60.4' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' /><path class='bond-6 atom-6 atom-7' d='M 234.4,58.2 L 268.4,47.9' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' /><path class='bond-7 atom-6 atom-8' d='M 234.4,58.2 L 224.6,100.0' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' /><path class='bond-8 atom-8 atom-9' d='M 224.6,100.0 L 249.0,122.7' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' /><path class='bond-9 atom-8 atom-10' d='M 224.6,100.0 L 183.6,112.4' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' /><path class='bond-9 atom-8 atom-10' d='M 219.2,94.9 L 185.3,105.2' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' /><path class='bond-10 atom-10 atom-3' d='M 183.6,112.4 L 152.3,83.1' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' /><path d='M 78.1,66.8 L 79.9,66.2 L 81.5,67.7' style='fill:none;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-miterlimit:10;stroke-opacity:1;' /><path d='M 109.7,94.1 L 111.2,95.6 L 113.3,94.9' style='fill:none;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-miterlimit:10;stroke-opacity:1;' /><path d='M 161.5,43.4 L 162.0,41.3 L 164.1,40.7' style='fill:none;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-miterlimit:10;stroke-opacity:1;' /><path d='M 201.0,29.5 L 203.1,28.9 L 204.6,30.3' style='fill:none;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-miterlimit:10;stroke-opacity:1;' /><path d='M 185.6,111.8 L 183.6,112.4 L 182.0,111.0' style='fill:none;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-miterlimit:10;stroke-opacity:1;' /><path class='atom-0' d='M 15.4 72.6L 17.0 72.6L 17.0 77.8L 23.3 77.8L 23.3 72.6L 24.9 72.6L 24.9 84.8L 23.3 84.8L 23.3 79.2L 17.0 79.2L 17.0 84.8L 15.4 84.8L 15.4 72.6' fill='#000000'/><path class='atom-0' d='M 27.3 84.4Q 27.6 83.6, 28.3 83.2Q 29.0 82.7, 30.0 82.7Q 31.2 82.7, 31.8 83.4Q 32.5 84.1, 32.5 85.2Q 32.5 86.4, 31.6 87.5Q 30.8 88.6, 29.0 89.9L 32.7 89.9L 32.7 90.9L 27.3 90.9L 27.3 90.1Q 28.8 89.0, 29.6 88.2Q 30.5 87.4, 31.0 86.7Q 31.4 86.0, 31.4 85.3Q 31.4 84.5, 31.0 84.1Q 30.6 83.6, 30.0 83.6Q 29.3 83.6, 28.9 83.9Q 28.5 84.2, 28.1 84.7L 27.3 84.4' fill='#000000'/><path class='atom-0' d='M 36.2 72.6L 40.1 79.1Q 40.5 79.7, 41.2 80.8Q 41.8 82.0, 41.8 82.1L 41.8 72.6L 43.5 72.6L 43.5 84.8L 41.8 84.8L 37.5 77.7Q 37.0 76.9, 36.5 76.0Q 36.0 75.0, 35.8 74.7L 35.8 84.8L 34.2 84.8L 34.2 72.6L 36.2 72.6' fill='#000000'/><path class='atom-7' d='M 269.8 45.8Q 269.8 42.9, 271.3 41.2Q 272.7 39.6, 275.4 39.6Q 278.1 39.6, 279.6 41.2Q 281.0 42.9, 281.0 45.8Q 281.0 48.7, 279.5 50.4Q 278.1 52.1, 275.4 52.1Q 272.7 52.1, 271.3 50.4Q 269.8 48.7, 269.8 45.8M 275.4 50.7Q 277.3 50.7, 278.3 49.5Q 279.3 48.2, 279.3 45.8Q 279.3 43.4, 278.3 42.2Q 277.3 41.0, 275.4 41.0Q 273.6 41.0, 272.6 42.2Q 271.6 43.4, 271.6 45.8Q 271.6 48.2, 272.6 49.5Q 273.6 50.7, 275.4 50.7' fill='#000000'/><path class='atom-7' d='M 282.9 39.7L 284.5 39.7L 284.5 44.9L 290.8 44.9L 290.8 39.7L 292.4 39.7L 292.4 51.9L 290.8 51.9L 290.8 46.3L 284.5 46.3L 284.5 51.9L 282.9 51.9L 282.9 39.7' fill='#000000'/><path class='atom-9' d='M 250.4 129.3Q 250.4 126.4, 251.8 124.8Q 253.3 123.2, 256.0 123.2Q 258.6 123.2, 260.1 124.8Q 261.5 126.4, 261.5 129.3Q 261.5 132.3, 260.1 134.0Q 258.6 135.6, 256.0 135.6Q 253.3 135.6, 251.8 134.0Q 250.4 132.3, 250.4 129.3M 256.0 134.3Q 257.8 134.3, 258.8 133.0Q 259.8 131.8, 259.8 129.3Q 259.8 126.9, 258.8 125.7Q 257.8 124.5, 256.0 124.5Q 254.1 124.5, 253.1 125.7Q 252.1 126.9, 252.1 129.3Q 252.1 131.8, 253.1 133.0Q 254.1 134.3, 256.0 134.3' fill='#000000'/><path class='atom-9' d='M 263.4 123.3L 265.1 123.3L 265.1 128.5L 271.3 128.5L 271.3 123.3L 272.9 123.3L 272.9 135.4L 271.3 135.4L 271.3 129.8L 265.1 129.8L 265.1 135.4L 263.4 135.4L 263.4 123.3' fill='#000000'/></svg>",
        "serotonin": "<?xml version='1.0' encoding='iso-8859-1'?><svg version='1.1' baseProfile='full'              xmlns='http://www.w3.org/2000/svg'                      xmlns:rdkit='http://www.rdkit.org/xml'                      xmlns:xlink='http://www.w3.org/1999/xlink'                  xml:space='preserve'width='305px' height='205px' viewBox='0 0 305 205'><!-- END OF HEADER --><path class='bond-0 atom-0 atom-1' d='M 44.7,60.0 L 69.1,84.2' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' /><path class='bond-1 atom-1 atom-2' d='M 69.1,84.2 L 110.4,73.0' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' /><path class='bond-2 atom-2 atom-3' d='M 110.4,73.0 L 140.8,103.2' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' /><path class='bond-3 atom-3 atom-4' d='M 140.8,103.2 L 134.2,145.5' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' /><path class='bond-3 atom-3 atom-4' d='M 146.4,108.8 L 141.3,141.9' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' /><path class='bond-4 atom-4 atom-5' d='M 134.2,145.5 L 166.5,161.8' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' /><path class='bond-5 atom-5 atom-6' d='M 178.5,158.8 L 202.7,134.4' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' /><path class='bond-6 atom-6 atom-7' d='M 202.7,134.4 L 245.5,136.5' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' /><path class='bond-6 atom-6 atom-7' d='M 206.7,128.2 L 242.1,129.9' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' /><path class='bond-7 atom-7 atom-8' d='M 245.5,136.5 L 268.7,100.5' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' /><path class='bond-8 atom-8 atom-9' d='M 268.7,100.5 L 249.1,62.4' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' /><path class='bond-8 atom-8 atom-9' d='M 261.3,100.2 L 245.1,68.6' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' /><path class='bond-9 atom-9 atom-10' d='M 249.1,62.4 L 267.3,34.1' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' /><path class='bond-10 atom-9 atom-11' d='M 249.1,62.4 L 206.3,60.3' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' /><path class='bond-11 atom-11 atom-12' d='M 206.3,60.3 L 183.1,96.3' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' /><path class='bond-11 atom-11 atom-12' d='M 209.7,66.9 L 190.5,96.7' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' /><path class='bond-12 atom-12 atom-3' d='M 183.1,96.3 L 140.8,103.2' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' /><path class='bond-13 atom-12 atom-6' d='M 183.1,96.3 L 202.7,134.4' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' /><path d='M 67.8,83.0 L 69.1,84.2 L 71.1,83.6' style='fill:none;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-miterlimit:10;stroke-opacity:1;' /><path d='M 108.3,73.5 L 110.4,73.0 L 111.9,74.5' style='fill:none;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-miterlimit:10;stroke-opacity:1;' /><path d='M 134.6,143.4 L 134.2,145.5 L 135.9,146.3' style='fill:none;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-miterlimit:10;stroke-opacity:1;' /><path d='M 243.3,136.4 L 245.5,136.5 L 246.6,134.7' style='fill:none;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-miterlimit:10;stroke-opacity:1;' /><path d='M 267.5,102.3 L 268.7,100.5 L 267.7,98.6' style='fill:none;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-miterlimit:10;stroke-opacity:1;' /><path d='M 208.4,60.4 L 206.3,60.3 L 205.1,62.1' style='fill:none;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-miterlimit:10;stroke-opacity:1;' /><path class='atom-0' d='M 15.3 47.9L 16.9 47.9L 16.9 53.1L 23.1 53.1L 23.1 47.9L 24.7 47.9L 24.7 60.1L 23.1 60.1L 23.1 54.5L 16.9 54.5L 16.9 60.1L 15.3 60.1L 15.3 47.9' fill='#000000'/><path class='atom-0' d='M 27.1 59.7Q 27.4 58.9, 28.1 58.5Q 28.8 58.0, 29.8 58.0Q 31.0 58.0, 31.7 58.7Q 32.4 59.4, 32.4 60.5Q 32.4 61.7, 31.5 62.8Q 30.6 63.9, 28.8 65.2L 32.5 65.2L 32.5 66.1L 27.1 66.1L 27.1 65.4Q 28.6 64.3, 29.5 63.5Q 30.4 62.7, 30.8 62.0Q 31.2 61.3, 31.2 60.6Q 31.2 59.8, 30.8 59.4Q 30.5 59.0, 29.8 59.0Q 29.1 59.0, 28.7 59.2Q 28.3 59.5, 28.0 60.1L 27.1 59.7' fill='#000000'/><path class='atom-0' d='M 36.0 47.9L 40.0 54.4Q 40.4 55.0, 41.0 56.2Q 41.6 57.3, 41.7 57.4L 41.7 47.9L 43.3 47.9L 43.3 60.1L 41.6 60.1L 37.3 53.1Q 36.8 52.2, 36.3 51.3Q 35.8 50.3, 35.6 50.1L 35.6 60.1L 34.1 60.1L 34.1 47.9L 36.0 47.9' fill='#000000'/><path class='atom-5' d='M 169.8 158.8L 173.8 165.2Q 174.2 165.8, 174.8 167.0Q 175.4 168.1, 175.5 168.2L 175.5 158.8L 177.1 158.8L 177.1 170.9L 175.4 170.9L 171.1 163.9Q 170.7 163.0, 170.1 162.1Q 169.6 161.2, 169.5 160.9L 169.5 170.9L 167.9 170.9L 167.9 158.8L 169.8 158.8' fill='#000000'/><path class='atom-5' d='M 167.7 172.1L 169.4 172.1L 169.4 177.3L 175.6 177.3L 175.6 172.1L 177.2 172.1L 177.2 184.2L 175.6 184.2L 175.6 178.6L 169.4 178.6L 169.4 184.2L 167.7 184.2L 167.7 172.1' fill='#000000'/><path class='atom-10' d='M 266.7 26.4Q 266.7 23.5, 268.2 21.9Q 269.6 20.3, 272.3 20.3Q 275.0 20.3, 276.4 21.9Q 277.9 23.5, 277.9 26.4Q 277.9 29.4, 276.4 31.1Q 275.0 32.7, 272.3 32.7Q 269.6 32.7, 268.2 31.1Q 266.7 29.4, 266.7 26.4M 272.3 31.4Q 274.1 31.4, 275.1 30.1Q 276.2 28.9, 276.2 26.4Q 276.2 24.1, 275.1 22.9Q 274.1 21.6, 272.3 21.6Q 270.4 21.6, 269.4 22.8Q 268.4 24.0, 268.4 26.4Q 268.4 28.9, 269.4 30.1Q 270.4 31.4, 272.3 31.4' fill='#000000'/><path class='atom-10' d='M 279.8 20.4L 281.4 20.4L 281.4 25.6L 287.6 25.6L 287.6 20.4L 289.3 20.4L 289.3 32.5L 287.6 32.5L 287.6 26.9L 281.4 26.9L 281.4 32.5L 279.8 32.5L 279.8 20.4' fill='#000000'/></svg>",
        "caffeine": "<?xml version='1.0' encoding='iso-8859-1'?><svg version='1.1' baseProfile='full'              xmlns='http://www.w3.org/2000/svg'                      xmlns:rdkit='http://www.rdkit.org/xml'                      xmlns:xlink='http://www.w3.org/1999/xlink'                  xml:space='preserve'width='266px' height='228px' viewBox='0 0 266 228'><!-- END OF HEADER --><path class='bond-0 atom-0 atom-1' d='M 252.3,155.5 L 227.7,124.8' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' /><path class='bond-1 atom-1 atom-2' d='M 223.0,107.1 L 234.0,67.0' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' /><path class='bond-2 atom-2 atom-3' d='M 234.0,67.0 L 198.7,43.9' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' /><path class='bond-2 atom-2 atom-3' d='M 225.2,70.4 L 194.5,50.3' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' /><path class='bond-3 atom-3 atom-4' d='M 184.5,44.9 L 152.0,70.9' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' /><path class='bond-4 atom-4 atom-5' d='M 152.0,70.9 L 170.0,118.3' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' /><path class='bond-4 atom-4 atom-5' d='M 161.1,73.4 L 175.1,110.5' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' /><path class='bond-5 atom-5 atom-6' d='M 170.0,118.3 L 137.9,157.6' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' /><path class='bond-6 atom-6 atom-7' d='M 133.5,156.9 L 148.9,197.4' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' /><path class='bond-6 atom-6 atom-7' d='M 140.6,154.2 L 156.0,194.7' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' /><path class='bond-7 atom-6 atom-8' d='M 137.9,157.6 L 95.0,150.6' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' /><path class='bond-8 atom-8 atom-9' d='M 80.7,158.1 L 55.7,188.6' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' /><path class='bond-9 atom-8 atom-10' d='M 84.5,140.5 L 69.9,102.0' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' /><path class='bond-10 atom-10 atom-11' d='M 72.7,98.6 L 28.8,91.4' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' /><path class='bond-10 atom-10 atom-11' d='M 71.5,106.1 L 27.5,98.9' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' /><path class='bond-11 atom-10 atom-12' d='M 69.9,102.0 L 94.9,71.4' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' /><path class='bond-12 atom-12 atom-13' d='M 98.7,53.9 L 84.1,15.3' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' /><path class='bond-13 atom-5 atom-1' d='M 170.0,118.3 L 213.5,116.3' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' /><path class='bond-14 atom-12 atom-4' d='M 109.1,63.9 L 152.0,70.9' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' /><path d='M 233.4,69.1 L 234.0,67.0 L 232.2,65.9' style='fill:none;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-miterlimit:10;stroke-opacity:1;' /><path d='M 139.5,155.6 L 137.9,157.6 L 135.7,157.2' style='fill:none;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-miterlimit:10;stroke-opacity:1;' /><path d='M 70.6,103.9 L 69.9,102.0 L 71.2,100.4' style='fill:none;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-miterlimit:10;stroke-opacity:1;' /><path class='atom-1' d='M 217.4 108.8L 222.1 116.4Q 222.6 117.1, 223.3 118.5Q 224.1 119.8, 224.1 119.9L 224.1 108.8L 226.0 108.8L 226.0 123.1L 224.1 123.1L 219.0 114.8Q 218.4 113.8, 217.8 112.7Q 217.2 111.6, 217.0 111.3L 217.0 123.1L 215.2 123.1L 215.2 108.8L 217.4 108.8' fill='#000000'/><path class='atom-3' d='M 188.4 32.0L 193.1 39.7Q 193.6 40.4, 194.3 41.8Q 195.1 43.1, 195.1 43.2L 195.1 32.0L 197.0 32.0L 197.0 46.4L 195.1 46.4L 190.0 38.1Q 189.4 37.1, 188.8 36.0Q 188.2 34.9, 188.0 34.5L 188.0 46.4L 186.1 46.4L 186.1 32.0L 188.4 32.0' fill='#000000'/><path class='atom-7' d='M 149.2 205.0Q 149.2 201.6, 150.9 199.7Q 152.6 197.7, 155.8 197.7Q 159.0 197.7, 160.7 199.7Q 162.4 201.6, 162.4 205.0Q 162.4 208.5, 160.7 210.5Q 158.9 212.5, 155.8 212.5Q 152.6 212.5, 150.9 210.5Q 149.2 208.5, 149.2 205.0M 155.8 210.8Q 158.0 210.8, 159.2 209.4Q 160.4 207.9, 160.4 205.0Q 160.4 202.2, 159.2 200.8Q 158.0 199.4, 155.8 199.4Q 153.6 199.4, 152.4 200.8Q 151.2 202.2, 151.2 205.0Q 151.2 207.9, 152.4 209.4Q 153.6 210.8, 155.8 210.8' fill='#000000'/><path class='atom-8' d='M 84.7 142.2L 89.4 149.8Q 89.8 150.6, 90.6 151.9Q 91.3 153.3, 91.4 153.4L 91.4 142.2L 93.3 142.2L 93.3 156.6L 91.3 156.6L 86.3 148.3Q 85.7 147.3, 85.1 146.2Q 84.4 145.1, 84.3 144.7L 84.3 156.6L 82.4 156.6L 82.4 142.2L 84.7 142.2' fill='#000000'/><path class='atom-11' d='M 13.3 93.8Q 13.3 90.4, 15.0 88.5Q 16.7 86.5, 19.9 86.5Q 23.1 86.5, 24.8 88.5Q 26.5 90.4, 26.5 93.8Q 26.5 97.3, 24.8 99.3Q 23.0 101.3, 19.9 101.3Q 16.7 101.3, 15.0 99.3Q 13.3 97.3, 13.3 93.8M 19.9 99.7Q 22.1 99.7, 23.3 98.2Q 24.5 96.7, 24.5 93.8Q 24.5 91.0, 23.3 89.6Q 22.1 88.2, 19.9 88.2Q 17.7 88.2, 16.5 89.6Q 15.3 91.0, 15.3 93.8Q 15.3 96.7, 16.5 98.2Q 17.7 99.7, 19.9 99.7' fill='#000000'/><path class='atom-12' d='M 98.8 55.6L 103.5 63.2Q 104.0 63.9, 104.8 65.3Q 105.5 66.6, 105.6 66.7L 105.6 55.6L 107.5 55.6L 107.5 69.9L 105.5 69.9L 100.4 61.6Q 99.9 60.6, 99.2 59.5Q 98.6 58.4, 98.4 58.1L 98.4 69.9L 96.6 69.9L 96.6 55.6L 98.8 55.6' fill='#000000'/></svg>",
        "methylphenidate": "<?xml version='1.0' encoding='iso-8859-1'?><svg version='1.1' baseProfile='full'              xmlns='http://www.w3.org/2000/svg'                      xmlns:rdkit='http://www.rdkit.org/xml'                      xmlns:xlink='http://www.w3.org/1999/xlink'                  xml:space='preserve'width='297px' height='266px' viewBox='0 0 297 266'><!-- END OF HEADER --><path class='bond-0 atom-0 atom-1' d='M 79.3,248.4 L 85.3,207.6' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' /><path class='bond-1 atom-1 atom-2' d='M 94.9,195.2 L 133.6,179.8' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' /><path class='bond-2 atom-2 atom-3' d='M 129.6,181.5 L 162.6,207.6' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' /><path class='bond-2 atom-2 atom-3' d='M 134.3,175.5 L 167.3,201.7' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' /><path class='bond-3 atom-2 atom-4' d='M 133.6,179.8 L 141.0,129.9' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' /><path class='bond-4 atom-4 atom-5' d='M 141.0,129.9 L 187.9,111.3' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' /><path class='bond-5 atom-5 atom-6' d='M 187.9,111.3 L 227.5,142.6' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' /><path class='bond-6 atom-6 atom-7' d='M 227.5,142.6 L 274.5,124.0' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' /><path class='bond-7 atom-7 atom-8' d='M 274.5,124.0 L 281.9,74.0' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' /><path class='bond-8 atom-8 atom-9' d='M 281.9,74.0 L 242.3,42.7' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' /><path class='bond-9 atom-9 atom-10' d='M 242.3,42.7 L 202.4,58.5' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' /><path class='bond-10 atom-4 atom-11' d='M 141.0,129.9 L 101.4,98.5' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' /><path class='bond-11 atom-11 atom-12' d='M 101.4,98.5 L 108.8,48.5' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' /><path class='bond-11 atom-11 atom-12' d='M 94.5,93.1 L 100.6,51.8' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' /><path class='bond-12 atom-12 atom-13' d='M 108.8,48.5 L 69.2,17.2' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' /><path class='bond-13 atom-13 atom-14' d='M 69.2,17.2 L 22.2,35.8' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' /><path class='bond-13 atom-13 atom-14' d='M 67.9,25.8 L 29.1,41.2' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' /><path class='bond-14 atom-14 atom-15' d='M 22.2,35.8 L 14.9,85.8' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' /><path class='bond-15 atom-15 atom-16' d='M 14.9,85.8 L 54.4,117.1' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' /><path class='bond-15 atom-15 atom-16' d='M 23.0,82.5 L 55.7,108.5' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' /><path class='bond-16 atom-10 atom-5' d='M 194.0,70.1 L 187.9,111.3' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' /><path class='bond-17 atom-16 atom-11' d='M 54.4,117.1 L 101.4,98.5' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' /><path d='M 131.7,180.6 L 133.6,179.8 L 134.0,177.3' style='fill:none;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-miterlimit:10;stroke-opacity:1;' /><path d='M 225.6,141.1 L 227.5,142.6 L 229.9,141.7' style='fill:none;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-miterlimit:10;stroke-opacity:1;' /><path d='M 272.1,124.9 L 274.5,124.0 L 274.9,121.5' style='fill:none;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-miterlimit:10;stroke-opacity:1;' /><path d='M 281.5,76.5 L 281.9,74.0 L 279.9,72.5' style='fill:none;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-miterlimit:10;stroke-opacity:1;' /><path d='M 244.2,44.3 L 242.3,42.7 L 240.3,43.5' style='fill:none;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-miterlimit:10;stroke-opacity:1;' /><path d='M 108.4,51.0 L 108.8,48.5 L 106.8,47.0' style='fill:none;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-miterlimit:10;stroke-opacity:1;' /><path d='M 71.1,18.8 L 69.2,17.2 L 66.8,18.1' style='fill:none;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-miterlimit:10;stroke-opacity:1;' /><path d='M 24.6,34.9 L 22.2,35.8 L 21.8,38.3' style='fill:none;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-miterlimit:10;stroke-opacity:1;' /><path d='M 15.2,83.3 L 14.9,85.8 L 16.8,87.3' style='fill:none;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-miterlimit:10;stroke-opacity:1;' /><path d='M 52.5,115.6 L 54.4,117.1 L 56.8,116.2' style='fill:none;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-miterlimit:10;stroke-opacity:1;' /><path class='atom-1' d='M 80.1 198.5Q 80.1 195.1, 81.8 193.1Q 83.5 191.2, 86.7 191.2Q 89.8 191.2, 91.5 193.1Q 93.2 195.1, 93.2 198.5Q 93.2 202.0, 91.5 203.9Q 89.8 205.9, 86.7 205.9Q 83.5 205.9, 81.8 203.9Q 80.1 202.0, 80.1 198.5M 86.7 204.3Q 88.9 204.3, 90.0 202.8Q 91.2 201.4, 91.2 198.5Q 91.2 195.7, 90.0 194.3Q 88.9 192.8, 86.7 192.8Q 84.5 192.8, 83.3 194.2Q 82.1 195.7, 82.1 198.5Q 82.1 201.4, 83.3 202.8Q 84.5 204.3, 86.7 204.3' fill='#000000'/><path class='atom-3' d='M 166.7 211.2Q 166.7 207.8, 168.3 205.9Q 170.0 204.0, 173.2 204.0Q 176.4 204.0, 178.1 205.9Q 179.8 207.8, 179.8 211.2Q 179.8 214.7, 178.1 216.7Q 176.3 218.7, 173.2 218.7Q 170.1 218.7, 168.3 216.7Q 166.7 214.7, 166.7 211.2M 173.2 217.0Q 175.4 217.0, 176.6 215.6Q 177.8 214.1, 177.8 211.2Q 177.8 208.4, 176.6 207.0Q 175.4 205.6, 173.2 205.6Q 171.0 205.6, 169.8 207.0Q 168.7 208.4, 168.7 211.2Q 168.7 214.1, 169.8 215.6Q 171.0 217.0, 173.2 217.0' fill='#000000'/><path class='atom-10' d='M 175.9 54.1L 177.8 54.1L 177.8 60.2L 185.1 60.2L 185.1 54.1L 187.1 54.1L 187.1 68.4L 185.1 68.4L 185.1 61.8L 177.8 61.8L 177.8 68.4L 175.9 68.4L 175.9 54.1' fill='#000000'/><path class='atom-10' d='M 192.1 54.1L 196.8 61.7Q 197.3 62.5, 198.0 63.8Q 198.8 65.2, 198.8 65.3L 198.8 54.1L 200.7 54.1L 200.7 68.4L 198.8 68.4L 193.7 60.2Q 193.2 59.2, 192.5 58.1Q 191.9 57.0, 191.7 56.6L 191.7 68.4L 189.9 68.4L 189.9 54.1L 192.1 54.1' fill='#000000'/></svg>",
        "atomoxetine": "<?xml version='1.0' encoding='iso-8859-1'?><svg version='1.1' baseProfile='full'              xmlns='http://www.w3.org/2000/svg'                      xmlns:rdkit='http://www.rdkit.org/xml'                      xmlns:xlink='http://www.w3.org/1999/xlink'                  xml:space='preserve'width='373px' height='223px' viewBox='0 0 373 223'><!-- END OF HEADER --><path class='bond-0 atom-0 atom-1' d='M 18.7,154.4 L 52.1,176.3' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' /><path class='bond-1 atom-1 atom-2' d='M 65.6,177.3 L 101.7,159.1' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' /><path class='bond-2 atom-2 atom-3' d='M 101.7,159.1 L 141.9,185.4' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' /><path class='bond-3 atom-3 atom-4' d='M 141.9,185.4 L 184.8,163.7' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' /><path class='bond-4 atom-4 atom-5' d='M 184.8,163.7 L 217.2,184.9' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' /><path class='bond-5 atom-5 atom-6' d='M 232.8,186.1 L 267.9,168.4' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' /><path class='bond-6 atom-6 atom-7' d='M 267.9,168.4 L 308.1,194.7' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' /><path class='bond-6 atom-6 atom-7' d='M 275.3,164.6 L 308.5,186.4' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' /><path class='bond-7 atom-7 atom-8' d='M 308.1,194.7 L 351.0,173.0' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' /><path class='bond-8 atom-8 atom-9' d='M 351.0,173.0 L 353.6,125.1' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' /><path class='bond-8 atom-8 atom-9' d='M 344.0,168.5 L 346.2,128.8' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' /><path class='bond-9 atom-9 atom-10' d='M 353.6,125.1 L 313.5,98.8' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' /><path class='bond-10 atom-10 atom-11' d='M 313.5,98.8 L 270.6,120.4' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' /><path class='bond-10 atom-10 atom-11' d='M 313.0,107.1 L 277.5,125.0' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' /><path class='bond-11 atom-4 atom-12' d='M 184.8,163.7 L 187.5,115.8' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' /><path class='bond-12 atom-12 atom-13' d='M 187.5,115.8 L 230.4,94.1' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' /><path class='bond-12 atom-12 atom-13' d='M 188.0,107.4 L 223.4,89.5' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' /><path class='bond-13 atom-13 atom-14' d='M 230.4,94.1 L 233.1,46.1' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' /><path class='bond-14 atom-14 atom-15' d='M 233.1,46.1 L 192.9,19.8' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' /><path class='bond-14 atom-14 atom-15' d='M 225.6,49.9 L 192.4,28.1' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' /><path class='bond-15 atom-15 atom-16' d='M 192.9,19.8 L 150.0,41.5' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' /><path class='bond-16 atom-16 atom-17' d='M 150.0,41.5 L 147.3,89.4' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' /><path class='bond-16 atom-16 atom-17' d='M 157.0,46.0 L 154.7,85.7' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' /><path class='bond-17 atom-17 atom-18' d='M 147.3,89.4 L 104.4,111.1' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' /><path class='bond-18 atom-11 atom-6' d='M 270.6,120.4 L 267.9,168.4' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' /><path class='bond-19 atom-17 atom-12' d='M 147.3,89.4 L 187.5,115.8' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' /><path d='M 99.9,160.0 L 101.7,159.1 L 103.7,160.4' style='fill:none;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-miterlimit:10;stroke-opacity:1;' /><path d='M 139.9,184.1 L 141.9,185.4 L 144.1,184.3' style='fill:none;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-miterlimit:10;stroke-opacity:1;' /><path d='M 306.1,193.4 L 308.1,194.7 L 310.2,193.6' style='fill:none;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-miterlimit:10;stroke-opacity:1;' /><path d='M 348.8,174.1 L 351.0,173.0 L 351.1,170.6' style='fill:none;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-miterlimit:10;stroke-opacity:1;' /><path d='M 353.5,127.5 L 353.6,125.1 L 351.6,123.7' style='fill:none;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-miterlimit:10;stroke-opacity:1;' /><path d='M 315.5,100.1 L 313.5,98.8 L 311.3,99.8' style='fill:none;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-miterlimit:10;stroke-opacity:1;' /><path d='M 272.7,119.3 L 270.6,120.4 L 270.4,122.8' style='fill:none;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-miterlimit:10;stroke-opacity:1;' /><path d='M 228.2,95.2 L 230.4,94.1 L 230.5,91.7' style='fill:none;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-miterlimit:10;stroke-opacity:1;' /><path d='M 232.9,48.5 L 233.1,46.1 L 231.1,44.8' style='fill:none;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-miterlimit:10;stroke-opacity:1;' /><path d='M 194.9,21.1 L 192.9,19.8 L 190.7,20.9' style='fill:none;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-miterlimit:10;stroke-opacity:1;' /><path d='M 152.1,40.4 L 150.0,41.5 L 149.9,43.9' style='fill:none;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-miterlimit:10;stroke-opacity:1;' /><path class='atom-1' d='M 55.8 173.9L 60.3 181.1Q 60.7 181.8, 61.4 183.1Q 62.2 184.4, 62.2 184.5L 62.2 173.9L 64.0 173.9L 64.0 187.5L 62.1 187.5L 57.4 179.6Q 56.8 178.7, 56.2 177.7Q 55.6 176.6, 55.5 176.3L 55.5 187.5L 53.7 187.5L 53.7 173.9L 55.8 173.9' fill='#000000'/><path class='atom-1' d='M 53.5 188.9L 55.4 188.9L 55.4 194.7L 62.3 194.7L 62.3 188.9L 64.2 188.9L 64.2 202.5L 62.3 202.5L 62.3 196.2L 55.4 196.2L 55.4 202.5L 53.5 202.5L 53.5 188.9' fill='#000000'/><path class='atom-5' d='M 218.8 190.1Q 218.8 186.8, 220.4 185.0Q 222.0 183.1, 225.0 183.1Q 228.0 183.1, 229.6 185.0Q 231.2 186.8, 231.2 190.1Q 231.2 193.4, 229.6 195.3Q 228.0 197.1, 225.0 197.1Q 222.0 197.1, 220.4 195.3Q 218.8 193.4, 218.8 190.1M 225.0 195.6Q 227.1 195.6, 228.2 194.2Q 229.3 192.8, 229.3 190.1Q 229.3 187.4, 228.2 186.0Q 227.1 184.7, 225.0 184.7Q 222.9 184.7, 221.8 186.0Q 220.7 187.4, 220.7 190.1Q 220.7 192.8, 221.8 194.2Q 222.9 195.6, 225.0 195.6' fill='#000000'/></svg>",
        "adrenaline": "<?xml version='1.0' encoding='iso-8859-1'?><svg version='1.1' baseProfile='full'              xmlns='http://www.w3.org/2000/svg'                      xmlns:rdkit='http://www.rdkit.org/xml'                      xmlns:xlink='http://www.w3.org/1999/xlink'                  xml:space='preserve'width='307px' height='236px' viewBox='0 0 307 236'><!-- END OF HEADER --><path class='bond-0 atom-0 atom-1' d='M 168.8,204.4 L 176.0,170.1' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' /><path class='bond-1 atom-1 atom-2' d='M 171.5,157.3 L 145.7,134.2' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' /><path class='bond-2 atom-2 atom-3' d='M 145.7,134.2 L 105.2,147.5' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' /><path class='bond-2 atom-2 atom-3' d='M 140.2,129.3 L 106.7,140.3' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' /><path class='bond-3 atom-3 atom-4' d='M 105.2,147.5 L 73.4,119.1' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' /><path class='bond-4 atom-4 atom-5' d='M 73.4,119.1 L 39.8,130.2' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' /><path class='bond-5 atom-4 atom-6' d='M 73.4,119.1 L 82.1,77.4' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' /><path class='bond-5 atom-4 atom-6' d='M 80.4,116.8 L 87.6,82.3' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' /><path class='bond-6 atom-6 atom-7' d='M 82.1,77.4 L 57.2,55.1' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' /><path class='bond-7 atom-6 atom-8' d='M 82.1,77.4 L 122.6,64.0' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' /><path class='bond-8 atom-8 atom-9' d='M 122.6,64.0 L 154.4,92.4' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' /><path class='bond-8 atom-8 atom-9' d='M 121.1,71.2 L 147.4,94.7' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' /><path class='bond-9 atom-9 atom-10' d='M 154.4,92.4 L 194.9,79.1' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' /><path class='bond-10 atom-10 atom-11' d='M 194.9,79.1 L 202.0,45.0' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' /><path class='bond-11 atom-10 atom-12' d='M 194.9,79.1 L 226.8,107.5' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' /><path class='bond-12 atom-12 atom-13' d='M 226.8,107.5 L 261.3,96.1' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' /><path class='bond-13 atom-9 atom-2' d='M 154.4,92.4 L 145.7,134.2' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' /><path d='M 107.2,146.9 L 105.2,147.5 L 103.6,146.1' style='fill:none;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-miterlimit:10;stroke-opacity:1;' /><path d='M 120.6,64.7 L 122.6,64.0 L 124.2,65.4' style='fill:none;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-miterlimit:10;stroke-opacity:1;' /><path d='M 225.2,106.1 L 226.8,107.5 L 228.5,106.9' style='fill:none;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-miterlimit:10;stroke-opacity:1;' /><path class='atom-1' d='M 174.9 156.6L 178.8 163.0Q 179.2 163.6, 179.8 164.8Q 180.5 165.9, 180.5 166.0L 180.5 156.6L 182.1 156.6L 182.1 168.7L 180.5 168.7L 176.2 161.7Q 175.7 160.9, 175.2 159.9Q 174.7 159.0, 174.5 158.7L 174.5 168.7L 172.9 168.7L 172.9 156.6L 174.9 156.6' fill='#000000'/><path class='atom-1' d='M 184.5 156.6L 186.1 156.6L 186.1 161.7L 192.3 161.7L 192.3 156.6L 193.9 156.6L 193.9 168.7L 192.3 168.7L 192.3 163.1L 186.1 163.1L 186.1 168.7L 184.5 168.7L 184.5 156.6' fill='#000000'/><path class='atom-5' d='M 15.4 126.5L 17.0 126.5L 17.0 131.6L 23.2 131.6L 23.2 126.5L 24.8 126.5L 24.8 138.6L 23.2 138.6L 23.2 133.0L 17.0 133.0L 17.0 138.6L 15.4 138.6L 15.4 126.5' fill='#000000'/><path class='atom-5' d='M 27.3 132.5Q 27.3 129.6, 28.7 128.0Q 30.2 126.4, 32.8 126.4Q 35.5 126.4, 36.9 128.0Q 38.4 129.6, 38.4 132.5Q 38.4 135.4, 36.9 137.1Q 35.5 138.8, 32.8 138.8Q 30.2 138.8, 28.7 137.1Q 27.3 135.5, 27.3 132.5M 32.8 137.4Q 34.7 137.4, 35.7 136.2Q 36.7 134.9, 36.7 132.5Q 36.7 130.1, 35.7 128.9Q 34.7 127.7, 32.8 127.7Q 31.0 127.7, 30.0 128.9Q 29.0 130.1, 29.0 132.5Q 29.0 134.9, 30.0 136.2Q 31.0 137.4, 32.8 137.4' fill='#000000'/><path class='atom-7' d='M 32.8 43.0L 34.4 43.0L 34.4 48.1L 40.6 48.1L 40.6 43.0L 42.2 43.0L 42.2 55.0L 40.6 55.0L 40.6 49.5L 34.4 49.5L 34.4 55.0L 32.8 55.0L 32.8 43.0' fill='#000000'/><path class='atom-7' d='M 44.7 49.0Q 44.7 46.1, 46.1 44.4Q 47.6 42.8, 50.2 42.8Q 52.9 42.8, 54.4 44.4Q 55.8 46.1, 55.8 49.0Q 55.8 51.9, 54.3 53.6Q 52.9 55.2, 50.2 55.2Q 47.6 55.2, 46.1 53.6Q 44.7 51.9, 44.7 49.0M 50.2 53.9Q 52.1 53.9, 53.1 52.6Q 54.1 51.4, 54.1 49.0Q 54.1 46.6, 53.1 45.4Q 52.1 44.2, 50.2 44.2Q 48.4 44.2, 47.4 45.4Q 46.4 46.6, 46.4 49.0Q 46.4 51.4, 47.4 52.6Q 48.4 53.9, 50.2 53.9' fill='#000000'/><path class='atom-11' d='M 198.1 37.3Q 198.1 34.4, 199.5 32.8Q 201.0 31.2, 203.6 31.2Q 206.3 31.2, 207.8 32.8Q 209.2 34.4, 209.2 37.3Q 209.2 40.3, 207.7 42.0Q 206.3 43.6, 203.6 43.6Q 201.0 43.6, 199.5 42.0Q 198.1 40.3, 198.1 37.3M 203.6 42.2Q 205.5 42.2, 206.5 41.0Q 207.5 39.8, 207.5 37.3Q 207.5 35.0, 206.5 33.8Q 205.5 32.6, 203.6 32.6Q 201.8 32.6, 200.8 33.8Q 199.8 35.0, 199.8 37.3Q 199.8 39.8, 200.8 41.0Q 201.8 42.2, 203.6 42.2' fill='#000000'/><path class='atom-11' d='M 211.1 31.3L 212.7 31.3L 212.7 36.5L 218.9 36.5L 218.9 31.3L 220.5 31.3L 220.5 43.4L 218.9 43.4L 218.9 37.8L 212.7 37.8L 212.7 43.4L 211.1 43.4L 211.1 31.3' fill='#000000'/><path class='atom-13' d='M 264.6 88.1L 268.6 94.5Q 269.0 95.2, 269.6 96.3Q 270.2 97.4, 270.3 97.5L 270.3 88.1L 271.9 88.1L 271.9 100.2L 270.2 100.2L 266.0 93.2Q 265.5 92.4, 264.9 91.4Q 264.4 90.5, 264.3 90.2L 264.3 100.2L 262.7 100.2L 262.7 88.1L 264.6 88.1' fill='#000000'/><path class='atom-13' d='M 274.2 88.1L 275.8 88.1L 275.8 93.3L 282.0 93.3L 282.0 88.1L 283.6 88.1L 283.6 100.2L 282.0 100.2L 282.0 94.6L 275.8 94.6L 275.8 100.2L 274.2 100.2L 274.2 88.1' fill='#000000'/><path class='atom-13' d='M 286.0 99.8Q 286.3 99.0, 287.0 98.6Q 287.7 98.2, 288.7 98.2Q 289.9 98.2, 290.5 98.8Q 291.2 99.5, 291.2 100.6Q 291.2 101.8, 290.3 102.9Q 289.5 104.0, 287.7 105.3L 291.3 105.3L 291.3 106.2L 286.0 106.2L 286.0 105.5Q 287.5 104.4, 288.3 103.6Q 289.2 102.9, 289.7 102.1Q 290.1 101.4, 290.1 100.7Q 290.1 99.9, 289.7 99.5Q 289.3 99.1, 288.7 99.1Q 288.0 99.1, 287.6 99.3Q 287.2 99.6, 286.9 100.2L 286.0 99.8' fill='#000000'/></svg>",
        "morphine": "<?xml version='1.0' encoding='iso-8859-1'?><svg version='1.1' baseProfile='full'              xmlns='http://www.w3.org/2000/svg'                      xmlns:rdkit='http://www.rdkit.org/xml'                      xmlns:xlink='http://www.w3.org/1999/xlink'                  xml:space='preserve'width='283px' height='254px' viewBox='0 0 283 254'><!-- END OF HEADER --><path class='bond-0 atom-0 atom-1' d='M 14.2,232.5 L 34.1,200.6' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' /><path class='bond-1 atom-1 atom-2' d='M 34.9,184.2 L 17.1,150.6' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' /><path class='bond-2 atom-2 atom-3' d='M 17.1,150.6 L 42.3,110.5' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' /><path class='bond-3 atom-3 atom-4' d='M 42.3,110.5 L 89.8,106.7' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' /><path class='bond-4 atom-4 atom-5' d='M 89.8,106.7 L 70.1,149.8' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' /><path class='bond-5 atom-5 atom-6' d='M 70.1,149.8 L 86.6,194.1' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' /><path class='bond-6 atom-6 atom-7' d='M 86.6,194.1 L 129.6,213.8' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' /><path class='bond-7 atom-7 atom-8' d='M 129.6,213.8 L 174.0,197.3' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' /><path class='bond-8 atom-8 atom-9' d='M 174.0,197.3 L 221.3,199.0' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' /><path class='bond-8 atom-8 atom-9' d='M 178.3,190.5 L 217.5,191.9' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' /><path class='bond-9 atom-9 atom-10' d='M 221.3,199.0 L 246.4,158.9' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' /><path class='bond-10 atom-10 atom-11' d='M 246.4,158.9 L 224.3,117.1' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' /><path class='bond-10 atom-10 atom-11' d='M 238.4,158.6 L 220.0,123.9' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' /><path class='bond-11 atom-11 atom-12' d='M 224.3,117.1 L 244.1,85.5' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' /><path class='bond-12 atom-11 atom-13' d='M 224.3,117.1 L 177.0,115.4' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' /><path class='bond-13 atom-13 atom-14' d='M 177.0,115.4 L 157.0,147.2' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' /><path class='bond-14 atom-13 atom-15' d='M 177.0,115.4 L 134.1,90.2' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' /><path class='bond-14 atom-13 atom-15' d='M 167.2,117.7 L 133.5,97.9' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' /><path class='bond-15 atom-15 atom-16' d='M 134.1,90.2 L 54.9,74.7' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' /><path class='bond-16 atom-16 atom-17' d='M 54.9,74.7 L 59.3,36.2' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' /><path class='bond-17 atom-16 atom-18' d='M 54.9,74.7 L 63.0,114.9' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' /><path class='bond-18 atom-18 atom-19' d='M 63.0,114.9 L 62.0,109.6' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' /><path class='bond-18 atom-18 atom-19' d='M 56.2,116.3 L 55.1,111.0' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' /><path class='bond-19 atom-6 atom-1' d='M 86.6,194.1 L 45.9,192.6' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' /><path class='bond-20 atom-15 atom-4' d='M 134.1,90.2 L 89.8,106.7' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' /><path class='bond-21 atom-16 atom-4' d='M 54.9,74.7 L 89.8,106.7' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' /><path class='bond-22 atom-19 atom-5' d='M 62.0,109.6 L 70.1,149.8' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' /><path class='bond-23 atom-14 atom-8' d='M 156.4,164.0 L 174.0,197.3' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' /><path d='M 18.0,152.3 L 17.1,150.6 L 18.4,148.6' style='fill:none;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-miterlimit:10;stroke-opacity:1;' /><path d='M 41.0,112.5 L 42.3,110.5 L 44.6,110.3' style='fill:none;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-miterlimit:10;stroke-opacity:1;' /><path d='M 127.5,212.8 L 129.6,213.8 L 131.9,213.0' style='fill:none;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-miterlimit:10;stroke-opacity:1;' /><path d='M 218.9,198.9 L 221.3,199.0 L 222.6,197.0' style='fill:none;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-miterlimit:10;stroke-opacity:1;' /><path d='M 245.2,160.9 L 246.4,158.9 L 245.3,156.8' style='fill:none;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-miterlimit:10;stroke-opacity:1;' /><path d='M 62.6,112.9 L 63.0,114.9 L 63.0,114.6' style='fill:none;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-miterlimit:10;stroke-opacity:1;' /><path d='M 62.0,109.8 L 62.0,109.6 L 62.4,111.6' style='fill:none;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-miterlimit:10;stroke-opacity:1;' /><path class='atom-1' d='M 36.3 185.7L 40.7 192.8Q 41.2 193.5, 41.9 194.8Q 42.6 196.0, 42.6 196.1L 42.6 185.7L 44.4 185.7L 44.4 199.1L 42.5 199.1L 37.8 191.3Q 37.3 190.4, 36.7 189.4Q 36.1 188.4, 35.9 188.0L 35.9 199.1L 34.2 199.1L 34.2 185.7L 36.3 185.7' fill='#000000'/><path class='atom-12' d='M 243.3 77.0Q 243.3 73.8, 244.9 72.0Q 246.4 70.2, 249.4 70.2Q 252.4 70.2, 254.0 72.0Q 255.6 73.8, 255.6 77.0Q 255.6 80.3, 254.0 82.1Q 252.4 84.0, 249.4 84.0Q 246.5 84.0, 244.9 82.1Q 243.3 80.3, 243.3 77.0M 249.4 82.4Q 251.5 82.4, 252.6 81.1Q 253.7 79.7, 253.7 77.0Q 253.7 74.4, 252.6 73.0Q 251.5 71.7, 249.4 71.7Q 247.4 71.7, 246.3 73.0Q 245.2 74.4, 245.2 77.0Q 245.2 79.7, 246.3 81.1Q 247.4 82.4, 249.4 82.4' fill='#000000'/><path class='atom-12' d='M 257.7 70.3L 259.5 70.3L 259.5 76.0L 266.3 76.0L 266.3 70.3L 268.1 70.3L 268.1 83.7L 266.3 83.7L 266.3 77.6L 259.5 77.6L 259.5 83.7L 257.7 83.7L 257.7 70.3' fill='#000000'/><path class='atom-14' d='M 145.7 155.5Q 145.7 152.3, 147.3 150.5Q 148.9 148.7, 151.8 148.7Q 154.8 148.7, 156.4 150.5Q 158.0 152.3, 158.0 155.5Q 158.0 158.8, 156.4 160.6Q 154.8 162.5, 151.8 162.5Q 148.9 162.5, 147.3 160.6Q 145.7 158.8, 145.7 155.5M 151.8 160.9Q 153.9 160.9, 155.0 159.6Q 156.1 158.2, 156.1 155.5Q 156.1 152.9, 155.0 151.5Q 153.9 150.2, 151.8 150.2Q 149.8 150.2, 148.7 151.5Q 147.6 152.9, 147.6 155.5Q 147.6 158.2, 148.7 159.6Q 149.8 160.9, 151.8 160.9' fill='#000000'/><path class='atom-17' d='M 54.2 27.7Q 54.2 24.5, 55.8 22.7Q 57.4 20.9, 60.3 20.9Q 63.3 20.9, 64.9 22.7Q 66.5 24.5, 66.5 27.7Q 66.5 31.0, 64.9 32.8Q 63.3 34.6, 60.3 34.6Q 57.4 34.6, 55.8 32.8Q 54.2 31.0, 54.2 27.7M 60.3 33.1Q 62.4 33.1, 63.5 31.8Q 64.6 30.4, 64.6 27.7Q 64.6 25.1, 63.5 23.7Q 62.4 22.4, 60.3 22.4Q 58.3 22.4, 57.2 23.7Q 56.1 25.0, 56.1 27.7Q 56.1 30.4, 57.2 31.8Q 58.3 33.1, 60.3 33.1' fill='#000000'/><path class='atom-17' d='M 68.6 21.0L 70.4 21.0L 70.4 26.7L 77.2 26.7L 77.2 21.0L 79.1 21.0L 79.1 34.4L 77.2 34.4L 77.2 28.2L 70.4 28.2L 70.4 34.4L 68.6 34.4L 68.6 21.0' fill='#000000'/></svg>"
    };

    // Export for Matter.js path (top of file)
    try { window.__CHEMTEX_MOLECULE_DATA__ = moleculeData; } catch (e) { }

    // ------------------------------
    // Matter.js physics path (preferred)
    // ------------------------------
    function shouldUseMatter() {
        try {
            // Opt-out toggle (useful if something breaks)
            var forced = localStorage.getItem('chemtex:popupPhysics');
            if (forced === 'legacy') return false;
            if (forced === 'matter') return true;
        } catch (e) { }

        return !!(window.Matter && window.MoleculeCollider);
    }

    function shouldDebugMatter() {
        try {
            return localStorage.getItem('chemtex:matterDebug') === '1';
        } catch (e) {
            return false;
        }
    }

    function getContainerSize() {
        // Prefer real measured size so we match CSS/layout.
        var rect = container.getBoundingClientRect();
        var w = Math.max(1, Math.round(rect.width || WIDTH));
        var h = Math.max(1, Math.round(rect.height || HEIGHT));
        return { width: w, height: h };
    }

    function createDebugEl(text) {
        var el = document.createElement('div');
        el.textContent = text;
        el.style.cssText = [
            'position:absolute',
            'left:6px',
            'top:6px',
            'z-index:9999',
            'font:12px/1.3 monospace',
            'color:#0b1',
            'background:rgba(0,0,0,0.55)',
            'padding:4px 6px',
            'border-radius:4px',
            'pointer-events:none'
        ].join(';');
        container.appendChild(el);
        return el;
    }

    function startMatterAnimation() {
        if (!shouldUseMatter()) return false;

        var Engine = window.Matter.Engine;
        var World = window.Matter.World;
        var Bodies = window.Matter.Bodies;
        var Body = window.Matter.Body;
        var Composite = window.Matter.Composite;
        var Runner = window.Matter.Runner;

        var size = getContainerSize();

        // Keep the legacy constants in sync-ish with the actual container.
        WIDTH = size.width;
        HEIGHT = size.height;

        var debug = shouldDebugMatter();
        var debugEl = debug ? createDebugEl('Matter.js: init…') : null;

        // If we already started once (popup can re-run scripts in dev reload), bail.
        if (container.__chemtexMatterRunning) return true;
        container.__chemtexMatterRunning = true;

        var engine = Engine.create({
            // In "space" we don't want bodies to fall asleep.
            enableSleeping: false
        });

        // No gravity: purely "floating" bounce physics.
        // (Matter default is { x: 0, y: 1 } which feels like gravity.)
        engine.world.gravity.x = 0;
        engine.world.gravity.y = 0;
        engine.world.gravity.scale = 0;

        // More iterations helps with stacked/complex shapes.
        engine.positionIterations = 10;
        engine.velocityIterations = 8;

        var world = engine.world;

        var wallThickness = 50;
        var walls = [
            Bodies.rectangle(WIDTH / 2, -wallThickness / 2, WIDTH + wallThickness * 2, wallThickness, { isStatic: true, restitution: 0.99, friction: 0, frictionStatic: 0, frictionAir: 0 }),
            Bodies.rectangle(WIDTH / 2, HEIGHT + wallThickness / 2, WIDTH + wallThickness * 2, wallThickness, { isStatic: true, restitution: 0.99, friction: 0, frictionStatic: 0, frictionAir: 0 }),
            Bodies.rectangle(-wallThickness / 2, HEIGHT / 2, wallThickness, HEIGHT + wallThickness * 2, { isStatic: true, restitution: 0.99, friction: 0, frictionStatic: 0, frictionAir: 0 }),
            Bodies.rectangle(WIDTH + wallThickness / 2, HEIGHT / 2, wallThickness, HEIGHT + wallThickness * 2, { isStatic: true, restitution: 0.99, friction: 0, frictionStatic: 0, frictionAir: 0 })
        ];
        World.add(world, walls);

        // Create bodies + DOM nodes.
        var data = window.__CHEMTEX_MOLECULE_DATA__ || moleculeData;
        var keys = Object.keys(data);
        if (!keys.length) return false;

        // Read settings from chrome.storage (popup UI sliders) if available, otherwise fall back.
        var count = targetMoleculeCount || 20;
        var scale = globalScale || 0.4;

        // These will be updated by storage get call
        var cursorGravityStrength = -0.9;
        var initialVelocity = 3.8;

        // Check for dark mode to set molecule stroke color
        var isDark = document.body.classList.contains('dark-mode');
        // We can listen for class changes on body if we want real-time updates without reload
        // but for now, we check at build time. 
        // Better: rebuildScene checks this.

        // Track cursor position for "cursor gravity".
        // IMPORTANT: if the cursor is over popup UI widgets (sliders etc), the container won't receive
        // pointermove events. So we listen at the document/window level (capture phase) and project
        // coordinates into the moleculeContainer's local space.
        var mouse = {
            x: WIDTH / 2,
            y: HEIGHT / 2,
            hasEverMoved: false,
            lastClientX: null,
            lastClientY: null,
            overPopup: false,
            lastInsideTs: 0
        };

        function updateMouseFromClient(clientX, clientY) {
            try {
                mouse.lastClientX = clientX;
                mouse.lastClientY = clientY;
                var r = container.getBoundingClientRect();
                mouse.x = clientX - r.left;
                mouse.y = clientY - r.top;
                mouse.hasEverMoved = true;
            } catch (e) { }
        }

        function isClientInsideContainerRect(clientX, clientY) {
            try {
                var r2 = container.getBoundingClientRect();
                return clientX >= r2.left && clientX <= r2.right && clientY >= r2.top && clientY <= r2.bottom;
            } catch (e) { }
            return false;
        }

        // Capture pointer movement anywhere in the popup.
        // Use capture=true so we still receive the event even if a slider stops propagation.
        document.addEventListener('pointermove', function (ev) {
            // Determine "over popup" by geometry instead of event targets.
            // (Targets can be unreliable with pointer-events:none layers.)
            mouse.overPopup = isClientInsideContainerRect(ev.clientX, ev.clientY);
            if (!mouse.overPopup) {
                // Don't let forces "stick" to an old cursor position.
                mouse.hasEverMoved = false;
                mouse.lastInsideTs = 0;
                return;
            }
            mouse.lastInsideTs = Date.now();
            updateMouseFromClient(ev.clientX, ev.clientY);
        }, { passive: true, capture: true });

        // Fallback for environments where pointer events behave oddly.
        document.addEventListener('mousemove', function (ev) {
            mouse.overPopup = isClientInsideContainerRect(ev.clientX, ev.clientY);
            if (!mouse.overPopup) {
                mouse.hasEverMoved = false;
                mouse.lastInsideTs = 0;
                return;
            }
            mouse.lastInsideTs = Date.now();
            updateMouseFromClient(ev.clientX, ev.clientY);
        }, { passive: true, capture: true });
        var bodies = [];
        var domNodes = [];

        function clearAllDom() {
            for (var k = 0; k < domNodes.length; k++) {
                try { domNodes[k].remove(); } catch (e) { }
            }
            domNodes = [];
        }

        function clearAllBodies() {
            try {
                for (var k2 = 0; k2 < bodies.length; k2++) {
                    World.remove(world, bodies[k2]);
                }
            } catch (e) { }
            bodies = [];
        }

        function rebuildScene() {
            // Re-read settings on rebuild.
            count = Math.max(0, Math.min(targetMoleculeCount || 20, keys.length));
            scale = globalScale || 0.4;

            clearAllBodies();
            clearAllDom();

            // Check dark mode
            var isDark = document.body.classList.contains('dark-mode');
            var strokeColor = isDark ? 'rgba(255, 255, 255, 0.15)' : '#000000';
            var strokeOpacity = isDark ? '0.15' : '1';

            if (debugEl) {
                debugEl.textContent = 'Matter.js: rebuilding… (' + count + ')';
            }

            // If count is 0, we just keep an empty scene.
            for (var i = 0; i < count; i++) {
                var key = keys[i % keys.length];
                var svgString = data[key];

                var div = document.createElement('div');
                div.className = 'floating-molecule';
                div.style.cssText = 'position:absolute; pointer-events:none; will-change: transform; transform-origin:center center;';

                // Hack to inject stroke color
                if (isDark) {
                    // Replace stroke:#000000 with white/transparent
                    svgString = svgString.replace(/stroke:#000000/g, "stroke:rgba(255,255,255,0.15)");
                    svgString = svgString.replace(/fill='#000000'/g, "fill='rgba(255,255,255,0.15)'");
                }

                div.innerHTML = svgString;
                container.appendChild(div);

                // Ensure SVG scales with globalScale.
                div.style.transform = 'translate(-9999px, -9999px) scale(' + scale + ')';

                // Build collider parts from SVG viewBox space.
                // We'll scale/offset into pixels when creating the Matter body.
                var parts;
                try {
                    parts = window.MoleculeCollider.buildCompoundFromSvgString(svgString, {
                        maxParts: 6,
                        samplesPerPath: 48,
                        simplifyEps: 0.6
                    });
                } catch (e) {
                    parts = null;
                    if (debugEl) debugEl.textContent = 'Matter.js: collider error (fallback)';
                }

                var fallbackW = 80;
                var fallbackH = 60;
                try {
                    var svgEl = div.querySelector('svg');
                    if (svgEl) {
                        fallbackW = parseFloat(svgEl.getAttribute('width') || fallbackW);
                        fallbackH = parseFloat(svgEl.getAttribute('height') || fallbackH);
                    }
                } catch (e) { }

                var x = 40 + (i * 37) % Math.max(1, (WIDTH - 80));
                var y = 40 + (i * 53) % Math.max(1, (HEIGHT - 120));

                var body;
                if (parts && parts.parts && parts.parts.length) {
                    // Convert SVG-viewBox-space polygons into screen-pixel-space vertices.
                    // Our DOM SVG is rendered with the same `scale`, and then we position the parent div
                    // via translate(x,y). So the physics vertices must match that pixel size.
                    var vb = parts.viewBox;
                    var svgToPx = scale;

                    // If the inline SVG has explicit width/height that differ from viewBox, respect that.
                    // This is common for the Kekule exports and is the most reliable visual sizing source.
                    try {
                        var svgEl2 = div.querySelector('svg');
                        if (svgEl2 && vb && vb.width && vb.height) {
                            var wAttr = parseFloat(svgEl2.getAttribute('width') || '');
                            var hAttr = parseFloat(svgEl2.getAttribute('height') || '');
                            if (isFinite(wAttr) && wAttr > 0 && isFinite(hAttr) && hAttr > 0) {
                                // Map SVG units -> rendered pixels before global scale.
                                var sx = wAttr / vb.width;
                                var sy = hAttr / vb.height;
                                // Use average to keep isotropic scaling.
                                svgToPx = scale * ((sx + sy) / 2);
                            }
                        }
                    } catch (e) { }

                    // For each centered part, we need to shift it by its original centroid relative to
                    // the SVG's overall center, then scale to pixels.
                    var centerX = (vb.minX + vb.width / 2);
                    var centerY = (vb.minY + vb.height / 2);

                    var verts = [];
                    for (var pi = 0; pi < parts.parts.length; pi++) {
                        var poly = parts.parts[pi];
                        var pc = (parts.partCenters && parts.partCenters[pi]) ? parts.partCenters[pi] : { x: centerX, y: centerY };
                        var ox = (pc.x - centerX);
                        var oy = (pc.y - centerY);
                        var scaled = poly.map(function (p) {
                            return { x: (p.x + ox) * svgToPx, y: (p.y + oy) * svgToPx };
                        });
                        verts.push(scaled);
                    }

                    // Matter.js expects the first element to be the parent body.
                    body = Body.create({
                        parts: verts,
                        restitution: 0.99,
                        friction: 0,
                        frictionStatic: 0,
                        frictionAir: 0,
                        density: 0.001
                    });

                    // --- Visual calibration ---
                    // Even with careful viewBox math, some exported SVGs have internal transforms/strokes
                    // that make the *rendered* bounds differ from viewBox bounds. That causes early
                    // "telepathic" collisions (physics body larger than appearance).
                    //
                    // We measure the actual SVG geometry bounds (in SVG units) and compare to Matter bounds
                    // (in px) to compute a correction factor.
                    try {
                        var svgCal = div.querySelector('svg');
                        if (svgCal && typeof svgCal.getBBox === 'function') {
                            var bb = svgCal.getBBox(); // in SVG user units
                            if (bb && isFinite(bb.width) && isFinite(bb.height) && bb.width > 0 && bb.height > 0) {
                                // Approximate rendered pixel size of the SVG shape.
                                var measW = bb.width * svgToPx;
                                var measH = bb.height * svgToPx;

                                var physW = body.bounds.max.x - body.bounds.min.x;
                                var physH = body.bounds.max.y - body.bounds.min.y;

                                if (isFinite(physW) && isFinite(physH) && physW > 1 && physH > 1) {
                                    // If physics body is bigger, scale it down toward the measured bounds.
                                    // We clamp to avoid crazy corrections.
                                    var kx = measW / physW;
                                    var ky = measH / physH;
                                    var k = (kx + ky) / 2;
                                    k = clamp(k, 0.6, 1.1);
                                    if (Math.abs(1 - k) > 0.02) {
                                        Body.scale(body, k, k);
                                    }
                                }
                            }
                        }
                    } catch (e) { }
                } else {
                    body = Bodies.rectangle(0, 0, fallbackW * scale, fallbackH * scale, {
                        restitution: 0.99,
                        friction: 0,
                        frictionStatic: 0,
                        frictionAir: 0,
                        density: 0.001
                    });
                }

                Body.setPosition(body, { x: x, y: y });
                var v0 = (typeof initialVelocity === 'number') ? initialVelocity : 3.8;
                var angle0 = Math.random() * Math.PI * 2;
                Body.setVelocity(body, { x: Math.cos(angle0) * v0, y: Math.sin(angle0) * v0 });
                Body.setAngularVelocity(body, (Math.random() - 0.5) * 0.08);

                // Ensure sleeping stays off per-body too.
                body.sleepThreshold = Infinity;

                bodies.push(body);
                domNodes.push(div);
            }

            if (bodies.length) {
                World.add(world, bodies);
            }
        }

        function readSettingsAndThen(cb) {
            function done() {
                // Apply the values we just fetched.
                count = Math.max(0, Math.min(targetMoleculeCount || 20, keys.length));
                scale = globalScale || 0.4;
                if (typeof cb === 'function') cb();
            }

            try {
                if (window.chrome && chrome.storage && chrome.storage.sync) {
                    chrome.storage.sync.get({
                        moleculeCount: targetMoleculeCount || 20,
                        moleculeScale: globalScale || 0.4,
                        cursorGravityStrength: -0.9,
                        initialVelocity: initialVelocity
                    }, function (s) {
                        try {
                            if (typeof s.moleculeCount === 'number') targetMoleculeCount = s.moleculeCount;
                            if (typeof s.moleculeScale === 'number') globalScale = s.moleculeScale;
                            if (typeof s.cursorGravityStrength === 'number') cursorGravityStrength = s.cursorGravityStrength;
                            if (typeof s.initialVelocity === 'number') initialVelocity = s.initialVelocity;
                        } catch (e2) { }
                        done();
                    });
                    return;
                }
            } catch (e) { }

            // No storage available; proceed with in-memory defaults.
            done();
        }

        // Initial build: wait for persisted slider settings so the first frame matches saved values.
        // Start animation INSIDE the callback to ensure settings are loaded before first tick
        readSettingsAndThen(function () {
            rebuildScene();
            // Start animation only after settings are loaded
            if (!animationId) {
                // Read settings first
                chrome.storage.sync.get({
                    moleculeCount: 20,
                    moleculeScale: 0.4,
                    cursorGravityStrength: -0.9,
                    initialVelocity: 3.8
                }, function (settings) {
                    targetMoleculeCount = settings.moleculeCount;
                    globalScale = settings.moleculeScale;
                    cursorGravityStrength = settings.cursorGravityStrength;
                    initialVelocity = settings.initialVelocity;

                    rebuildScene();
                    requestAnimationFrame(tick);
                });

                // Listen for storage changes to update physics live
                chrome.storage.onChanged.addListener(function (changes, area) {
                    if (area === 'sync') {
                        if (changes.cursorGravityStrength) cursorGravityStrength = changes.cursorGravityStrength.newValue;
                        if (changes.initialVelocity) initialVelocity = changes.initialVelocity.newValue;
                        if (changes.moleculeCount || changes.moleculeScale) {
                            targetMoleculeCount = changes.moleculeCount ? changes.moleculeCount.newValue : targetMoleculeCount;
                            globalScale = changes.moleculeScale ? changes.moleculeScale.newValue : globalScale;
                            rebuildScene();
                        }
                    }
                });

                return true;
            }

            // Expose update function for popup.js to call
            window.updateMoleculeAnimation = function () {
                chrome.storage.sync.get({
                    moleculeCount: 20,
                    moleculeScale: 0.4,
                    cursorGravityStrength: -0.9,
                    initialVelocity: 3.8
                }, function (settings) {
                    targetMoleculeCount = settings.moleculeCount;
                    globalScale = settings.moleculeScale;
                    cursorGravityStrength = settings.cursorGravityStrength;
                    initialVelocity = settings.initialVelocity;
                    rebuildScene();
                });
            };
        });

        // Optional: visualize colliders (rough outlines) using an overlay canvas.
        var overlay = null;
        var overlayCtx = null;
        if (debug) {
            overlay = document.createElement('canvas');
            overlay.width = WIDTH;
            overlay.height = HEIGHT;
            overlay.style.cssText = 'position:absolute;left:0;top:0;z-index:9998;pointer-events:none;';
            container.appendChild(overlay);
            overlayCtx = overlay.getContext('2d');
        }

        var runner = Runner.create({
            isFixed: true,
            delta: 1000 / 60
        });
        Runner.run(runner, engine);

        function clamp(v, min, max) {
            return Math.max(min, Math.min(max, v));
        }

        function enforceMinSpeed(b, minSpeed) {
            // Prevent numerical solver losses from making bodies appear to "cool down".
            // We preserve direction; only boost magnitude when it's too small.
            var vx = b.velocity.x;
            var vy = b.velocity.y;
            var speedSq = vx * vx + vy * vy;
            var minSq = minSpeed * minSpeed;
            if (!isFinite(speedSq)) return;

            if (speedSq < 1e-6) {
                // Random nudge.
                var angle = Math.random() * Math.PI * 2;
                Body.setVelocity(b, { x: Math.cos(angle) * minSpeed, y: Math.sin(angle) * minSpeed });
                return;
            }
            if (speedSq < minSq) {
                var speed = Math.sqrt(speedSq);
                var s = minSpeed / speed;
                Body.setVelocity(b, { x: vx * s, y: vy * s });
            }
        }

        function tick() {
            // Cursor gravity: attraction force toward the cursor.
            // We intentionally use a shaped response so very small slider values are gentle.
            // (User expectation: 0.1 should be *barely* noticeable, not a tractor beam.)
            // Only apply when pointer is over the popup (or very shortly after, as a safety fallback).
            // The fallback covers rare cases where pointerout isn't delivered during popup close.
            // Only apply cursor forces if we have a *recent* in-popup cursor update.
            // (Prevents attraction to a stale/last-known cursor spot after leaving the popup.)
            var recentlyInside = (Date.now() - (mouse.lastInsideTs || 0)) < 80;
            if (cursorGravityStrength && mouse.hasEverMoved && mouse.overPopup && recentlyInside) {
                // Map slider -> physics strength.
                //  - 0.0 = off
                //  - 0.1 stays tiny
                //  - -4..4 supported (negative = repulsion)
                var raw = cursorGravityStrength;
                var sign = raw < 0 ? -1 : 1;
                var strength = Math.abs(raw);
                strength = Math.pow(strength / 4, 2.2) * 4; // shaped curve for 0..4 range

                for (var cg = 0; cg < bodies.length; cg++) {
                    var bb = bodies[cg];
                    var dx = mouse.x - bb.position.x;
                    var dy = mouse.y - bb.position.y;

                    // Distance falloff: strong near the cursor, weaker far away.
                    var distSq = dx * dx + dy * dy + 80 * 80;
                    var inv = 1 / Math.sqrt(distSq);

                    // Force magnitude is scaled by mass so large compounds respond similarly to small ones.
                    // Reduced constant (previous was too strong even at 0.1).
                    var base = 0.00055 * strength * bb.mass;
                    var fx = dx * inv * base * sign;
                    var fy = dy * inv * base * sign;

                    // Safety cap to prevent "snap" when extremely close to the cursor.
                    var maxStepForce = 0.0012 * bb.mass; // per physics step
                    fx = clamp(fx, -maxStepForce, maxStepForce);
                    fy = clamp(fy, -maxStepForce, maxStepForce);

                    Body.applyForce(bb, bb.position, { x: fx, y: fy });
                }
            }

            // Render DOM transforms from physics bodies.
            for (var j = 0; j < bodies.length; j++) {
                var b = bodies[j];
                var node = domNodes[j];

                // Velocity clamp to reduce tunneling for fast movers.
                // (Do clamp first, then enforce minimum speed so the floor still holds.)
                Body.setVelocity(b, {
                    x: clamp(b.velocity.x, -MAX_SPEED, MAX_SPEED),
                    y: clamp(b.velocity.y, -MAX_SPEED, MAX_SPEED)
                });

                // In space mode, keep a visible baseline kinetic energy.
                // You can think of this like a tiny perpetual micro-thruster.
                enforceMinSpeed(b, 0.65);

                // Translate so the SVG is centered on the body.
                var tx = b.position.x;
                var ty = b.position.y;
                var rot = b.angle;
                node.style.transform = 'translate(' + tx + 'px,' + ty + 'px) translate(-50%,-50%) rotate(' + rot + 'rad) scale(' + scale + ')';
            }

            if (overlayCtx && overlay) {
                overlayCtx.clearRect(0, 0, overlay.width, overlay.height);
                overlayCtx.save();
                overlayCtx.strokeStyle = 'rgba(0,255,80,0.55)';
                overlayCtx.lineWidth = 1;
                var allBodies = Composite.allBodies(world);
                for (var bi = 0; bi < allBodies.length; bi++) {
                    var bb = allBodies[bi];
                    if (bb.isStatic) continue;
                    var verts = bb.vertices;
                    if (!verts || !verts.length) continue;
                    overlayCtx.beginPath();
                    overlayCtx.moveTo(verts[0].x, verts[0].y);
                    for (var vi = 1; vi < verts.length; vi++) overlayCtx.lineTo(verts[vi].x, verts[vi].y);
                    overlayCtx.closePath();
                    overlayCtx.stroke();
                }
                overlayCtx.restore();
            }

            if (debugEl) {
                debugEl.textContent = 'Matter.js: running | count=' + (targetMoleculeCount || 0) +
                    ' | scale=' + (globalScale || 0) +
                    ' | v0=' + (initialVelocity || 0) +
                    ' | cursorG=' + (cursorGravityStrength || 0) +
                    ' | ptr=' + (mouse.hasEverMoved ? 'on' : 'off') +
                    ' | over=' + (mouse.overPopup ? 'yes' : 'no') +
                    ' | bodies=' + bodies.length;
            }

            animationId = requestAnimationFrame(tick);
        }

        // Let popup.js trigger updates from sliders. We implement a rebuild rather than no-op.
        window.updateMoleculeAnimation = function () {
            try {
                if (window.chrome && chrome.storage && chrome.storage.sync) {
                    chrome.storage.sync.get({ moleculeCount: targetMoleculeCount || 20, moleculeScale: globalScale || 0.4, cursorGravityStrength: -0.9, initialVelocity: initialVelocity }, function (s) {
                        if (typeof s.moleculeCount === 'number') targetMoleculeCount = s.moleculeCount;
                        if (typeof s.moleculeScale === 'number') globalScale = s.moleculeScale;
                        if (typeof s.cursorGravityStrength === 'number') cursorGravityStrength = s.cursorGravityStrength;
                        if (typeof s.initialVelocity === 'number') initialVelocity = s.initialVelocity;
                        rebuildScene();
                    });
                } else {
                    rebuildScene();
                }
            } catch (e) {
                rebuildScene();
            }
        };

        // Animation is started inside readSettingsAndThen callback, not here
        return true;

    }

    // If the Matter.js engine is available, start it and bypass the legacy physics.
    // (Legacy remains as a fallback for safety.)
    if (startMatterAnimation()) {
        return;
    }

    // Global state
    var molecules = [];
    var animationId = null;

    // Global settings - must match popup.js defaults
    var globalScale = 0.4;
    var targetMoleculeCount = 20;

    // Molecule class with improved physics
    function Molecule(type, svgString, x, y) {
        this.type = type;
        this.element = document.createElement('div');
        this.element.className = 'floating-molecule';
        this.element.style.cssText = 'position: absolute; pointer-events: none; will-change: transform; transform-origin: center center;';
        this.element.innerHTML = svgString;
        container.appendChild(this.element);

        // Extract dimensions from SVG
        var svg = this.element.querySelector('svg');
        var width = 60;
        var height = 60;

        if (svg) {
            var wAttr = svg.getAttribute('width');
            var hAttr = svg.getAttribute('height');
            if (wAttr) width = parseFloat(wAttr);
            if (hAttr) height = parseFloat(hAttr);

            // Apply global scale
            width *= globalScale;
            height *= globalScale;

            svg.style.width = width + 'px';
            svg.style.height = height + 'px';
        }

        this.width = width;
        this.height = height;
        this.x = x;
        this.y = y;

        // Calculate collision radius - use smaller of width/height and shrink to fit actual content
        // SVGs typically have ~20-30% whitespace around the molecule
        var CONTENT_SHRINK = 0.45; // Molecules typically fill about 45-50% of SVG bounding box
        this.radius = Math.min(width, height) * 0.5 * CONTENT_SHRINK;

        // Initialize velocity with random direction
        var angle = Math.random() * Math.PI * 2;
        var speed = 0.5 + Math.random() * 0.6;
        this.vx = Math.cos(angle) * speed;
        this.vy = Math.sin(angle) * speed;

        this.rotation = Math.random() * 360;
        this.rotationSpeed = (Math.random() - 0.5) * 0.6;
        this.opacity = 0.5 + Math.random() * 0.3;
        this.element.style.opacity = this.opacity;
        this.mass = Math.max(1, (this.width * this.height) / 1000);

        this.updatePosition();
    }

    Molecule.prototype.updatePosition = function () {
        this.element.style.transform = 'translate(' + this.x + 'px, ' + this.y + 'px) rotate(' + this.rotation + 'deg)';
    };

    // Clamp velocity to prevent tunneling
    Molecule.prototype.clampVelocity = function () {
        var speed = Math.sqrt(this.vx * this.vx + this.vy * this.vy);
        if (speed > MAX_SPEED) {
            var scale = MAX_SPEED / speed;
            this.vx *= scale;
            this.vy *= scale;
        }
    };

    // Update position for a single sub-step
    Molecule.prototype.updatePhysics = function (dt) {
        this.x += this.vx * dt;
        this.y += this.vy * dt;
        this.rotation += this.rotationSpeed * dt;

        // Wall collisions - use the visual bounds (full width/height)
        var margin = 2;

        if (this.x < margin) {
            this.x = margin;
            this.vx = Math.abs(this.vx) * WALL_RESTITUTION;
            this.rotationSpeed *= -0.5;
        } else if (this.x > WIDTH - this.width - margin) {
            this.x = WIDTH - this.width - margin;
            this.vx = -Math.abs(this.vx) * WALL_RESTITUTION;
            this.rotationSpeed *= -0.5;
        }

        if (this.y < margin) {
            this.y = margin;
            this.vy = Math.abs(this.vy) * WALL_RESTITUTION;
            this.rotationSpeed *= -0.5;
        } else if (this.y > HEIGHT - this.height - margin) {
            this.y = HEIGHT - this.height - margin;
            this.vy = -Math.abs(this.vy) * WALL_RESTITUTION;
            this.rotationSpeed *= -0.5;
        }
    };

    // Get center position
    Molecule.prototype.getCenterX = function () { return this.x + this.width / 2; };
    Molecule.prototype.getCenterY = function () { return this.y + this.height / 2; };

    // Calculate distance between molecule centers
    Molecule.prototype.distanceTo = function (other) {
        var dx = this.getCenterX() - other.getCenterX();
        var dy = this.getCenterY() - other.getCenterY();
        return Math.sqrt(dx * dx + dy * dy);
    };

    // Circle-based collision detection - uses the tighter radius
    Molecule.prototype.collidesWith = function (other) {
        var dist = this.distanceTo(other);
        var minDist = this.radius + other.radius;
        return dist < minDist;
    };

    // Circle-based collision handling
    Molecule.prototype.handleCollision = function (other) {
        var dx = this.getCenterX() - other.getCenterX();
        var dy = this.getCenterY() - other.getCenterY();
        var dist = Math.sqrt(dx * dx + dy * dy);

        // Avoid division by zero
        if (dist < 0.001) {
            dx = 1;
            dy = 0;
            dist = 1;
        }

        // Collision normal (unit vector from other to this)
        var nx = dx / dist;
        var ny = dy / dist;

        // Calculate minimum separation distance
        var minDist = this.radius + other.radius;
        var overlap = minDist - dist;

        // Separate the molecules so they no longer overlap
        if (overlap > 0) {
            var separation = overlap * 0.5 + 0.5; // Add small buffer
            this.x += nx * separation;
            this.y += ny * separation;
            other.x -= nx * separation;
            other.y -= ny * separation;
        }

        // Calculate relative velocity along collision normal
        var dvx = this.vx - other.vx;
        var dvy = this.vy - other.vy;
        var dvn = dvx * nx + dvy * ny;

        // Only resolve if molecules are approaching
        if (dvn > 0) return;

        // Calculate impulse using conservation of momentum
        var totalMass = this.mass + other.mass;
        var impulse = -(1 + RESTITUTION) * dvn / totalMass;

        // Apply impulse to velocities
        this.vx += impulse * other.mass * nx;
        this.vy += impulse * other.mass * ny;
        other.vx -= impulse * this.mass * nx;
        other.vy -= impulse * this.mass * ny;

        // Add rotation effect from collision
        this.rotationSpeed += (Math.random() - 0.5) * 0.5;
        other.rotationSpeed += (Math.random() - 0.5) * 0.5;

        // Clamp velocities to prevent runaway speeds
        this.clampVelocity();
        other.clampVelocity();
    };

    // Main animation loop with sub-stepping
    function animate() {
        // Sub-stepping loop to prevent tunneling
        var dt = 1.0 / SUB_STEPS;

        for (var step = 0; step < SUB_STEPS; step++) {
            // Update physics for each molecule
            for (var i = 0; i < molecules.length; i++) {
                molecules[i].updatePhysics(dt);
            }

            // Check and resolve collisions
            for (var j = 0; j < molecules.length; j++) {
                for (var k = j + 1; k < molecules.length; k++) {
                    if (molecules[j].collidesWith(molecules[k])) {
                        molecules[j].handleCollision(molecules[k]);
                    }
                }
            }
        }

        // Update visual positions
        for (var m = 0; m < molecules.length; m++) {
            molecules[m].updatePosition();
        }

        animationId = requestAnimationFrame(animate);
    }

    function startAnimation() {
        if (animationId) cancelAnimationFrame(animationId);

        container.innerHTML = '';
        molecules = [];

        chrome.storage.sync.get({
            moleculeCount: 20,
            moleculeScale: 0.6,
            cursorGravityStrength: -0.9,
            initialVelocity: 3.0
        }, function (settings) {
            targetMoleculeCount = settings.moleculeCount;
            globalScale = settings.moleculeScale;
            // Update global physics vars if they exist (defined in previous versions of this file or assumed)
            // Since this file seems to have multiple versions of logic, ensure we set vars used by updatePhysic
            if (typeof cursorGravityStrength !== 'undefined') cursorGravityStrength = settings.cursorGravityStrength;
            if (typeof MAX_SPEED !== 'undefined') MAX_SPEED = settings.initialVelocity;

            // Always use mol2chemfig SVGs
            var svgs = moleculeData;
            var types = Object.keys(svgs);

            // Shuffle molecule types
            types.sort(function () { return Math.random() - 0.5; });

            // Create molecules
            var count = Math.min(targetMoleculeCount, 40); // Cap at 40 for performance

            var attempts = 0;
            var maxAttempts = count * 15;

            while (molecules.length < count && attempts < maxAttempts) {
                attempts++;
                var type = types[molecules.length % types.length];
                var svgString = svgs[type];

                if (!svgString) continue;

                // Calculate temp size for positioning
                var tempDiv = document.createElement('div');
                tempDiv.innerHTML = svgString;
                var tempSvg = tempDiv.querySelector('svg');
                var w = 60, h = 60;
                if (tempSvg) {
                    var wAttr = tempSvg.getAttribute('width');
                    var hAttr = tempSvg.getAttribute('height');
                    if (wAttr) w = parseFloat(wAttr) * globalScale;
                    if (hAttr) h = parseFloat(hAttr) * globalScale;
                }

                var x = 10 + Math.random() * (WIDTH - w - 20);
                var y = 10 + Math.random() * (HEIGHT - h - 20);

                // Calculate the collision radius for this potential molecule
                var CONTENT_SHRINK = 0.45;
                var newRadius = Math.min(w, h) * 0.5 * CONTENT_SHRINK;
                var newCenterX = x + w / 2;
                var newCenterY = y + h / 2;

                // Check for overlap with existing molecules using circle collision
                var valid = true;
                for (var j = 0; j < molecules.length; j++) {
                    var other = molecules[j];
                    var dx = newCenterX - other.getCenterX();
                    var dy = newCenterY - other.getCenterY();
                    var dist = Math.sqrt(dx * dx + dy * dy);
                    var minDist = newRadius + other.radius + 5; // Small buffer for initial spacing

                    if (dist < minDist) {
                        valid = false;
                        break;
                    }
                }

                if (valid) {
                    molecules.push(new Molecule(type, svgString, x, y));
                }
            }

            if (molecules.length > 0) {
                animate();
            }
        });
    }

    // Expose update function
    window.updateMoleculeAnimation = startAnimation;

    // Initial start
    startAnimation();

})();
