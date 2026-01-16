// DIAGNOSTIC SCRIPT - Paste this in browser console to find the red background source

console.log('%c🔍 DIAGNOSTIC: Finding Red Background Source', 'background: red; color: white; font-size: 16px; padding: 10px;');

// Find all molecule images
const moleculeImages = document.querySelectorAll('img[data-molecule-viewer], .molecule-diagram, .molecule-container');

console.log(`Found ${moleculeImages.length} molecule elements`);

moleculeImages.forEach((el, index) => {
    const computed = window.getComputedStyle(el);
    const parent = el.parentElement;
    const parentComputed = parent ? window.getComputedStyle(parent) : null;

    console.log(`\n=== Element ${index + 1} ===`);
    console.log('Element:', el);
    console.log('Tag:', el.tagName);
    console.log('Classes:', el.className);

    // Check background
    if (computed.backgroundColor && computed.backgroundColor !== 'rgba(0, 0, 0, 0)' && computed.backgroundColor !== 'transparent') {
        console.log('%c❌ BACKGROUND COLOR:', 'color: red; font-weight: bold;', computed.backgroundColor);
    } else {
        console.log('%c✅ Background: transparent', 'color: green;');
    }

    // Check filter
    if (computed.filter && computed.filter !== 'none') {
        console.log('%c❌ FILTER:', 'color: red; font-weight: bold;', computed.filter);
    } else {
        console.log('%c✅ Filter: none', 'color: green;');
    }

    // Check parent
    if (parent) {
        console.log('Parent:', parent.tagName, parent.className);
        if (parentComputed.backgroundColor && parentComputed.backgroundColor !== 'rgba(0, 0, 0, 0)' && parentComputed.backgroundColor !== 'transparent') {
            console.log('%c❌ PARENT BACKGROUND:', 'color: red; font-weight: bold;', parentComputed.backgroundColor);
        }
        if (parentComputed.filter && parentComputed.filter !== 'none') {
            console.log('%c❌ PARENT FILTER:', 'color: red; font-weight: bold;', parentComputed.filter);
        }
    }

    // Check inline styles
    if (el.style.backgroundColor) {
        console.log('%c⚠️ INLINE background:', 'color: orange;', el.style.backgroundColor);
    }
    if (el.style.filter) {
        console.log('%c⚠️ INLINE filter:', 'color: orange;', el.style.filter);
    }
    if (parent && parent.style.backgroundColor) {
        console.log('%c⚠️ PARENT INLINE background:', 'color: orange;', parent.style.backgroundColor);
    }
    if (parent && parent.style.filter) {
        console.log('%c⚠️ PARENT INLINE filter:', 'color: orange;', parent.style.filter);
    }
});

console.log('\n%c🔍 Check the output above to find what is applying the red background!', 'background: blue; color: white; font-size: 14px; padding: 8px;');
