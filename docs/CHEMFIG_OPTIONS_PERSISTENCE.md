# Chemfig Options Persistence - Implementation Summary

## Problem Statement
When users apply mol2chemfig options in the UI (test_m2cf_full.html), the options don't persist:
1. Options reset to defaults after hitting "search and convert" again
2. Cache links don't show after applying options
3. Users must manually reapply options for each conversion

## Solution Implemented

### 1. localStorage Persistence
**File**: `test_m2cf_full.html`

Added persistent storage for chemfig options using browser's localStorage:

```javascript
// State variable to track saved options
let savedOptions = {
    selections: [],
    angle: 0,
    indentation: 4,
    h2: 'keep'
};
```

### 2. Save Options on Apply
When users click "Apply Options", the settings are saved to localStorage:

```javascript
// In applyMoleculeOptions() function
savedOptions = {
    selections: selected,
    angle: angle,
    indentation: indentation,
    h2: h2
};
localStorage.setItem('chemfig_options', JSON.stringify(savedOptions));
console.log('💾 Saved chemfig options to localStorage:', savedOptions);
```

### 3. Load Saved Options on Page Load
Added `loadSavedChemfigOptions()` function that runs on page load:

```javascript
function loadSavedChemfigOptions() {
    const savedOptionsStr = localStorage.getItem('chemfig_options');
    if (savedOptionsStr) {
        savedOptions = JSON.parse(savedOptionsStr);

        // Apply saved checkbox selections
        OPTIONS.forEach((opt, idx) => {
            const checkbox = document.getElementById(`mol_opt_${idx}`);
            if (checkbox && savedOptions.selections.includes(opt.value)) {
                checkbox.checked = true;
            }
        });

        // Apply saved numeric values
        document.getElementById('moleculeAngle').value = savedOptions.angle;
        document.getElementById('moleculeIndent').value = savedOptions.indentation;
        document.getElementById('moleculeH2').value = savedOptions.h2;

        updateActiveOptions();
    }
}
```

Called during initialization:
```javascript
document.addEventListener('DOMContentLoaded', () => {
    checkBackendStatus();
    initializeComposers();
    populateOptions();
    hookOptionCheckboxes();
    loadSettings();
    loadSavedChemfigOptions(); // Load saved chemfig options
});
```

### 4. Auto-Apply Options on Conversion
Modified `submitMolecule()` to automatically apply saved options after generating chemfig:

```javascript
// Auto-apply saved options if they exist
if (savedOptions.selections.length > 0) {
    console.log('🔄 Auto-applying saved options after molecule generation');
    btn.innerHTML = '<span class="loading">⏳</span> Applying saved options...';
    setTimeout(async () => {
        await applyMoleculeOptions();
        btn.disabled = false;
        btn.textContent = 'Generate Chemfig';
    }, 500);
    return;
}
```

### 5. Display Cache Links
Enhanced the apply options response handler to show cache links in the output:

```javascript
// Display chemfig output with cache links if available
let outputText = result.chemfig;

// Add cache links section
if (result.pdflink || result.svglink) {
    outputText += '\n\n--- Cache Links ---\n';
    if (result.svglink) {
        const svgCacheUrl = result.svglink.startsWith('http') ? result.svglink : `${API_BASE.replace('/m2cf', '')}${result.svglink}`;
        outputText += `SVG: ${svgCacheUrl}\n`;
    }
    if (result.pdflink) {
        const pdfCacheUrl = result.pdflink.startsWith('http') ? result.pdflink : `${API_BASE.replace('/m2cf', '')}${result.pdflink}`;
        outputText += `PDF: ${pdfCacheUrl}\n`;
    }
}

document.getElementById('moleculeOutput').textContent = outputText;
```

### 6. User Feedback
Added success banner after applying options:

```javascript
// Show success message
const successDiv = document.createElement('div');
successDiv.className = 'success-banner';
successDiv.textContent = '✅ Options applied and saved! They will be used for future conversions.';
successDiv.style.marginBottom = '15px';
const errorDiv = document.getElementById('moleculeError');
errorDiv.innerHTML = '';
errorDiv.appendChild(successDiv);
errorDiv.style.display = 'block';
setTimeout(() => {
    errorDiv.style.display = 'none';
}, 4000);
```

## Backend Support

The mol2chemfig backend (port 8000) already returns the necessary data:

### /m2cf/submit endpoint
Returns base molecule with no options:
- `chemfig`: LaTeX code
- `svglink`: SVG data URI or path
- `pdflink`: PDF data URI or path
- `chem_data`: Original input
- `chem_format`: Input format (smiles/mol)

### /m2cf/apply endpoint
Returns molecule with applied options:
- `chemfig`: LaTeX code with options applied
- `svglink`: SVG with options applied
- `pdflink`: PDF with options applied

Example request:
```bash
curl -X POST http://localhost:8000/m2cf/apply \
  -H "Content-Type: application/json" \
  -d '{
    "chem_data":"CCO",
    "chem_format":"smiles",
    "selections":["-o","-c"],
    "angle":0,
    "indentation":4,
    "h2":"keep"
  }'
```

## User Experience Flow

1. **First Visit**:
   - User searches for a molecule (e.g., "aspirin")
   - Molecule is generated with default settings

2. **Configuring Options**:
   - User selects desired options (aromatic circles, show carbons, etc.)
   - User clicks "Apply Options"
   - Options are saved to localStorage
   - Cache links are displayed
   - Success message confirms options are saved

3. **Subsequent Conversions**:
   - User searches for a new molecule (e.g., "caffeine")
   - Saved options are automatically applied
   - No need to click "Apply Options" again
   - User sees molecule with their preferred settings

4. **Persistence Across Sessions**:
   - Options remain saved even after closing browser
   - Next time user opens the page, their preferences are restored

## Testing

To test the implementation:

1. **Start the backend**:
   ```bash
   cd C:\Users\Kapil\Personal\PROJECTS\Chemparser
   docker-compose up -d
   ```

2. **Open the UI**:
   - Navigate to `http://localhost:9000/test_m2cf_full.html`
   - Or open `test_m2cf_full.html` directly in browser

3. **Test scenario**:
   - Search for "aspirin"
   - Check "aromatic circles" and "show carbons" options
   - Click "Apply Options"
   - Verify success message appears
   - Verify cache links are displayed in output
   - Search for "caffeine"
   - Verify options are automatically applied
   - Refresh the page
   - Verify checkboxes remain checked
   - Search for another molecule
   - Verify options are automatically applied

## Files Modified

1. **test_m2cf_full.html** (Lines affected):
   - Lines 799-804: Added savedOptions state variable
   - Line 813: Added loadSavedChemfigOptions() call
   - Lines 1340-1348: Added localStorage save in applyMoleculeOptions()
   - Lines 1387-1442: Enhanced response handling with cache links and success message
   - Lines 1297-1307: Added auto-apply logic in submitMolecule()
   - Lines 1594-1640: Added loadSavedChemfigOptions() function

## Benefits

1. **Better UX**: Users don't need to reapply options for each molecule
2. **Persistence**: Settings survive page refreshes and browser restarts
3. **Transparency**: Cache links are clearly displayed
4. **Feedback**: Success messages confirm when options are saved
5. **Consistency**: All molecules use the same preferred settings

## Notes

- Options are stored per-browser (localStorage is browser-specific)
- Clearing browser data will reset saved options
- Options can be reset by unchecking all options and clicking "Apply Options"
- The backend returns data URIs for immediate display, but cache links can be extracted from the response
- Options apply to molecule conversions, not to the layered view (which uses instant toggling)
