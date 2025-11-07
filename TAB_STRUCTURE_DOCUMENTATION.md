# Tab Structure Documentation

## Page-Level Tabs vs Internal Tabs

The MoleculeViewer now has **two levels of tabs**:

### Level 1: Page-Level Tabs (Top-Most)
These switch between completely different applications running side-by-side.

```
┌────────────────────────────────────────────────────────────┐
│  📊 MoleculeViewer  │  🧬 Mol2ChemFig (LaTeX)             │  ← PAGE TABS
├────────────────────────────────────────────────────────────┤
│                                                            │
│  MoleculeViewer Page Content                              │
│  ┌─────────────────────────────┐                          │
│  │  SMILES  │  Name  ← INPUT TABS │                          │
│  │          │                 │                          │
│  │  [Input Field]          │                          │
│  │  [Examples]             │                          │
│  │  [Output]               │                          │
│  └─────────────────────────────┘                          │
│                                                            │
└────────────────────────────────────────────────────────────┘
```

### Level 2: Internal Tabs (Within MoleculeViewer Page)
These switch between input methods for MoleculeViewer only.

```
┌────────────────────────────────────────────────────────────┐
│  SMILES  │  Name                    ← INTERNAL INPUT TABS  │
├────────────────────────────────────────────────────────────┤
│                                                            │
│  Input Method 1: SMILES                                   │
│  ┌──────────────────────────────────────────────────┐    │
│  │ Enter SMILES String:                             │    │
│  │ [c1ccccc1...........................]             │    │
│  │ [Convert to SVG]                                 │    │
│  └──────────────────────────────────────────────────┘    │
│                                                            │
│  Input Method 2: Chemical Name (hidden until clicked)    │
│  ┌──────────────────────────────────────────────────┐    │
│  │ Enter Chemical Name:                             │    │
│  │ [benzene.........................]                │    │
│  │ [Lookup & Convert]                               │    │
│  └──────────────────────────────────────────────────┘    │
│                                                            │
└────────────────────────────────────────────────────────────┘
```

## HTML Structure

### Page-Level Tab System
```html
<!-- Top-level page tabs -->
<div class="page-tabs">
    <button class="page-tab-button active" onclick="switchPageTab('moleculeviewer-page')">
        📊 MoleculeViewer
    </button>
    <button class="page-tab-button" onclick="switchPageTab('mol2chemfig-page')">
        🧬 Mol2ChemFig (LaTeX)
    </button>
</div>

<!-- MoleculeViewer Page -->
<div id="moleculeviewer-page" class="page-content active">
    <!-- All MoleculeViewer content here -->
    <div class="main-content">
        <!-- Original tabs, input, output -->
    </div>
</div>

<!-- Mol2ChemFig Page -->
<div id="mol2chemfig-page" class="page-content">
    <!-- All Mol2ChemFig content here -->
    <div class="m2cf-container">
        <!-- Mol2ChemFig UI -->
    </div>
</div>
```

### MoleculeViewer Internal Tabs (Within Its Page)
```html
<div class="tabs">
    <button class="tab-button active" onclick="switchTab('smiles-tab')">SMILES</button>
    <button class="tab-button" onclick="switchTab('nomenclature-tab')">Name</button>
</div>

<!-- SMILES Input Tab -->
<div id="smiles-tab" class="tab-content active">
    <textarea id="smiles-input">...</textarea>
    <button onclick="convertSMILES()">Convert to SVG</button>
</div>

<!-- Nomenclature Input Tab -->
<div id="nomenclature-tab" class="tab-content">
    <input type="text" id="nomenclature-input">
    <button onclick="convertNomenclature()">Lookup & Convert</button>
</div>
```

## CSS Classes

### Page-Level Tab Styling
```css
.page-tabs {
    display: flex;
    gap: 0;
    margin-bottom: 30px;
    border-bottom: 3px solid #e0e0e0;
}

.page-tab-button {
    padding: 15px 30px;
    background: none;
    border: none;
    cursor: pointer;
    color: #999;
    font-weight: 700;
    font-size: 1.1em;
    border-bottom: 3px solid transparent;
    transition: all 0.3s;
}

.page-tab-button.active {
    color: #667eea;
    border-bottom-color: #667eea;
}

.page-content {
    display: none;
}

.page-content.active {
    display: block;
}
```

### MoleculeViewer Internal Tab Styling (Unchanged)
```css
.tabs {
    display: flex;
    gap: 0;
    margin-bottom: 20px;
    border-bottom: 2px solid #e0e0e0;
}

.tab-button {
    flex: 1;
    padding: 12px;
    background: none;
    border: none;
    cursor: pointer;
    color: #999;
    font-weight: 600;
    border-bottom: 3px solid transparent;
    margin: 0;
    border-radius: 0;
    transition: all 0.3s;
}

.tab-button.active {
    color: #667eea;
    border-bottom-color: #667eea;
}

.tab-content {
    display: none;
}

.tab-content.active {
    display: block;
}
```

## JavaScript Functions

### Page-Level Tab Switching
```javascript
function switchPageTab(pageId) {
    // Hide all page content
    document.querySelectorAll('.page-content').forEach(page => {
        page.classList.remove('active');
    });
    document.querySelectorAll('.page-tab-button').forEach(btn => {
        btn.classList.remove('active');
    });
    
    // Show selected page
    document.getElementById(pageId).classList.add('active');
    event.target.classList.add('active');
}
```

### MoleculeViewer Internal Tab Switching
```javascript
function switchTab(tabId) {
    // Hide all tabs
    document.querySelectorAll('.tab-content').forEach(tab => {
        tab.classList.remove('active');
    });
    document.querySelectorAll('.tab-button').forEach(btn => {
        btn.classList.remove('active');
    });
    
    // Show selected tab
    document.getElementById(tabId).classList.add('active');
    event.target.classList.add('active');
}
```

## Tab Switching Flow

### User Clicks "🧬 Mol2ChemFig" Page Tab
1. **Event:** `onclick="switchPageTab('mol2chemfig-page')"`
2. **Function:** `switchPageTab()` executes
3. **DOM Changes:**
   - Remove `.active` from all `.page-content` elements
   - Remove `.active` from all `.page-tab-button` elements
   - Add `.active` to `#mol2chemfig-page` (makes it visible)
   - Add `.active` to clicked button (highlights it)
4. **Result:** Mol2ChemFig page is now visible, MoleculeViewer page is hidden

### User Clicks "Name" Tab Within MoleculeViewer
1. **Event:** `onclick="switchTab('nomenclature-tab')"`
2. **Function:** `switchTab()` executes
3. **DOM Changes:**
   - Remove `.active` from all `.tab-content` elements
   - Remove `.active` from all `.tab-button` elements
   - Add `.active` to `#nomenclature-tab` (makes it visible)
   - Add `.active` to clicked button (highlights it)
4. **Result:** Chemical name input is now visible, SMILES input is hidden
5. **Note:** This only affects MoleculeViewer tab, Mol2ChemFig tab is unaffected

## State Management

### Per-Page State
Each page maintains its own state:

**MoleculeViewer Page:**
- `document.getElementById('smiles-input').value`
- `document.getElementById('nomenclature-input').value`
- `document.getElementById('width-input').value`
- Form checkbox states (fancy-bonds, aromatic-circles, etc.)

**Mol2ChemFig Page:**
- `document.getElementById('m2cf-smiles-input').value`
- `currentM2CFMolecule` object (stores SMILES, chemfig, options)
- Form checkbox states (m2cf-aromatic, m2cf-carbon, etc.)

### Persistence Across Tab Switches
When switching between pages:
- All input values are preserved (browser doesn't reset them)
- State variables remain in memory
- User can switch back and see exactly what they left

Example:
1. User enters `c1ccccc1` in MoleculeViewer SMILES field
2. Clicks "🧬 Mol2ChemFig" tab
3. MoleculeViewer SMILES field is hidden but still contains `c1ccccc1`
4. User generates a molecule in Mol2ChemFig
5. Clicks "📊 MoleculeViewer" tab
6. MoleculeViewer SMILES field still shows `c1ccccc1` ✓

## CSS Cascading

### Page-Level Tabs (Larger, More Prominent)
- Padding: `15px 30px` (more space)
- Font size: `1.1em` (larger)
- Font weight: `700` (bold)
- Border: `3px` thick

### Internal Tabs (Smaller, Within Content)
- Padding: `12px` (compact)
- Font size: `default` (smaller)
- Font weight: `600` (semi-bold)
- Border: `3px` thick

This visual hierarchy helps users understand that page-level tabs are more important.

## Responsive Design

### On Desktop (> 768px)
- Both tab levels display normally
- Page-level tabs take full width
- Internal tabs take full width
- Side-by-side grid layouts work

### On Mobile (< 768px)
- Tab styling remains the same
- Text size on buttons might wrap
- Grid layouts stack vertically
- All functionality still works

---

## Summary

| Aspect | Page-Level Tabs | Internal Tabs |
|--------|---|---|
| **Purpose** | Switch between applications | Switch input methods |
| **Affected** | Entire page/UI changes | Only input section |
| **Function** | `switchPageTab()` | `switchTab()` |
| **Button Class** | `.page-tab-button` | `.tab-button` |
| **Container Class** | `.page-content` | `.tab-content` |
| **Active State** | `.page-content.active` | `.tab-content.active` |
| **Backend** | Both (RDKit or Docker) | MoleculeViewer only |
| **Independence** | Completely separate | Linked within page |

Both tab systems work together to create a professional, organized interface for molecular visualization!
