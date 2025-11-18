# 🎨 VISUAL GUIDE - What You See

## Browser View (localhost:5000)

### MoleculeViewer Tab (Default - What You See When You First Open)

```
╔════════════════════════════════════════════════════════════════╗
║  🧪 MoleculeViewer  🧬 Mol2ChemFig (LaTeX)                    ║  ← PAGE TABS
║  SMILES to SVG Converter with Nomenclature Lookup & M2CF...    ║
╠════════════════════════════════════════════════════════════════╣
║                                                                ║
║  Input                              Examples                  ║
║  ┌────────────────┐                 ┌────────────────────┐   ║
║  │ SMILES │ Name  │ ← INTERNAL TAB  │ Benzene Button    │   ║
║  ├────────────────┤                 │ Acetic Acid       │   ║
║  │ [SMILES input] │                 │ Ibuprofen         │   ║
║  │ [Convert →]    │                 │ + more...          │   ║
║  │                │                 └────────────────────┘   ║
║  │ Width: 600     │                                          ║
║  │ Height: 500    │                                          ║
║  │                │                                          ║
║  │ 🎨 Options     │                                          ║
║  │ ☑ Fancy Bonds  │                                          ║
║  │ ☑ Aromatic...  │                                          ║
║  └────────────────┘                                          ║
║                                                                ║
║  Molecule Visualization                                       ║
║  ┌────────────────────────────────────────────────────────┐   ║
║  │  [Enter SMILES and click Convert]                     │   ║
║  └────────────────────────────────────────────────────────┘   ║
║                                                                ║
║  Molecular Information                                        ║
║  ┌─────────────┬──────────────┐                              ║
║  │ Mol Weight  │ Formula      │                              ║
║  │ Atoms       │ Bonds        │                              ║
║  └─────────────┴──────────────┘                              ║
║                                                                ║
╚════════════════════════════════════════════════════════════════╝
```

---

## Mol2ChemFig Tab (Click the 🧬 Tab)

### After Clicking "🧬 Mol2ChemFig (LaTeX)" Button

```
╔════════════════════════════════════════════════════════════════╗
║  🧪 MoleculeViewer  🧬 Mol2ChemFig (LaTeX)                    ║  ← NOW THIS TAB IS HIGHLIGHTED
║  SMILES to SVG Converter with Nomenclature Lookup & M2CF...    ║
╠════════════════════════════════════════════════════════════════╣
║                                                                ║
║  🧬 Mol2ChemFig Generator                                      ║
║  ┌───────────────────────────────────────────────────────┐   ║
║  │ Enter SMILES or Chemical Structure:                   │   ║
║  │ ┌──────────────────────────────────────────────────┐  │   ║
║  │ │ [c1ccccc1.....................]                 │  │   ║
║  │ └──────────────────────────────────────────────────┘  │   ║
║  │                                                       │   ║
║  │ Rendering Options:                                  │   ║
║  │ ┌────────────┐ ┌────────────┐ ┌────────────┐        │   ║
║  │ │ ☑ 🔵       │ │ ☐ 🅲       │ │ ☐ 🅼      │        │   ║
║  │ │ Aromatic   │ │ Show       │ │ Show      │        │   ║
║  │ │ Circles    │ │ Carbons    │ │ Methyls   │        │   ║
║  │ └────────────┘ └────────────┘ └────────────┘        │   ║
║  │ ┌────────────┐ ┌────────────┐ ┌────────────┐        │   ║
║  │ │ ☐ 🔢       │ │ ☑ ✨       │ │ ☐ 📦      │        │   ║
║  │ │ Atom       │ │ Fancy      │ │ Compact   │        │   ║
║  │ │ Numbers    │ │ Bonds      │ │ View      │        │   ║
║  │ └────────────┘ └────────────┘ └────────────┘        │   ║
║  │                                                       │   ║
║  │ [🎨 Generate SVG]  [📚 Example: Benzene]            │   ║
║  └───────────────────────────────────────────────────────┘   ║
║                                                                ║
║  Molecule Rendering (LaTeX ChemFig Format)                    ║
║  ┌───────────────────────────────────────────────────────┐   ║
║  │                                                       │   ║
║  │          benzene structure with circles              │   ║
║  │              (SVG preview)                           │   ║
║  │                                                       │   ║
║  │          [⬇️ Download SVG]                           │   ║
║  │          💾 Downloads as .svg file                   │   ║
║  │                                                       │   ║
║  └───────────────────────────────────────────────────────┘   ║
║                                                                ║
║  ChemFig Code (for LaTeX)                                     ║
║  ┌───────────────────────────────────────────────────────┐   ║
║  │ \chemfig{*6(=(-...)(-...)=(-...)(-...)=(-...)(-...))} │   ║
║  │                                                       │   ║
║  │ Use with \chemfig{...} in your LaTeX document       │   ║
║  └───────────────────────────────────────────────────────┘   ║
║                                                                ║
╚════════════════════════════════════════════════════════════════╝
```

---

## Tab Switching Animation (What Happens)

### Timeline of Events

```
CLICK "🧬 Mol2ChemFig" TAB
         ↓
    switchPageTab('mol2chemfig-page')
         ↓
    Hide .moleculeviewer-page (display: none)
    Show .mol2chemfig-page (display: block)
    Highlight 🧬 button (color: #667eea, border-bottom: 3px)
    Unhighlight 📊 button (color: #999, border-bottom: none)
         ↓
    [INSTANT - CSS change, no reload]
         ↓
    User sees new interface
```

---

## What's Different - Side by Side

### MoleculeViewer Tab

```
INPUT FORM
┌──────────────────────────────┐
│ SMILES | Name                │  ← Choose input method
├──────────────────────────────┤
│ Enter SMILES String:         │
│ [c1ccccc1 or CCO...]         │
│                              │
│ OR                           │
│                              │
│ Enter Chemical Name:         │
│ [benzene, aspirin, etc...]   │
│                              │
│ Width: 600px                 │
│ Height: 500px                │
│                              │
│ 🎨 Visualization Options     │
│ ☑ Fancy Bonds                │
│ ☑ Aromatic Circles           │
│ ☐ Show Carbons               │
│ ☐ Show Methyls               │
│ ☐ Atom Numbers               │
│ ☐ Recalc Coordinates         │
│                              │
│ Hydrogens: [keep ▼]          │
│ Flip: [No Flip ▼]            │
│ Rotate: [0° ▼]               │
│                              │
│ [Convert to SVG]             │
└──────────────────────────────┘

OUTPUT
┌──────────────────────────────┐
│ SVG rendered by RDKit        │
│                              │
│ [Download SVG]               │
│ Cache URL (24 hours)         │
│                              │
│ Molecular Weight: 78.11      │
│ Formula: C6H6                │
│ Atoms: 12   Bonds: 6         │
└──────────────────────────────┘
```

### Mol2ChemFig Tab

```
INPUT FORM
┌──────────────────────────────┐
│ Rendering Options            │
├──────────────────────────────┤
│ Enter SMILES:                │
│ [c1ccccc1.................]   │
│                              │
│ Rendering Options (all at once)
│ ☑ 🔵 Aromatic Circles (-o)  │
│ ☐ 🅲 Show Carbons (-c)      │
│ ☐ 🅼 Show Methyls (-m)      │
│ ☐ 🔢 Atom Numbers (-n)      │
│ ☑ ✨ Fancy Bonds (-f)        │
│ ☐ 📦 Compact View (-z)      │
│                              │
│ [🎨 Generate SVG]            │
│ [📚 Example: Benzene]        │
└──────────────────────────────┘

OUTPUT
┌──────────────────────────────┐
│ SVG rendered by LaTeX        │
│ (with ChemFig code applied)  │
│                              │
│ [⬇️ Download SVG]            │
│                              │
│ ChemFig Code:                │
│ \chemfig{*6(...)}            │
│ (for LaTeX documents)        │
└──────────────────────────────┘
```

---

## Color & Design Scheme

### MoleculeViewer Tab
- **Header:** Purple gradient (`#667eea` to `#764ba2`)
- **Background:** White cards with shadows
- **Active Tab:** Purple text + border
- **Buttons:** Purple gradient with hover effect

### Mol2ChemFig Tab (Same Design Language)
- **Header:** Same purple gradient
- **Background:** White cards with shadows
- **Options:** Rounded buttons with emojis
- **Buttons:** Purple primary, gray secondary

---

## State Example

### Scenario: User Tests Both Backends

```
USER ACTION                          SCREEN STATE
────────────────────────────────────────────────────────────
1. Opens http://localhost:5000/      MoleculeViewer tab (active)
                                     [SMILES field empty]
                                     [Mol2ChemFig tab hidden]

2. Types "c1ccccc1" in SMILES         [SMILES field: "c1ccccc1"]
   (still on MoleculeViewer tab)      [Output: empty]

3. Clicks "Convert to SVG"            [Output: RDKit benzene SVG]
                                     [Info: RDKit properties]

4. Clicks "🧬 Mol2ChemFig" tab        MoleculeViewer tab (hidden)
                                     Mol2ChemFig tab (now active)
                                     [SMILES field in m2cf: empty]
                                     [Output: empty]

5. Types "c1ccccc1" in mol2chemfig    [m2cf SMILES field: "c1ccccc1"]
   SMILES field                       [MoleculeViewer preserved]

6. Checks "Aromatic Circles"          [m2cf options: aromatic checked]

7. Clicks "Generate SVG"              [Output: LaTeX benzene SVG]
                                     [ChemFig code: shown]

8. Clicks "📊 MoleculeViewer" tab     MoleculeViewer tab (now active)
                                     [SMILES field: still "c1ccccc1"]
                                     [Output: still RDKit benzene]
                                     [m2cf tab hidden but state saved]

9. Changes SMILES to "CCO"            [SMILES field: "CCO"]

10. Converts again                    [Output: RDKit ethanol SVG]

11. Clicks "🧬 Mol2ChemFig" tab       [m2cf SMILES field: still "c1ccccc1"]
                                     [Options: aromatic still checked]
                                     [Previous output: still there]
```

**Key observation:** State is completely independent and preserved!

---

## Error Handling

### If Docker Backend is Down

```
User clicks "Generate SVG" in Mol2ChemFig tab
         ↓
generateM2CF() runs
         ↓
fetch() to localhost:8000 fails
         ↓
catch (error)
         ↓
showM2CFMessage('Error: Connection refused', 'error')
         ↓
User sees red error box: "Error: Connection refused"
```

### If SMILES is Invalid

```
User enters: "invalid_smiles_string"
         ↓
User clicks "Generate SVG"
         ↓
Backend returns: { error: "Invalid SMILES" }
         ↓
if (data.error) showM2CFMessage(`Error: ${data.error}`, 'error')
         ↓
User sees red error box: "Error: Invalid SMILES"
```

---

## Download File Example

### User Downloads SVG from Mol2ChemFig

```
User clicks "[⬇️ Download SVG]"
         ↓
downloadM2CFSVG() called
         ↓
base64 SVG decoded
         ↓
Browser file dialog opens
         ↓
File saved as: "molecule_1731XXX.svg"
         ↓
File appears in Downloads folder
         ↓
User can open in:
├─ Browser (for viewing)
├─ Inkscape (for editing)
├─ LaTeX document (as image)
└─ Other SVG viewers
```

---

## Summary: Visual Differences

| Aspect | MoleculeViewer | Mol2ChemFig |
|--------|---|---|
| **Layout** | Grid: Input+Examples | Stack: Input, Output |
| **Input** | 2 tab choices | Single SMILES |
| **Options** | Dropdown menus + checkboxes | 6 checkboxes only |
| **Output** | SVG + Properties | SVG + ChemFig code |
| **Color of active buttons** | Purple gradient | Purple gradient |
| **Information shown** | Molecular weight, formula, etc. | ChemFig code for LaTeX |

---

This is what you'll see when browsing to `http://localhost:5000/`!
