# 🎨 Visual Guide: How Text Replacement Works

---

## 📝 What You Type in ChatGPT

```
"Let's analyze the structure of chem:acetone and compare it to chem:benzene"
```

---

## 🔄 What the Extension Does (Behind the Scenes)

### Step 1: Detection
```
Extension scans for pattern: chem:XXXXX
Found: "chem:acetone"
Found: "chem:benzene"
```

### Step 2: API Calls
```
Call 1: GET https://YOUR-APP.herokuapp.com/img/nomenclature?nomenclature=acetone
Call 2: GET https://YOUR-APP.herokuapp.com/img/nomenclature?nomenclature=benzene
```

### Step 3: Backend Response
```json
{
  "success": true,
  "cache_url": "https://YOUR-APP.herokuapp.com/cache/smiles_CC(=O)C_72b4c84c.svg",
  "svg": "<svg width='300' height='200'>...acetone structure...</svg>",
  "expires_in_hours": 24
}
```

---

## 🎯 What You See in ChatGPT (After Replacement)

```
"Let's analyze the structure of 

[ACETONE MOLECULE IMAGE]
📍 https://YOUR-APP.herokuapp.com/cache/smiles_CC(=O)C_72b4c84c.svg
⏱️ Expires in 24 hours

and compare it to 

[BENZENE MOLECULE IMAGE]
📍 https://YOUR-APP.herokuapp.com/cache/smiles_c1ccccc1_a3b9f21e.svg
⏱️ Expires in 24 hours
"
```

### Visual Breakdown

**BEFORE (raw text):**
```
chem:acetone
```

**AFTER (rendered):**
```
┌─────────────────────────────┐
│                             │
│     [Molecule Drawing]      │
│                             │
│         H₃C─CO─CH₃          │
│                             │
└─────────────────────────────┘
📍 Download: https://your-app.herokuapp.com/cache/smiles_CC(=O)C_hash.svg
```

---

## 🌍 Worldwide Link Sharing

### Scenario: You share the link with a colleague

**You copy from ChatGPT:**
```
https://mol2chemfig-kapil.herokuapp.com/cache/smiles_CC(=O)C_72b4c84c77f80f62f6005fe6ec837e72.svg
```

**Your colleague in Japan clicks it:**
```
✅ Browser downloads: acetone.svg
✅ Can edit in Inkscape/Adobe
✅ Works for 24 hours
✅ After 24 hours: 404 Not Found
```

---

## 📊 Supported Input Formats

### Format 1: Chemical Names (Nomenclature)
```
Input:  chem:acetone
Output: [Acetone structure SVG]

Input:  chem:benzene  
Output: [Benzene ring SVG]

Input:  chem:aspirin
Output: [Aspirin structure SVG]
```

### Format 2: SMILES Strings
```
Input:  chem:CCO
Output: [Ethanol structure SVG]

Input:  chem:CC(=O)C
Output: [Acetone structure SVG]

Input:  chem:c1ccccc1
Output: [Benzene ring SVG]
```

---

## 🎯 Real Example Walkthrough

### Example 1: Simple Molecule
```
ChatGPT Prompt:
"Show me chem:ethanol"

Extension Detects:
- Pattern: chem:ethanol
- Type: nomenclature (no special chars)

Extension Calls:
GET https://your-app.herokuapp.com/img/nomenclature?nomenclature=ethanol

Backend Converts:
ethanol → CCO (SMILES) → SVG

Extension Replaces Text With:
[ETHANOL SVG IMAGE]
📍 https://your-app.herokuapp.com/cache/smiles_CCO_abc123.svg
```

### Example 2: SMILES Format
```
ChatGPT Prompt:
"Analyze chem:CC(=O)C"

Extension Detects:
- Pattern: chem:CC(=O)C
- Type: SMILES (contains = and parentheses)

Extension Calls:
GET https://your-app.herokuapp.com/img/smiles?smiles=CC(=O)C

Backend Converts:
CC(=O)C (SMILES) → SVG

Extension Replaces Text With:
[ACETONE SVG IMAGE]
📍 https://your-app.herokuapp.com/cache/smiles_CC(=O)C_def456.svg
```

---

## 🔧 Interactive Features

### Feature 1: Click to Copy Link
```
User clicks on 📍 icon
→ Link copied to clipboard
→ Toast notification: "✅ Copied to clipboard!"
```

### Feature 2: Click to Download
```
User clicks on image
→ Browser downloads .svg file
→ File named: smiles_CCO_hash.svg
```

### Feature 3: Rotation Controls
```
User right-clicks image
→ Context menu: "Rotate 90°"
→ Image rotates without re-fetching
```

---

## 📸 Before & After Comparison

### BEFORE Extension
```
ChatGPT Output:
"The structure chem:benzene contains aromatic rings"

(Raw text, no rendering)
```

### AFTER Extension
```
ChatGPT Output:
"The structure [BENZENE RING IMAGE with 6-carbon aromatic ring] contains aromatic rings"

📍 Download: https://your-app.herokuapp.com/cache/smiles_c1ccccc1_hash.svg
(Clickable SVG image with worldwide download link)
```

---

## 🎨 Visual Rendering Details

### Image Container
```css
┌────────────────────────────────────┐
│  Background: #f9f9f9              │
│  Border: 1px solid #e0e0e0        │
│  Border-radius: 6px               │
│  Padding: 8px                     │
│                                    │
│  ┌──────────────────────────┐    │
│  │                          │    │
│  │   [Molecule SVG Image]   │    │
│  │   Max-width: 300px       │    │
│  │   Max-height: 200px      │    │
│  │                          │    │
│  └──────────────────────────┘    │
│                                    │
│  📍 https://your-app.../cache/... │
│  ⏱️ Expires in 24 hours          │
│                                    │
└────────────────────────────────────┘
```

---

## 🌟 Live Demo Scenarios

### Scenario 1: Chemistry Homework
```
Student asks ChatGPT:
"What is the structure of chem:glucose?"

ChatGPT sees:
"What is the structure of chem:glucose?"

Extension transforms to:
"What is the structure of [GLUCOSE STRUCTURE IMAGE]?"

Student can:
✅ See the structure visually
✅ Download the SVG for their report
✅ Share the link with classmates
```

### Scenario 2: Research Discussion
```
Researcher says:
"Compare chem:aspirin and chem:ibuprofen"

Extension shows:
[ASPIRIN IMAGE] vs [IBUPROFEN IMAGE]

Both with downloadable links:
📍 Link 1: aspirin SVG
📍 Link 2: ibuprofen SVG
```

### Scenario 3: Multiple Molecules
```
Input:
"Analyze chem:CCO, chem:acetone, and chem:benzene"

Output:
"Analyze [ETHANOL], [ACETONE], and [BENZENE]"

3 images rendered
3 download links provided
All work worldwide for 24 hours
```

---

## 🔄 Text Replacement Process

### Step-by-Step Transformation

**Original Text:**
```
"The compound chem:acetone is a solvent"
```

**Step 1: Scan (Extension)**
```
Found: "chem:acetone" at position 14
Type: nomenclature (plain text)
```

**Step 2: Fetch (API Call)**
```
GET /img/nomenclature?nomenclature=acetone
Response: {svg: "...", cache_url: "..."}
```

**Step 3: Replace (DOM Manipulation)**
```html
Before:
"The compound chem:acetone is a solvent"

After:
"The compound <div class="molecule-container">
  <img src="data:image/svg+xml;base64,...">
  <a href="https://.../cache/...">📍 Download</a>
</div> is a solvent"
```

**Step 4: Render (User Sees)**
```
"The compound [ACETONE IMAGE + LINK] is a solvent"
```

---

## 🎯 Summary

### What `chem:` Does
1. **Detects** chemical names or SMILES
2. **Converts** to SVG via your backend
3. **Replaces** text with image + link
4. **Provides** worldwide download link

### What You Get
- ✅ Instant molecule visualization
- ✅ Downloadable SVG files
- ✅ Worldwide shareable links
- ✅ 24-hour cache validity
- ✅ Just like CodeCogs, but better!

---

## 🚀 Ready to Use

After deploying to Heroku:

1. Type `chem:benzene` anywhere
2. See instant molecule rendering
3. Click link to download SVG
4. Share link worldwide
5. Enjoy your chemistry visualization! 🧪

---

**Created:** November 4, 2025  
**System:** MoleculeViewer + Chrome Extension  
**Status:** Production Ready
