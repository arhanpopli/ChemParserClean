# 🌍 GitHub Pages + Heroku Deployment - LIVE NOW! 🚀

## ✨ What You Get

✅ **Worldwide accessible cache links** for 24 hours  
✅ **GitHub Pages** frontend hosting (FREE)  
✅ **Heroku/Railway/Vercel** backend (FREE tier available)  
✅ **Automatic CI/CD** via GitHub Actions  
✅ **Global URL sharing** - anyone can download SVGs from your link  

Example worldwide-accessible cache URL:
```
https://your-heroku-app.herokuapp.com/cache/smiles_CCO_72b4c84c77f80f62f6005fe6ec837e72.svg
```

---

## 🚀 Deploy in 5 Minutes

### STEP 1: Create GitHub Repo
```
https://github.com/new
Name: Mol2chemfig
Make it PUBLIC ✅
Create Repository
```

### STEP 2: Push Code
```bash
cd c:\Users\Kapil\Personal\PROJECTS\Mol2chemfig
git config user.name "Your Name"
git config user.email "your@email.com"
git add .
git commit -m "🧬 MoleculeViewer - worldwide accessible caching"
git remote add origin https://github.com/YOUR_USERNAME/Mol2chemfig.git
git branch -M main
git push -u origin main
```

### STEP 3: Deploy Backend
Choose ONE:

**Option A: Heroku (Recommended)**
```bash
# Sign up: https://heroku.com
# Install CLI from https://devcenter.heroku.com/articles/heroku-cli

heroku login
heroku create mol2chemfig-kapil
git push heroku main
heroku config:set PUBLIC_BASE_URL=https://mol2chemfig-kapil.herokuapp.com
```

**Option B: Railway (Super Easy)**
```
1. Go to https://railway.app
2. Click "New Project"
3. Select "Deploy from GitHub repo"
4. Select your Mol2chemfig repo
5. Add env var: PUBLIC_BASE_URL=https://[your-railway-url].up.railway.app
6. Deploy!
```

### STEP 4: Enable GitHub Pages
In your repo:
1. Settings → Pages
2. Source: "Deploy from a branch"
3. Branch: gh-pages
4. Save

### STEP 5: Update & Push Production Config
```bash
# Edit MoleculeViewer/.env
PUBLIC_BASE_URL=https://mol2chemfig-kapil.herokuapp.com

# Push changes
git add .
git commit -m "Update production URL"
git push origin main
```

---

## ✅ Test It Now!

After deployment:
1. Visit: `https://YOUR_USERNAME.github.io/Mol2chemfig`
2. Enter: `CCO` (ethanol)
3. Click: "SAVE AS SVG"
4. See: Worldwide accessible link! 🌍

---

## 🎯 What Happens

```
┌─────────────────────────────────────────────────────┐
│                  YOUR BROWSER                       │
│  https://your-username.github.io/Mol2chemfig       │
└────────────────────┬────────────────────────────────┘
                     │ 
                     │ Fetch SVG data
                     ↓
┌─────────────────────────────────────────────────────┐
│              HEROKU BACKEND (API)                   │
│  https://mol2chemfig-kapil.herokuapp.com            │
│  ├─ /img/smiles?smiles=CCO                         │
│  ├─ /img/nomenclature?nomenclature=acetone         │
│  └─ /cache/smiles_CCO_hash.svg ← WORLDWIDE LINK!  │
└─────────────────────────────────────────────────────┘
```

**Result:** Anyone with the cache link can download the SVG for 24 hours! 🎉

---

## 🔗 GitHub Pages URL
After enabling GitHub Pages, your site is live at:
```
https://YOUR_GITHUB_USERNAME.github.io/Mol2chemfig
```

Replace `YOUR_GITHUB_USERNAME` with your actual GitHub username!

---

## 📝 Notes

- **Cache expires in 24 hours** - after which the link returns 404
- **Worldwide accessible** - works from any IP address
- **Free tier** - GitHub Pages + Heroku free tier = $0/month
- **Custom domain** - can add later for professional look

## Next: Production URL Sharing

Once live, share links like:
```
Download my molecule: https://mol2chemfig-kapil.herokuapp.com/cache/smiles_CCO_72b4c84c77f80f62f6005fe6ec837e72.svg
```

Anyone can download it for 24 hours! 🌍
4. ✅ Structures render with your chosen settings

---

## 3-Minute Setup

### 1️⃣ Start Server (Terminal 1)
```bash
cd MoleculeViewer
python run_server.py
```

### 2️⃣ Load Extension (Chrome)
- `chrome://extensions/` 
- Developer mode: ON
- Load unpacked → select `chem-extension/`

### 3️⃣ Open Popup
- Click extension icon
- Select "🧪 MoleculeViewer (Best)"
- See new "🧪 MoleculeViewer Options" section!

### 4️⃣ Test
- Go to ChatGPT
- Use: `\chemfig{C-C-C}`
- Renders from localhost:5000 ✓
- Press F12 → check console for `🔬 Using MoleculeViewer server`

---

## 8 Rendering Options

Toggle these in the popup (only show when MoleculeViewer selected):

| Option | Toggle/Dropdown | What It Does |
|--------|---|---|
| 🔘 Show Carbon Atoms | Toggle | Shows C labels |
| 🔘 Show Methyl Groups | Toggle | Shows CH₃ labels |
| 🔘 Aromatic Circles | Toggle | Circles in benzene rings |
| 🔘 Fancy Bonds | Toggle | Better bond lines |
| 🔘 Atom Numbers | Toggle | Number atoms 1,2,3... |
| 🔘 Flip Horizontal | Toggle | Mirror left-right |
| 🔘 Flip Vertical | Toggle | Mirror up-down |
| 📝 Hydrogen Display | Dropdown | keep / add / delete |

---

## Where Files Are

```
chem-extension/
├── content.js          ← Settings, rendering logic, image loading
├── popup.html          ← New MoleculeViewer Options section
└── popup.js            ← Event handlers for new options

MoleculeViewer/
├── app/api.py          ← Server with /api/smiles-to-svg endpoint
└── run_server.py       ← Start server here

Documentation/
├── EXTENSION_FIX_SUMMARY.md       ← What was fixed
├── EXTENSION_UPDATE_LOG.md        ← Technical details
├── INTEGRATION_COMPLETE.md        ← Full summary
└── TESTING_CHECKLIST.md           ← 40+ tests
```

---

## Troubleshooting

### ❓ Still showing CodeCogs URL?
- [ ] Is server running? (`python run_server.py`)
- [ ] Is popup showing "🧪 MoleculeViewer (Best)" selected?
- [ ] Reload page (F5)
- [ ] Check console (F12) for errors

### ❓ MoleculeViewer Options don't show?
- [ ] Select "🧪 MoleculeViewer (Best)" in dropdown
- [ ] Options should appear below
- [ ] If not, check for errors in console

### ❓ Options don't seem to work?
- [ ] Toggle option (e.g., "Show Carbon Atoms" ON)
- [ ] See success message: "Show carbons enabled..."
- [ ] **Reload page with F5** (options apply on reload)
- [ ] Check console: should see your options logged

### ❓ Server connection error?
- [ ] Terminal: `python run_server.py`
- [ ] Browser: Check server is at `localhost:5000`
- [ ] F12 Console: Look for `❌ MoleculeViewer fetch failed`
- [ ] Falls back to CodeCogs automatically

---

## What's Different

```
BEFORE (Issue):
  Extension popup → Select MoleculeViewer
  ❌ Option appears in dropdown
  ❌ But still uses CodeCogs
  ❌ No UI to configure rendering
  ❌ Settings lost on reload

NOW (Fixed):
  Extension popup → Select MoleculeViewer  
  ✅ New "MoleculeViewer Options" section appears
  ✅ Uses localhost:5000/api/smiles-to-svg
  ✅ 8 rendering options to configure
  ✅ Settings saved & persist between sessions
  ✅ Options sent in POST request
  ✅ Structures render with your preferences
```

---

## Key Code Changes

### Settings Loading (Fixed)
```javascript
// Was losing defaults
chrome.storage.sync.get(settings, (result) => {
  settings = result;  // ❌ Lost defaults
});

// Now keeps defaults + merges stored
chrome.storage.sync.get(null, (result) => {
  settings = { ...settings, ...result };  // ✅ Works!
});
```

### Rendering (Improved)
```javascript
// Was simple GET string
return `http://localhost:5000/api/render-smiles?smiles=${smiles}`;

// Now full POST with all options
return {
  isMoleculeViewer: true,
  smiles: smiles,
  options: {
    show_carbons: settings.showCarbons,
    show_methyls: settings.showMethyls,
    // ... all 8 options ...
  }
};
```

### Image Loading (New)
```javascript
// New function for async MoleculeViewer rendering
function loadMoleculeViewerImage(img) {
  fetch('http://localhost:5000/api/smiles-to-svg', {
    method: 'POST',
    body: JSON.stringify({
      smiles: moleculeData.smiles,
      options: moleculeData.options  // ← Your settings!
    })
  })
  .then(response => response.json())
  .then(data => {
    // Render SVG with your chosen options
    img.src = URL.createObjectURL(new Blob([data.svg]));
  });
}
```

---

## Test Scenarios

### Quick Test (2 min)
```
1. Start server: python run_server.py
2. Load extension
3. Go to ChatGPT
4. Type: \chemfig{C-C-C}
5. See structure render ✓
```

### Full Test (10 min)
```
1. Open popup → select MoleculeViewer
2. Toggle "Show Carbon Atoms" ON
3. See success message
4. Go to ChatGPT
5. Type: \chemfig{C-C-C}
6. Reload page
7. Structure shows C labels ✓
```

### Comprehensive Test
See: `TESTING_CHECKLIST.md` (40+ test cases)

---

## Console Commands (Developer)

**See all logs:**
```javascript
chemRendererLogs()  // Returns table of all logs
```

**Check settings:**
```javascript
console.log(settings)  // Current settings
```

**Check rendering options:**
```javascript
chrome.storage.sync.get(null, (result) => console.log(result))
```

---

## FAQ

**Q: Do I need to edit any Python code?**
A: No! The server endpoint already exists and handles all options.

**Q: Do I need to restart the browser?**
A: No, just reload the webpage (F5) for settings to apply.

**Q: Can I use other rendering engines?**
A: Yes! Dropdown still has CodeCogs, LaTeX.Online, QuickLaTeX options.

**Q: What if server is not running?**
A: Falls back to CodeCogs automatically (no errors).

**Q: Do settings persist?**
A: Yes! Uses `chrome.storage.sync` (syncs across browsers too).

**Q: Can I reset to defaults?**
A: Clear Chrome storage: DevTools → Application → Clear.

---

## Performance

- ⚡ Max 3 structures load simultaneously (prevents lag)
- 🚀 Lazy-loading: only renders visible structures
- 🔄 Fast reload with caching
- 📱 Works on mobile too

---

## Next Steps

1. ✅ Server running: `python run_server.py`
2. ✅ Extension loaded unpacked
3. ✅ Select MoleculeViewer in popup
4. ✅ Test on ChatGPT or webpage
5. 📝 Report any issues
6. 🎉 Use the new rendering options!

---

## Still Have Questions?

See these documents:
- **What changed?** → `EXTENSION_UPDATE_LOG.md`
- **How does it work?** → `EXTENSION_FIX_SUMMARY.md`
- **Test everything:** → `TESTING_CHECKLIST.md`
- **Full details?** → `INTEGRATION_COMPLETE.md`

---

**Status: ✅ READY TO USE**

Extension now properly uses MoleculeViewer with all rendering options! 🧪✨
