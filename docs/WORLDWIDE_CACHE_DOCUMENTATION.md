# 🌍 MoleculeViewer - Worldwide Accessible Cache Links Documentation

**Document Created:** November 4, 2025  
**Status:** Ready for Production Deployment  
**Privacy:** Code stays private, links are worldwide accessible

---

## 📋 Table of Contents

1. [Overview](#overview)
2. [How It Works](#how-it-works)
3. [Worldwide Accessibility](#worldwide-accessibility)
4. [Current Local Setup](#current-local-setup)
5. [Production Deployment Guide](#production-deployment-guide)
6. [Extension Integration](#extension-integration)
7. [Testing & Verification](#testing--verification)
8. [FAQ & Troubleshooting](#faq--troubleshooting)

---

## 🎯 Overview

### What is MoleculeViewer?

A **worldwide-accessible molecule caching service** that:
- ✅ Converts SMILES strings & chemical names to SVG images
- ✅ Caches SVGs for 24 hours with unique worldwide-accessible links
- ✅ Works exactly like CodeCogs (but self-hosted)
- ✅ Integrates with Chrome Extension for ChatGPT/web content
- ✅ Keeps your source code **PRIVATE** while serving public links

### Key Features

| Feature | Details |
|---------|---------|
| **SMILES Support** | `CCO`, `c1ccccc1`, complex molecules ✅ |
| **Nomenclature Support** | `acetone`, `benzene`, `aspirin` ✅ |
| **Cache Duration** | 24 hours then auto-delete |
| **Cache Size** | ~1-5 KB per SVG |
| **Worldwide Access** | Yes (just like CodeCogs!) |
| **CodeCogs Alternative** | Better: no rate limits, self-hosted |
| **Privacy** | Your code is private, links are public |

---

## 🔄 How It Works

### Architecture Diagram

```
┌─────────────────────────────────────────────────────┐
│                  YOUR BROWSER                       │
│  (Chrome Extension detects chem: prefix)            │
└──────────────────────┬──────────────────────────────┘
                       │
                       │ Sends SMILES/nomenclature
                       ↓
┌─────────────────────────────────────────────────────┐
│         FLASK BACKEND (Heroku/Railway)              │
│  - Parses chemistry via RDKit                       │
│  - Generates SVG image                             │
│  - Caches SVG with MD5 hash filename               │
│  - Returns JSON with cache link                    │
└──────────────────────┬──────────────────────────────┘
                       │
                       │ Returns: {
                       │   "success": true,
                       │   "cache_url": "https://your-backend/cache/smiles_CCO_hash.svg",
                       │   "svg": "<svg>...</svg>",
                       │   "expires_in_hours": 24
                       │ }
                       ↓
┌─────────────────────────────────────────────────────┐
│         WORLDWIDE ACCESSIBLE LINK                   │
│  https://your-backend.herokuapp.com/cache/         │
│  smiles_CCO_72b4c84c77f80f62f6005fe6ec837e72.svg  │
│                                                     │
│  Anyone, anywhere can download this SVG!          │
└─────────────────────────────────────────────────────┘
```

### Data Flow

1. **User types `chem:CCO`** in ChatGPT or webpage
2. **Chrome Extension** detects `chem:` prefix
3. **Extension sends** SMILES to your backend API
4. **Backend generates** SVG using RDKit
5. **SVG is cached** in `/svg-cache/` directory
6. **Unique URL returned** to extension
7. **Extension displays** molecule image
8. **Anyone with link** can download SVG for 24 hours

---

## 🌍 Worldwide Accessibility

### Just Like CodeCogs!

#### CodeCogs (External Service - Rate Limited)
```
https://latex.codecogs.com/svg.image?${formula}
❌ Rate limited
❌ Depends on external service
❌ Can't customize rendering
```

#### MoleculeViewer (Your Service - No Limits)
```
https://your-backend.herokuapp.com/cache/smiles_CCO_hash.svg
✅ No rate limits
✅ Self-hosted & controlled
✅ Full rendering customization
✅ 24-hour cache for performance
```

### How Anyone Can Access Your Links

**Scenario:** You send a link to a friend in Japan
```
https://mol2chemfig-kapil.herokuapp.com/cache/smiles_CCO_72b4c84c77f80f62f6005fe6ec837e72.svg
```

**Friend in Tokyo opens link:**
- ✅ Downloads the SVG successfully
- ✅ Can edit in Inkscape/Adobe
- ✅ Works for 24 hours
- ✅ After 24 hours: link expires (404)

**Why it works worldwide:**
- Heroku serves from global CDN
- DNS resolves worldwide
- CORS enabled for cross-domain requests
- No geographic restrictions

---

## 💻 Current Local Setup

### What You Have Right Now

**Flask Server Running Locally:**
```
http://192.168.1.4:5000
```

**Local Endpoints:**
```
GET /img/smiles?smiles=CCO
GET /img/nomenclature?nomenclature=acetone
GET /cache/{filename}
```

**Example Response:**
```json
{
  "success": true,
  "smiles": "CCO",
  "cache_url": "http://192.168.1.4:5000/cache/smiles_CCO_72b4c84c77f80f62f6005fe6ec837e72.svg",
  "expires_in_hours": 24,
  "svg": "<svg>...</svg>"
}
```

### Files You Have

```
MoleculeViewer/
├── app/
│   ├── api.py              ← Flask backend with caching
│   ├── chemistry.py        ← RDKit integration
│   └── config.py
├── templates/
│   └── index.html          ← Web UI with download button
├── .env                    ← Configuration
├── requirements.txt        ← Python dependencies
├── Procfile                ← Heroku deployment config
└── run.py                  ← Entry point
```

**Configuration:**
```env
# Current (Local):
PUBLIC_BASE_URL=http://192.168.1.4:5000

# Future (Production):
PUBLIC_BASE_URL=https://mol2chemfig-kapil.herokuapp.com
```

---

## 🚀 Production Deployment Guide

### When You're Ready to Deploy

#### Step 1: Create Private GitHub Repo
```
https://github.com/new
Name: Mol2chemfig
Visibility: PRIVATE ✅
```

#### Step 2: Push Code
```powershell
cd c:\Users\Kapil\Personal\PROJECTS\Mol2chemfig
git remote add origin https://github.com/YOUR_USERNAME/Mol2chemfig.git
git branch -M main
git push -u origin main
```

#### Step 3: Deploy to Heroku
```powershell
heroku login
heroku create mol2chemfig-kapil
cd MoleculeViewer
git push heroku main

# Set production URL
heroku config:set PUBLIC_BASE_URL=https://mol2chemfig-kapil.herokuapp.com
```

#### Step 4: Update .env
```env
PUBLIC_BASE_URL=https://mol2chemfig-kapil.herokuapp.com
```

#### Step 5: Live!
```
Frontend: https://your-username.github.io/Mol2chemfig (or GitHub Pages)
Backend: https://mol2chemfig-kapil.herokuapp.com
Cache Links: https://mol2chemfig-kapil.herokuapp.com/cache/smiles_CCO_hash.svg
```

---

## 🔌 Extension Integration

### How Chrome Extension Uses It

**File:** `chem-extension/content.js`

```javascript
// Detect chem: prefix
if (text.match(/chem:(\w+)/)) {
    const query = RegExp.$1;  // Extract "acetone" from "chem:acetone"
    
    // Fetch from your backend
    fetch('http://192.168.1.4:5000/img/nomenclature?nomenclature=' + query)
        .then(r => r.json())
        .then(data => {
            // Display SVG image
            showImage(data.svg);
            
            // Show cache link
            showLink(data.cache_url);
        })
}
```

### Extension Workflow

```
ChatGPT: "Look at this chem:aspirin"
    ↓
Extension detects: chem:aspirin
    ↓
Extension sends: GET /img/nomenclature?nomenclature=aspirin
    ↓
Backend returns: 
{
  "success": true,
  "cache_url": "https://your-backend/cache/smiles_CC(=O)Oc1ccccc1C(=O)O_hash.svg",
  "svg": "<svg>...</svg>"
}
    ↓
Extension displays: [SVG Image] + [Download Link]
```

---

## 🧪 Testing & Verification

### Test Locally (What You Should Do Now)

#### Test 1: Simple SMILES
```powershell
# Open browser or use curl
http://192.168.1.4:5000/img/smiles?smiles=CCO

# Should return JSON with cache_url
```

#### Test 2: Nomenclature
```powershell
http://192.168.1.4:5000/img/nomenclature?nomenclature=acetone

# Should convert acetone → SMILES → SVG
```

#### Test 3: Download SVG
```powershell
# Visit the cache_url from response
http://192.168.1.4:5000/cache/smiles_CCO_hash.svg

# Should download .svg file
```

#### Test 4: Extension Integration
1. Open ChatGPT
2. Type: `chem:benzene`
3. Extension should render molecule
4. Should show download link
5. Click link → SVG downloads

### Test Production (After Deployment)

#### Test 1: Worldwide Access
```powershell
# From any device, anywhere:
https://mol2chemfig-kapil.herokuapp.com/cache/smiles_CCO_hash.svg

# Should work worldwide ✅
```

#### Test 2: Verify Cache
```powershell
# Generate same molecule twice
GET /img/smiles?smiles=CCO (first time)
GET /img/smiles?smiles=CCO (second time)

# Should return same cache_url both times (using hash deduplication)
```

#### Test 3: Time Limit
```powershell
# Generate molecule
GET /img/smiles?smiles=CCO

# Wait 24 hours...
GET /cache/smiles_CCO_hash.svg
# Should return 404 (expired)
```

---

## ❓ FAQ & Troubleshooting

### Q: Will my code be visible?
**A:** No! Only your compiled Flask API endpoints are public. Source code stays private on GitHub.

### Q: Can people misuse the links?
**A:** They expire in 24 hours. Each link is unique (MD5 hash). No authentication needed but hard to guess/brute-force.

### Q: What if I want to add authentication?
**A:** Easy! Add API key requirement to endpoints later:
```python
@app.route('/img/smiles')
def img_smiles():
    api_key = request.args.get('api_key')
    if api_key != SECRET_KEY:
        return {"error": "Unauthorized"}, 401
    # ... rest of code
```

### Q: How much storage?
**A:** 
- Each SVG: ~1-5 KB
- Cache: 24 hours only
- Example: 1000 molecules/day = ~5 MB (auto-deleted daily)

### Q: Cost?
**A:** 
- **Heroku:** Free tier (0.5 dyno) = always free or ~$7/month for always-on
- **GitHub:** Free for private repos
- **Total:** Free to ~$7/month

### Q: What if I need to clear cache?
**A:**
```powershell
# SSH into Heroku
heroku run bash
# Remove old files
rm svg-cache/*

# Or redeploy to reset
git push heroku main
```

### Q: Can I use Railway instead of Heroku?
**A:** Yes! Railway is actually easier:
1. Go to https://railway.app
2. Connect GitHub repo
3. Add environment variable: `PUBLIC_BASE_URL=https://your-railway-url.up.railway.app`
4. Done! ✅

### Q: What's the difference from CodeCogs?
**A:**

| Feature | CodeCogs | MoleculeViewer |
|---------|----------|----------------|
| Rate Limit | 50/day (free) | Unlimited |
| Cost | Free (limited) | Free (self-hosted) |
| Customization | No | Yes |
| Privacy | External service | Your control |
| Setup | Instant | 5 minutes |
| Caching | CodeCogs | Yours |

---

## 📦 To Deploy Later

When you're ready, just run:

```powershell
# 1. Create private GitHub repo
# 2. Push code
git remote add origin https://github.com/YOUR_USERNAME/Mol2chemfig.git
git push -u origin main

# 3. Deploy backend
heroku create mol2chemfig-kapil
git push heroku main

# 4. Set environment
heroku config:set PUBLIC_BASE_URL=https://mol2chemfig-kapil.herokuapp.com

# 5. Live!
```

Your worldwide-accessible links will be ready! 🌍

---

## 🎯 Current Status

✅ **LOCAL:** Working at http://192.168.1.4:5000  
✅ **EXTENSION:** Can access local links  
✅ **CACHING:** 24-hour cache system implemented  
✅ **DOWNLOAD:** SVG download button working  
⏳ **PRODUCTION:** Ready to deploy (awaiting your decision)

---

## 📞 Next Steps

1. **Test locally** with extension (see Testing section)
2. **Verify worldwide access** from multiple devices
3. **When ready:** Follow Production Deployment Guide
4. **Optional:** Add authentication/API keys

**Save this document for later reference!** 📄

---

**Created:** November 4, 2025  
**System:** MoleculeViewer v1.0  
**Author:** Development Team  
**Status:** Production Ready
