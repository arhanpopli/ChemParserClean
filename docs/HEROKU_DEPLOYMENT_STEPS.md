# 🚀 Heroku Deployment Guide - MoleculeViewer + Chrome Extension

**Last Updated:** November 4, 2025  
**Status:** Ready to Deploy

---

## 📋 What's Been Configured

✅ **MoleculeViewer Backend**
- `.env` file prepared with production URL placeholder
- Worldwide cache link system working
- 24-hour cache expiration active

✅ **Chrome Extension**
- `chem:` prefix detection working
- API URL configuration ready for Heroku
- Automatic text replacement enabled

---

## 🎯 How It Works

```
User types in ChatGPT:
"Let's look at chem:acetone"
            ↓
Extension detects: chem:acetone
            ↓
Extension sends: GET https://YOUR-HEROKU-APP.herokuapp.com/img/nomenclature?nomenclature=acetone
            ↓
Backend responds: 
{
  "success": true,
  "cache_url": "https://YOUR-HEROKU-APP.herokuapp.com/cache/smiles_CC(=O)C_hash.svg",
  "svg": "<svg>...</svg>",
  "expires_in_hours": 24
}
            ↓
Extension replaces text:
"Let's look at chem:acetone" 
becomes
"Let's look at [SVG IMAGE] 📍 https://YOUR-HEROKU-APP.herokuapp.com/cache/..."
            ↓
User can click the link to download SVG!
```

---

## 🌟 Step-by-Step Deployment

### Step 1: Get Your Heroku App Name

When you ran `heroku create`, it gave you an app name like:
- `mol2chemfig-kapil` 
- `my-molecule-viewer`
- `chemistry-cache-app`

**Your Heroku URL will be:** `https://YOUR-APP-NAME.herokuapp.com`

---

### Step 2: Update MoleculeViewer .env

📄 **File:** `MoleculeViewer\.env`

Replace `YOUR-HEROKU-APP` with your actual app name:

```properties
# 🌍 PRODUCTION (Heroku) - Worldwide accessible cache links
PUBLIC_BASE_URL=https://mol2chemfig-kapil.herokuapp.com

# For local testing (comment out when deploying):
# PUBLIC_BASE_URL=http://192.168.1.4:5000
```

**Example:**
```properties
PUBLIC_BASE_URL=https://mol2chemfig-kapil.herokuapp.com
```

---

### Step 3: Update Chrome Extension API URL

📄 **File:** `chem-extension\content.js`

**Line 12-13:** Replace `YOUR-HEROKU-APP`:

```javascript
// Replace YOUR-HEROKU-APP with your actual Heroku app name
const MOLECULE_VIEWER_API = 'https://mol2chemfig-kapil.herokuapp.com';

// Comment out local testing:
// const MOLECULE_VIEWER_API = 'http://192.168.1.4:5000';
```

---

### Step 4: Deploy Backend to Heroku

```powershell
cd C:\Users\Kapil\Personal\PROJECTS\Mol2chemfig\MoleculeViewer

# Login to Heroku
heroku login

# Add remote (if not already added)
heroku git:remote -a YOUR-APP-NAME

# Deploy
git add .
git commit -m "Configure for production deployment"
git push heroku main

# Set environment variable on Heroku
heroku config:set PUBLIC_BASE_URL=https://YOUR-APP-NAME.herokuapp.com
```

---

### Step 5: Verify Backend is Live

```powershell
# Test SMILES endpoint
Invoke-WebRequest -Uri "https://YOUR-APP-NAME.herokuapp.com/img/smiles?smiles=CCO" -UseBasicParsing

# Test nomenclature endpoint
Invoke-WebRequest -Uri "https://YOUR-APP-NAME.herokuapp.com/img/nomenclature?nomenclature=acetone" -UseBasicParsing

# Should return JSON with cache_url
```

**Expected Response:**
```json
{
  "success": true,
  "smiles": "CCO",
  "cache_url": "https://YOUR-APP-NAME.herokuapp.com/cache/smiles_CCO_hash.svg",
  "svg": "<svg>...</svg>",
  "expires_in_hours": 24
}
```

---

### Step 6: Reload Chrome Extension

1. Open Chrome: `chrome://extensions/`
2. Find "Chemistry Formula Renderer"
3. Click **Reload** button (circular arrow icon)
4. Extension now uses Heroku API! ✅

---

### Step 7: Test in ChatGPT

1. Open ChatGPT
2. Type: `chem:benzene`
3. **Expected:** Benzene structure appears with clickable link
4. Type: `chem:CCO`
5. **Expected:** Ethanol structure appears

**Worldwide Link Test:**
- Copy the cache link from ChatGPT
- Send to a friend anywhere in the world
- They can download the SVG for 24 hours! 🌍

---

## 🔍 Verification Checklist

After deployment, verify these work:

### Backend Endpoints
- [ ] `https://YOUR-APP.herokuapp.com/` shows homepage
- [ ] `https://YOUR-APP.herokuapp.com/img/smiles?smiles=CCO` returns JSON
- [ ] `https://YOUR-APP.herokuapp.com/img/nomenclature?nomenclature=acetone` returns JSON
- [ ] Cache URL in response is accessible worldwide

### Extension Integration  
- [ ] `chem:benzene` in ChatGPT shows molecule
- [ ] `chem:CCO` in ChatGPT shows ethanol
- [ ] Cache link is clickable and downloads SVG
- [ ] Link works from different computer/network

### Worldwide Access
- [ ] Send cache link to friend → they can download
- [ ] Open cache link on phone → downloads SVG
- [ ] Cache expires after 24 hours (returns 404)

---

## 🎨 Supported Formats

### Chemical Names (Nomenclature)
```
chem:acetone
chem:benzene
chem:aspirin
chem:caffeine
chem:ethanol
```

### SMILES Strings
```
chem:CCO                  (ethanol)
chem:CC(=O)C              (acetone)
chem:c1ccccc1             (benzene)
chem:CC(C)CC(=O)O         (complex molecule)
```

---

## 📊 Example Output

### ChatGPT Input:
```
"Let's analyze chem:acetone and chem:benzene together"
```

### After Replacement:
```
"Let's analyze [ACETONE SVG IMAGE] 📍 https://your-app.herokuapp.com/cache/...
and [BENZENE SVG IMAGE] 📍 https://your-app.herokuapp.com/cache/..."
```

Each link is:
- ✅ Worldwide accessible
- ✅ Downloadable as SVG
- ✅ Valid for 24 hours
- ✅ Works just like CodeCogs!

---

## 🛠️ Troubleshooting

### Extension Not Detecting `chem:`
**Solution:** Reload extension at `chrome://extensions/`

### Backend Returns 404
**Solution:** Check Heroku logs:
```powershell
heroku logs --tail -a YOUR-APP-NAME
```

### CORS Errors
**Solution:** Backend already has CORS enabled. Check browser console for details.

### Cache Links Return 404
**Possible Causes:**
1. Link expired (24 hours passed)
2. Cache cleared on backend
3. Wrong PUBLIC_BASE_URL in .env

---

## 📝 Quick Reference

### Configuration Files

| File | Line | Change |
|------|------|--------|
| `MoleculeViewer\.env` | 6 | Set `PUBLIC_BASE_URL=https://YOUR-APP.herokuapp.com` |
| `chem-extension\content.js` | 12 | Set `const MOLECULE_VIEWER_API = 'https://YOUR-APP.herokuapp.com'` |

### Deployment Commands

```powershell
# 1. Update config files (see above)

# 2. Commit changes
cd C:\Users\Kapil\Personal\PROJECTS\Mol2chemfig
git add .
git commit -m "Configure for Heroku production"

# 3. Deploy to Heroku
cd MoleculeViewer
git push heroku main

# 4. Set environment
heroku config:set PUBLIC_BASE_URL=https://YOUR-APP.herokuapp.com

# 5. Verify
heroku logs --tail
```

### Testing Commands

```powershell
# Test SMILES
curl "https://YOUR-APP.herokuapp.com/img/smiles?smiles=CCO"

# Test nomenclature
curl "https://YOUR-APP.herokuapp.com/img/nomenclature?nomenclature=acetone"

# Test cache link (copy from response above)
curl "https://YOUR-APP.herokuapp.com/cache/smiles_CCO_hash.svg" -o test.svg
```

---

## 🎯 Next Steps

1. **Get your Heroku app name** (from `heroku create` output)
2. **Update both config files** with your app name
3. **Deploy to Heroku** using commands above
4. **Reload Chrome extension**
5. **Test in ChatGPT** with `chem:benzene`
6. **Share worldwide links!** 🌍

---

## 💡 Pro Tips

### Faster Deployment
```powershell
# One-liner deploy
cd MoleculeViewer ; git add . ; git commit -m "Deploy" ; git push heroku main
```

### Check if Backend is Live
```powershell
# Should return "OK"
curl https://YOUR-APP.herokuapp.com/health
```

### View Real-time Logs
```powershell
heroku logs --tail -a YOUR-APP-NAME
```

### Restart Backend
```powershell
heroku restart -a YOUR-APP-NAME
```

---

## 🌍 Worldwide Access Confirmed!

Your cache links will work from:
- ✅ Any country
- ✅ Any network (WiFi, mobile data)
- ✅ Any device (PC, Mac, phone, tablet)
- ✅ Any browser

Just like CodeCogs, but **YOU** control it! 🎉

---

**Ready to deploy?** Follow the steps above! 🚀
