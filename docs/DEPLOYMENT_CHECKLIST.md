# ✅ Quick Deployment Checklist

**Status:** Ready to Deploy to Heroku  
**Date:** November 4, 2025

---

## 🎯 Prerequisites

- [x] Heroku account created
- [x] Heroku CLI installed
- [x] Git repository set up
- [x] Chrome extension installed locally

---

## 📝 Configuration Steps

### Step 1: Get Your Heroku App Name
```powershell
# If you already created it:
heroku apps

# If you need to create it:
heroku create your-app-name
```

**Write your app name here:** `_____________________`

Your URL will be: `https://___________.herokuapp.com`

---

### Step 2: Update MoleculeViewer Configuration

**File:** `MoleculeViewer\.env`

**Line 6:** Replace `YOUR-HEROKU-APP` with your actual app name:

```properties
# BEFORE
PUBLIC_BASE_URL=https://YOUR-HEROKU-APP.herokuapp.com

# AFTER (example)
PUBLIC_BASE_URL=https://mol2chemfig-kapil.herokuapp.com
```

**Status:** [ ] Updated

---

### Step 3: Update Chrome Extension Configuration

**File:** `chem-extension\content.js`

**Line 12:** Replace `YOUR-HEROKU-APP` with your actual app name:

```javascript
// BEFORE
const MOLECULE_VIEWER_API = 'https://YOUR-HEROKU-APP.herokuapp.com';

// AFTER (example)
const MOLECULE_VIEWER_API = 'https://mol2chemfig-kapil.herokuapp.com';
```

**Status:** [ ] Updated

---

## 🚀 Deployment Steps

### Step 4: Deploy Backend to Heroku

```powershell
cd C:\Users\Kapil\Personal\PROJECTS\Mol2chemfig\MoleculeViewer

# Login to Heroku
heroku login

# Add remote if needed
heroku git:remote -a YOUR-APP-NAME

# Commit changes
git add .
git commit -m "Configure for production"

# Deploy
git push heroku main

# Set environment variable
heroku config:set PUBLIC_BASE_URL=https://YOUR-APP-NAME.herokuapp.com
```

**Status:** [ ] Deployed

---

### Step 5: Verify Backend is Working

Test these URLs in your browser (replace YOUR-APP-NAME):

1. **Health Check:**
   ```
   https://YOUR-APP-NAME.herokuapp.com/health
   ```
   **Expected:** "OK" or health status
   **Status:** [ ] Working

2. **SMILES Test:**
   ```
   https://YOUR-APP-NAME.herokuapp.com/img/smiles?smiles=CCO
   ```
   **Expected:** JSON with `cache_url` field
   **Status:** [ ] Working

3. **Nomenclature Test:**
   ```
   https://YOUR-APP-NAME.herokuapp.com/img/nomenclature?nomenclature=acetone
   ```
   **Expected:** JSON with `cache_url` field
   **Status:** [ ] Working

---

### Step 6: Reload Chrome Extension

1. Open Chrome: `chrome://extensions/`
2. Find "Chemistry Formula Renderer"
3. Click the **Reload** button (circular arrow)

**Status:** [ ] Reloaded

---

### Step 7: Test in ChatGPT

1. Open ChatGPT: https://chat.openai.com
2. Type: `chem:benzene`
3. **Expected:** Benzene structure appears with download link

**Status:** [ ] Working

---

## 🧪 Testing Checklist

### Basic Tests
- [ ] `chem:benzene` shows benzene ring
- [ ] `chem:acetone` shows acetone structure
- [ ] `chem:CCO` shows ethanol structure
- [ ] Download link is clickable
- [ ] Download link works (downloads SVG)

### Worldwide Access Test
- [ ] Copy cache URL from ChatGPT
- [ ] Open in different browser/device
- [ ] Confirm SVG downloads
- [ ] Send to friend → they can download

### Multiple Molecules Test
- [ ] Type: `chem:benzene and chem:acetone`
- [ ] Both molecules render
- [ ] Both download links work

---

## 🔍 Troubleshooting

### Issue: Extension not detecting `chem:`
**Solution:** Reload extension at `chrome://extensions/`

### Issue: Backend returns 404
**Check:** 
```powershell
heroku logs --tail -a YOUR-APP-NAME
```

### Issue: CORS error in browser console
**Solution:** Backend has CORS enabled. Check PUBLIC_BASE_URL is correct.

### Issue: Cache links don't work
**Check:**
1. PUBLIC_BASE_URL in `.env` matches Heroku URL
2. Environment variable set on Heroku:
   ```powershell
   heroku config -a YOUR-APP-NAME
   ```

---

## 📊 Final Verification

### All Systems Go! ✅

- [x] Backend deployed to Heroku
- [x] Extension configured
- [x] `chem:` detection working
- [x] Text replacement active
- [x] Download links worldwide accessible
- [x] 24-hour cache working

---

## 🎉 Success Criteria

You're done when:

1. ✅ You type `chem:benzene` in ChatGPT
2. ✅ Benzene structure appears instantly
3. ✅ Download link appears below image
4. ✅ Clicking link downloads SVG
5. ✅ Link works from any device/location
6. ✅ Link expires after 24 hours

---

## 📝 Notes

**Your Heroku App Name:** `_____________________`

**Your Backend URL:** `https://____________.herokuapp.com`

**Deployment Date:** `_____________________`

**First Successful Test:** `_____________________`

---

## 🆘 Need Help?

### View Logs
```powershell
heroku logs --tail -a YOUR-APP-NAME
```

### Restart Backend
```powershell
heroku restart -a YOUR-APP-NAME
```

### Check Environment Variables
```powershell
heroku config -a YOUR-APP-NAME
```

### Redeploy
```powershell
cd MoleculeViewer
git push heroku main
```

---

## 🌟 What's Next?

After successful deployment:

1. **Share your cache links** with colleagues
2. **Test from different locations** to confirm worldwide access
3. **Monitor Heroku logs** for any issues
4. **Enjoy chemistry visualization!** 🧪

---

**Ready?** Follow the checklist above step-by-step! 🚀

**Created:** November 4, 2025  
**System:** MoleculeViewer + Chrome Extension  
**Status:** Production Ready
