# 🎯 QUICK START - Deploy in 5 Steps

---

## Step 1️⃣: Get Your Heroku App Name

```powershell
heroku apps
```

**Your app name:** `_________________`

---

## Step 2️⃣: Update 2 Files

### File 1: `MoleculeViewer\.env` (line 6)
```properties
PUBLIC_BASE_URL=https://YOUR-APP-NAME.herokuapp.com
```

### File 2: `chem-extension\content.js` (line 12)
```javascript
const MOLECULE_VIEWER_API = 'https://YOUR-APP-NAME.herokuapp.com';
```

---

## Step 3️⃣: Deploy to Heroku

```powershell
cd C:\Users\Kapil\Personal\PROJECTS\Mol2chemfig\MoleculeViewer
heroku git:remote -a YOUR-APP-NAME
git add .
git commit -m "Production deployment"
git push heroku main
heroku config:set PUBLIC_BASE_URL=https://YOUR-APP-NAME.herokuapp.com
```

---

## Step 4️⃣: Reload Extension

1. Open: `chrome://extensions/`
2. Find: "Chemistry Formula Renderer"
3. Click: Reload button (🔄)

---

## Step 5️⃣: Test It!

Open ChatGPT and type:
```
chem:benzene
```

**Expected:** Benzene ring appears with download link! ✅

---

## 🎉 You're Done!

Now you can use:
- `chem:acetone` → Shows acetone structure
- `chem:CCO` → Shows ethanol structure
- `chem:aspirin` → Shows aspirin structure

All with worldwide downloadable links! 🌍

---

## 📚 Need More Help?

Read these guides:
1. **DEPLOYMENT_CHECKLIST.md** - Full checklist
2. **HEROKU_DEPLOYMENT_STEPS.md** - Detailed guide
3. **VISUAL_GUIDE_TEXT_REPLACEMENT.md** - How it works
4. **READY_FOR_DEPLOYMENT.md** - Complete summary

---

**Created:** November 4, 2025  
**Status:** Ready to Deploy 🚀
