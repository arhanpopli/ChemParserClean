# GitHub Pages + Flask Backend Setup Guide

## Quick Start - Deploy in 5 Steps

### Step 1: Create GitHub Repository
1. Go to https://github.com/new
2. Name it: `Mol2chemfig` (or your preferred name)
3. Make it **PUBLIC**
4. Click "Create repository"

### Step 2: Push Code to GitHub
```bash
cd Mol2chemfig
git add .
git commit -m "Initial MoleculeViewer setup"
git remote add origin https://github.com/YOUR_USERNAME/Mol2chemfig.git
git branch -M main
git push -u origin main
```

### Step 3: Enable GitHub Pages
1. Go to your repo → Settings → Pages
2. Select "Deploy from a branch"
3. Choose branch: `gh-pages`
4. Click Save

### Step 4: Deploy Flask Backend
Choose ONE option:

#### Option A: Free Heroku (Recommended)
1. Create account at https://heroku.com
2. Install Heroku CLI
3. Run:
```bash
heroku login
heroku create your-app-name
git push heroku main
```
4. Set environment variable:
```bash
heroku config:set PUBLIC_BASE_URL=https://your-app-name.herokuapp.com
```

#### Option B: Railway.app (Even Easier)
1. Go to https://railway.app
2. Connect GitHub repo
3. Deploy with one click
4. Set `PUBLIC_BASE_URL` environment variable

#### Option C: Vercel (Fastest)
1. Go to https://vercel.com
2. Import your GitHub repo
3. Deploy
4. Set environment variables

### Step 5: Update Configuration
After deployment, update `.env` with your backend URL:
```env
PUBLIC_BASE_URL=https://your-heroku-app.herokuapp.com
# or
PUBLIC_BASE_URL=https://your-railway-app.up.railway.app
# or
PUBLIC_BASE_URL=https://your-vercel-app.vercel.app
```

## Worldwide Access URLs

Once deployed, you'll get:
- **Frontend (GitHub Pages)**: `https://USERNAME.github.io/Mol2chemfig`
- **Backend (Heroku/Railway/Vercel)**: `https://your-backend-domain.com`
- **Cache Links**: `https://your-backend-domain.com/cache/smiles_CCO_hash.svg`

All cache links will be **worldwide accessible** for 24 hours!

## Local Testing

Before deploying, test locally:
```bash
cd MoleculeViewer
python run.py
# Open http://192.168.1.4:5000
```

## Architecture

```
┌─────────────────────────────────────┐
│   GitHub Pages (Frontend)           │
│   HTML + JavaScript UI              │
│   https://username.github.io/mol2   │
└────────────┬────────────────────────┘
             │ Fetches SVG data
             ↓
┌─────────────────────────────────────┐
│   Heroku/Railway (Backend)          │
│   Flask API Server                  │
│   https://your-backend.herokuapp.com │
│   ├─ /img/smiles                    │
│   ├─ /img/nomenclature              │
│   └─ /cache/{filename}              │
└─────────────────────────────────────┘
```

## Support

For issues, check:
1. GitHub Actions tab → Workflows (for deployment errors)
2. Heroku/Railway dashboard (for backend errors)
3. Browser console (F12) for frontend errors

Need help? File an issue on GitHub!
