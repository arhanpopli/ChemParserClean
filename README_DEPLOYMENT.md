# 🧬 MoleculeViewer

**Convert SMILES & Chemical Names to SVG - Worldwide Accessible Caching**

![Status](https://img.shields.io/badge/Status-Production%20Ready-brightgreen)
![License](https://img.shields.io/badge/License-MIT-blue)
![Python](https://img.shields.io/badge/Python-3.8+-blue)

---

## 🌍 Features

- ✅ **SMILES to SVG** - Convert chemical formulas to beautiful SVG renderings
- ✅ **Nomenclature Lookup** - Convert chemical names to SMILES to SVG
- ✅ **24-Hour Caching** - Generated images cached with unique URLs
- ✅ **Worldwide Accessible** - Share cache links globally
- ✅ **Chrome Extension** - Right-click to render molecules from `chem:` prefix
- ✅ **GitHub Pages + Heroku** - Free tier deployment included
- ✅ **Real-time Rendering** - Powered by RDKit

---

## 🚀 Quick Deploy (GitHub Pages + Heroku)

See **[QUICK_START.md](./QUICK_START.md)** for 5-minute deployment!

### TL;DR:
```bash
# 1. Create GitHub repo
https://github.com/new

# 2. Push code
git push origin main

# 3. Deploy backend
heroku create your-app-name
git push heroku main

# 4. Enable GitHub Pages in Settings

# 5. Live at: https://username.github.io/Mol2chemfig ✨
```

---

## 💻 Local Development

### Requirements
- Python 3.8+
- RDKit
- Flask

### Setup
```bash
cd MoleculeViewer
pip install -r requirements.txt
python run.py
```

Visit: http://192.168.1.4:5000 or http://localhost:5000

---

## 📡 API Endpoints

### Generate SVG from SMILES
```bash
GET /img/smiles?smiles=CCO&width=300&height=200
```

Response:
```json
{
  "success": true,
  "smiles": "CCO",
  "cache_url": "http://192.168.1.4:5000/cache/smiles_CCO_hash.svg",
  "expires_in_hours": 24,
  "svg": "<svg>...</svg>"
}
```

### Generate SVG from Nomenclature
```bash
GET /img/nomenclature?nomenclature=acetone&width=300&height=200
```

### Download Cached SVG
```bash
GET /cache/smiles_CCO_hash.svg
```

---

## 🔗 Worldwide Accessible Links

Once deployed, share cached SVG links like:
```
https://your-heroku-app.herokuapp.com/cache/smiles_CCO_72b4c84c77f80f62f6005fe6ec837e72.svg
```

Anyone can download the SVG for **24 hours** from anywhere in the world! 🌍

---

## 📁 Project Structure

```
Mol2chemfig/
├── MoleculeViewer/
│   ├── app/
│   │   ├── api.py              # Flask API endpoints
│   │   ├── chemistry.py        # RDKit integration
│   │   └── config.py           # Configuration
│   ├── templates/
│   │   └── index.html          # Web UI
│   ├── svg-cache/              # 24-hour cache storage
│   ├── .env                    # Configuration (localhost)
│   ├── .env.production         # Configuration (production)
│   ├── requirements.txt        # Python dependencies
│   ├── Procfile                # Heroku deployment config
│   └── run.py                  # Entry point
├── chem-extension/             # Chrome extension
├── .github/workflows/
│   └── deploy.yml              # GitHub Actions CI/CD
├── QUICK_START.md              # Deployment guide
└── DEPLOYMENT_GUIDE.md         # Detailed setup
```

---

## 🧪 Testing

### Local API Test
```bash
curl "http://192.168.1.4:5000/img/smiles?smiles=CCO"
```

### Chrome Extension
1. `chrome://extensions`
2. Enable Developer mode
3. Load unpacked folder: `chem-extension/`
4. Use in ChatGPT: `chem:acetone`

---

## 📚 Documentation

- [Quick Start Deployment](./QUICK_START.md) - Deploy in 5 minutes
- [Detailed Deployment Guide](./DEPLOYMENT_GUIDE.md) - Step-by-step setup
- [Architecture Overview](./ARCHITECTURE.md) - System design

---

## 🔐 Security

- ✅ CORS enabled for cross-domain requests
- ✅ Environment variables for sensitive config
- ✅ 24-hour cache expiry for cleanup
- ⚠️ Add authentication for production (see docs)

---

## 🚀 Production Deployment

### Option 1: Heroku (Recommended)
```bash
heroku create your-app-name
git push heroku main
heroku config:set PUBLIC_BASE_URL=https://your-app-name.herokuapp.com
```

### Option 2: Railway
```
https://railway.app → New Project → Deploy from GitHub
```

### Option 3: Vercel
```
https://vercel.com → Import Repository → Deploy
```

---

## 📊 Performance

- **Render time:** ~50-100ms per molecule
- **Cache hit:** 1-5ms (file lookup + download)
- **Cache size:** ~1-5KB per SVG
- **Concurrent users:** ~100+ on free tier

---

## 🐛 Troubleshooting

**Cache links not accessible?**
- Check `PUBLIC_BASE_URL` matches your backend domain
- Verify port is open (5000, 8080, etc.)

**Chemistry rendering fails?**
- Check RDKit installation: `python -c "from rdkit import Chem; print('✅ RDKit OK')"`
- Verify SMILES format is valid

**GitHub Pages not updating?**
- Check GitHub Actions tab for workflow status
- Clear browser cache (Ctrl+Shift+Delete)

---

## 🤝 Contributing

Fork, commit, and submit PRs! 

Areas for contribution:
- [ ] Database backend (replace file cache)
- [ ] User authentication
- [ ] Advanced rendering options
- [ ] Mobile app version
- [ ] 3D molecule viewer

---

## 📄 License

MIT License - feel free to use in personal/commercial projects!

---

## 📞 Support

- Issues: GitHub Issues tab
- Docs: See [DEPLOYMENT_GUIDE.md](./DEPLOYMENT_GUIDE.md)
- Email: your-email@example.com

---

## 🎯 Roadmap

- [ ] User authentication system
- [ ] Database storage (vs file-based)
- [ ] Advanced molecule search
- [ ] Batch SVG generation
- [ ] REST API with API keys
- [ ] WebAssembly rendering (client-side)
- [ ] Mobile app (React Native)
- [ ] 3D interactive viewer

---

**Made with ❤️ by Kapil** | [Star on GitHub ⭐](https://github.com/yourusername/Mol2chemfig)
