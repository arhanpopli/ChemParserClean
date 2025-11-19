# ChemParser - Vercel Deployment Guide

This guide will help you deploy ChemParser to Vercel.

## Prerequisites

1. A Vercel account (sign up at https://vercel.com)
2. Git repository with your code
3. Vercel CLI (optional, for local testing)

## Quick Deploy

### Option 1: Deploy via Vercel Dashboard (Recommended)

1. **Push your code to GitHub:**
   ```bash
   git add .
   git commit -m "Add Vercel deployment configuration"
   git push origin claude/fix-404-error-01BfWFJyL7NKcyQda2LWEgfU
   ```

2. **Go to Vercel Dashboard:**
   - Visit https://vercel.com/new
   - Click "Import Project"
   - Select your GitHub repository

3. **Configure the project:**
   - **Framework Preset:** Other
   - **Root Directory:** `./` (leave as default)
   - **Build Command:** Leave empty or use `echo "No build required"`
   - **Output Directory:** Leave empty
   - **Install Command:** `npm install`

4. **Environment Variables (Optional):**
   - Add any environment variables if needed
   - For example: `NODE_ENV=production`

5. **Deploy:**
   - Click "Deploy"
   - Wait for deployment to complete (usually 1-2 minutes)

### Option 2: Deploy via Vercel CLI

1. **Install Vercel CLI:**
   ```bash
   npm install -g vercel
   ```

2. **Login to Vercel:**
   ```bash
   vercel login
   ```

3. **Deploy:**
   ```bash
   cd /path/to/ChemParserClean
   vercel
   ```

4. **Follow the prompts:**
   - Set up and deploy: Yes
   - Which scope: Select your account
   - Link to existing project: No
   - Project name: chemparser (or your choice)
   - Directory: `./` (just press Enter)
   - Override settings: No

5. **Production deployment:**
   ```bash
   vercel --prod
   ```

## Project Structure

After setup, your Vercel deployment will have:

```
ChemParserClean/
├── api/                      # Serverless functions
│   ├── health.js            # Health check endpoint
│   └── img/
│       ├── smiles.js        # SMILES to SVG conversion
│       └── nomenclature.js  # Name to SVG conversion
├── index.html               # Landing page
├── vercel.json              # Vercel configuration
├── package.json             # Node dependencies
├── requirements.txt         # Python dependencies
└── .vercelignore           # Files to exclude from deployment
```

## API Endpoints

Once deployed, your API will be available at:

- **Base URL:** `https://your-project.vercel.app`
- **Health Check:** `GET /api/health`
- **SMILES to SVG:** `GET /api/img/smiles?smiles=CCO&width=300&height=200`
- **Nomenclature to SVG:** `GET /api/img/nomenclature?nomenclature=ethanol&width=300&height=200`

## Testing Your Deployment

1. **Test the health endpoint:**
   ```bash
   curl https://your-project.vercel.app/api/health
   ```

2. **Test SMILES conversion:**
   ```bash
   curl "https://your-project.vercel.app/api/img/smiles?smiles=c1ccccc1"
   ```

3. **Visit the landing page:**
   ```
   https://your-project.vercel.app
   ```

## Troubleshooting

### 404 Error - NOT_FOUND

**Problem:** Getting 404 errors when accessing endpoints.

**Solution:**
- Make sure you've pushed all files to the repository
- Check that `vercel.json` is in the root directory
- Verify that the `api/` folder contains the serverless functions
- Try redeploying: `vercel --prod --force`

### Function Timeout

**Problem:** API requests timing out after 10 seconds.

**Solution:**
- Vercel free tier has a 10-second limit for serverless functions
- Upgrade to Pro plan for 60-second limit
- Optimize your Python code for faster execution

### Missing Dependencies

**Problem:** Errors about missing npm or Python packages.

**Solution:**
- Ensure `package.json` lists all Node.js dependencies
- Ensure `requirements.txt` lists all Python dependencies
- Check Vercel build logs for specific errors

### Python/RDKit Issues

**Problem:** RDKit not available or Python errors.

**Solution:**
- RDKit is challenging to deploy on Vercel due to size limits
- The current implementation includes a fallback that shows SMILES text
- For full RDKit support, consider:
  - Using a separate service (Railway, Render, or Heroku)
  - Using Docker-based deployment
  - Integrating with external chemistry APIs

## Environment Variables

Set these in Vercel Dashboard > Settings > Environment Variables:

- `NODE_ENV=production` (optional)
- Add any custom API keys or configuration

## Custom Domain

1. Go to Vercel Dashboard > Your Project > Settings > Domains
2. Add your custom domain
3. Follow DNS configuration instructions
4. Wait for SSL certificate provisioning (automatic)

## Monitoring

- **View Logs:** Vercel Dashboard > Your Project > Logs
- **Analytics:** Vercel Dashboard > Your Project > Analytics
- **Function Stats:** View execution time and memory usage

## Limits (Free Tier)

- **Bandwidth:** 100GB/month
- **Function Execution:** 100GB-hours/month
- **Function Duration:** 10 seconds max
- **Deployments:** Unlimited
- **Team Members:** 1 (just you)

## Upgrading to Pro

For production use, consider Vercel Pro ($20/month):
- 60-second function timeout
- 1TB bandwidth
- Priority support
- Team collaboration

## Alternative Deployment Options

If Vercel doesn't meet your needs:

1. **Railway.app** - Better for Python/RDKit
2. **Render.com** - Free tier with Docker support
3. **Heroku** - Classic PaaS (no free tier)
4. **Fly.io** - Good for Docker deployments
5. **AWS Lambda** - More complex but powerful

## Support

- **Vercel Docs:** https://vercel.com/docs
- **Community:** https://github.com/vercel/vercel/discussions
- **Status:** https://www.vercel-status.com/

## Next Steps

1. **Test all endpoints** to ensure they work correctly
2. **Set up monitoring** to track usage and errors
3. **Configure custom domain** if needed
4. **Add authentication** for production use
5. **Implement rate limiting** to prevent abuse

## Important Notes

- The current Vercel setup provides basic functionality
- Full RDKit integration may require additional services
- Consider hybrid approach: Vercel for frontend, separate service for heavy chemistry processing
- Keep an eye on function execution time (10s limit on free tier)

---

**Deployment Status:** Ready for Vercel deployment
**Last Updated:** 2025-11-19
**Version:** 1.0.0
