# Start Search Server Batch File

## Quick Start

Simply double-click:
```
start-search-server.bat
```

The server will start on **port 8001** and display:
```
✓ Starting server...

Datasets loaded successfully.
Search API Server running at http://localhost:8001/
```

## What It Does

1. **Kills any existing node processes** - Cleans up old instances
2. **Waits 2 seconds** - Ensures clean shutdown
3. **Starts the search server** - On port 8001

## How to Use

### Option 1: Double-click (Easiest)
1. Navigate to `C:\Users\Kapil\Personal\STUFF\Chemparser\`
2. Double-click `start-search-server.bat`
3. A command window opens and runs the server
4. Keep the window open while using the extension

### Option 2: Command Line
```powershell
cd C:\Users\Kapil\Personal\STUFF\Chemparser
.\start-search-server.bat
```

### Option 3: Add to Startup
1. Press `Win + R`
2. Type `shell:startup`
3. Create a shortcut to `start-search-server.bat` there
4. Server will start automatically when you log in

## Test the Server

Open PowerShell and run:
```powershell
Invoke-RestMethod -Uri "http://localhost:8001/search?q=caffeine"
```

You should see:
```
name      cid embed_url
----      --- ---------
Caffeine 2519 https://embed.molview.org/v1/?cid=2519
```

## Verify It's Working

The Chrome extension will:
- ✅ Load RCSB protein images
- ✅ Search for molecules
- ✅ Display 3D viewers
- ✅ Work without the PHP server on port 8000

## Stopping the Server

Simply close the command window or press `Ctrl+C`

## Advanced: Keep Server Running

To make the server restart automatically if it crashes, use a more advanced batch file:

```batch
:start
echo Starting search server...
node Molview\molview\search-server.js
echo Server crashed! Restarting in 5 seconds...
timeout /t 5 /nobreak
goto start
```

## Troubleshooting

### "Node not found" error
- Make sure Node.js is installed: `node --version`
- The batch file must be in `C:\Users\Kapil\Personal\STUFF\Chemparser\`

### "Address already in use" error
- Another process is using port 8001
- The batch file tries to kill it automatically
- If it still fails, manually kill it:
  ```powershell
  Get-Process -Name "node" | Stop-Process -Force
  ```

### Port 8001 not accessible
- Check if Windows Firewall is blocking it
- Try: `netstat -ano | findstr ":8001"`

## Files

- `start-search-server.bat` - Main batch file (this one)
- `Molview/molview/search-server.js` - The actual server
- `Molview/molview/src/datasets/` - Search data files

## Integration

This batch file is used by:
- Chrome extension (localhost:8001 searches)
- Test HTML files
- Any local application that needs chemical compound data

Make sure it's running whenever you use the extension!
