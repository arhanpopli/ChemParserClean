/**
 * ChemParser Launcher Control Server
 * Node.js backend that allows the HTML launcher to start/stop servers
 * Port: 3000
 */

const express = require('express');
const cors = require('cors');
const { exec, spawn } = require('child_process');
const path = require('path');
const fs = require('fs');

const app = express();
const PORT = 3000;

app.use(cors());
app.use(express.json());

// Track running processes
const runningProcesses = {};

// Serve the launcher HTML
app.use(express.static(__dirname));

// Get status of all servers
app.get('/api/status', async (req, res) => {
  const status = {};
  
  // Check each server
  const servers = [
    { id: 'pubchem', port: 5002, url: 'http://localhost:5002/health' },
    { id: 'moleculeviewer', port: 5000, url: 'http://localhost:5000/health' },
    { id: 'mol2chemfig', port: 5001, url: 'http://localhost:5001/health' }
  ];

  for (const server of servers) {
    try {
      const response = await fetch(server.url, { signal: AbortSignal.timeout(2000) });
      status[server.id] = {
        running: response.ok,
        port: server.port
      };
    } catch {
      status[server.id] = {
        running: false,
        port: server.port
      };
    }
  }

  res.json(status);
});

// Start PubChem server
app.post('/api/start/pubchem', (req, res) => {
  console.log('🚀 Starting PubChem Server...');
  
  const serverPath = path.join(__dirname, 'MoleculeViewer', 'pubchem');
  const scriptPath = path.join(serverPath, 'server.js');

  if (runningProcesses.pubchem) {
    return res.json({ success: false, message: 'PubChem server already running' });
  }

  // Kill any existing node processes on port 5002
  exec('taskkill /F /IM node.exe /FI "WINDOWTITLE eq PubChem*" 2>nul', () => {
    setTimeout(() => {
      const process = spawn('node', [scriptPath], {
        cwd: serverPath,
        detached: false,
        windowsHide: false
      });

      runningProcesses.pubchem = process;

      process.stdout.on('data', (data) => {
        console.log(`[PubChem] ${data.toString().trim()}`);
      });

      process.stderr.on('data', (data) => {
        console.error(`[PubChem Error] ${data.toString().trim()}`);
      });

      process.on('exit', (code) => {
        console.log(`[PubChem] Process exited with code ${code}`);
        delete runningProcesses.pubchem;
      });

      res.json({ success: true, message: 'PubChem server started on port 5002' });
    }, 1000);
  });
});

// Stop PubChem server
app.post('/api/stop/pubchem', (req, res) => {
  console.log('⏹️ Stopping PubChem Server...');
  
  if (runningProcesses.pubchem) {
    runningProcesses.pubchem.kill();
    delete runningProcesses.pubchem;
  }

  // Also kill any node processes on port 5002
  exec('netstat -ano | findstr :5002', (err, stdout) => {
    if (stdout) {
      const lines = stdout.split('\n');
      lines.forEach(line => {
        const match = line.match(/LISTENING\s+(\d+)/);
        if (match) {
          const pid = match[1];
          exec(`taskkill /F /PID ${pid}`, () => {
            console.log(`Killed process ${pid}`);
          });
        }
      });
    }
  });

  res.json({ success: true, message: 'PubChem server stopped' });
});

// Start MoleculeViewer server
app.post('/api/start/moleculeviewer', (req, res) => {
  console.log('🚀 Starting MoleculeViewer Server...');
  
  const serverPath = path.join(__dirname, 'MoleculeViewer');
  const scriptPath = path.join(serverPath, 'server.js');

  if (runningProcesses.moleculeviewer) {
    return res.json({ success: false, message: 'MoleculeViewer already running' });
  }

  const process = spawn('node', [scriptPath], {
    cwd: serverPath,
    detached: false
  });

  runningProcesses.moleculeviewer = process;

  process.stdout.on('data', (data) => {
    console.log(`[MoleculeViewer] ${data.toString().trim()}`);
  });

  res.json({ success: true, message: 'MoleculeViewer started on port 5000' });
});

// Stop MoleculeViewer
app.post('/api/stop/moleculeviewer', (req, res) => {
  if (runningProcesses.moleculeviewer) {
    runningProcesses.moleculeviewer.kill();
    delete runningProcesses.moleculeviewer;
  }

  res.json({ success: true, message: 'MoleculeViewer stopped' });
});

// Start ALL servers
app.post('/api/start/all', (req, res) => {
  console.log('🚀 Starting ALL Servers...');
  
  // Start PubChem
  const pubchemPath = path.join(__dirname, 'MoleculeViewer', 'pubchem');
  const pubchemScript = path.join(pubchemPath, 'server.js');
  
  if (!runningProcesses.pubchem) {
    const pubchemProcess = spawn('node', [pubchemScript], {
      cwd: pubchemPath,
      detached: false
    });

    runningProcesses.pubchem = pubchemProcess;

    pubchemProcess.stdout.on('data', (data) => {
      console.log(`[PubChem] ${data.toString().trim()}`);
    });
  }

  res.json({ success: true, message: 'All servers started' });
});

// Stop ALL servers
app.post('/api/stop/all', (req, res) => {
  console.log('⏹️ Stopping ALL Servers...');
  
  Object.keys(runningProcesses).forEach(key => {
    if (runningProcesses[key]) {
      runningProcesses[key].kill();
      delete runningProcesses[key];
    }
  });

  // Force kill any remaining node/python processes
  exec('taskkill /F /IM node.exe 2>nul', () => {
    console.log('Killed all node processes');
  });

  res.json({ success: true, message: 'All servers stopped' });
});

// View logs
app.get('/api/logs/:serverId', (req, res) => {
  const { serverId } = req.params;
  
  // Return placeholder logs (in real implementation, read from log files)
  res.json({
    logs: [
      { time: new Date().toISOString(), level: 'info', message: `${serverId} server running` }
    ]
  });
});

// Health check
app.get('/health', (req, res) => {
  res.json({
    status: 'ok',
    service: 'ChemParser Launcher Control',
    port: PORT,
    runningServers: Object.keys(runningProcesses)
  });
});

// Start the control server
app.listen(PORT, () => {
  console.log('\n' + '='.repeat(70));
  console.log('🎮 CHEMPARSER LAUNCHER CONTROL SERVER');
  console.log('='.repeat(70));
  console.log(`✓ Control Panel: http://localhost:${PORT}/launcher.html`);
  console.log(`✓ API Endpoint:  http://localhost:${PORT}/api/`);
  console.log('='.repeat(70));
  console.log('\nThis server allows the HTML launcher to start/stop other servers.');
  console.log('Keep this window open!\n');
});

// Cleanup on exit
process.on('SIGINT', () => {
  console.log('\n🛑 Shutting down launcher control server...');
  
  Object.keys(runningProcesses).forEach(key => {
    if (runningProcesses[key]) {
      runningProcesses[key].kill();
    }
  });

  process.exit(0);
});
