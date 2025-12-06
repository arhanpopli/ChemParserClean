const { spawn } = require('child_process');

const child = spawn('node', ['Molview/molview/search-server.js'], {
    cwd: 'C:\\Users\\Kapil\\Personal\\STUFF\\Chemparser',
    stdio: 'inherit'
});

child.on('error', (error) => {
    console.error(`Error: ${error.message}`);
});

child.on('exit', (code, signal) => {
    console.log(`Child process exited with code ${code} and signal ${signal}`);
});

// Keep the process alive
setInterval(() => {
    // Do nothing, just keep alive
}, 10000);
