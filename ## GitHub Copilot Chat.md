## GitHub Copilot Chat

- Extension Version: 0.31.5 (prod)
- VS Code: vscode/1.104.0
- OS: Windows

## Network

User Settings:
```json
  "github.copilot.advanced.debug.useElectronFetcher": true,
  "github.copilot.advanced.debug.useNodeFetcher": false,
  "github.copilot.advanced.debug.useNodeFetchFetcher": true
```

Connecting to https://api.github.com:
- DNS ipv4 Lookup: 20.207.73.85 (14 ms)
- DNS ipv6 Lookup: Error (902 ms): getaddrinfo ENOTFOUND api.github.com
- Proxy URL: None (1 ms)
- Electron fetch (configured): HTTP 200 (125 ms)
- Node.js https: HTTP 200 (157 ms)
- Node.js fetch: HTTP 200 (397 ms)

Connecting to https://api.githubcopilot.com/_ping:
- DNS ipv4 Lookup: 140.82.114.21 (11 ms)
- DNS ipv6 Lookup: Error (12 ms): getaddrinfo ENOTFOUND api.githubcopilot.com
- Proxy URL: None (27 ms)
- Electron fetch (configured): HTTP 200 (859 ms)
- Node.js https: HTTP 200 (871 ms)
- Node.js fetch: HTTP 200 (885 ms)

## Documentation

In corporate networks: [Troubleshooting firewall settings for GitHub Copilot](https://docs.github.com/en/copilot/troubleshooting-github-copilot/troubleshooting-firewall-settings-for-github-copilot).