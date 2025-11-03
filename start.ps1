# Start the compose stack from this folder
Set-Location -Path "$PSScriptRoot"

# Ensure .env is loaded by docker compose
docker compose up -d
