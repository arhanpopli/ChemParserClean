# ✅ Mol2chemfig Docker Setup - COMPLETE

## Status: RUNNING ✅

Both Docker containers are now running successfully!

---

## 🚀 Access Points

### Frontend (Web UI)
- **URL**: http://localhost:8080
- **Port**: 8080 (maps to container 9000)
- **Container**: m2cf_frontend
- **Status**: ✅ Up and running

### Backend (API Server)
- **URL**: http://localhost:8000
- **Port**: 8000 (maps to container 5000)
- **Container**: m2cf_backend
- **Status**: ✅ Up and running

---

## 📦 Docker Images

- **Frontend**: `pychemist/m2cf_web_frontend:latest` (440MB)
- **Backend**: `pychemist/m2cf_web_backend:latest` (3.96GB)

Both images have been successfully pulled and are running.

---

## 🔧 Container Management

### View running containers:
```powershell
docker ps
```

### View all containers:
```powershell
docker ps -a
```

### Check container logs:
```powershell
docker logs m2cf_backend
docker logs m2cf_frontend
```

### Stop all services:
```powershell
docker-compose down
```

### Start services again:
```powershell
cd C:\Users\Kapil\Personal\PROJECTS\Mol2chemfig
docker-compose up -d
```

---

## 📋 What's Running

### m2cf_backend (Container)
- **Image**: pychemist/m2cf_web_backend:latest
- **Port**: 0.0.0.0:8000->5000/tcp
- **Status**: Up 31+ seconds
- **Function**: Python Flask API server for Mol2chemfig

### m2cf_frontend (Container)
- **Image**: pychemist/m2cf_web_frontend:latest
- **Port**: 0.0.0.0:8080->9000/tcp
- **Status**: Up 31+ seconds
- **Function**: Quasar web UI for Mol2chemfig

---

## 🛠️ Launcher Script

Created: `C:\Users\Kapil\Personal\PROJECTS\Mol2chemfig\start-docker.bat`

This batch file will:
1. Verify Docker is installed
2. Check Docker daemon is running
3. Start all services with docker-compose

Simply double-click to launch!

---

## 📊 Configuration (.env)

```
BACKEND_HOST_PORT=8000
BACKEND_CONTAINER_PORT=5000
FRONTEND_HOST_PORT=8080
FRONTEND_CONTAINER_PORT=9000
API_URL=http://localhost:8000
```

---

## ✅ Verification

To verify services are running:

```powershell
# Check containers
docker ps

# Check backend is responding
docker logs m2cf_backend -n 5

# Check frontend is responding
docker logs m2cf_frontend -n 5
```

---

## 🎯 Next Steps

1. **Access Frontend**: Open http://localhost:8080 in your browser
2. **Test API**: Backend is available at http://localhost:8000
3. **View Logs**: Use `docker logs <container_name>` to troubleshoot
4. **Stop Services**: Run `docker-compose down` when done

---

## 📝 Notes

- Both containers are set to restart unless stopped (`restart: unless-stopped`)
- Backend uses volume mapping for data persistence
- Frontend is served via Nginx with Quasar framework
- Debugger is enabled on the backend (visible in logs)

---

**Status**: ✅ All systems operational!

Services will continue running until you stop them with `docker-compose down`.
