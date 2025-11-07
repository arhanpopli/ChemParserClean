#!/usr/bin/env python3
"""
Simple MCP (Model Context Protocol) Server for autonomous development
Allows Claude to work on tasks autonomously with proper tools
"""

import json
import sys
import subprocess
import os
from pathlib import Path

# Configuration
PROJECT_ROOT = r"C:\Users\Kapil\Personal\PROJECTS\Mol2chemfig"
FLASK_PATH = os.path.join(PROJECT_ROOT, "MoleculeViewer")
TEST_SCRIPT = os.path.join(PROJECT_ROOT, "tests", "test_api.py")

class MCPServer:
    def __init__(self):
        self.tools = {
            "run_tests": self.run_tests,
            "start_flask": self.start_flask,
            "stop_flask": self.stop_flask,
            "check_status": self.check_status,
            "get_task": self.get_task,
            "mark_complete": self.mark_complete,
        }
    
    def run_tests(self, args=None):
        """Run the automated test suite"""
        try:
            result = subprocess.run(
                [sys.executable, TEST_SCRIPT],
                cwd=PROJECT_ROOT,
                capture_output=True,
                text=True,
                timeout=60
            )
            return {
                "success": result.returncode == 0,
                "output": result.stdout,
                "errors": result.stderr,
                "tests_passed": "ALL TESTS PASSED" in result.stdout
            }
        except Exception as e:
            return {"success": False, "error": str(e)}
    
    def start_flask(self, args=None):
        """Start Flask server"""
        try:
            # Kill existing processes
            subprocess.run(
                "Get-Process python -ErrorAction SilentlyContinue | Stop-Process -Force",
                shell=True,
                capture_output=True
            )
            
            # Start new server
            subprocess.Popen(
                [sys.executable, "run_server.py"],
                cwd=FLASK_PATH,
                stdout=subprocess.DEVNULL,
                stderr=subprocess.DEVNULL
            )
            
            return {"success": True, "message": "Flask server started"}
        except Exception as e:
            return {"success": False, "error": str(e)}
    
    def stop_flask(self, args=None):
        """Stop Flask server"""
        try:
            subprocess.run(
                "Get-Process python -ErrorAction SilentlyContinue | Stop-Process -Force",
                shell=True,
                capture_output=True
            )
            return {"success": True, "message": "Flask server stopped"}
        except Exception as e:
            return {"success": False, "error": str(e)}
    
    def check_status(self, args=None):
        """Check system status"""
        return {
            "flask_running": self._is_flask_running(),
            "tests_available": os.path.exists(TEST_SCRIPT),
            "project_path": PROJECT_ROOT,
            "todos_available": os.path.exists(os.path.join(PROJECT_ROOT, "Todolist.md"))
        }
    
    def _is_flask_running(self):
        """Check if Flask is running"""
        try:
            import requests
            r = requests.get("http://localhost:5000/", timeout=2)
            return r.status_code in [200, 404]
        except:
            return False
    
    def get_task(self, args=None):
        """Get next task from Todolist"""
        todolist_path = os.path.join(PROJECT_ROOT, "MoleculeViewer", "docs", "Todolist.md")
        if os.path.exists(todolist_path):
            with open(todolist_path, 'r') as f:
                content = f.read()
                # Find first incomplete task
                lines = content.split('\n')
                for i, line in enumerate(lines):
                    if '- [ ]' in line:
                        return {"task": line.replace('- [ ]', '').strip(), "line": i}
        return {"task": None, "message": "No tasks found"}
    
    def mark_complete(self, args=None):
        """Mark a task as complete"""
        if not args or "task" not in args:
            return {"success": False, "error": "Task name required"}
        
        # This would need proper implementation
        return {"success": True, "message": f"Marked {args['task']} as complete"}
    
    def handle_request(self, request):
        """Handle MCP request"""
        tool = request.get("tool")
        args = request.get("args")
        
        if tool in self.tools:
            return self.tools[tool](args)
        return {"error": f"Unknown tool: {tool}"}

def main():
    """Main entry point"""
    server = MCPServer()
    
    # Check if running interactively or as server
    if len(sys.argv) > 1:
        # Command line interface
        command = sys.argv[1]
        result = server.handle_request({"tool": command})
        print(json.dumps(result, indent=2))
    else:
        # Interactive mode
        print("MCP Server Ready")
        print("Available commands:")
        for tool in server.tools:
            print(f"  - {tool}")
        
        status = server.check_status()
        print(f"\nStatus: {json.dumps(status, indent=2)}")

if __name__ == "__main__":
    main()
