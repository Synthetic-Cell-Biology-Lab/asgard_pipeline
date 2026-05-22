#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
BACKEND_DIR="$SCRIPT_DIR/backend"
FRONTEND_DIR="$SCRIPT_DIR/frontend"

BACKEND_PORT=8000
FRONTEND_PORT=5173

# ── Cleanup on exit ───────────────────────────────────────────────────────────
cleanup() {
    echo ""
    echo "Shutting down..."
    [[ -n "${BACKEND_PID:-}"  ]] && kill "$BACKEND_PID"  2>/dev/null
    [[ -n "${FRONTEND_PID:-}" ]] && kill "$FRONTEND_PID" 2>/dev/null
    wait 2>/dev/null
    echo "Done."
}
trap cleanup EXIT INT TERM

# ── Checks ────────────────────────────────────────────────────────────────────
if ! command -v uvicorn &>/dev/null; then
    echo "❌ uvicorn not found. Activate your conda/venv environment first."
    exit 1
fi

if ! command -v npm &>/dev/null; then
    echo "❌ npm not found."
    exit 1
fi

if [ ! -f "$BACKEND_DIR/app.py" ]; then
    echo "❌ Backend not found at $BACKEND_DIR/app.py"
    exit 1
fi

if [ ! -f "$FRONTEND_DIR/package.json" ]; then
    echo "❌ Frontend not found at $FRONTEND_DIR/package.json"
    exit 1
fi

# ── Backend ───────────────────────────────────────────────────────────────────
echo "🚀 Starting backend  →  http://localhost:${BACKEND_PORT}"
cd "$BACKEND_DIR"
uvicorn app:app --reload --port "$BACKEND_PORT" &
BACKEND_PID=$!

# ── Frontend ──────────────────────────────────────────────────────────────────
echo "🖥  Starting frontend →  http://localhost:${FRONTEND_PORT}"
cd "$FRONTEND_DIR"
npm run dev -- --port "$FRONTEND_PORT" &
FRONTEND_PID=$!

# ── Wait for backend to be ready ──────────────────────────────────────────────
echo "⏳ Waiting for backend..."
for i in $(seq 1 20); do
    if curl -sf "http://localhost:${BACKEND_PORT}/configs" &>/dev/null; then
        echo "✅ Backend ready."
        break
    fi
    sleep 0.5
done

echo ""
echo "=============================================="
echo "  Asgard Pipeline GUI"
echo "  Frontend : http://localhost:${FRONTEND_PORT}"
echo "  API      : http://localhost:${BACKEND_PORT}"
echo "  Press Ctrl+C to stop."
echo "=============================================="
echo ""

wait