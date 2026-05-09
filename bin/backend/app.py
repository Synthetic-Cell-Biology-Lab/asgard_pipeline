from fastapi import FastAPI, HTTPException
from fastapi.middleware.cors import CORSMiddleware
from fastapi.responses import StreamingResponse
from pydantic import BaseModel
from pathlib import Path
from typing import Dict, Any, List, Optional
import yaml
import subprocess
import json
import uuid
import threading
import time
import mimetypes

app = FastAPI()

app.add_middleware(
    CORSMiddleware,
    allow_origins=["*"],
    allow_credentials=True,
    allow_methods=["*"],
    allow_headers=["*"],
)

BASE_DIR = Path(__file__).resolve().parents[2]
CONFIGS_DIR = BASE_DIR / "processes"
DATABASE_DIR = BASE_DIR / "database"
RUN_PIPELINE = BASE_DIR / "bin" / "run_pipeline.sh"

CONFIGS_DIR.mkdir(parents=True, exist_ok=True)
DATABASE_DIR.mkdir(parents=True, exist_ok=True)

RUNS: Dict[str, Dict[str, Any]] = {}


class SetDirRequest(BaseModel):
    dir: str


class ConfigFileInfo(BaseModel):
    name: str
    size: int
    mtime: float


class ConfigListResponse(BaseModel):
    configs: List[ConfigFileInfo]
    dir: str


class SaveConfigRequest(BaseModel):
    filename: str
    params: Dict[str, Any]


class ConfigResponse(BaseModel):
    name: str
    raw: str
    parsed: Dict[str, Any]


class RunCreateRequest(BaseModel):
    config: str


def _safe_resolve(base: Path, relative_path: str) -> Path:
    candidate = (base / relative_path).resolve()
    if not str(candidate).startswith(str(base.resolve())):
        raise HTTPException(status_code=400, detail="Path escapes base directory")
    return candidate


@app.post("/set-dir")
def set_dir(req: SetDirRequest):
    global CONFIGS_DIR
    expanded = Path(req.dir).expanduser().resolve()
    expanded.mkdir(parents=True, exist_ok=True)
    CONFIGS_DIR = expanded
    return {"dir": str(CONFIGS_DIR)}


@app.get("/configs", response_model=ConfigListResponse)
def list_configs():
    configs = []
    for file in CONFIGS_DIR.iterdir():
        if file.suffix in [".yaml", ".yml"]:
            stat = file.stat()
            configs.append(ConfigFileInfo(name=file.name, size=stat.st_size, mtime=stat.st_mtime))
    configs.sort(key=lambda x: x.mtime, reverse=True)
    return {"configs": configs, "dir": str(CONFIGS_DIR)}


@app.get("/configs/{filename}", response_model=ConfigResponse)
def get_config(filename: str):
    safe_name = Path(filename).name
    full_path = CONFIGS_DIR / safe_name
    if not full_path.exists():
        raise HTTPException(status_code=404, detail="Not found")
    raw = full_path.read_text()
    parsed = yaml.safe_load(raw) or {}
    return {"name": safe_name, "raw": raw, "parsed": parsed}


@app.post("/configs")
def save_config(req: SaveConfigRequest):
    safe_name = Path(req.filename).name
    if not safe_name.endswith((".yaml", ".yml")):
        raise HTTPException(status_code=400, detail="Filename must end in .yaml or .yml")
    full_path = CONFIGS_DIR / safe_name
    full_path.write_text(yaml.dump(req.params, sort_keys=False))
    return {"saved": safe_name}


@app.get("/browse")
def browse(path: str = ""):
    target = _safe_resolve(BASE_DIR, path or "")
    if not target.exists():
        raise HTTPException(404, "Path not found")
    if not target.is_dir():
        raise HTTPException(400, "Not a directory")

    entries = []
    for item in sorted(target.iterdir(), key=lambda x: (not x.is_dir(), x.name.lower())):
        stat = item.stat()
        entries.append(
            {
                "name": item.name,
                "path": str(item.relative_to(BASE_DIR)),
                "type": "directory" if item.is_dir() else "file",
                "size": stat.st_size if item.is_file() else None,
                "mtime": stat.st_mtime,
                "extension": item.suffix.lower() if item.is_file() else None,
            }
        )

    rel = "" if target == BASE_DIR else str(target.relative_to(BASE_DIR))
    parent = None if target == BASE_DIR else str(target.parent.relative_to(BASE_DIR))
    return {"current_path": rel, "parent_path": parent, "entries": entries}


@app.get("/file")
def get_file(path: str):
    target = _safe_resolve(BASE_DIR, path)
    if not target.exists() or not target.is_file():
        raise HTTPException(404, "File not found")

    mime, _ = mimetypes.guess_type(str(target))
    mime = mime or "application/octet-stream"
    is_text = mime.startswith("text/") or target.suffix.lower() in {".yaml", ".yml", ".log", ".nwk", ".fasta", ".fa", ".aln", ".afa", ".tsv", ".csv"}

    meta = {
        "name": target.name,
        "path": str(target.relative_to(BASE_DIR)),
        "size": target.stat().st_size,
        "mtime": target.stat().st_mtime,
        "mime": mime,
        "extension": target.suffix.lower(),
    }

    if is_text:
        return {"kind": "text", "metadata": meta, "content": target.read_text(errors="replace")}

    return {"kind": "binary", "metadata": meta, "download_url": f"/download?path={path}"}


@app.get("/download")
def download(path: str):
    target = _safe_resolve(BASE_DIR, path)
    if not target.exists() or not target.is_file():
        raise HTTPException(404, "File not found")
    return StreamingResponse(open(target, "rb"), media_type=mimetypes.guess_type(str(target))[0] or "application/octet-stream")


def _collect_logs(run_id: str, process: subprocess.Popen):
    run = RUNS[run_id]
    for line in process.stdout:
        run["logs"].append({"ts": time.time(), "type": "stdout", "line": line.rstrip("\n")})
    process.wait()
    run["status"] = "succeeded" if process.returncode == 0 else "failed"
    run["returncode"] = process.returncode


@app.post("/runs")
def create_run(req: RunCreateRequest):
    safe_name = Path(req.config).name
    config_path = CONFIGS_DIR / safe_name
    if not config_path.exists():
        raise HTTPException(404, "Config file not found")

    cmd = ["bash", str(RUN_PIPELINE), str(config_path)]
    proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True, bufsize=1)

    run_id = str(uuid.uuid4())
    RUNS[run_id] = {
        "id": run_id,
        "status": "running",
        "config": safe_name,
        "started_at": time.time(),
        "pid": proc.pid,
        "process": proc,
        "logs": [{"ts": time.time(), "type": "info", "line": " ".join(cmd)}],
    }

    threading.Thread(target=_collect_logs, args=(run_id, proc), daemon=True).start()
    return {k: v for k, v in RUNS[run_id].items() if k not in {"process", "logs"}}


@app.get("/runs/{run_id}")
def run_status(run_id: str):
    run = RUNS.get(run_id)
    if not run:
        raise HTTPException(404, "Run not found")
    return {k: v for k, v in run.items() if k not in {"process", "logs"}}


@app.get("/runs/{run_id}/logs")
def run_logs(run_id: str, offset: int = 0):
    run = RUNS.get(run_id)
    if not run:
        raise HTTPException(404, "Run not found")
    return {"run_id": run_id, "offset": offset, "next_offset": len(run["logs"]), "logs": run["logs"][offset:]}


@app.get("/runs/{run_id}/stream")
def run_stream(run_id: str, offset: int = 0):
    if run_id not in RUNS:
        raise HTTPException(404, "Run not found")

    def event_stream():
        idx = offset
        while True:
            run = RUNS.get(run_id)
            if not run:
                break
            logs = run["logs"]
            while idx < len(logs):
                yield f"data: {json.dumps(logs[idx])}\n\n"
                idx += 1
            if run["status"] != "running" and idx >= len(logs):
                yield f"data: {json.dumps({'type':'exit','line':run['status'],'ts':time.time()})}\n\n"
                break
            time.sleep(0.25)

    return StreamingResponse(event_stream(), media_type="text/event-stream")
