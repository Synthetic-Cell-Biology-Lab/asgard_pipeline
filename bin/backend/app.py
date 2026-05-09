from fastapi import FastAPI, HTTPException
from fastapi.middleware.cors import CORSMiddleware
from fastapi.responses import StreamingResponse
from pydantic import BaseModel
from pathlib import Path
from typing import Dict, Any, List
import yaml
import subprocess
import os
import json

app = FastAPI()

# ── CORS ─────────────────────────────────────────────────────

app.add_middleware(
    CORSMiddleware,
    allow_origins=["*"],
    allow_credentials=True,
    allow_methods=["*"],
    allow_headers=["*"],
)

# ── CONFIG ───────────────────────────────────────────────────
BASE_DIR = Path(__file__).resolve().parents[2]

CONFIGS_DIR = BASE_DIR / "processes"

DATABASE = BASE_DIR/ "database"


CONFIGS_DIR.mkdir(exist_ok=True)

RUN_PIPELINE = (
    Path(__file__).resolve().parents[2]
    / "bin"
    / "run_pipeline.sh"
)


cmd = [
    "bash",
    str(RUN_PIPELINE),
    "{configFile}"
]

CONFIGS_DIR.mkdir(parents=True, exist_ok=True)

# ── PYDANTIC MODELS ──────────────────────────────────────────

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


# ── ROUTES ───────────────────────────────────────────────────

@app.post("/set-dir")
def set_dir(req: SetDirRequest):

    global CONFIGS_DIR

    expanded = Path(req.dir).expanduser()

    expanded.mkdir(parents=True, exist_ok=True)

    CONFIGS_DIR = expanded

    return {"dir": str(CONFIGS_DIR)}


@app.get("/configs", response_model=ConfigListResponse)
def list_configs():

    configs = []

    for file in CONFIGS_DIR.iterdir():

        if file.suffix in [".yaml", ".yml"]:

            stat = file.stat()

            configs.append(
                ConfigFileInfo(
                    name=file.name,
                    size=stat.st_size,
                    mtime=stat.st_mtime
                )
            )

    configs.sort(key=lambda x: x.mtime, reverse=True)

    return {
        "configs": configs,
        "dir": str(CONFIGS_DIR)
    }


@app.get("/configs/{filename}", response_model=ConfigResponse)
def get_config(filename: str):

    safe_name = Path(filename).name

    full_path = CONFIGS_DIR / safe_name

    if not full_path.exists():
        raise HTTPException(status_code=404, detail="Not found")

    raw = full_path.read_text()

    parsed = yaml.safe_load(raw)

    return {
        "name": safe_name,
        "raw": raw,
        "parsed": parsed
    }


@app.post("/configs")
def save_config(req: SaveConfigRequest):

    safe_name = Path(req.filename).name

    if not safe_name.endswith((".yaml", ".yml")):
        raise HTTPException(
            status_code=400,
            detail="Filename must end in .yaml or .yml"
        )

    full_path = CONFIGS_DIR / safe_name

    yaml_str = yaml.dump(
        req.params,
        sort_keys=False
    )

    full_path.write_text(yaml_str)

    return {"saved": safe_name}


@app.get("/run")
def run_pipeline(config: str):

    safe_name = Path(config).name

    config_path = CONFIGS_DIR / safe_name

    if not config_path.exists():
        raise HTTPException(
            status_code=404,
            detail="Config file not found"
        )

    cmd = [
        "bash",
        str(RUN_PIPELINE),
        str(config_path)
    ]

    def event_stream():

        yield (
            f"data: "
            f"{json.dumps({'type': 'info', 'data': ' '.join(cmd)})}\n\n"
        )

        proc = subprocess.Popen(
            cmd,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            text=True,
            bufsize=1
        )

        for line in proc.stdout:

            yield (
                f"data: "
                f"{json.dumps({'type': 'stdout', 'data': line})}\n\n"
            )

        proc.wait()

        if proc.returncode == 0:
            msg = "✓ Pipeline finished successfully."
        else:
            msg = f"✗ Exited with code {proc.returncode}"

        yield (
            f"data: "
            f"{json.dumps({'type': 'exit', 'data': msg})}\n\n"
        )

    return StreamingResponse(
        event_stream(),
        media_type="text/event-stream"
    )


@app.get("/browse")
def browse(path: str = DATABASE):

    target = (BASE_DIR / path).resolve()

    if not target.exists():
        raise HTTPException(404, "Path not found")

    if not target.is_dir():
        raise HTTPException(400, "Not a directory")

    entries = []

    for item in sorted(target.iterdir()):

        entries.append({
            "name": item.name,
            "type": "directory" if item.is_dir() else "file",
            "size": item.stat().st_size if item.is_file() else None
        })

    return {
        "current_path": str(target.relative_to(BASE_DIR)),
        "entries": entries
    }