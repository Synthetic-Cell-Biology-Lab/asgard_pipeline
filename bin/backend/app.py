from fastapi import FastAPI, HTTPException
from fastapi.middleware.cors import CORSMiddleware
from fastapi.responses import StreamingResponse
from pydantic import BaseModel
from pathlib import Path
from typing import Dict, Any, List
import yaml
import subprocess
import json
import uuid
import threading
import time
import mimetypes
import re

app = FastAPI()

app.add_middleware(
    CORSMiddleware,
    allow_origins=["*"],
    allow_credentials=True,
    allow_methods=["*"],
    allow_headers=["*"],
)

# app.py is at asgard_pipeline/bin/backend/app.py
# parents[0] = backend/
# parents[1] = bin/
# parents[2] = asgard_pipeline/  ← BASE_DIR
BASE_DIR      = Path(__file__).resolve().parents[2]
DATABASE_DIR  = BASE_DIR / "database"
RUN_PIPELINE  = BASE_DIR / "bin" / "run_pipeline.sh"
DATABASE_DIR.mkdir(parents=True, exist_ok=True)

_state = {"configs_dir": BASE_DIR / "processes"}

# Templates live at asgard_pipeline/templates/configs/
TEMPLATES_DIR = BASE_DIR / "templates" / "configs"
TEMPLATES_DIR.mkdir(parents=True, exist_ok=True)

def get_configs_dir() -> Path:
    return _state["configs_dir"]

def set_configs_dir(path: Path):
    path.mkdir(parents=True, exist_ok=True)
    _state["configs_dir"] = path

set_configs_dir(get_configs_dir())

RUNS: Dict[str, Dict[str, Any]] = {}


# ── Pydantic models ───────────────────────────────────────────────────────────

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


# ── Helpers ───────────────────────────────────────────────────────────────────

def _safe_resolve(base: Path, relative_path: str) -> Path:
    candidate = (base / relative_path).resolve()
    if not str(candidate).startswith(str(base.resolve())):
        raise HTTPException(status_code=400, detail="Path escapes base directory")
    return candidate


def _collect_logs(run_id: str, process: subprocess.Popen):
    run = RUNS[run_id]
    for line in process.stdout:
        run["logs"].append({"ts": time.time(), "type": "stdout", "line": line.rstrip("\n")})
    process.wait()
    run["status"] = "succeeded" if process.returncode == 0 else "failed"
    run["returncode"] = process.returncode


def _parse_dot(dot_str: str) -> dict:
    """Parse snakemake --dag dot output into {nodes, edges}."""
    nodes = {}
    edges = []
    for m in re.finditer(r'(\d+)\[label\s*=\s*"([^"]+)"', dot_str):
        node_id = m.group(1)
        label = m.group(2).split('\\n')[0].strip()
        nodes[node_id] = {"id": node_id, "label": label, "rule": label, "status": "pending"}
    for m in re.finditer(r'(\d+)\s*->\s*(\d+)', dot_str):
        edges.append({"source": m.group(1), "target": m.group(2)})
    return {"nodes": list(nodes.values()), "edges": edges}


def _annotate_dag(dag: dict, logs: list) -> dict:
    """Colour DAG nodes by status based on log lines."""
    completed, running, failed = set(), set(), set()
    job_to_rule: Dict[str, str] = {}

    rule_pattern = r'([A-Za-z0-9_.-]+)'

    for entry in logs:
        line = (entry.get("line", "") or "").strip()
        if not line:
            continue

        # Rule started: "rule <name>:"
        m = re.search(rf'^rule {rule_pattern}:', line)
        if m:
            running.add(m.group(1))

        # Capture job id to rule mapping from execution lines.
        # Examples:
        # "jobid: 3"
        # "rule align_sequences:"
        m_jobid = re.search(r'jobid:\s*(\d+)', line)
        if m_jobid:
            current_jobid = m_jobid.group(1)
            # if same line contains rule
            m_rule_inline = re.search(rf'rule {rule_pattern}:', line)
            if m_rule_inline:
                rule_name = m_rule_inline.group(1)
                job_to_rule[current_jobid] = rule_name
                running.add(rule_name)

        # Failures: "Error in rule <name>:"
        m_fail = re.search(rf'Error in rule {rule_pattern}:', line)
        if m_fail:
            rule_name = m_fail.group(1)
            failed.add(rule_name)
            running.discard(rule_name)

        # Completed job: "Finished job 3."
        m_finished = re.search(r'Finished job\s+(\d+)', line)
        if m_finished:
            jobid = m_finished.group(1)
            rule_name = job_to_rule.get(jobid)
            if rule_name:
                completed.add(rule_name)
                running.discard(rule_name)

        # Pipeline-level completion fallback
        if "done" in line and "steps" in line:
            for rule_name in list(running):
                if rule_name not in failed:
                    completed.add(rule_name)
            running.clear()

    debug = {
        "completed_rules": sorted(completed),
        "running_rules": sorted(running),
        "failed_rules": sorted(failed),
    }

    unmatched = []
    for node in dag["nodes"]:
        rule = node.get("rule", node.get("label", "")).strip()
        if rule in failed:
            node["status"] = "error"
        elif rule in completed:
            node["status"] = "done"
        elif rule in running:
            node["status"] = "running"
        else:
            node["status"] = "pending"
            if rule and rule not in completed and rule not in failed:
                unmatched.append(rule)

    debug["unmatched_nodes"] = sorted(set(unmatched))
    dag["debug"] = debug
    return dag


def _list_yaml_files(directory: Path) -> List[ConfigFileInfo]:
    """Return all .yaml/.yml files in a directory, sorted newest first."""
    items = []
    for file in directory.iterdir():
        if file.suffix in (".yaml", ".yml"):
            stat = file.stat()
            items.append(ConfigFileInfo(name=file.name, size=stat.st_size, mtime=stat.st_mtime))
    items.sort(key=lambda x: x.mtime, reverse=True)
    return items


# ── Config endpoints ──────────────────────────────────────────────────────────

@app.post("/set-dir")
def set_dir(req: SetDirRequest):
    expanded = Path(req.dir).expanduser().resolve()
    set_configs_dir(expanded)
    return {"dir": str(get_configs_dir())}


@app.get("/configs", response_model=ConfigListResponse)
def list_configs():
    configs_dir = get_configs_dir()
    return {"configs": _list_yaml_files(configs_dir), "dir": str(configs_dir)}


@app.get("/configs/{filename}", response_model=ConfigResponse)
def get_config(filename: str):
    safe_name = Path(filename).name
    full_path = get_configs_dir() / safe_name
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
    full_path = get_configs_dir() / safe_name
    full_path.write_text(yaml.dump(req.params, sort_keys=False))
    return {"saved": safe_name}


# ── Template endpoints ────────────────────────────────────────────────────────

@app.get("/templates", response_model=ConfigListResponse)
def list_templates():
    return {"configs": _list_yaml_files(TEMPLATES_DIR), "dir": str(TEMPLATES_DIR)}


@app.get("/templates/{filename}", response_model=ConfigResponse)
def get_template(filename: str):
    safe_name = Path(filename).name
    full_path = TEMPLATES_DIR / safe_name
    if not full_path.exists():
        raise HTTPException(status_code=404, detail="Template not found")
    raw = full_path.read_text()
    parsed = yaml.safe_load(raw) or {}
    return {"name": safe_name, "raw": raw, "parsed": parsed}


# ── File browser endpoints ────────────────────────────────────────────────────

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
        entries.append({
            "name": item.name,
            "path": str(item.relative_to(BASE_DIR)),
            "type": "directory" if item.is_dir() else "file",
            "size": stat.st_size if item.is_file() else None,
            "mtime": stat.st_mtime,
            "extension": item.suffix.lower() if item.is_file() else None,
        })

    rel    = "" if target == BASE_DIR else str(target.relative_to(BASE_DIR))
    parent = None if target == BASE_DIR else str(target.parent.relative_to(BASE_DIR))
    return {"current_path": rel, "parent_path": parent, "entries": entries}


@app.get("/file")
def get_file(path: str):
    target = _safe_resolve(BASE_DIR, path)
    if not target.exists() or not target.is_file():
        raise HTTPException(404, "File not found")

    mime, _ = mimetypes.guess_type(str(target))
    mime = mime or "application/octet-stream"
    is_text = mime.startswith("text/") or target.suffix.lower() in {
        ".yaml", ".yml", ".log", ".nwk", ".fasta", ".fa", ".aln", ".afa", ".tsv", ".csv"
    }

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
    return StreamingResponse(
        open(target, "rb"),
        media_type=mimetypes.guess_type(str(target))[0] or "application/octet-stream"
    )


# ── Run endpoints ─────────────────────────────────────────────────────────────

@app.post("/runs")
def create_run(req: RunCreateRequest):
    if not RUN_PIPELINE.exists():
        raise HTTPException(500, f"run_pipeline.sh not found at {RUN_PIPELINE}")

    safe_name   = Path(req.config).name
    config_path = get_configs_dir() / safe_name
    if not config_path.exists():
        raise HTTPException(404, f"Config file not found: {config_path}")

    cmd  = ["bash", str(RUN_PIPELINE), str(config_path)]
    proc = subprocess.Popen(
        cmd,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        text=True,
        bufsize=1,
    )

    run_id = str(uuid.uuid4())
    RUNS[run_id] = {
        "id":          run_id,
        "status":      "running",
        "config":      safe_name,
        "config_path": str(config_path),
        "started_at":  time.time(),
        "pid":         proc.pid,
        "process":     proc,
        "logs":        [{"ts": time.time(), "type": "info", "line": " ".join(cmd)}],
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
    return {
        "run_id":      run_id,
        "offset":      offset,
        "next_offset": len(run["logs"]),
        "logs":        run["logs"][offset:],
    }


@app.get("/runs/{run_id}/stream")
def run_stream(run_id: str, offset: int = 0):
    if run_id not in RUNS:
        raise HTTPException(404, "Run not found")

    def event_stream():
        idx = offset
        while True:
            run  = RUNS.get(run_id)
            if not run:
                break
            logs = run["logs"]
            while idx < len(logs):
                yield f"data: {json.dumps(logs[idx])}\n\n"
                idx += 1
            if run["status"] != "running" and idx >= len(logs):
                yield f"data: {json.dumps({'type': 'exit', 'line': run['status'], 'ts': time.time()})}\n\n"
                break
            time.sleep(0.25)

    return StreamingResponse(event_stream(), media_type="text/event-stream")


# ── DAG endpoint ──────────────────────────────────────────────────────────────

def _load_config(config_path: str) -> dict:
    with open(config_path) as f:
        return yaml.safe_load(f) or {}


@app.get("/runs/{run_id}/dag")
def run_dag(run_id: str):
    run = RUNS.get(run_id)

    if not run:
        raise HTTPException(404, "Run not found")

    config_path = run.get("config_path")

    if not config_path:
        raise HTTPException(400, "No config path recorded for this run")

    try:
        # ------------------------------------------------------------------
        # Load config
        # ------------------------------------------------------------------

        cfg = _load_config(config_path)

        pipeline = cfg.get("pipeline")
        protein  = cfg.get("protein_name")
        run_name = cfg.get("run_id")

        if not pipeline:
            raise HTTPException(400, "Missing 'pipeline' in config")

        if not protein:
            raise HTTPException(400, "Missing 'protein' in config")

        if not run_name:
            raise HTTPException(400, "Missing 'run_id' in config")

        # ------------------------------------------------------------------
        # Resolve DAG path
        # ------------------------------------------------------------------

        result_dir = Path(
            cfg.get(
                "parent_dir",
                BASE_DIR
                / "database"
                / "protein_sets"
                / protein
                / run_name
            )
        )

        metadata_dir = result_dir / "metadata"

        dag_dot = (
            metadata_dir
            / f"{pipeline}_{run_name}_dag.dot"
        )

        if not dag_dot.exists():
            raise HTTPException(
                404,
                f"DAG file not found: {dag_dot}"
            )

        # ------------------------------------------------------------------
        # Load and parse DAG
        # ------------------------------------------------------------------

        dot_text = dag_dot.read_text()

        if not dot_text.strip():
            raise HTTPException(
                500,
                f"DAG file is empty: {dag_dot}"
            )

        dag = _parse_dot(dot_text)

        # ------------------------------------------------------------------
        # Annotate DAG with runtime statuses
        # ------------------------------------------------------------------

        dag = _annotate_dag(
            dag,
            run.get("logs", [])
        )

    except subprocess.TimeoutExpired:
        raise HTTPException(504, "DAG generation timed out")
        return {
            "run_id": run_id,
            "pipeline": pipeline,
            "protein": protein,
            "dag_path": str(dag_dot),
            "nodes": dag["nodes"],
            "edges": dag["edges"],
        }

    except HTTPException:
        raise

    except Exception as e:
        raise HTTPException(
            500,
            f"Failed to load DAG: {str(e)}"
        )
