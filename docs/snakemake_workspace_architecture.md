# Snakemake Workspace Architecture & Implementation Plan

This document defines the development plan for evolving ASGARD into a local computational biology workspace built on top of Snakemake.

## Core Goal

Build a system that can:

- browse filesystem outputs
- edit configs with structured forms
- launch Snakemake pipelines through wrappers
- stream logs live
- preview scientific artifacts
- track run history and provenance

## Architecture

```text
Frontend (React + Vite)
    ↓
FastAPI Backend
    ↓
Filesystem + Snakemake wrappers + Database
```

## Non-Negotiable Principles

1. **Filesystem-first foundation** before adding DB complexity.
2. **Never execute Snakemake directly** from API handlers; always call shell wrappers (e.g. `bin/run_pipeline.sh`).
3. **Protein-set-aware UX** should become first-class navigation, not just raw directory browsing.
4. **Artifact-aware previews** (tree/sequence/table/log/image/etc.) are required.

## Recommended Development Order

1. Filesystem browsing
2. File previews
3. Config editing
4. Pipeline launching
5. Live log streaming
6. Run persistence/reconnect
7. Protein-set workspace UI
8. Scientific viewers
9. DAG visualization
10. Database layer hardening

## Phase Plan

### Phase 1 — Core Filesystem Explorer (first priority)

**Backend**
- Add `GET /browse?path=`
- Validate and normalize paths against `BASE_DIR`
- Prevent traversal outside allowed root
- Return entries with directory/file type, size, extension, and modified time

**Frontend**
- Build explorer primitives:
  - `FileExplorer`
  - `FileRow`
  - `Breadcrumbs`
  - `PreviewPane`
- Add folder navigation, breadcrumbs, up navigation, file-type icons, and selected-file highlighting

### Phase 2 — File Preview System

**Backend**
- Add `GET /file?path=`
- Return text payload + metadata + MIME for text-like files
- For binary/image files return a download/asset URL

**Frontend preview routing**
- `.nwk/.treefile/.contree` → tree preview
- `.fasta/.aln/.afa` → sequence preview
- `.png/.svg` → image preview
- `.log` → log preview
- `.yaml/.yml` → config preview/editor
- `.tsv/.csv` → table preview

### Phase 3 — Config Management

**Backend**
- Extend existing config endpoints to return richer metadata:
  - tags
  - pipeline type
  - protein set
  - timestamps

**Frontend**
- Add `ConfigSidebar`, `ConfigEditor`, `ConfigSearch`
- Move from plain-text YAML editing to structured editing for:
  - strings
  - booleans
  - numbers
  - lists
  - nested dictionaries

### Phase 4 — Pipeline Execution API

Replace direct run streaming endpoint shape with run resource model:

- `POST /runs` (start run)
- `GET /runs/{id}` (status/metadata)
- `GET /runs/{id}/logs` (retrievable buffered logs)

Run object should include at minimum: `id`, `status`, `config`, `started_at`, and `pid`.

### Phase 5 — Live Log Streaming

- Add `GET /runs/{id}/stream` with SSE
- Merge stdout/stderr, add timestamps, support auto-scroll and downloadable logs
- Browser disconnect must not terminate run

### Phase 6 — Protein-Set Workspace UX

Promote biological organization in navigation:

- Sidebar root: protein sets (`Actin`, `Tubulin`, `FtsZ`, etc.)
- Workspace sections by artifact class (`Trees`, `Alignments`, `HMMER`, `Synteny`, `Structures`, `Logs`)

### Phase 7 — Scientific Artifact Viewers

- Tree viewer (PhyloCanvas/jsPhyloSVG or equivalent)
- Alignment viewer with coloring, consensus, scrolling
- Table viewer with virtualized rendering + sort/filter
- Synteny visualization later (e.g., clinker/SVG pan+zoom)

### Phase 8 — DAG Visualization

- Generate DAG via Snakemake
- Convert DOT output to SVG
- Render zoomable graph with rule status overlays

### Phase 9 — Database Layer

After filesystem + run model stabilize, add SQLite for:

- runs
- configs
- artifacts
- protein sets

### Phase 10 — Advanced Features

Later extensions:

- auth + multi-user labs
- queueing/HPC/SLURM submission
- scientific overlays and advanced analytics

## Suggested Project Structure

### Frontend

```text
src/
├── components/
│   ├── explorer/
│   ├── preview/
│   ├── runs/
│   ├── configs/
│   ├── workspace/
│   └── common/
├── hooks/
├── services/
│   ├── api.js
│   ├── filesystem.js
│   ├── runs.js
│   └── configs.js
├── pages/
├── App.jsx
└── App.css
```

### Backend

```text
backend/
├── app.py
├── routes/
│   ├── browse.py
│   ├── configs.py
│   ├── runs.py
│   └── previews.py
├── services/
│   ├── filesystem_service.py
│   ├── snakemake_service.py
│   ├── preview_service.py
│   └── run_service.py
├── models/
└── database/
```

## Product Framing

This should be treated as **a computational biology workspace built on Snakemake**, not merely a GUI wrapper around Snakemake commands.
