# 🧬 ASGARD Pipeline

A modular, reproducible bioinformatics workflow framework built using **Snakemake**, designed for structured protein set extraction, curation, and downstream analysis.

---

## 🚀 Overview

**ASGARD Pipeline** is a multi-pipeline framework where each biological workflow lives in its own modular directory while sharing:

* A unified executor (`run_pipeline.sh`)
* Config-driven execution
* Structured logging
* Conda-based reproducibility
* Human-in-the-loop review gates
* Clear separation between infrastructure and biological logic

The system is designed to scale as more pipelines are added.

---

## 🏗 Repository Structure

```
asgard_pipeline/
│
├── bin/
│   ├── run_pipeline.sh
│   ├── pipelines/
│   │   ├── protein_pipeline/
│   │   │   ├── Snakefile
│   │   │   └── rules/
│   │   └── ...
│   └── envs/
│
├── processes/
│   ├── protein_pipeline_IPS_ftsz.yaml
│   └── ...
│
├── database/
├── logs/
└── README.md
```

### Key Components

| Component         | Purpose                                        |
| ----------------- | ---------------------------------------------- |
| `run_pipeline.sh` | Unified executor for all pipelines             |
| `pipelines/`      | Modular workflow definitions                   |
| `processes/`      | Config files for individual runs               |
| `envs/`           | Conda environment YAML files                   |
| `logs/`           | Structured execution logs                      |
| `.snakemake/`     | Auto-generated workflow state and environments |

---

## ▶️ Running a Pipeline

Pipelines are executed via:

```bash
bash bin/run_pipeline.sh processes/<config.yaml>
```

The configuration file determines:

* Which pipeline to run
* Input files
* Target protein
* Run ID
* Log directory
* Execution cores

---

## 📄 Configuration File Structure

Example:

```yaml
pipeline: protein_pipeline
protein_name: ftsz
run_id: ftsz_test_run

log_dir: logs/
cores: 24

database: path/to/interpro.parquet
protein_file: path/to/protein.csv
fasta_file: path/to/all_sequences.fasta

search_string: ftsz
rstring: null
```

### Important Notes

* Exactly one of `search_string` or `rstring` must be set.
* `run_id` determines log continuity across reruns.
* If `run_id` is unchanged, logs append to the same file.

---

## 📜 Logging System

Each run generates a structured log:

```
<log_dir>/<pipeline>_<protein>_<run_id>.log
```

The log contains:

* Run metadata
* Config snapshot
* Execution timestamps
* Full Snakemake output
* Rule-level execution details

If a pipeline includes manual checkpoints, rerunning with the same `run_id` appends to the same log.

---

## 🧪 Conda Environments

* Managed automatically via `--use-conda`
* Stored under `.snakemake/conda/`
* Hash-based isolation ensures reproducibility
* Uses modern Conda solver (libmamba backend)

Environments are recreated only if:

* The YAML changes
* Dependencies change
* Channels change

---

## 🔍 Manual Review Gates

Some pipelines implement human-in-the-loop review steps.

Typical behavior:

1. Pipeline generates `.unr.fasta`
2. Execution pauses
3. User manually curates file
4. Save curated file as `.rev.fasta`
5. Create marker:

```bash
touch REVIEW_DONE.flag
```

6. Re-run pipeline to continue

This ensures controlled biological validation.

---

## 🛡 Design Principles

* Modular workflow separation
* Config-driven execution
* File-state driven logic
* Reproducible environments
* Structured logging
* Human validation checkpoints
* Scalable pipeline architecture

---


## 🗺 Snakemake Workspace Roadmap

A full architecture and phased implementation plan for the FastAPI + React workspace system is now documented in:

- `docs/snakemake_workspace_architecture.md`

This roadmap defines the development order, API evolution, preview strategy, run model, and long-term product framing for evolving ASGARD into a computational biology workspace.

---
## 📈 Future Extensions

Planned expansions may include:

* Alignment pipelines
* Phylogeny pipelines
* Annotation workflows
* Automated report generation
* Cluster/HPC profile integration
* Config schema validation
* Pipeline registry system

---

## 🧠 Philosophy

ASGARD Pipeline is designed not just to run workflows, but to:

* Ensure reproducibility
* Enforce structured execution
* Maintain traceability
* Enable clean expansion
* Support research-grade bioinformatics

---

## 📌 Pipeline-Specific Documentation

Each pipeline should contain its own `README.md` describing:

* Biological objective
* Input/output structure
* Special rules
* Manual steps
* Required config parameters

This keeps infrastructure documentation separate from biological workflow documentation.

---

## 🏁 Summary

ASGARD Pipeline provides a structured, extensible, and reproducible framework for bioinformatics workflows, designed for long-term maintainability and scalability.

