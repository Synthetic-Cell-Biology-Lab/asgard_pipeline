#!/usr/bin/env bash

set -euo pipefail

CONFIG_FILE=${1:-}

# -------------------------------
# Resolve paths
# -------------------------------

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
BASE_DIR="$(dirname "$SCRIPT_DIR")"
MANIFEST="$BASE_DIR/processes/run_manifest.tsv"

cd "$BASE_DIR"

source .env

cd "$BASE_DIR"

# -------------------------------
# Validate config
# -------------------------------

if [ -z "$CONFIG_FILE" ]; then
  echo "Usage: ./run_pipeline.sh <config.yaml>"
  exit 1
fi

if [ ! -f "$CONFIG_FILE" ]; then
  echo "❌ Config file not found: $CONFIG_FILE"
  exit 1
fi


read_config_values () {
python3 - "$@" <<'EOF'
import sys, yaml

config_file = sys.argv[1]
keys = sys.argv[2:]

with open(config_file) as f:
    cfg = yaml.safe_load(f) or {}

for k in keys:
    value = cfg
    for part in k.split("."):
        if isinstance(value, dict):
            value = value.get(part)
        else:
            value = None
            break
    print(value if value is not None else "")
EOF
}

# -------------------------------
# Read config values
# -------------------------------

mapfile -t VALUES < <(read_config_values "$CONFIG_FILE" run.pipeline run.log_dir run.reason run.cores run.id run.parent_dir run.protein_name)

PIPELINE="${VALUES[0]}"
LOG_DIR="${VALUES[1]}"
RUN_REASON="${VALUES[2]}"
CORES="${VALUES[3]}"
RUN_ID="${VALUES[4]}"
PARENT_DIR="${VALUES[5]}"
PROTEIN="${VALUES[6]}"

WORKING_DIR="$PWD"

# -------------------------------
# Defaults
# -------------------------------

CORES=${CORES:-1}
LOG_DIR=${LOG_DIR:-logs}

if [ -z "$RUN_ID" ]; then
    RUN_ID=$(date +"%Y%m%d_%H%M%S")
fi

if [ -z "$PIPELINE" ]; then
  echo "❌ 'pipeline' not defined in config."
  exit 1
fi

# -------------------------------
# Snakefile path
# -------------------------------

SNAKEFILE="${WORKING_DIR}/bin/pipelines/${PIPELINE}/Snakefile"

if [ ! -f "$SNAKEFILE" ]; then
  echo "❌ Snakefile not found for pipeline: $PIPELINE"
  echo "Expected at: $SNAKEFILE"
  exit 1
fi

# -------------------------------
# Output structure
# -------------------------------

RESULT_DIR=${PARENT_DIR:-"$BASE_DIR/database/protein_sets/${PROTEIN}/${RUN_ID}"}

META_DIR="${RESULT_DIR}/metadata"

mkdir -p "$RESULT_DIR"
mkdir -p "$META_DIR"
mkdir -p "$LOG_DIR"

# -------------------------------
# Log + DAG paths
# -------------------------------

LOG_FILE="${META_DIR}/${PIPELINE}_${RUN_ID}.log"

DAG_DOT="${META_DIR}/${PIPELINE}_${RUN_ID}_dag.dot"
DAG_SVG="${META_DIR}/${PIPELINE}_${RUN_ID}_dag.svg"

# -------------------------------
# Manifest
# -------------------------------

if [ ! -f "$MANIFEST" ]; then
    printf "timestamp\trun_id\tpipeline\tprotein\tconfig\tfinal_config\tresult_dir\tcores\tstatus\tgit_commit\treason\n" > "$MANIFEST"
fi


# -------------------------------
# Preserve config
# -------------------------------

cp "$CONFIG_FILE" "${META_DIR}/config_used.yaml"

# -------------------------------
# Log header
# -------------------------------

{
echo "=============================================="
echo "🚀 Pipeline Execution Log"
echo "Pipeline: $PIPELINE"
echo "Config: $CONFIG_FILE"
echo "Start Time: $(date)"
echo "Cores: $CORES"
echo "Run ID: $RUN_ID"
echo "Result Directory: $RESULT_DIR"
echo "----------------------------------------------"
echo "Reason:"
echo "$RUN_REASON"
echo "=============================================="
echo ""
} >> "$LOG_FILE"

{
echo "=============================================="
echo "Config File"
echo "=============================================="
echo ""
cat "$CONFIG_FILE"
echo ""
echo "=============================================="
echo ""
} >> "$LOG_FILE"

# -------------------------------
# Conda setup
# -------------------------------

source /root/miniconda3/etc/profile.d/conda.sh

export CONDA_EXE="/root/miniconda3/bin/conda"
export PATH="/root/miniconda3/condabin:$PATH"

conda activate snakemake

# -------------------------------
# Generate DAG
# -------------------------------

echo "🧩 Generating DAG..."

snakemake \
  --snakefile "$SNAKEFILE" \
  --configfile "$CONFIG_FILE" \
  --config pipeline_log="$LOG_FILE" \
  --use-conda \
  --dag > "$DAG_DOT"

dot -Tsvg "$DAG_DOT" > "$DAG_SVG"

echo "✅ DAG saved to:"
echo "   DOT: $DAG_DOT"
echo "   SVG: $DAG_SVG"

# -------------------------------
# Run Snakemake
# -------------------------------

STATUS="SUCCESS"

if ! snakemake \
    --snakefile "$SNAKEFILE" \
    --configfile "$CONFIG_FILE" \
    --config pipeline_log="$LOG_FILE" \
    --use-conda \
    --rerun-incomplete \
    --cores "$CORES" \
    --printshellcmds \
    2>&1 | tee -a "$LOG_FILE"
then
    STATUS="FAILED"
fi

# -------------------------------
# update manifest
# -------------------------------

printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" \
    "$(date --iso-8601=seconds)" \
    "$RUN_ID" \
    "$PIPELINE" \
    "$PROTEIN" \
    "$CONFIG_FILE" \
    "${META_DIR}/config_used.yaml" \
    "$RESULT_DIR" \
    "$CORES" \
    "$STATUS" \
    "$RUN_REASON" \
>> "$MANIFEST"

if [ "$STATUS" = "FAILED" ]; then
    exit 1
fi



# -------------------------------
# Footer
# -------------------------------


{
echo ""
echo "=============================================="
echo "✅ Finished at $(date)"
echo "=============================================="
} | tee -a "$LOG_FILE"

echo ""
echo "📁 Result directory: $RESULT_DIR"
echo "📝 Log file: $LOG_FILE"
echo "🧩 DAG SVG: $DAG_SVG"