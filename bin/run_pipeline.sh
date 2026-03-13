#!/usr/bin/env bash

set -euo pipefail

CONFIG_FILE=${1:-}

if [ -z "$CONFIG_FILE" ]; then
  echo "Usage: ./run_pipeline.sh <config.yaml>"
  exit 1
fi

if [ ! -f "$CONFIG_FILE" ]; then
  echo "Config file not found: $CONFIG_FILE"
  exit 1
fi

# -------------------------------
# YAML reader helper
# -------------------------------

read_config_value () {
python - <<EOF
import yaml
with open("$CONFIG_FILE") as f:
    cfg = yaml.safe_load(f) or {}
value = cfg.get("$1", "")
if value is None:
    value = ""
print(value)
EOF
}

PIPELINE=$(read_config_value pipeline)
PROTEIN=$(read_config_value protein_name)
LOG_DIR=$(read_config_value log_dir)
RUN_REASON=$(read_config_value reason)
CORES=$(read_config_value cores)
WORKING_DIR="$PWD"


# Defaults
CORES=${CORES:-1}
LOG_DIR=${LOG_DIR:-logs}

if [ -z "$PIPELINE" ]; then
  echo "❌ 'pipeline' not defined in config."
  exit 1
fi

if [ -z "$PROTEIN" ]; then
  echo "❌ 'protein_name' not defined in config."
  exit 1
fi

SNAKEFILE="${WORKING_DIR}/bin/pipelines/${PIPELINE}/Snakefile"

if [ ! -f "$SNAKEFILE" ]; then
  echo "❌ Snakefile not found for pipeline: $PIPELINE"
  echo "Expected at: $SNAKEFILE"
  exit 1
fi

# -------------------------------
# Create log filename
# -------------------------------

mkdir -p "$LOG_DIR"

RUN_ID=$(read_config_value run_id)

# Fallback to timestamp only if run_id not provided
if [ -z "$RUN_ID" ]; then
    RUN_ID=$(date +"%Y%m%d_%H%M%S")
fi

LOG_FILE="${LOG_DIR}/${PIPELINE}_${PROTEIN}_${RUN_ID}.log"

# -------------------------------
# Log header
# -------------------------------

{
echo "=============================================="
echo "🚀 Pipeline Execution Log"
echo "Pipeline: $PIPELINE"
echo "Protein: $PROTEIN"
echo "Config: $CONFIG_FILE"
echo "Start Time: $(date)"
echo "Cores: $CORES"
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
# Run Snakemake
# -------------------------------

snakemake \
  --snakefile "$SNAKEFILE" \
  --configfile "$CONFIG_FILE" \
  --config pipeline_log="$LOG_FILE" \
  --use-conda \
  --cores "$CORES" \
  --printshellcmds \
  --verbose \
  --conda-frontend mamba \
  2>&1 | tee -a "$LOG_FILE"


# -------------------------------
# Footer
# -------------------------------

{
echo ""
echo "=============================================="
echo "✅ Finished at $(date)"
echo "=============================================="
} | tee -a "$LOG_FILE"

echo "Log saved to: $LOG_FILE"
