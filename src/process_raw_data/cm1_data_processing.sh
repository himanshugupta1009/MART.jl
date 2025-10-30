#!/bin/bash

set -euo pipefail

# -------------------------------
# CONFIGURATION (EDIT THESE)
# -------------------------------
SRC_ENDPOINT_ID="bf383737-9573-4010-9be0-9c454dc66239"           # Remote collection ID
DST_ENDPOINT_ID="c028965f-6e5d-11f0-bc79-0e9afee528db"            # Your Globus Connect Personal (local) collection ID
DST_DIR="/home/himanshu/Documents/Research/MART.jl/dataset/"         # Destination path as seen by the *destination endpoint*
OUTPUT_DIR="/media/storage/himanshu_storage/MART/Processed_CM1/"    # Where your Julia script writes processed files (POSIX path)
PROCESSOR="/home/himanshu/Documents/Research/MART.jl/src/process_raw_data/cm1_globus.jl"  # Path to your Julia script
FILE_LIST="/home/himanshu/Documents/Research/MART.jl/src/process_raw_data/cm1_file_names.txt"                           # One source-path per line (relative to SRC endpoint root)
LOG_FILE="/home/himanshu/Documents/Research/MART.jl/src/process_raw_data/processing_log.txt"

# Make sure local dirs exist (POSIX)
mkdir -p "$OUTPUT_DIR" "$DST_DIR"

log() { echo "[$(date '+%Y-%m-%d %H:%M:%S')] $*" | tee -a "$LOG_FILE"; }

# -------------------------------
# VALIDATION
# -------------------------------
command -v globus >/dev/null || { echo "ERROR: globus CLI not found"; exit 1; }
command -v /home/himanshu/julia-1.10.2/bin/julia >/dev/null || { echo "ERROR: julia not found"; exit 1; }

# -------------------------------
# PROCESS ONE FILE
# -------------------------------
process_one() {
  local src_path="$1"   # e.g., /remote/dir/cm1out_000019.nc
  local filename
  filename="$(basename "$src_path")"
  local dst_path="${DST_DIR}/${filename}"  # where the destination endpoint will place the file

  log "START: $src_path"
  # 1) Submit transfer and capture TASK_ID
  TASK_ID="$(
    globus transfer \
      "${SRC_ENDPOINT_ID}:${src_path}" \
      "${DST_ENDPOINT_ID}:${dst_path}" \
      --sync-level checksum \
      --format UNIX \
      --jmespath 'task_id'
  )" || { log "ERROR: Failed to submit transfer for $src_path"; return 1; }

  log "Submitted transfer (task_id=$TASK_ID). Waiting..."

  # 2) Wait for the task to finish and check final status
  globus task wait "$TASK_ID" >/dev/null || true
  STATUS="$(globus task show "$TASK_ID" --format UNIX --jmespath 'status' 2>/dev/null || echo FAILED)"

  if [[ "$STATUS" != "SUCCEEDED" ]]; then
    log "ERROR: Transfer failed for $src_path (task $TASK_ID, status=$STATUS)"
    return 1
  fi

  # 3) Sanity-check local file exists (POSIX path matches endpoint path here)
  if [[ ! -f "$dst_path" ]]; then
    log "ERROR: Transfer said SUCCEEDED but file not found at $dst_path"
    return 1
  fi

  # 4) Process with Julia
  log "Processing with Julia: $filename"
  cmd="/home/himanshu/julia-1.10.2/bin/julia \"$PROCESSOR\" \"$dst_path\" \"$OUTPUT_DIR\""
  echo "Running command: $cmd"
  if ! eval "$cmd" >>"$LOG_FILE" 2>&1; then
    log "ERROR: Julia processing failed for $filename"
    return 1
  fi

  # 5) Delete raw file
  rm -f "$dst_path"
  log "DONE: $filename (raw deleted)"
}

# -------------------------------
# MAIN LOOP (one file at a time)
# -------------------------------
echo "Looking for files in folder"
if [[ ! -f "$FILE_LIST" ]]; then
  echo "ERROR: FILE_LIST not found: $FILE_LIST" >&2
  exit 1
fi
echo "Looking for files in folder"

i=0
total=$(wc -l < "$FILE_LIST" | tr -d ' ')
echo "Total files are $total"

while IFS= read -r src_path || [[ -n "$src_path" ]]; do
     src_path="${src_path%%$'\r'}"  # remove trailing \r
     echo "DEBUG: src_path='$src_path'"

     echo "DEBUG: File $i of $total"

     [[ -z "$src_path" ]] && { echo "DEBUG: skip empty line"; continue; }

     echo "DEBUG: About to process '$src_path'"
     if ! process_one "$src_path"; then
         echo "DEBUG: process_one failed for '$src_path'"
     fi
     i=$((i+1))
done < "$FILE_LIST"

echo "All files processed. Check the log for details: $LOG_FILE"
