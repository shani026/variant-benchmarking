#!/usr/bin/env bash
set -euo pipefail
project_dir="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd)"
build="${1:-}"
case "$build" in hg19|hg38) ;; *) echo 'Usage: code/run-nfcore.sh hg19|hg38 [additional Nextflow arguments]' >&2; exit 2 ;; esac
shift
cd "$project_dir"
exec nextflow run nf-core/variantbenchmarking -r 1.5.0 -profile docker \
    -params-file "code/config/params_${build}.yaml" \
    -c code/config/mac_local.config -work-dir "work/${build}" "$@"
