#!/usr/bin/env bash
set -euo pipefail

OUT="analysis/hg38/tables/hg38_svtype_counts.csv"
mkdir -p analysis/hg38/tables

echo "source,caller,svtype,count" > "$OUT"

count_svtypes () {
    local source="$1"
    local caller="$2"
    local vcf="$3"

    echo "Counting $source / $caller: $vcf" >&2

    bcftools query -f '%INFO/SVTYPE\n' "$vcf" \
      | awk -v source="$source" -v caller="$caller" '
          $1 != "." && $1 != "" {
              count[$1]++
          }
          END {
              for (sv in count) {
                  print source "," caller "," sv "," count[sv]
              }
          }' >> "$OUT"
}

count_svtypes "truth"   "HG002_truth" "Truth_vcfs/hg38/HG002_GRCh38_v5.0q_stvar.vcf.gz"
count_svtypes "caller"  "delly"       "caller_vcfs/HG002_delly_lr.hg38.renamed.vcf.gz"
count_svtypes "caller"  "dragen"      "caller_vcfs/HG002_dragen.sv.hg38.renamed.vcf.gz"
count_svtypes "caller"  "sniffles"    "caller_vcfs/HG002_sniffles.hg38.renamed.vcf.gz"

echo "Wrote $OUT"
