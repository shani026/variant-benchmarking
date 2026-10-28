from pathlib import Path
from collections import Counter
import subprocess
import shutil
import json
import datetime
import math
import re
import html
import base64
import mimetypes

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

# ============================================================
# hg19-only project paths
# ============================================================

ROOT = Path(__file__).resolve().parents[2]

OUT = ROOT / "reports" / "hg19"
TABLES = OUT / "tables"
PLOTS = OUT / "plots"
HTML = OUT / "html"
LOGS = OUT / "logs"
SELF = OUT / "self_survivor"
PLAIN = SELF / "plain_vcfs"

for d in [TABLES, PLOTS, HTML, LOGS, SELF, PLAIN]:
    d.mkdir(parents=True, exist_ok=True)

STRUCT = ROOT / "benchmark-runs/nfcore/hg19/structural"
PIPELINE_INFO = ROOT / "benchmark-runs/nfcore/hg19/pipeline_info"

# ------------------------------------------------------------
# hg19 summary metric tables
# ------------------------------------------------------------

SUMMARY_TABLES = {
    "Truvari": STRUCT / "summary/tables/truvari/truvari.summary.csv",
    "SVbenchmark": STRUCT / "summary/tables/svbenchmark/svbenchmark.summary.csv",
    "RTGtools": STRUCT / "summary/tables/rtgtools/rtgtools.summary.csv",
}

# ------------------------------------------------------------
# hg19 summary comparison CSVs
# ------------------------------------------------------------

COMPARISON_CSVS = {
    "Truvari": {
        "FN": STRUCT / "summary/comparisons/truvari/truvari.FN.csv",
        "FP": STRUCT / "summary/comparisons/truvari/truvari.FP.csv",
        "TP_base": STRUCT / "summary/comparisons/truvari/truvari.TP_base.csv",
        "TP_comp": STRUCT / "summary/comparisons/truvari/truvari.TP_comp.csv",
    },
    "SVbenchmark": {
        "FN": STRUCT / "summary/comparisons/svbenchmark/svbenchmark.FN.csv",
        "FP": STRUCT / "summary/comparisons/svbenchmark/svbenchmark.FP.csv",
        "TP_base": STRUCT / "summary/comparisons/svbenchmark/svbenchmark.TP_base.csv",
        "TP_comp": STRUCT / "summary/comparisons/svbenchmark/svbenchmark.TP_comp.csv",
    },
    "RTGtools": {
        "FN": STRUCT / "summary/comparisons/rtgtools/rtgtools.FN.csv",
        "FP": STRUCT / "summary/comparisons/rtgtools/rtgtools.FP.csv",
        "TP_base": STRUCT / "summary/comparisons/rtgtools/rtgtools.TP_base.csv",
        "TP_comp": STRUCT / "summary/comparisons/rtgtools/rtgtools.TP_comp.csv",
    },
}

# ------------------------------------------------------------
# hg19 raw truth and caller VCFs
# ------------------------------------------------------------

RAW_VCFS = {
    "Truth": ROOT / "Truth_vcfs/hg19/HG002_GRCh37_v5.0q_stvar.vcf.gz",
    "Delly": ROOT / "caller_vcfs/DELLY.hg19.norm_recall_HG002.merged__filt_DUPtoINS.renamed.vcf.gz",
    "DRAGEN": ROOT / "caller_vcfs/HG002_dragen.sv.hg19.renamed.vcf.gz",
    "Sniffles": ROOT / "caller_vcfs/HG002_sniffles.hg19.renamed.vcf.gz",
}

# ------------------------------------------------------------
# hg19 prepared VCFs from nf-core preprocessing
# ------------------------------------------------------------

PREPARED_VCFS = {
    "Truth": STRUCT / "HG002/preprocess/HG002_GRCh37_v5.0q_stvar.reheader.filtermissing.vcf.gz",
    "Delly": STRUCT / "HG002_delly_hg19/preprocess/DELLY.hg19.norm_recall_HG002.merged__filt_DUPtoINS.renamed.reheader.filtermissing.vcf.gz",
    "DRAGEN": STRUCT / "HG002_dragen_hg19/preprocess/HG002_dragen.sv.hg19.renamed.reheader.filtermissing.vcf.gz",
    "Sniffles": STRUCT / "HG002_sniffles_hg19/preprocess/HG002_sniffles.hg19.renamed.reheader.filtermissing.vcf.gz",
}

# ------------------------------------------------------------
# hg19 bcftools and nf-core SURVIVOR stats
# ------------------------------------------------------------

BCFTOOLS_STATS = {
    "Truth": STRUCT / "HG002/stats/bcftools/HG002.bcftools_stats.txt",
    "Delly": STRUCT / "HG002_delly_hg19/stats/bcftools/HG002_delly_hg19.delly.bcftools_stats.txt",
    "DRAGEN": STRUCT / "HG002_dragen_hg19/stats/bcftools/HG002_dragen_hg19.dragen.bcftools_stats.txt",
    "Sniffles": STRUCT / "HG002_sniffles_hg19/stats/bcftools/HG002_sniffles_hg19.sniffles.bcftools_stats.txt",
}

NFCORE_SURVIVOR_STATS = {
    "Truth": STRUCT / "HG002/stats/survivor/HG002_mqc.stats",
    "Delly": STRUCT / "HG002_delly_hg19/stats/survivor/HG002_delly_hg19.delly_mqc.stats",
    "DRAGEN": STRUCT / "HG002_dragen_hg19/stats/survivor/HG002_dragen_hg19.dragen_mqc.stats",
    "Sniffles": STRUCT / "HG002_sniffles_hg19/stats/survivor/HG002_sniffles_hg19.sniffles_mqc.stats",
}

# ------------------------------------------------------------
# hg19 benchmark event VCFs
# ------------------------------------------------------------

BENCHMARK_VCFS = {
    "Truvari": {
        "Delly": {
            "FN": STRUCT / "HG002_delly_hg19/benchmarks/truvari/HG002_delly_hg19.HG002.delly.fn.vcf.gz",
            "FP": STRUCT / "HG002_delly_hg19/benchmarks/truvari/HG002_delly_hg19.HG002.delly.fp.vcf.gz",
            "TP_base": STRUCT / "HG002_delly_hg19/benchmarks/truvari/HG002_delly_hg19.HG002.delly.tp-base.vcf.gz",
            "TP_comp": STRUCT / "HG002_delly_hg19/benchmarks/truvari/HG002_delly_hg19.HG002.delly.tp-comp.vcf.gz",
        },
        "DRAGEN": {
            "FN": STRUCT / "HG002_dragen_hg19/benchmarks/truvari/HG002_dragen_hg19.HG002.dragen.fn.vcf.gz",
            "FP": STRUCT / "HG002_dragen_hg19/benchmarks/truvari/HG002_dragen_hg19.HG002.dragen.fp.vcf.gz",
            "TP_base": STRUCT / "HG002_dragen_hg19/benchmarks/truvari/HG002_dragen_hg19.HG002.dragen.tp-base.vcf.gz",
            "TP_comp": STRUCT / "HG002_dragen_hg19/benchmarks/truvari/HG002_dragen_hg19.HG002.dragen.tp-comp.vcf.gz",
        },
        "Sniffles": {
            "FN": STRUCT / "HG002_sniffles_hg19/benchmarks/truvari/HG002_sniffles_hg19.HG002.sniffles.fn.vcf.gz",
            "FP": STRUCT / "HG002_sniffles_hg19/benchmarks/truvari/HG002_sniffles_hg19.HG002.sniffles.fp.vcf.gz",
            "TP_base": STRUCT / "HG002_sniffles_hg19/benchmarks/truvari/HG002_sniffles_hg19.HG002.sniffles.tp-base.vcf.gz",
            "TP_comp": STRUCT / "HG002_sniffles_hg19/benchmarks/truvari/HG002_sniffles_hg19.HG002.sniffles.tp-comp.vcf.gz",
        },
    },
    "SVbenchmark": {
        "Delly": {
            "FN": STRUCT / "HG002_delly_hg19/benchmarks/svanalyzer/HG002_delly_hg19.HG002.delly.falsenegatives.vcf.gz",
            "FP": STRUCT / "HG002_delly_hg19/benchmarks/svanalyzer/HG002_delly_hg19.HG002.delly.falsepositives.vcf.gz",
        },
        "DRAGEN": {
            "FN": STRUCT / "HG002_dragen_hg19/benchmarks/svanalyzer/HG002_dragen_hg19.HG002.dragen.falsenegatives.vcf.gz",
            "FP": STRUCT / "HG002_dragen_hg19/benchmarks/svanalyzer/HG002_dragen_hg19.HG002.dragen.falsepositives.vcf.gz",
        },
        "Sniffles": {
            "FN": STRUCT / "HG002_sniffles_hg19/benchmarks/svanalyzer/HG002_sniffles_hg19.HG002.sniffles.falsenegatives.vcf.gz",
            "FP": STRUCT / "HG002_sniffles_hg19/benchmarks/svanalyzer/HG002_sniffles_hg19.HG002.sniffles.falsepositives.vcf.gz",
        },
    },
    "RTGtools": {
        "Delly": {
            "FN": STRUCT / "HG002_delly_hg19/benchmarks/rtgtools/HG002_delly_hg19.HG002.delly.fn.vcf.gz",
            "FP": STRUCT / "HG002_delly_hg19/benchmarks/rtgtools/HG002_delly_hg19.HG002.delly.fp.vcf.gz",
            "TP_base": STRUCT / "HG002_delly_hg19/benchmarks/rtgtools/HG002_delly_hg19.HG002.delly.tp-baseline.vcf.gz",
            "TP_comp": STRUCT / "HG002_delly_hg19/benchmarks/rtgtools/HG002_delly_hg19.HG002.delly.tp.vcf.gz",
        },
        "DRAGEN": {
            "FN": STRUCT / "HG002_dragen_hg19/benchmarks/rtgtools/HG002_dragen_hg19.HG002.dragen.fn.vcf.gz",
            "FP": STRUCT / "HG002_dragen_hg19/benchmarks/rtgtools/HG002_dragen_hg19.HG002.dragen.fp.vcf.gz",
            "TP_base": STRUCT / "HG002_dragen_hg19/benchmarks/rtgtools/HG002_dragen_hg19.HG002.dragen.tp-baseline.vcf.gz",
            "TP_comp": STRUCT / "HG002_dragen_hg19/benchmarks/rtgtools/HG002_dragen_hg19.HG002.dragen.tp.vcf.gz",
        },
        "Sniffles": {
            "FN": STRUCT / "HG002_sniffles_hg19/benchmarks/rtgtools/HG002_sniffles_hg19.HG002.sniffles.fn.vcf.gz",
            "FP": STRUCT / "HG002_sniffles_hg19/benchmarks/rtgtools/HG002_sniffles_hg19.HG002.sniffles.fp.vcf.gz",
            "TP_base": STRUCT / "HG002_sniffles_hg19/benchmarks/rtgtools/HG002_sniffles_hg19.HG002.sniffles.tp-baseline.vcf.gz",
            "TP_comp": STRUCT / "HG002_sniffles_hg19/benchmarks/rtgtools/HG002_sniffles_hg19.HG002.sniffles.tp.vcf.gz",
        },
    },
}

MULTIQC = STRUCT / "multiqc"
MULTIQC_PLOTS = MULTIQC / "multiqc_plots/png"
EXISTING_SUMMARY_PLOTS = STRUCT / "summary/plots"

CALLER_ORDER = ["Delly", "DRAGEN", "Sniffles"]
METHOD_ORDER = ["Truvari", "SVbenchmark", "RTGtools"]
STATUS_ORDER = ["TP_base", "TP_comp", "FP", "FN"]
SVLEN_ORDER = [
    "<50 bp",
    "50-99 bp",
    "100-499 bp",
    "500-999 bp",
    "1-4.9 kb",
    "5-9.9 kb",
    "10-49.9 kb",
    ">=50 kb",
    "NA/BND/UNKNOWN",
]

# ============================================================
# Helpers
# ============================================================

def run_cmd(cmd):
    p = subprocess.run(cmd, text=True, capture_output=True)
    if p.returncode != 0:
        return None, p.stderr.strip(), p.returncode
    return p.stdout, "", p.returncode

def savefig(name):
    plt.tight_layout()
    plt.savefig(PLOTS / name, dpi=220, bbox_inches="tight")
    plt.close()

def write_csv(df, name):
    out = TABLES / name
    df.to_csv(out, index=False)
    return out

def normalize_caller(x):
    s = str(x).lower()
    if "delly" in s:
        return "Delly"
    if "dragen" in s:
        return "DRAGEN"
    if "sniffles" in s:
        return "Sniffles"
    return str(x)

def normalize_method(x):
    s = str(x).lower()
    if "truvari" in s:
        return "Truvari"
    if "svbenchmark" in s or "svanalyzer" in s:
        return "SVbenchmark"
    if "rtg" in s:
        return "RTGtools"
    return str(x)

def parse_length(pos, end, svlen):
    if svlen and svlen not in [".", "NA", ""]:
        raw = str(svlen).split(",")[0]
        try:
            return abs(int(float(raw)))
        except Exception:
            pass
    if end and end not in [".", "NA", ""]:
        try:
            return abs(int(float(end)) - int(float(pos)))
        except Exception:
            pass
    return None

def svlen_bin(length):
    if length is None:
        return "NA/BND/UNKNOWN"
    try:
        v = abs(float(length))
    except Exception:
        return "NA/BND/UNKNOWN"
    if math.isnan(v):
        return "NA/BND/UNKNOWN"
    if v < 50:
        return "<50 bp"
    if v < 100:
        return "50-99 bp"
    if v < 500:
        return "100-499 bp"
    if v < 1000:
        return "500-999 bp"
    if v < 5000:
        return "1-4.9 kb"
    if v < 10000:
        return "5-9.9 kb"
    if v < 50000:
        return "10-49.9 kb"
    return ">=50 kb"

def slug(x):
    return re.sub(r"[^A-Za-z0-9]+", "_", str(x)).strip("_")

def df_to_html(df, max_rows=150):
    if df is None or df.empty:
        return "<p>No rows found.</p>"
    d = df.copy()
    note = ""
    if len(d) > max_rows:
        d = d.head(max_rows)
        note = f"<p class='small-note'>Showing first {max_rows} rows. Full table is saved in the tables folder.</p>"
    return note + d.to_html(index=False, classes="tbl", border=0, escape=True, float_format=lambda x: f"{x:.6g}")

def img_data_uri(path):
    path = Path(path)
    if not path.exists():
        return ""
    mime = mimetypes.guess_type(path)[0] or "image/png"
    encoded = base64.b64encode(path.read_bytes()).decode("ascii")
    return f"data:{mime};base64,{encoded}"

def plot_card(filename, title):
    path = PLOTS / filename
    if not path.exists():
        return ""
    return f"""
    <div class="plot-card">
      <h3>{html.escape(title)}</h3>
      <img src="{img_data_uri(path)}" alt="{html.escape(title)}">
    </div>
    """

def plot_section(title, items):
    cards = []
    for filename, caption in items:
        c = plot_card(filename, caption)
        if c:
            cards.append(c)
    if not cards:
        return ""
    return f"""
    <section>
      <h2>{html.escape(title)}</h2>
      <div class="plot-grid">
        {''.join(cards)}
      </div>
    </section>
    """

def plot_grouped_bar(df, index, columns, values, title, ylabel, filename, ylim=None, logy=False, stacked=False):
    if df is None or df.empty:
        return
    pivot = df.pivot_table(index=index, columns=columns, values=values, aggfunc="sum", fill_value=0, observed=False)
    if pivot.empty:
        return
    ax = pivot.plot(kind="bar", figsize=(10.5, 5.7), stacked=stacked, logy=logy)
    ax.set_title(title)
    ax.set_xlabel(index)
    ax.set_ylabel(ylabel)
    if ylim:
        ax.set_ylim(*ylim)
    ax.legend(title=columns, bbox_to_anchor=(1.02, 1), loc="upper left")
    plt.xticks(rotation=35, ha="right")
    savefig(filename)

def plot_single_bar(df, x, y, title, ylabel, filename, logy=False, ylim=None):
    if df is None or df.empty:
        return
    ax = df.set_index(x)[y].plot(kind="bar", figsize=(10, 5), logy=logy)
    ax.set_title(title)
    ax.set_xlabel(x)
    ax.set_ylabel(ylabel)
    if ylim:
        ax.set_ylim(*ylim)
    plt.xticks(rotation=35, ha="right")
    savefig(filename)

def plot_heatmap(matrix_df, title, filename, value_label):
    if matrix_df is None or matrix_df.empty:
        return
    data = matrix_df.values.astype(float)
    fig, ax = plt.subplots(figsize=(7.5, 4.8))
    im = ax.imshow(data, aspect="auto")
    ax.set_xticks(np.arange(len(matrix_df.columns)))
    ax.set_xticklabels(matrix_df.columns)
    ax.set_yticks(np.arange(len(matrix_df.index)))
    ax.set_yticklabels(matrix_df.index)
    for i in range(data.shape[0]):
        for j in range(data.shape[1]):
            ax.text(j, i, f"{data[i, j]:.4f}", ha="center", va="center", fontsize=8)
    ax.set_title(title)
    fig.colorbar(im, ax=ax, label=value_label)
    savefig(filename)

def plot_scatter_metrics(metrics, x, y, filename):
    if metrics is None or metrics.empty:
        return
    fig, ax = plt.subplots(figsize=(8, 5.5))
    for _, r in metrics.iterrows():
        ax.scatter(r[x], r[y], s=80)
        ax.annotate(
            f"{r['caller']} · {r['method']}",
            (r[x], r[y]),
            xytext=(5, 5),
            textcoords="offset points",
            fontsize=8,
        )
    ax.set_title(f"hg19 {y} vs {x}")
    ax.set_xlabel(x)
    ax.set_ylabel(y)
    ax.grid(True, alpha=0.3)
    savefig(filename)

# ============================================================
# 1. Load benchmark summary metrics
# ============================================================

def load_summary_table(path, method):
    df = pd.read_csv(path)
    df["method"] = method
    df = df.rename(columns={
        "Caller": "caller",
        "Tool": "tool",
        "File": "source_file",
        "Threshold": "threshold",
        "TP_base": "tp_base",
        "TP_comp": "tp_comp",
        "FP": "fp",
        "FN": "fn",
        "Precision": "precision",
        "Recall": "recall",
        "F1": "f1",
        "F1 Score": "f1",
    })

    for c in ["tool", "source_file", "threshold"]:
        if c not in df.columns:
            df[c] = ""
    for c in ["tp_base", "tp_comp", "fp", "fn", "precision", "recall", "f1"]:
        if c not in df.columns:
            df[c] = 0

    df["caller"] = df["caller"].apply(normalize_caller)
    df["method"] = df["method"].apply(normalize_method)

    for c in ["tp_base", "tp_comp", "fp", "fn"]:
        df[c] = pd.to_numeric(df[c], errors="coerce").fillna(0).astype(int)
    for c in ["precision", "recall", "f1"]:
        df[c] = pd.to_numeric(df[c], errors="coerce").fillna(0.0)

    return df[[
        "method", "caller", "tool", "source_file", "threshold",
        "tp_base", "tp_comp", "fp", "fn", "precision", "recall", "f1"
    ]]

metric_frames = []
missing_summary_rows = []

for method, path in SUMMARY_TABLES.items():
    if path.exists():
        metric_frames.append(load_summary_table(path, method))
    else:
        missing_summary_rows.append({"method": method, "expected_file": str(path), "status": "missing"})

metrics = pd.concat(metric_frames, ignore_index=True) if metric_frames else pd.DataFrame()

if not metrics.empty:
    metrics["method"] = pd.Categorical(metrics["method"], categories=METHOD_ORDER, ordered=True)
    metrics["caller"] = pd.Categorical(metrics["caller"], categories=CALLER_ORDER, ordered=True)
    metrics = metrics.sort_values(["method", "caller"])
    metrics["caller_tool"] = metrics["caller"].astype(str) + " · " + metrics["method"].astype(str)

write_csv(metrics, "benchmark_summary_metrics_all_tools.csv")
write_csv(pd.DataFrame(missing_summary_rows), "missing_summary_tables.csv")

mean_by_caller = pd.DataFrame()
if not metrics.empty:
    mean_by_caller = (
        metrics
        .groupby("caller", observed=False)[["tp_base", "tp_comp", "fp", "fn", "precision", "recall", "f1"]]
        .mean(numeric_only=True)
        .reset_index()
    )
write_csv(mean_by_caller, "mean_metrics_by_caller_across_tools.csv")

method_delta_rows = []
if not metrics.empty:
    for caller in CALLER_ORDER:
        sub = metrics[metrics["caller"].astype(str) == caller]
        for a, b in [("SVbenchmark", "Truvari"), ("RTGtools", "Truvari"), ("RTGtools", "SVbenchmark")]:
            aa = sub[sub["method"].astype(str) == a]
            bb = sub[sub["method"].astype(str) == b]
            if not aa.empty and not bb.empty:
                method_delta_rows.append({
                    "caller": caller,
                    "comparison": f"{a} minus {b}",
                    "delta_precision": float(aa["precision"].iloc[0]) - float(bb["precision"].iloc[0]),
                    "delta_recall": float(aa["recall"].iloc[0]) - float(bb["recall"].iloc[0]),
                    "delta_f1": float(aa["f1"].iloc[0]) - float(bb["f1"].iloc[0]),
                    "delta_tp_base": int(aa["tp_base"].iloc[0]) - int(bb["tp_base"].iloc[0]),
                    "delta_tp_comp": int(aa["tp_comp"].iloc[0]) - int(bb["tp_comp"].iloc[0]),
                    "delta_fp": int(aa["fp"].iloc[0]) - int(bb["fp"].iloc[0]),
                    "delta_fn": int(aa["fn"].iloc[0]) - int(bb["fn"].iloc[0]),
                })

method_delta_df = pd.DataFrame(method_delta_rows)
write_csv(method_delta_df, "benchmark_tool_metric_differences.csv")

# ============================================================
# 2. Parse bcftools stats
# ============================================================

def parse_bcftools_stats(path, source):
    rows = []
    if not path.exists():
        return rows

    wanted = {
        "number of records": "records",
        "number of SNPs": "snps",
        "number of MNPs": "mnps",
        "number of indels": "indels",
        "number of others": "others",
        "number of multiallelic sites": "multiallelic_sites",
        "number of multiallelic SNP sites": "multiallelic_snp_sites",
        "number of transitions": "transitions",
        "number of transversions": "transversions",
    }

    for line in path.read_text(errors="replace").splitlines():
        if not line.startswith("SN"):
            continue
        parts = line.split("\t")
        if len(parts) < 4:
            continue
        label = parts[2].strip().rstrip(":")
        value = parts[3].strip()
        if label in wanted:
            try:
                value_num = int(float(value))
            except Exception:
                value_num = value
            rows.append({
                "source": source,
                "metric": wanted[label],
                "value": value_num,
            })
    return rows

bcftools_rows = []
for source, path in BCFTOOLS_STATS.items():
    bcftools_rows.extend(parse_bcftools_stats(path, source))

bcftools_stats_df = pd.DataFrame(bcftools_rows)
write_csv(bcftools_stats_df, "bcftools_stats_summary_parsed.csv")

bcftools_wide = pd.DataFrame()
if not bcftools_stats_df.empty:
    bcftools_wide = (
        bcftools_stats_df
        .pivot_table(index="source", columns="metric", values="value", aggfunc="first", fill_value=0)
        .reset_index()
    )
write_csv(bcftools_wide, "bcftools_stats_summary_wide.csv")

# ============================================================
# 3. VCF summary parser using bcftools query
# ============================================================

def stream_vcf_summary(vcf_path, source, stage, method="", caller="", status=""):
    vcf_path = Path(vcf_path)

    result = {
        "source": source,
        "stage": stage,
        "method": method,
        "caller": caller,
        "status": status,
        "vcf": str(vcf_path),
        "exists": vcf_path.exists(),
        "total_records": 0,
        "filter_pass_records": 0,
        "filter_nonpass_records": 0,
        "query_error": "",
        "filter_counts": Counter(),
        "svtype_counts": Counter(),
        "svlen_bin_counts": Counter(),
        "chrom_counts": Counter(),
    }

    if not vcf_path.exists():
        result["query_error"] = "missing file"
        return result

    if shutil.which("bcftools") is None:
        result["query_error"] = "bcftools not found in PATH"
        return result

    fmt = "%CHROM\\t%POS\\t%INFO/END\\t%INFO/SVTYPE\\t%INFO/SVLEN\\t%FILTER\\n"
    cmd = ["bcftools", "query", "-u", "-f", fmt, str(vcf_path)]

    p = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True, bufsize=1)

    try:
        for line in p.stdout:
            parts = line.rstrip("\n").split("\t")
            if len(parts) != 6:
                continue

            chrom, pos, end, svtype, svlen, filt = parts
            svtype = svtype if svtype and svtype != "." else "UNKNOWN"
            filt = filt if filt and filt != "." else "UNKNOWN"
            length = parse_length(pos, end, svlen)
            bin_label = svlen_bin(length)

            result["total_records"] += 1

            if filt == "PASS":
                result["filter_pass_records"] += 1
            else:
                result["filter_nonpass_records"] += 1

            result["filter_counts"][filt] += 1
            result["svtype_counts"][svtype] += 1
            result["svlen_bin_counts"][bin_label] += 1
            result["chrom_counts"][chrom] += 1

        stderr = p.stderr.read()
        rc = p.wait()
        if rc != 0:
            result["query_error"] = stderr.strip()
    except Exception as e:
        result["query_error"] = str(e)
        try:
            p.kill()
        except Exception:
            pass

    return result

def summarize_vcf_map(vcf_map, stage):
    summaries = []
    for source, path in vcf_map.items():
        summaries.append(stream_vcf_summary(path, source=source, stage=stage))
    return summaries

input_summaries = summarize_vcf_map(RAW_VCFS, "raw") + summarize_vcf_map(PREPARED_VCFS, "prepared")

input_totals_rows = []
input_filter_rows = []
input_svtype_rows = []
input_svlen_rows = []
input_chrom_rows = []

for s in input_summaries:
    input_totals_rows.append({
        "source": s["source"],
        "stage": s["stage"],
        "exists": s["exists"],
        "total_records": s["total_records"],
        "filter_pass_records": s["filter_pass_records"],
        "filter_nonpass_records": s["filter_nonpass_records"],
        "query_error": s["query_error"],
    })
    for k, v in s["filter_counts"].items():
        input_filter_rows.append({"source": s["source"], "stage": s["stage"], "filter": k, "count": v})
    for k, v in s["svtype_counts"].items():
        input_svtype_rows.append({"source": s["source"], "stage": s["stage"], "svtype": k, "count": v})
    for k, v in s["svlen_bin_counts"].items():
        input_svlen_rows.append({"source": s["source"], "stage": s["stage"], "svlen_bin": k, "count": v})
    for k, v in s["chrom_counts"].items():
        input_chrom_rows.append({"source": s["source"], "stage": s["stage"], "chrom": k, "count": v})

input_totals_df = pd.DataFrame(input_totals_rows)
input_filter_df = pd.DataFrame(input_filter_rows)
input_svtype_df = pd.DataFrame(input_svtype_rows)
input_svlen_df = pd.DataFrame(input_svlen_rows)
input_chrom_df = pd.DataFrame(input_chrom_rows)

write_csv(input_totals_df, "input_vcf_total_record_counts_raw_and_prepared.csv")
write_csv(input_filter_df, "input_vcf_filter_counts_raw_and_prepared.csv")
write_csv(input_svtype_df, "input_vcf_svtype_counts_raw_and_prepared.csv")
write_csv(input_svlen_df, "input_vcf_svlen_bin_counts_raw_and_prepared.csv")
write_csv(input_chrom_df, "input_vcf_chromosome_counts_raw_and_prepared.csv")

# ============================================================
# 4. Benchmark event VCF summaries
# ============================================================

event_summaries = []

for method, caller_dict in BENCHMARK_VCFS.items():
    for caller, status_dict in caller_dict.items():
        for status, path in status_dict.items():
            event_summaries.append(
                stream_vcf_summary(
                    path,
                    source=f"{method}_{caller}_{status}",
                    stage="benchmark_event",
                    method=method,
                    caller=caller,
                    status=status,
                )
            )

event_totals_rows = []
event_filter_rows = []
event_svtype_rows = []
event_svlen_rows = []
event_chrom_rows = []

for s in event_summaries:
    event_totals_rows.append({
        "method": s["method"],
        "caller": s["caller"],
        "status": s["status"],
        "exists": s["exists"],
        "total_records": s["total_records"],
        "filter_pass_records": s["filter_pass_records"],
        "filter_nonpass_records": s["filter_nonpass_records"],
        "query_error": s["query_error"],
    })
    for k, v in s["filter_counts"].items():
        event_filter_rows.append({
            "method": s["method"],
            "caller": s["caller"],
            "status": s["status"],
            "filter": k,
            "count": v,
        })
    for k, v in s["svtype_counts"].items():
        event_svtype_rows.append({
            "method": s["method"],
            "caller": s["caller"],
            "status": s["status"],
            "svtype": k,
            "count": v,
        })
    for k, v in s["svlen_bin_counts"].items():
        event_svlen_rows.append({
            "method": s["method"],
            "caller": s["caller"],
            "status": s["status"],
            "svlen_bin": k,
            "count": v,
        })
    for k, v in s["chrom_counts"].items():
        event_chrom_rows.append({
            "method": s["method"],
            "caller": s["caller"],
            "status": s["status"],
            "chrom": k,
            "count": v,
        })

event_totals_df = pd.DataFrame(event_totals_rows)
event_filter_df = pd.DataFrame(event_filter_rows)
event_svtype_df = pd.DataFrame(event_svtype_rows)
event_svlen_df = pd.DataFrame(event_svlen_rows)
event_chrom_df = pd.DataFrame(event_chrom_rows)

write_csv(event_totals_df, "benchmark_event_vcf_total_counts.csv")
write_csv(event_filter_df, "benchmark_event_vcf_filter_counts.csv")
write_csv(event_svtype_df, "benchmark_event_vcf_svtype_counts.csv")
write_csv(event_svlen_df, "benchmark_event_vcf_svlen_bin_counts.csv")
write_csv(event_chrom_df, "benchmark_event_vcf_chromosome_counts.csv")

# ============================================================
# 5. Summary comparison CSV counts
# ============================================================

comparison_csv_rows = []

for method, status_map in COMPARISON_CSVS.items():
    for status, path in status_map.items():
        row = {
            "method": method,
            "status": status,
            "exists": path.exists(),
            "rows": 0,
            "columns": "",
            "parse_error": "",
        }
        if path.exists():
            try:
                with open(path, "r", errors="replace") as fh:
                    header = fh.readline().strip()
                    row["columns"] = header
                    row["rows"] = sum(1 for _ in fh)
            except Exception as e:
                row["parse_error"] = str(e)
        comparison_csv_rows.append(row)

comparison_csv_df = pd.DataFrame(comparison_csv_rows)
write_csv(comparison_csv_df, "summary_comparison_csv_file_counts.csv")

# ============================================================
# 6. nf-core SURVIVOR file index
# ============================================================

survivor_file_rows = []
for source, path in NFCORE_SURVIVOR_STATS.items():
    survivor_file_rows.append({
        "source": source,
        "exists": path.exists(),
        "size_bytes": path.stat().st_size if path.exists() else 0,
    })

survivor_files_df = pd.DataFrame(survivor_file_rows)
write_csv(survivor_files_df, "nfcore_survivor_stats_files_index.csv")

# ============================================================
# 7. Generate input / SURVIVOR-style plots
# ============================================================

if not input_totals_df.empty:
    input_totals_df["source_stage"] = input_totals_df["source"] + " · " + input_totals_df["stage"]

    plot_single_bar(
        input_totals_df,
        x="source_stage",
        y="total_records",
        title="hg19 input VCF total records: raw and prepared",
        ylabel="Record count",
        filename="input_total_records_raw_prepared_log.png",
        logy=True,
    )

    pass_plot = input_totals_df.melt(
        id_vars=["source_stage"],
        value_vars=["filter_pass_records", "filter_nonpass_records"],
        var_name="filter_group",
        value_name="count",
    )
    plot_grouped_bar(
        pass_plot,
        index="source_stage",
        columns="filter_group",
        values="count",
        title="hg19 input VCF FILTER=PASS and non-PASS records",
        ylabel="Record count",
        filename="input_filter_pass_nonpass_records_log.png",
        logy=True,
        stacked=True,
    )

if not bcftools_wide.empty:
    class_cols = [c for c in ["snps", "indels", "others", "mnps"] if c in bcftools_wide.columns]
    if class_cols:
        bc_plot = bcftools_wide[["source"] + class_cols].melt(
            id_vars=["source"],
            value_vars=class_cols,
            var_name="variant_class",
            value_name="count",
        )
        plot_grouped_bar(
            bc_plot,
            index="source",
            columns="variant_class",
            values="count",
            title="hg19 bcftools variant classes from prepared VCFs",
            ylabel="Count",
            filename="bcftools_variant_classes_prepared_log.png",
            logy=True,
        )

prepared_svtype = input_svtype_df[input_svtype_df["stage"] == "prepared"].copy()
if not prepared_svtype.empty:
    plot_grouped_bar(
        prepared_svtype,
        index="svtype",
        columns="source",
        values="count",
        title="hg19 SURVIVOR-style SVTYPE counts from prepared VCFs",
        ylabel="Count",
        filename="survivor_style_svtype_counts_prepared_log.png",
        logy=True,
    )

    frac = prepared_svtype.copy()
    frac["fraction"] = frac.groupby("source")["count"].transform(lambda x: x / x.sum() if x.sum() else 0)
    write_csv(frac, "input_vcf_svtype_fractions_prepared.csv")
    plot_grouped_bar(
        frac,
        index="source",
        columns="svtype",
        values="fraction",
        title="hg19 SURVIVOR-style SVTYPE composition from prepared VCFs",
        ylabel="Fraction",
        filename="survivor_style_svtype_fraction_prepared.png",
        stacked=True,
    )

prepared_svlen = input_svlen_df[input_svlen_df["stage"] == "prepared"].copy()
if not prepared_svlen.empty:
    prepared_svlen["svlen_bin"] = pd.Categorical(prepared_svlen["svlen_bin"], categories=SVLEN_ORDER, ordered=True)
    prepared_svlen = prepared_svlen.sort_values("svlen_bin")

    plot_grouped_bar(
        prepared_svlen,
        index="svlen_bin",
        columns="source",
        values="count",
        title="hg19 SURVIVOR-style SVLEN bin counts from prepared VCFs",
        ylabel="Count",
        filename="survivor_style_svlen_bin_counts_prepared_log.png",
        logy=True,
    )

    frac = prepared_svlen.copy()
    frac["fraction"] = frac.groupby("source")["count"].transform(lambda x: x / x.sum() if x.sum() else 0)
    write_csv(frac, "input_vcf_svlen_fractions_prepared.csv")
    plot_grouped_bar(
        frac,
        index="source",
        columns="svlen_bin",
        values="fraction",
        title="hg19 SURVIVOR-style SVLEN composition from prepared VCFs",
        ylabel="Fraction",
        filename="survivor_style_svlen_fraction_prepared.png",
        stacked=True,
    )

prepared_chrom = input_chrom_df[input_chrom_df["stage"] == "prepared"].copy()
if not prepared_chrom.empty:
    chrom_order = sorted(prepared_chrom["chrom"].unique(), key=lambda x: (len(str(x)), str(x)))
    prepared_chrom["chrom"] = pd.Categorical(prepared_chrom["chrom"], categories=chrom_order, ordered=True)
    plot_grouped_bar(
        prepared_chrom,
        index="chrom",
        columns="source",
        values="count",
        title="hg19 prepared VCF record counts by chromosome",
        ylabel="Count",
        filename="prepared_vcf_chromosome_counts_log.png",
        logy=True,
    )

# ============================================================
# 8. Benchmark metric plots
# ============================================================

if not metrics.empty:
    for metric, label, filename in [
        ("precision", "Precision", "benchmark_precision_by_caller_tool.png"),
        ("recall", "Recall", "benchmark_recall_by_caller_tool.png"),
        ("f1", "F1", "benchmark_f1_by_caller_tool.png"),
    ]:
        plot_grouped_bar(
            metrics,
            index="caller",
            columns="method",
            values=metric,
            title=f"hg19 {label} by caller and benchmarking tool",
            ylabel=label,
            filename=filename,
            ylim=(0, 1.05) if metric == "precision" else None,
        )

    max_recall = float(metrics["recall"].max()) if len(metrics) else 0
    max_f1 = float(metrics["f1"].max()) if len(metrics) else 0

    plot_grouped_bar(
        metrics,
        index="caller",
        columns="method",
        values="recall",
        title="hg19 recall by caller and benchmarking tool, zoomed scale",
        ylabel="Recall",
        filename="benchmark_recall_by_caller_tool_zoomed.png",
        ylim=(0, max(0.01, max_recall * 1.35)),
    )

    plot_grouped_bar(
        metrics,
        index="caller",
        columns="method",
        values="f1",
        title="hg19 F1 by caller and benchmarking tool, zoomed scale",
        ylabel="F1",
        filename="benchmark_f1_by_caller_tool_zoomed.png",
        ylim=(0, max(0.02, max_f1 * 1.35)),
    )

    for metric, label, filename in [
        ("tp_base", "TP-base", "benchmark_tp_base_by_caller_tool_log.png"),
        ("tp_comp", "TP-comp", "benchmark_tp_comp_by_caller_tool_log.png"),
        ("fp", "FP", "benchmark_fp_by_caller_tool_log.png"),
        ("fn", "FN", "benchmark_fn_by_caller_tool_log.png"),
    ]:
        plot_grouped_bar(
            metrics,
            index="caller",
            columns="method",
            values=metric,
            title=f"hg19 {label} by caller and benchmarking tool",
            ylabel="Count",
            filename=filename,
            logy=True,
        )

    combo = metrics[["caller_tool", "precision", "recall", "f1"]].copy()
    combo_long = combo.melt(id_vars=["caller_tool"], var_name="metric", value_name="value")
    plot_grouped_bar(
        combo_long,
        index="caller_tool",
        columns="metric",
        values="value",
        title="hg19 precision, recall, and F1 together",
        ylabel="Metric value",
        filename="benchmark_precision_recall_f1_together.png",
        ylim=(0, 1.05),
    )

    for method in METHOD_ORDER:
        sub = metrics[metrics["method"].astype(str) == method].copy()
        if sub.empty:
            continue
        long = sub[["caller", "precision", "recall", "f1"]].melt(id_vars=["caller"], var_name="metric", value_name="value")
        plot_grouped_bar(
            long,
            index="caller",
            columns="metric",
            values="value",
            title=f"hg19 {method}: precision, recall, and F1 by caller",
            ylabel="Metric value",
            filename=f"{method.lower()}_precision_recall_f1_by_caller.png",
            ylim=(0, 1.05),
        )

    for caller in CALLER_ORDER:
        sub = metrics[metrics["caller"].astype(str) == caller].copy()
        if sub.empty:
            continue
        long = sub[["method", "precision", "recall", "f1"]].melt(id_vars=["method"], var_name="metric", value_name="value")
        plot_grouped_bar(
            long,
            index="method",
            columns="metric",
            values="value",
            title=f"hg19 {caller}: precision, recall, and F1 by benchmarking tool",
            ylabel="Metric value",
            filename=f"{caller.lower()}_precision_recall_f1_by_tool.png",
            ylim=(0, 1.05),
        )

    for metric, label in [("precision", "Precision"), ("recall", "Recall"), ("f1", "F1")]:
        mat = metrics.pivot_table(index="caller", columns="method", values=metric, aggfunc="mean", fill_value=0, observed=False)
        plot_heatmap(mat, f"hg19 {label} heatmap: caller × benchmarking tool", f"heatmap_{metric}_caller_by_tool.png", label)

    plot_scatter_metrics(metrics, "precision", "recall", "scatter_precision_vs_recall.png")
    plot_scatter_metrics(metrics, "precision", "f1", "scatter_precision_vs_f1.png")
    plot_scatter_metrics(metrics, "recall", "f1", "scatter_recall_vs_f1.png")

    if not method_delta_df.empty:
        for comparison in method_delta_df["comparison"].unique():
            sub = method_delta_df[method_delta_df["comparison"] == comparison].copy()
            long = sub.melt(
                id_vars=["caller"],
                value_vars=["delta_precision", "delta_recall", "delta_f1"],
                var_name="delta_metric",
                value_name="difference",
            )
            fname = re.sub(r"[^A-Za-z0-9]+", "_", comparison.lower()).strip("_")
            plot_grouped_bar(
                long,
                index="caller",
                columns="delta_metric",
                values="difference",
                title=f"hg19 metric differences: {comparison}",
                ylabel="Difference",
                filename=f"method_difference_{fname}.png",
            )

# ============================================================
# 9. Benchmark event plots
# ============================================================

if not event_totals_df.empty:
    for status in STATUS_ORDER:
        sub = event_totals_df[event_totals_df["status"] == status].copy()
        if sub.empty:
            continue
        plot_grouped_bar(
            sub,
            index="caller",
            columns="method",
            values="total_records",
            title=f"hg19 {status} event VCF records by caller and method",
            ylabel="Record count",
            filename=f"event_{status.lower()}_counts_by_caller_method_log.png",
            logy=True,
        )

    plot_grouped_bar(
        event_totals_df,
        index="status",
        columns="method",
        values="total_records",
        title="hg19 benchmark event VCF record counts by status and method",
        ylabel="Record count",
        filename="event_vcf_counts_by_status_method_log.png",
        logy=True,
    )

if not event_svtype_df.empty:
    for status in STATUS_ORDER:
        sub = event_svtype_df[event_svtype_df["status"] == status].copy()
        if sub.empty:
            continue
        plot_grouped_bar(
            sub,
            index="svtype",
            columns="caller",
            values="count",
            title=f"hg19 {status} SVTYPE counts by caller",
            ylabel="Count",
            filename=f"event_{status.lower()}_svtype_counts_by_caller_log.png",
            logy=True,
        )
        plot_grouped_bar(
            sub,
            index="svtype",
            columns="method",
            values="count",
            title=f"hg19 {status} SVTYPE counts by benchmarking method",
            ylabel="Count",
            filename=f"event_{status.lower()}_svtype_counts_by_method_log.png",
            logy=True,
        )

if not event_svlen_df.empty:
    event_svlen_df["svlen_bin"] = pd.Categorical(event_svlen_df["svlen_bin"], categories=SVLEN_ORDER, ordered=True)
    event_svlen_df = event_svlen_df.sort_values("svlen_bin")

    for status in STATUS_ORDER:
        sub = event_svlen_df[event_svlen_df["status"] == status].copy()
        if sub.empty:
            continue
        plot_grouped_bar(
            sub,
            index="svlen_bin",
            columns="caller",
            values="count",
            title=f"hg19 {status} SVLEN bin counts by caller",
            ylabel="Count",
            filename=f"event_{status.lower()}_svlen_counts_by_caller_log.png",
            logy=True,
        )
        plot_grouped_bar(
            sub,
            index="svlen_bin",
            columns="method",
            values="count",
            title=f"hg19 {status} SVLEN bin counts by benchmarking method",
            ylabel="Count",
            filename=f"event_{status.lower()}_svlen_counts_by_method_log.png",
            logy=True,
        )

if not event_chrom_df.empty:
    for status in STATUS_ORDER:
        sub = event_chrom_df[event_chrom_df["status"] == status].copy()
        if sub.empty:
            continue
        plot_grouped_bar(
            sub,
            index="chrom",
            columns="caller",
            values="count",
            title=f"hg19 {status} chromosome counts by caller",
            ylabel="Count",
            filename=f"event_{status.lower()}_chromosome_counts_by_caller_log.png",
            logy=True,
        )

# ============================================================
# 10. Existing nf-core / MultiQC plots
# ============================================================

selected_existing_names = [
    "variant_calling_summary.png",
    "survivor-cnt.png",
    "survivor-pct.png",
    "bcftools-stats-subtypes-cnt.png",
    "bcftools-stats-subtypes-pct.png",
    "bcftools_stats_indel-lengths-log.png",
    "truvari-bench-pre-rec-plot.png",
    "upset_truvari_tp_fn_mqc.png",
    "upset_truvari_tp_fp_mqc.png",
    "upset_svbenchmark_tp_fn_mqc.png",
    "upset_svbenchmark_tp_fp_mqc.png",
    "upset_rtgtools_tp_fp_mqc.png",
]

existing_plot_rows = []
selected_existing = []

for base in [MULTIQC_PLOTS, EXISTING_SUMMARY_PLOTS]:
    if not base.exists():
        continue
    for p in sorted(base.rglob("*.png")):
        existing_plot_rows.append({
            "file": str(p),
            "name": p.name,
        })
        if p.name in selected_existing_names:
            dest = PLOTS / f"existing_{p.name}"
            try:
                dest.write_bytes(p.read_bytes())
                selected_existing.append(dest.name)
            except Exception:
                pass

existing_plots_df = pd.DataFrame(existing_plot_rows)
write_csv(existing_plots_df, "existing_nfcore_plot_files.csv")

# ============================================================
# 11. Self-computed SURVIVOR stats through Docker
# ============================================================

DOCKER_IMAGE = "informationsea/survivor"

def parse_survivor_table_text(text, source, origin):
    rows = []
    lines = [l.strip() for l in str(text).splitlines() if l.strip()]
    for i in range(len(lines) - 1):
        header = lines[i].split()
        vals = lines[i + 1].split()
        if len(header) == len(vals) and any(h in header for h in ["Tot", "DEL", "DUP", "INS", "INV", "TRA"]):
            for h, v in zip(header, vals):
                try:
                    rows.append({
                        "source": source,
                        "metric": h,
                        "value": float(v),
                        "origin": origin,
                    })
                except Exception:
                    pass
            break
    return rows

def parse_survivor_file(path, source, origin):
    rows = []
    path = Path(path)
    if not path.exists():
        return rows
    text = path.read_text(errors="replace")
    rows.extend(parse_survivor_table_text(text, source, origin))

    for line in text.splitlines():
        raw = line.strip()
        if not raw or raw.startswith("#"):
            continue
        parts = raw.replace(":", " ").split()
        nums = []
        labels = []
        for p in parts:
            try:
                nums.append(float(p.replace(",", "")))
            except Exception:
                labels.append(p)
        if nums and labels:
            rows.append({
                "source": source,
                "metric": " ".join(labels),
                "value": nums[-1],
                "origin": origin,
            })
    return rows

# Remove old self-survivor outputs except plain_vcfs directory
for p in SELF.glob("*"):
    if p.is_file():
        p.unlink()

self_command_rows = []
self_file_rows = []
self_rows = []

# Decompress prepared VCFs to plain VCF
for source, gzvcf in PREPARED_VCFS.items():
    plain_vcf = PLAIN / f"{slug(source)}.prepared.vcf"

    if not gzvcf.exists():
        self_command_rows.append({
            "source": source,
            "step": "decompress",
            "command": "",
            "return_code": "",
            "status": "input VCF missing",
            "stdout": "",
            "stderr": str(gzvcf),
        })
        continue

    if shutil.which("bcftools") is None:
        self_command_rows.append({
            "source": source,
            "step": "decompress",
            "command": "bcftools",
            "return_code": "",
            "status": "bcftools not found",
            "stdout": "",
            "stderr": "bcftools not found in PATH",
        })
        continue

    cmd = ["bcftools", "view", "-o", str(plain_vcf), str(gzvcf)]
    stdout, stderr, rc = run_cmd(cmd)

    self_command_rows.append({
        "source": source,
        "step": "decompress",
        "command": " ".join(cmd),
        "return_code": rc,
        "status": "success" if rc == 0 and plain_vcf.exists() else "failed",
        "stdout": (stdout or "")[-1000:],
        "stderr": (stderr or "")[-1000:],
    })

# Run SURVIVOR stats with Docker
docker_ok = shutil.which("docker") is not None

for source in PREPARED_VCFS:
    plain_vcf = PLAIN / f"{slug(source)}.prepared.vcf"
    out_prefix = SELF / f"{slug(source)}_self_survivor"

    if not docker_ok:
        self_command_rows.append({
            "source": source,
            "step": "SURVIVOR stats",
            "command": "docker",
            "return_code": "",
            "status": "docker not found",
            "stdout": "",
            "stderr": "Docker not found in PATH",
        })
        continue

    if not plain_vcf.exists():
        self_command_rows.append({
            "source": source,
            "step": "SURVIVOR stats",
            "command": "",
            "return_code": "",
            "status": "plain VCF missing",
            "stdout": "",
            "stderr": str(plain_vcf),
        })
        continue

    container_vcf = "/data/" + str(plain_vcf.relative_to(ROOT))
    container_out = "/data/" + str(out_prefix.relative_to(ROOT))

    cmd = [
        "docker", "run", "--rm",
        "--platform", "linux/amd64",
        "-v", f"{ROOT}:/data",
        DOCKER_IMAGE,
        "SURVIVOR", "stats",
        container_vcf,
        "-1", "-1", "-1",
        container_out,
    ]

    stdout, stderr, rc = run_cmd(cmd)
    produced = sorted(SELF.glob(f"{out_prefix.name}*"))

    self_command_rows.append({
        "source": source,
        "step": "SURVIVOR stats",
        "command": " ".join(cmd),
        "return_code": rc,
        "status": "success" if rc == 0 and produced else "failed",
        "stdout": (stdout or "")[-2000:],
        "stderr": (stderr or "")[-2000:],
    })

    if stdout:
        self_rows.extend(parse_survivor_table_text(stdout, source, "self_SURVIVOR_stdout"))

    for f in produced:
        if f.is_file():
            self_file_rows.append({
                "source": source,
                "filename": f.name,
                "size_bytes": f.stat().st_size,
            })
            self_rows.extend(parse_survivor_file(f, source, "self_SURVIVOR_file"))

nfcore_rows = []
for source, f in NFCORE_SURVIVOR_STATS.items():
    nfcore_rows.extend(parse_survivor_file(f, source, "nfcore_SURVIVOR_file"))

self_commands_df = pd.DataFrame(self_command_rows)
self_files_df = pd.DataFrame(self_file_rows)
self_survivor_df = pd.DataFrame(self_rows)
nfcore_survivor_df = pd.DataFrame(nfcore_rows)

write_csv(self_commands_df, "self_survivor_command_log.csv")
write_csv(self_files_df, "self_survivor_output_files.csv")
write_csv(self_survivor_df, "self_survivor_parsed_stats.csv")
write_csv(nfcore_survivor_df, "nfcore_survivor_parsed_stats.csv")

survivor_diff_df = pd.DataFrame()
if not self_survivor_df.empty and not nfcore_survivor_df.empty:
    a = self_survivor_df.groupby(["source", "metric"], as_index=False)["value"].max()
    b = nfcore_survivor_df.groupby(["source", "metric"], as_index=False)["value"].max()
    merged = a.merge(b, on=["source", "metric"], suffixes=("_self", "_nfcore"), how="inner")
    if not merged.empty:
        merged["difference_self_minus_nfcore"] = merged["value_self"] - merged["value_nfcore"]
        survivor_diff_df = merged

write_csv(survivor_diff_df, "self_survivor_vs_nfcore_survivor_differences.csv")

# Self SURVIVOR plots
if not self_commands_df.empty:
    status_counts = self_commands_df["status"].value_counts().reset_index()
    status_counts.columns = ["status", "count"]
    plot_single_bar(
        status_counts,
        x="status",
        y="count",
        title="hg19 self-run SURVIVOR command status",
        ylabel="Count",
        filename="self_survivor_command_status.png",
    )

if not self_survivor_df.empty:
    keep = ["Tot", "DEL", "DUP", "INS", "INV", "TRA"]
    plot_df = self_survivor_df[self_survivor_df["metric"].isin(keep)].copy()
    if not plot_df.empty:
        plot_df = plot_df.groupby(["source", "metric"], as_index=False)["value"].max()
        plot_grouped_bar(
            plot_df,
            index="metric",
            columns="source",
            values="value",
            title="hg19 self-computed SURVIVOR stats",
            ylabel="Count",
            filename="self_survivor_metrics_log.png",
            logy=True,
        )

if not survivor_diff_df.empty:
    d = survivor_diff_df.copy()
    d["abs_difference"] = d["difference_self_minus_nfcore"].abs()
    d = d.sort_values("abs_difference", ascending=False).head(30)
    d["source_metric"] = d["source"] + " · " + d["metric"]
    plot_single_bar(
        d,
        x="source_metric",
        y="difference_self_minus_nfcore",
        title="hg19 self SURVIVOR minus nf-core SURVIVOR differences",
        ylabel="Difference",
        filename="self_survivor_vs_nfcore_differences.png",
    )

# ============================================================
# 12. Build embedded HTML report
# ============================================================

def metric_cards_html():
    if metrics.empty:
        return "<p>No benchmark metrics found.</p>"

    cards = []
    for metric, label in [("precision", "Highest precision"), ("recall", "Highest recall"), ("f1", "Highest F1")]:
        idx = metrics[metric].astype(float).idxmax()
        r = metrics.loc[idx]
        cards.append(f"""
        <div class="card">
          <div class="value">{float(r[metric]):.6f}</div>
          <div class="label">{label}</div>
          <div class="sub">{html.escape(str(r['caller']))} · {html.escape(str(r['method']))}</div>
        </div>
        """)

    for metric, label in [("tp_base", "Largest TP-base"), ("fp", "Largest FP"), ("fn", "Largest FN")]:
        idx = metrics[metric].astype(float).idxmax()
        r = metrics.loc[idx]
        cards.append(f"""
        <div class="card">
          <div class="value">{int(r[metric]):,}</div>
          <div class="label">{label}</div>
          <div class="sub">{html.escape(str(r['caller']))} · {html.escape(str(r['method']))}</div>
        </div>
        """)

    return "<div class='cards'>" + "".join(cards) + "</div>"

def input_cards_html():
    if input_totals_df.empty:
        return "<p>No input VCF totals found.</p>"
    cards = []
    prepared = input_totals_df[input_totals_df["stage"] == "prepared"].copy()
    for _, r in prepared.iterrows():
        cards.append(f"""
        <div class="card">
          <div class="value">{int(r['total_records']):,}</div>
          <div class="label">{html.escape(str(r['source']))} prepared records</div>
          <div class="sub">FILTER=PASS records: {int(r['filter_pass_records']):,}</div>
        </div>
        """)
    return "<div class='cards'>" + "".join(cards) + "</div>"

benchmark_plot_items = [
    ("benchmark_precision_by_caller_tool.png", "Precision by caller and benchmarking tool"),
    ("benchmark_recall_by_caller_tool.png", "Recall by caller and benchmarking tool"),
    ("benchmark_recall_by_caller_tool_zoomed.png", "Recall by caller and benchmarking tool, zoomed"),
    ("benchmark_f1_by_caller_tool.png", "F1 by caller and benchmarking tool"),
    ("benchmark_f1_by_caller_tool_zoomed.png", "F1 by caller and benchmarking tool, zoomed"),
    ("benchmark_precision_recall_f1_together.png", "Precision, recall, and F1 together"),
    ("benchmark_tp_base_by_caller_tool_log.png", "TP-base counts by caller and tool"),
    ("benchmark_tp_comp_by_caller_tool_log.png", "TP-comp counts by caller and tool"),
    ("benchmark_fp_by_caller_tool_log.png", "FP counts by caller and tool"),
    ("benchmark_fn_by_caller_tool_log.png", "FN counts by caller and tool"),
]

input_plot_items = [
    ("input_total_records_raw_prepared_log.png", "Input VCF total records, raw and prepared"),
    ("input_filter_pass_nonpass_records_log.png", "Input VCF FILTER=PASS and non-PASS records"),
    ("bcftools_variant_classes_prepared_log.png", "bcftools variant classes from prepared VCFs"),
    ("survivor_style_svtype_counts_prepared_log.png", "SURVIVOR-style SVTYPE counts from prepared VCFs"),
    ("survivor_style_svtype_fraction_prepared.png", "SURVIVOR-style SVTYPE composition from prepared VCFs"),
    ("survivor_style_svlen_bin_counts_prepared_log.png", "SURVIVOR-style SVLEN bin counts from prepared VCFs"),
    ("survivor_style_svlen_fraction_prepared.png", "SURVIVOR-style SVLEN composition from prepared VCFs"),
    ("prepared_vcf_chromosome_counts_log.png", "Prepared VCF record counts by chromosome"),
]

heatmap_scatter_items = [
    ("heatmap_precision_caller_by_tool.png", "Precision heatmap"),
    ("heatmap_recall_caller_by_tool.png", "Recall heatmap"),
    ("heatmap_f1_caller_by_tool.png", "F1 heatmap"),
    ("scatter_precision_vs_recall.png", "Precision vs recall"),
    ("scatter_precision_vs_f1.png", "Precision vs F1"),
    ("scatter_recall_vs_f1.png", "Recall vs F1"),
]

event_plot_items = [
    ("event_vcf_counts_by_status_method_log.png", "Benchmark event VCF counts by status and method"),
    ("event_fn_counts_by_caller_method_log.png", "FN records by caller and method"),
    ("event_fp_counts_by_caller_method_log.png", "FP records by caller and method"),
    ("event_tp_base_counts_by_caller_method_log.png", "TP-base records by caller and method"),
    ("event_tp_comp_counts_by_caller_method_log.png", "TP-comp records by caller and method"),
    ("event_fn_svtype_counts_by_caller_log.png", "FN SVTYPE counts by caller"),
    ("event_fn_svtype_counts_by_method_log.png", "FN SVTYPE counts by method"),
    ("event_fn_svlen_counts_by_caller_log.png", "FN SVLEN bin counts by caller"),
    ("event_fn_svlen_counts_by_method_log.png", "FN SVLEN bin counts by method"),
    ("event_fn_chromosome_counts_by_caller_log.png", "FN chromosome counts by caller"),
    ("event_fp_svtype_counts_by_caller_log.png", "FP SVTYPE counts by caller"),
    ("event_fp_svlen_counts_by_caller_log.png", "FP SVLEN bin counts by caller"),
    ("event_tp_base_svtype_counts_by_caller_log.png", "TP-base SVTYPE counts by caller"),
    ("event_tp_base_svlen_counts_by_caller_log.png", "TP-base SVLEN bin counts by caller"),
]

per_tool_items = []
for method in METHOD_ORDER:
    per_tool_items.append((f"{method.lower()}_precision_recall_f1_by_caller.png", f"{method}: precision, recall, and F1 by caller"))

per_caller_items = []
for caller in CALLER_ORDER:
    per_caller_items.append((f"{caller.lower()}_precision_recall_f1_by_tool.png", f"{caller}: precision, recall, and F1 by tool"))

method_difference_items = []
for p in sorted(PLOTS.glob("method_difference_*.png")):
    method_difference_items.append((p.name, p.stem.replace("_", " ")))

existing_plot_items = [(f, f.replace("existing_", "").replace(".png", "").replace("_", " ")) for f in selected_existing]

self_survivor_items = [
    ("self_survivor_command_status.png", "Self-run SURVIVOR command status"),
    ("self_survivor_metrics_log.png", "Self-computed SURVIVOR statistics"),
    ("self_survivor_vs_nfcore_differences.png", "Self SURVIVOR vs nf-core SURVIVOR differences"),
]

html_text = f"""<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="UTF-8">
<title>hg19 SV benchmarking results report</title>
<style>
:root {{
  --bg:#0d1117;
  --surface:#161b22;
  --surface2:#1f2937;
  --border:#30363d;
  --text:#e6edf3;
  --subtext:#c9d1d9;
  --muted:#8b949e;
  --accent:#58a6ff;
}}
* {{ box-sizing:border-box; }}
body {{
  margin:0;
  background:var(--bg);
  color:var(--text);
  font-family:Arial, sans-serif;
  font-size:14px;
  line-height:1.55;
}}
.header {{
  padding:34px 48px 28px;
  border-bottom:1px solid var(--border);
  background:linear-gradient(135deg,#0d1117 0%,#161b22 65%,#1c2a3a 100%);
}}
.header h1 {{
  margin:0;
  color:var(--accent);
  font-size:30px;
}}
.header .meta {{
  color:var(--muted);
  margin-top:8px;
}}
.badge {{
  display:inline-block;
  border:1px solid rgba(88,166,255,.35);
  background:rgba(88,166,255,.10);
  color:var(--accent);
  padding:2px 8px;
  border-radius:4px;
  font-size:11px;
  text-transform:uppercase;
  font-weight:700;
  margin-right:6px;
}}
.container {{
  max-width:1450px;
  margin:0 auto;
  padding:34px 48px 80px;
}}
section {{
  margin-bottom:46px;
}}
h2 {{
  color:var(--subtext);
  border-bottom:1px solid var(--border);
  padding-bottom:8px;
  margin-top:0;
}}
h3 {{
  color:var(--subtext);
}}
.note {{
  background:rgba(88,166,255,.07);
  border-left:4px solid var(--accent);
  padding:12px 16px;
  margin-bottom:18px;
  color:var(--subtext);
}}
.cards {{
  display:grid;
  grid-template-columns:repeat(auto-fit,minmax(180px,1fr));
  gap:12px;
  margin-bottom:20px;
}}
.card {{
  background:var(--surface);
  border:1px solid var(--border);
  border-radius:8px;
  padding:15px;
}}
.value {{
  font-size:25px;
  font-weight:800;
  color:var(--accent);
}}
.label {{
  color:var(--muted);
  font-size:11px;
  text-transform:uppercase;
  letter-spacing:1px;
}}
.sub {{
  color:var(--subtext);
  font-size:12px;
  margin-top:4px;
}}
.plot-grid {{
  display:grid;
  grid-template-columns:repeat(auto-fit,minmax(520px,1fr));
  gap:20px;
}}
.plot-card {{
  background:var(--surface);
  border:1px solid var(--border);
  border-radius:8px;
  padding:15px;
}}
.plot-card h3 {{
  margin:0 0 10px;
  font-size:15px;
}}
.plot-card img {{
  width:100%;
  background:white;
  border:1px solid var(--border);
  border-radius:4px;
}}
.tbl {{
  width:100%;
  border-collapse:collapse;
  background:var(--surface);
  margin:14px 0 24px;
}}
.tbl th {{
  background:var(--surface2);
  color:var(--accent);
  text-align:left;
  padding:8px 10px;
  font-size:11px;
  text-transform:uppercase;
}}
.tbl td {{
  border-top:1px solid var(--border);
  padding:7px 10px;
  color:var(--subtext);
  font-size:13px;
  vertical-align:top;
}}
code {{
  background:var(--surface2);
  color:var(--subtext);
  padding:2px 5px;
  border-radius:3px;
}}
.small-note {{
  color:var(--muted);
  font-size:12px;
}}
</style>
</head>
<body>

<div class="header">
  <span class="badge">hg19</span>
  <span class="badge">GRCh37</span>
  <span class="badge">GIAB HG002</span>
  <span class="badge">Truvari</span>
  <span class="badge">SVbenchmark</span>
  <span class="badge">RTGtools</span>
  <h1>hg19 SV benchmarking results report</h1>
  <div class="meta">
    Generated {datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")} from the hg19 result files.
  </div>
</div>

<div class="container">

<section>
  <h2>1. Results overview</h2>
  <div class="note">
    This hg19 report summarizes benchmark metrics, input VCF statistics, SURVIVOR-style VCF summaries, benchmark event VCFs, and self-computed SURVIVOR statistics.
  </div>
  {metric_cards_html()}
</section>

<section>
  <h2>2. Input record overview</h2>
  <div class="note">
    These cards show the number of records in the prepared truth and caller VCFs after nf-core preprocessing.
    FILTER=PASS records means records where the VCF FILTER column is exactly PASS.
  </div>
  {input_cards_html()}
</section>

{plot_section("3. Input VCF composition and SURVIVOR-style plots", input_plot_items)}

<section>
  <h2>4. bcftools stats tables</h2>
  <h3>Parsed bcftools statistics</h3>
  {df_to_html(bcftools_stats_df, max_rows=100)}
  <h3>bcftools statistics, wide format</h3>
  {df_to_html(bcftools_wide, max_rows=100)}
</section>

<section>
  <h2>5. Benchmark summary metrics</h2>
  {df_to_html(metrics, max_rows=100)}
</section>

{plot_section("6. Benchmark metric plots", benchmark_plot_items)}
{plot_section("7. Heatmaps and scatter plots", heatmap_scatter_items)}
{plot_section("8. Per-tool metric plots", per_tool_items)}
{plot_section("9. Per-caller metric plots", per_caller_items)}
{plot_section("10. Benchmarking-tool difference plots", method_difference_items)}

<section>
  <h2>11. Benchmarking-tool metric differences</h2>
  {df_to_html(method_delta_df, max_rows=200)}
</section>

<section>
  <h2>12. Benchmark event VCF total counts</h2>
  {df_to_html(event_totals_df, max_rows=200)}
</section>

{plot_section("13. TP / FP / FN event plots", event_plot_items)}

<section>
  <h2>14. Event SVTYPE, SVLEN, chromosome tables</h2>
  <h3>Event SVTYPE counts</h3>
  {df_to_html(event_svtype_df, max_rows=200)}
  <h3>Event SVLEN bin counts</h3>
  {df_to_html(event_svlen_df, max_rows=200)}
  <h3>Event chromosome counts</h3>
  {df_to_html(event_chrom_df, max_rows=200)}
</section>

<section>
  <h2>15. Summary comparison CSV counts</h2>
  {df_to_html(comparison_csv_df, max_rows=100)}
</section>

<section>
  <h2>16. Input VCF tables</h2>
  <h3>Total counts</h3>
  {df_to_html(input_totals_df, max_rows=100)}
  <h3>SVTYPE counts</h3>
  {df_to_html(input_svtype_df, max_rows=200)}
  <h3>SVLEN bin counts</h3>
  {df_to_html(input_svlen_df, max_rows=200)}
  <h3>FILTER counts</h3>
  {df_to_html(input_filter_df, max_rows=200)}
</section>

<section>
  <h2>17. nf-core SURVIVOR stats files</h2>
  {df_to_html(survivor_files_df, max_rows=100)}
</section>

{plot_section("18. Existing nf-core / MultiQC plots included for reference", existing_plot_items)}

<section>
  <h2>19. Self-computed SURVIVOR statistics</h2>
  <div class="note">
    This section decompresses the prepared hg19 VCFs to plain VCF, runs <code>SURVIVOR stats</code> through Docker using <code>-1 -1 -1</code> thresholds, parses the self-computed outputs, and compares them with the nf-core SURVIVOR statistics files.
  </div>

  <h3>Self-run SURVIVOR command log</h3>
  {df_to_html(self_commands_df, max_rows=80)}

  <h3>Self-run SURVIVOR output files</h3>
  {df_to_html(self_files_df, max_rows=80)}

  {plot_section("Self-computed SURVIVOR plots", self_survivor_items)}

  <h3>Parsed self-computed SURVIVOR statistics</h3>
  {df_to_html(self_survivor_df, max_rows=100)}

  <h3>Parsed nf-core SURVIVOR statistics</h3>
  {df_to_html(nfcore_survivor_df, max_rows=100)}

  <h3>Self-computed SURVIVOR vs nf-core SURVIVOR differences</h3>
  {df_to_html(survivor_diff_df, max_rows=100)}
</section>

</div>
</body>
</html>
"""

REPORT = HTML / "hg19_sv_benchmarking_full_results_report_WITH_SELF_SURVIVOR.html"
REPORT.write_text(html_text)

summary_txt = LOGS / "hg19_report_generation_summary.txt"
summary_txt.write_text(
    "\n".join([
        f"Generated: {datetime.datetime.now().isoformat()}",
        f"Report: {REPORT}",
        f"Tables: {TABLES}",
        f"Plots: {PLOTS}",
        f"Metric rows: {len(metrics)}",
        f"Input total rows: {len(input_totals_df)}",
        f"Benchmark event rows: {len(event_totals_df)}",
        f"Self SURVIVOR rows: {len(self_survivor_df)}",
    ])
)

print()
print("DONE")
print("Final hg19 report:")
print(REPORT)
print()
print("Tables:")
print(TABLES)
print()
print("Plots:")
print(PLOTS)
print()
print("Log summary:")
print(summary_txt)
