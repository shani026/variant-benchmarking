from pathlib import Path
import subprocess
import datetime
import shutil
import textwrap
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[2]

OUT = ROOT / "analysis" / "hg38_complete"
TABLES = OUT / "tables"
PLOTS = OUT / "plots"
HTML = OUT / "html"
LOGS = OUT / "logs"

for d in [TABLES, PLOTS, HTML, LOGS]:
    d.mkdir(parents=True, exist_ok=True)

TRUVARI_CSV = ROOT / "benchmark-runs/nfcore/hg38/structural/summary/tables/truvari/truvari.summary.csv"
SVBENCHMARK_CSV = ROOT / "benchmark-runs/nfcore/hg38/structural/summary/tables/svbenchmark/svbenchmark.summary.csv"

VCFS = {
    "Truth": ROOT / "Truth_vcfs/hg38/HG002_GRCh38_v5.0q_stvar.vcf.gz",
    "Delly": ROOT / "caller_vcfs/HG002_delly_lr.hg38.renamed.vcf.gz",
    "DRAGEN": ROOT / "caller_vcfs/HG002_dragen.sv.hg38.renamed.vcf.gz",
    "Sniffles": ROOT / "caller_vcfs/HG002_sniffles.hg38.renamed.vcf.gz",
}

CALLER_ORDER = ["delly", "dragen", "sniffles"]
CALLER_LABEL = {
    "delly": "Delly",
    "dragen": "DRAGEN",
    "sniffles": "Sniffles",
}

METHOD_ORDER = ["Truvari", "SVbenchmark"]

def run_cmd(cmd):
    p = subprocess.run(cmd, text=True, capture_output=True)
    if p.returncode != 0:
        return None, p.stderr.strip()
    return p.stdout, None

def require_file(path):
    if not path.exists():
        raise FileNotFoundError(f"Missing required file: {path}")

for p in [TRUVARI_CSV, SVBENCHMARK_CSV, *VCFS.values()]:
    require_file(p)

if shutil.which("bcftools") is None:
    raise RuntimeError("bcftools is not available in PATH. Install/activate bcftools first.")

# -----------------------------
# 1. Load benchmark metrics
# -----------------------------

def load_metric_csv(path, method):
    df = pd.read_csv(path)
    df["method"] = method
    df = df.rename(columns={
        "Caller": "caller",
        "TP_base": "tp_base",
        "TP_comp": "tp_comp",
        "FP": "fp",
        "FN": "fn",
        "Precision": "precision",
        "Recall": "recall",
        "F1": "f1",
        "Tool": "tool",
        "File": "source_file",
    })

    keep = ["method", "tool", "source_file", "caller", "tp_base", "tp_comp", "fp", "fn", "precision", "recall", "f1"]
    df = df[keep]

    for col in ["tp_base", "tp_comp", "fp", "fn"]:
        df[col] = pd.to_numeric(df[col], errors="coerce").fillna(0).astype(int)

    for col in ["precision", "recall", "f1"]:
        df[col] = pd.to_numeric(df[col], errors="coerce").fillna(0.0)

    return df

metrics = pd.concat([
    load_metric_csv(TRUVARI_CSV, "Truvari"),
    load_metric_csv(SVBENCHMARK_CSV, "SVbenchmark"),
], ignore_index=True)

metrics["caller_label"] = metrics["caller"].map(CALLER_LABEL).fillna(metrics["caller"])
metrics["method"] = pd.Categorical(metrics["method"], categories=METHOD_ORDER, ordered=True)
metrics["caller_label"] = pd.Categorical(metrics["caller_label"], categories=["Delly", "DRAGEN", "Sniffles"], ordered=True)
metrics = metrics.sort_values(["method", "caller_label"])

metrics.to_csv(TABLES / "hg38_truvari_svbenchmark_metrics.csv", index=False)

mean_by_caller = (
    metrics
    .groupby(["caller", "caller_label"], observed=True)[["precision", "recall", "f1", "tp_base", "tp_comp", "fp", "fn"]]
    .mean(numeric_only=True)
    .reset_index()
)

mean_by_caller.to_csv(TABLES / "hg38_mean_metrics_by_caller.csv", index=False)

# Delta: SVbenchmark - Truvari
tru = metrics[metrics["method"] == "Truvari"].set_index("caller")
svb = metrics[metrics["method"] == "SVbenchmark"].set_index("caller")

delta_rows = []
for caller in CALLER_ORDER:
    if caller in tru.index and caller in svb.index:
        delta_rows.append({
            "caller": caller,
            "caller_label": CALLER_LABEL[caller],
            "delta_precision_svbenchmark_minus_truvari": svb.loc[caller, "precision"] - tru.loc[caller, "precision"],
            "delta_recall_svbenchmark_minus_truvari": svb.loc[caller, "recall"] - tru.loc[caller, "recall"],
            "delta_f1_svbenchmark_minus_truvari": svb.loc[caller, "f1"] - tru.loc[caller, "f1"],
        })

delta_df = pd.DataFrame(delta_rows)
delta_df.to_csv(TABLES / "hg38_method_discrepancy_svbenchmark_minus_truvari.csv", index=False)

# -----------------------------
# 2. Parse actual VCFs for SURVIVOR-style statistics
# -----------------------------

def query_vcf_records(vcf):
    fmt = "%CHROM\\t%POS\\t%INFO/END\\t%INFO/SVTYPE\\t%INFO/SVLEN\\t%FILTER\\n"
    stdout, err = run_cmd(["bcftools", "query", "-f", fmt, str(vcf)])

    if err:
        raise RuntimeError(f"bcftools failed for {vcf}: {err}")

    rows = []

    for line in stdout.splitlines():
        parts = line.split("\t")
        if len(parts) != 6:
            continue

        chrom, pos, end, svtype, svlen, filt = parts

        svtype = svtype if svtype and svtype != "." else "UNKNOWN"
        filt = filt if filt and filt != "." else "UNKNOWN"

        length_value = None

        if svlen and svlen != ".":
            raw = svlen.split(",")[0]
            try:
                length_value = abs(int(float(raw)))
            except Exception:
                length_value = None

        if length_value is None and end and end != ".":
            try:
                length_value = abs(int(end) - int(pos))
            except Exception:
                length_value = None

        if length_value is None:
            length_bin = "NA/BND"
        elif length_value < 50:
            length_bin = "<50 bp"
        elif length_value < 100:
            length_bin = "50-99 bp"
        elif length_value < 500:
            length_bin = "100-499 bp"
        elif length_value < 1000:
            length_bin = "500-999 bp"
        elif length_value < 5000:
            length_bin = "1-4.9 kb"
        elif length_value < 10000:
            length_bin = "5-9.9 kb"
        elif length_value < 50000:
            length_bin = "10-49.9 kb"
        else:
            length_bin = ">=50 kb"

        rows.append({
            "chrom": chrom,
            "pos": pos,
            "end": end,
            "svtype": svtype,
            "svlen": length_value,
            "svlen_bin": length_bin,
            "filter": filt,
        })

    return pd.DataFrame(rows)

vcf_records = {}

for label, vcf in VCFS.items():
    df = query_vcf_records(vcf)
    df["source"] = label
    df["vcf"] = str(vcf)
    vcf_records[label] = df

all_vcf_records = pd.concat(vcf_records.values(), ignore_index=True)
all_vcf_records.to_csv(TABLES / "hg38_actual_vcf_records_svtype_svlen_filter.csv", index=False)

svtype_counts = (
    all_vcf_records
    .groupby(["source", "svtype"], observed=True)
    .size()
    .reset_index(name="count")
)
svtype_counts.to_csv(TABLES / "hg38_survivor_style_svtype_counts.csv", index=False)

svtype_fraction = svtype_counts.copy()
svtype_fraction["fraction"] = svtype_fraction.groupby("source")["count"].transform(lambda x: x / x.sum())
svtype_fraction.to_csv(TABLES / "hg38_survivor_style_svtype_fractions.csv", index=False)

svlen_counts = (
    all_vcf_records
    .groupby(["source", "svlen_bin"], observed=True)
    .size()
    .reset_index(name="count")
)
svlen_counts.to_csv(TABLES / "hg38_survivor_style_svlen_bin_counts.csv", index=False)

filter_counts = (
    all_vcf_records
    .groupby(["source", "filter"], observed=True)
    .size()
    .reset_index(name="count")
)
filter_counts.to_csv(TABLES / "hg38_filter_counts.csv", index=False)

vcf_totals = (
    all_vcf_records
    .groupby("source", observed=True)
    .agg(
        total_records=("svtype", "size"),
        pass_records=("filter", lambda x: int((x == "PASS").sum())),
        nonpass_records=("filter", lambda x: int((x != "PASS").sum())),
        svtypes=("svtype", lambda x: ", ".join(sorted(set(x)))),
    )
    .reset_index()
)
vcf_totals.to_csv(TABLES / "hg38_vcf_totals.csv", index=False)

# Collect nf-core SURVIVOR stats files if present
survivor_stat_paths = []
for pattern in ["*_mqc.stats", "*.stats", "*SURVIVOR*", "*survivor*"]:
    survivor_stat_paths.extend((ROOT / "benchmark-runs/nfcore/hg38/structural").rglob(pattern))

survivor_stat_paths = sorted(set([p for p in survivor_stat_paths if p.is_file()]))

with open(LOGS / "hg38_survivor_stats_files_found.txt", "w") as fh:
    for p in survivor_stat_paths:
        fh.write(str(p) + "\n")

# Optional local SURVIVOR availability check
with open(LOGS / "hg38_survivor_binary_check.txt", "w") as fh:
    if shutil.which("SURVIVOR"):
        fh.write(f"SURVIVOR found: {shutil.which('SURVIVOR')}\n")
    else:
        fh.write("SURVIVOR binary not found in PATH. Report uses nf-core SURVIVOR outputs if present and independent bcftools-derived SURVIVOR-style VCF statistics.\n")

# -----------------------------
# 3. Plot helpers
# -----------------------------

def savefig(name):
    plt.tight_layout()
    plt.savefig(PLOTS / name, dpi=220, bbox_inches="tight")
    plt.close()

def label_bars(ax, percent=False, integer=False, small=False):
    for container in ax.containers:
        labels = []
        for v in container.datavalues:
            if pd.isna(v):
                labels.append("")
            elif percent:
                labels.append(f"{v:.3f}")
            elif integer:
                labels.append(f"{int(v):,}")
            else:
                labels.append(f"{v:.3f}")
        ax.bar_label(container, labels=labels, fontsize=7 if small else 8, padding=2)

def plot_metric_by_method(metric, ylabel, filename, ylim=None):
    pivot = metrics.pivot_table(index="caller_label", columns="method", values=metric, aggfunc="mean", observed=True)
    ax = pivot.plot(kind="bar", figsize=(8.5, 5))
    ax.set_title(f"hg38 {ylabel} by caller and benchmarking tool")
    ax.set_xlabel("Caller")
    ax.set_ylabel(ylabel)
    if ylim:
        ax.set_ylim(*ylim)
    ax.legend(title="Benchmarking tool")
    plt.xticks(rotation=0)
    label_bars(ax, percent=metric in ["precision", "recall", "f1"], small=True)
    savefig(filename)

def plot_method_specific(method):
    df = metrics[metrics["method"] == method].copy()
    df = df.set_index("caller_label")[["precision", "recall", "f1"]]
    ax = df.plot(kind="bar", figsize=(8, 5))
    ax.set_title(f"{method}: precision, recall and F1 by caller")
    ax.set_xlabel("Caller")
    ax.set_ylabel("Metric value")
    ax.set_ylim(0, max(1.05, df.max().max() * 1.15))
    ax.legend(title="Metric")
    plt.xticks(rotation=0)
    label_bars(ax, percent=True, small=True)
    savefig(f"hg38_{method.lower()}_precision_recall_f1_by_caller.png")

def plot_caller_specific(caller):
    label = CALLER_LABEL[caller]
    df = metrics[metrics["caller"] == caller].copy()
    df = df.set_index("method")[["precision", "recall", "f1"]]
    ax = df.plot(kind="bar", figsize=(8, 5))
    ax.set_title(f"{label}: precision, recall and F1 by benchmarking tool")
    ax.set_xlabel("Benchmarking tool")
    ax.set_ylabel("Metric value")
    ax.set_ylim(0, max(1.05, df.max().max() * 1.15))
    ax.legend(title="Metric")
    plt.xticks(rotation=0)
    label_bars(ax, percent=True, small=True)
    savefig(f"hg38_{caller}_precision_recall_f1_by_tool.png")

def plot_scatter(x, y, filename):
    fig, ax = plt.subplots(figsize=(7, 5))

    for _, r in metrics.iterrows():
        ax.scatter(r[x], r[y], s=80)
        ax.annotate(
            f"{r['caller_label']} · {r['method']}",
            (r[x], r[y]),
            xytext=(5, 5),
            textcoords="offset points",
            fontsize=8,
        )

    ax.set_title(f"hg38 {y.upper()} vs {x.upper()}")
    ax.set_xlabel(x.capitalize())
    ax.set_ylabel(y.capitalize())
    ax.grid(True, alpha=0.3)
    savefig(filename)

def plot_table_counts(df, index, columns, values, title, ylabel, filename, stacked=False, normalize=False):
    pivot = df.pivot_table(index=index, columns=columns, values=values, aggfunc="sum", fill_value=0, observed=True)

    if normalize:
        pivot = pivot.div(pivot.sum(axis=1), axis=0).fillna(0)

    ax = pivot.plot(kind="bar", stacked=stacked, figsize=(10.5, 6))
    ax.set_title(title)
    ax.set_xlabel(index)
    ax.set_ylabel(ylabel)
    ax.legend(title=columns, bbox_to_anchor=(1.02, 1), loc="upper left")
    plt.xticks(rotation=0)
    savefig(filename)

# -----------------------------
# 4. Benchmark metric plots
# -----------------------------

plot_metric_by_method("recall", "Recall", "hg38_recall_by_caller_and_tool.png", ylim=(0, max(0.01, metrics["recall"].max() * 1.35)))
plot_metric_by_method("precision", "Precision", "hg38_precision_by_caller_and_tool.png", ylim=(0, 1.05))
plot_metric_by_method("f1", "F1", "hg38_f1_by_caller_and_tool.png", ylim=(0, max(0.02, metrics["f1"].max() * 1.35)))
plot_metric_by_method("fn", "False negatives", "hg38_false_negatives_by_caller_and_tool.png")

# Together comparison: all metrics side-by-side per caller+tool
combo = metrics.copy()
combo["caller_tool"] = combo["caller_label"].astype(str) + " · " + combo["method"].astype(str)
combo_plot = combo.set_index("caller_tool")[["precision", "recall", "f1"]]
ax = combo_plot.plot(kind="bar", figsize=(11, 5.5))
ax.set_title("hg38 precision, recall and F1 together")
ax.set_xlabel("Caller · benchmarking tool")
ax.set_ylabel("Metric value")
ax.set_ylim(0, 1.05)
ax.legend(title="Metric")
plt.xticks(rotation=35, ha="right")
label_bars(ax, percent=True, small=True)
savefig("hg38_precision_recall_f1_together.png")

for method in METHOD_ORDER:
    plot_method_specific(method)

for caller in CALLER_ORDER:
    plot_caller_specific(caller)

plot_scatter("precision", "recall", "hg38_precision_vs_recall.png")
plot_scatter("precision", "f1", "hg38_precision_vs_f1.png")
plot_scatter("recall", "f1", "hg38_recall_vs_f1.png")

# Method discrepancy plots
if not delta_df.empty:
    delta_plot = delta_df.set_index("caller_label")[
        [
            "delta_precision_svbenchmark_minus_truvari",
            "delta_recall_svbenchmark_minus_truvari",
            "delta_f1_svbenchmark_minus_truvari",
        ]
    ]
    ax = delta_plot.plot(kind="bar", figsize=(9, 5))
    ax.axhline(0, linewidth=1)
    ax.set_title("SVbenchmark minus Truvari metric difference")
    ax.set_xlabel("Caller")
    ax.set_ylabel("Difference")
    ax.legend(["Δ precision", "Δ recall", "Δ F1"], title="Metric difference")
    plt.xticks(rotation=0)
    label_bars(ax, small=True)
    savefig("hg38_svbenchmark_minus_truvari_differences.png")

# -----------------------------
# 5. SURVIVOR-style VCF plots
# -----------------------------

plot_table_counts(
    svtype_counts,
    index="svtype",
    columns="source",
    values="count",
    title="SURVIVOR-style SVTYPE counts from actual VCFs",
    ylabel="Record count",
    filename="hg38_survivor_svtype_counts_by_vcf.png",
    stacked=False,
)

plot_table_counts(
    svtype_fraction,
    index="source",
    columns="svtype",
    values="fraction",
    title="SURVIVOR-style SVTYPE composition fraction",
    ylabel="Fraction of records",
    filename="hg38_survivor_svtype_fraction_stacked.png",
    stacked=True,
    normalize=False,
)

svlen_order = ["<50 bp", "50-99 bp", "100-499 bp", "500-999 bp", "1-4.9 kb", "5-9.9 kb", "10-49.9 kb", ">=50 kb", "NA/BND"]
svlen_counts["svlen_bin"] = pd.Categorical(svlen_counts["svlen_bin"], categories=svlen_order, ordered=True)
svlen_counts = svlen_counts.sort_values("svlen_bin")

plot_table_counts(
    svlen_counts,
    index="svlen_bin",
    columns="source",
    values="count",
    title="SURVIVOR-style SV length bin counts from actual VCFs",
    ylabel="Record count",
    filename="hg38_survivor_svlen_bin_counts.png",
    stacked=False,
)

plot_table_counts(
    svlen_counts,
    index="source",
    columns="svlen_bin",
    values="count",
    title="SURVIVOR-style SV length composition",
    ylabel="Record count",
    filename="hg38_survivor_svlen_composition_stacked.png",
    stacked=True,
)

plot_table_counts(
    filter_counts,
    index="source",
    columns="filter",
    values="count",
    title="VCF FILTER status counts",
    ylabel="Record count",
    filename="hg38_filter_status_counts.png",
    stacked=True,
)

ax = vcf_totals.set_index("source")[["total_records", "pass_records", "nonpass_records"]].plot(kind="bar", figsize=(9, 5))
ax.set_title("Total, PASS and non-PASS records by VCF")
ax.set_xlabel("VCF source")
ax.set_ylabel("Record count")
plt.xticks(rotation=0)
label_bars(ax, integer=True, small=True)
savefig("hg38_total_pass_nonpass_records.png")

# -----------------------------
# 6. HTML helpers
# -----------------------------

def df_to_html(df, float_cols=True):
    return df.to_html(
        index=False,
        border=0,
        classes="tbl",
        float_format=(lambda x: f"{x:.6g}") if float_cols else None,
    )

def img(filename, caption):
    return f"""
    <div class="plot-card">
      <h3>{caption}</h3>
      <img src="../plots/{filename}" alt="{caption}">
    </div>
    """

def metric_cards():
    cards = ""
    for _, r in mean_by_caller.iterrows():
        caller = r["caller_label"]
        cards += f"""
        <div class="caller-block">
          <h3>{caller}</h3>
          <div class="cards">
            <div class="card"><div class="value good">{r['precision']:.4f}</div><div class="label">Mean precision</div></div>
            <div class="card"><div class="value bad">{r['recall']:.4f}</div><div class="label">Mean recall</div></div>
            <div class="card"><div class="value bad">{r['f1']:.4f}</div><div class="label">Mean F1</div></div>
            <div class="card"><div class="value">{int(round(r['tp_base'])):,}</div><div class="label">Mean TP-base</div></div>
            <div class="card"><div class="value warn">{int(round(r['fp'])):,}</div><div class="label">Mean FP</div></div>
            <div class="card"><div class="value bad">{int(round(r['fn'])):,}</div><div class="label">Mean FN</div></div>
          </div>
        </div>
        """
    return cards

def conclusion_text():
    mean_rank_recall = mean_by_caller.sort_values("recall", ascending=False)
    mean_rank_f1 = mean_by_caller.sort_values("f1", ascending=False)
    mean_rank_precision = mean_by_caller.sort_values("precision", ascending=False)

    best_recall = mean_rank_recall.iloc[0]
    best_f1 = mean_rank_f1.iloc[0]
    best_precision = mean_rank_precision.iloc[0]

    worst_recall = mean_by_caller.sort_values("recall", ascending=True).iloc[0]

    max_recall = mean_by_caller["recall"].max()
    low_recall_note = ""
    if max_recall < 0.05:
        low_recall_note = """
        <li><strong>Recall is globally very low.</strong> Even the best caller recovers less than 5% of the truth set under these benchmark settings. This means the main issue is not just caller ranking; it is the large false-negative burden and/or differences in SV representation, truth-set scope, SVTYPE composition, or benchmark compatibility.</li>
        """

    return f"""
    <ul>
      <li><strong>Best mean recall:</strong> {best_recall['caller_label']} ({best_recall['recall']:.4f}).</li>
      <li><strong>Best mean F1:</strong> {best_f1['caller_label']} ({best_f1['f1']:.4f}).</li>
      <li><strong>Best mean precision:</strong> {best_precision['caller_label']} ({best_precision['precision']:.4f}).</li>
      <li><strong>Lowest mean recall:</strong> {worst_recall['caller_label']} ({worst_recall['recall']:.4f}).</li>
      {low_recall_note}
      <li><strong>Precision and recall are decoupled.</strong> High precision with very low recall means the caller/tool combination is conservative: many calls that are matched are correct, but most truth variants are missed.</li>
      <li><strong>SVbenchmark and Truvari disagree slightly.</strong> Use the method-difference plot to see whether the benchmarking tool changes the ranking or only changes absolute metric values.</li>
      <li><strong>RTGtools was removed from this report</strong> because the hg38 RTGtools summary returned zero values for all callers, so it does not contribute useful comparative interpretation for this run.</li>
      <li><strong>SURVIVOR-style plots are included</strong> to check whether low recall is explained by VCF composition: missing SV types, heavy BND/NA events, very different SVLEN distributions, or filter differences between callers and truth.</li>
    </ul>
    """

def files_found_html():
    if not survivor_stat_paths:
        return "<p>No nf-core SURVIVOR stats files were found under benchmark-runs/nfcore/hg38/structural.</p>"

    items = "\n".join([f"<li><code>{p}</code></li>" for p in survivor_stat_paths])
    return f"<ul>{items}</ul>"

plot_sections = ""

plot_sections += img("hg38_recall_by_caller_and_tool.png", "Recall by caller and benchmarking tool")
plot_sections += img("hg38_precision_by_caller_and_tool.png", "Precision by caller and benchmarking tool")
plot_sections += img("hg38_f1_by_caller_and_tool.png", "F1 by caller and benchmarking tool")
plot_sections += img("hg38_precision_recall_f1_together.png", "Precision, recall and F1 together")
plot_sections += img("hg38_false_negatives_by_caller_and_tool.png", "False negatives by caller and benchmarking tool")

plot_sections += img("hg38_precision_vs_recall.png", "Precision vs recall")
plot_sections += img("hg38_precision_vs_f1.png", "Precision vs F1")
plot_sections += img("hg38_recall_vs_f1.png", "Recall vs F1")

if (PLOTS / "hg38_svbenchmark_minus_truvari_differences.png").exists():
    plot_sections += img("hg38_svbenchmark_minus_truvari_differences.png", "Benchmarking-tool discrepancy: SVbenchmark minus Truvari")

for method in METHOD_ORDER:
    plot_sections += img(f"hg38_{method.lower()}_precision_recall_f1_by_caller.png", f"{method}: precision, recall and F1 by caller")

for caller in CALLER_ORDER:
    plot_sections += img(f"hg38_{caller}_precision_recall_f1_by_tool.png", f"{CALLER_LABEL[caller]}: precision, recall and F1 by tool")

survivor_plot_sections = ""
survivor_plot_sections += img("hg38_survivor_svtype_counts_by_vcf.png", "SURVIVOR-style SVTYPE counts by VCF")
survivor_plot_sections += img("hg38_survivor_svtype_fraction_stacked.png", "SURVIVOR-style SVTYPE composition fractions")
survivor_plot_sections += img("hg38_survivor_svlen_bin_counts.png", "SURVIVOR-style SV length bin counts")
survivor_plot_sections += img("hg38_survivor_svlen_composition_stacked.png", "SURVIVOR-style SV length composition")
survivor_plot_sections += img("hg38_filter_status_counts.png", "VCF FILTER status counts")
survivor_plot_sections += img("hg38_total_pass_nonpass_records.png", "Total, PASS and non-PASS records")

# -----------------------------
# 7. Final HTML
# -----------------------------

html_text = f"""<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="UTF-8">
<title>hg38 SV benchmarking complete analysis</title>
<style>
:root {{
  --bg:#0d1117;
  --surface:#161b22;
  --border:#30363d;
  --accent:#58a6ff;
  --green:#3fb950;
  --warn:#d29922;
  --danger:#f85149;
  --muted:#8b949e;
  --text:#e6edf3;
  --subtext:#c9d1d9;
  --highlight:#1f2937;
}}
* {{ box-sizing:border-box; }}
body {{
  margin:0;
  background:var(--bg);
  color:var(--text);
  font-family:Arial, sans-serif;
  font-size:14px;
  line-height:1.6;
}}
.header {{
  padding:36px 48px 28px;
  border-bottom:1px solid var(--border);
  background:linear-gradient(135deg,#0d1117 0%,#161b22 60%,#1c2a3a 100%);
}}
.header h1 {{
  margin:0;
  color:var(--accent);
  font-size:30px;
}}
.meta {{
  color:var(--muted);
  margin-top:8px;
}}
.badge {{
  display:inline-block;
  color:var(--accent);
  border:1px solid rgba(88,166,255,.3);
  background:rgba(88,166,255,.1);
  padding:2px 8px;
  border-radius:4px;
  font-size:11px;
  font-weight:bold;
  text-transform:uppercase;
  margin-right:6px;
}}
.container {{
  max-width:1400px;
  margin:0 auto;
  padding:34px 48px 70px;
}}
.section {{
  margin-bottom:44px;
}}
.section-label {{
  color:var(--muted);
  text-transform:uppercase;
  letter-spacing:2px;
  border-bottom:1px solid var(--border);
  padding-bottom:7px;
  margin-bottom:18px;
  font-size:12px;
  font-weight:bold;
}}
.note {{
  background:rgba(88,166,255,.07);
  border-left:4px solid var(--accent);
  padding:14px 18px;
  margin-bottom:18px;
  color:var(--subtext);
}}
.warning {{
  background:rgba(248,81,73,.08);
  border-left:4px solid var(--danger);
  padding:14px 18px;
  margin-bottom:18px;
  color:var(--subtext);
}}
.cards {{
  display:grid;
  grid-template-columns:repeat(auto-fit,minmax(170px,1fr));
  gap:12px;
}}
.card {{
  background:var(--surface);
  border:1px solid var(--border);
  border-radius:8px;
  padding:16px;
}}
.value {{
  font-size:26px;
  font-weight:bold;
  color:var(--accent);
}}
.value.good {{ color:var(--green); }}
.value.warn {{ color:var(--warn); }}
.value.bad {{ color:var(--danger); }}
.label {{
  color:var(--muted);
  font-size:11px;
  text-transform:uppercase;
  letter-spacing:1px;
}}
.caller-block {{
  margin-bottom:24px;
}}
.caller-block h3 {{
  color:var(--subtext);
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
  padding:16px;
}}
.plot-card h3 {{
  margin-top:0;
  color:var(--subtext);
  font-size:15px;
}}
.plot-card img {{
  width:100%;
  max-width:100%;
  background:white;
  border:1px solid var(--border);
  border-radius:4px;
}}
.tbl {{
  width:100%;
  border-collapse:collapse;
  background:var(--surface);
  margin-bottom:24px;
}}
.tbl th {{
  background:var(--highlight);
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
}}
code {{
  background:var(--highlight);
  color:var(--subtext);
  padding:2px 5px;
  border-radius:3px;
}}
.file-list {{
  font-size:12px;
  color:var(--subtext);
}}
hr {{
  border:0;
  border-top:1px solid var(--border);
  margin:36px 0;
}}
</style>
</head>
<body>

<div class="header">
  <span class="badge">nf-core/variantbenchmarking</span>
  <span class="badge">hg38</span>
  <span class="badge">GIAB HG002</span>
  <span class="badge">Truvari + SVbenchmark</span>
  <h1>hg38 SV benchmarking complete analysis</h1>
  <div class="meta">
    Generated: {datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")} |
    Callers: Delly · DRAGEN · Sniffles |
    Truth: GIAB HG002 GRCh38 v5.0q |
    RTGtools removed because it returned zero metrics for this run.
  </div>
</div>

<div class="container">

<div class="section">
  <div class="section-label">Executive summary</div>
  <div class="warning">
    This report reads actual files from your local paths, not pasted values. It compares caller results against the GIAB truth set using Truvari and SVbenchmark, and includes SURVIVOR-style VCF composition plots to explain why recall looks low.
  </div>
  <div class="note">
    Input metric files:
    <br><code>{TRUVARI_CSV}</code>
    <br><code>{SVBENCHMARK_CSV}</code>
    <br><br>
    VCFs checked:
    <br><code>{VCFS["Truth"]}</code>
    <br><code>{VCFS["Delly"]}</code>
    <br><code>{VCFS["DRAGEN"]}</code>
    <br><code>{VCFS["Sniffles"]}</code>
  </div>
  {conclusion_text()}
</div>

<hr>

<div class="section">
  <div class="section-label">Per-caller mean metrics across benchmarking tools</div>
  {metric_cards()}
</div>

<hr>

<div class="section">
  <div class="section-label">Benchmarking metric plots</div>
  <div class="plot-grid">
    {plot_sections}
  </div>
</div>

<hr>

<div class="section">
  <div class="section-label">SURVIVOR-style VCF composition plots</div>
  <div class="note">
    These plots are generated directly from the actual truth/caller VCFs using <code>bcftools query</code>. They are intended to reproduce the kind of SV composition checks that help explain benchmarking discrepancies: SVTYPE distribution, SVLEN distribution, FILTER status, and total record counts.
  </div>
  <div class="plot-grid">
    {survivor_plot_sections}
  </div>
</div>

<hr>

<div class="section">
  <div class="section-label">Tables: precision, recall and F1</div>
  {df_to_html(metrics)}
</div>

<div class="section">
  <div class="section-label">Tables: mean metrics by caller</div>
  {df_to_html(mean_by_caller)}
</div>

<div class="section">
  <div class="section-label">Tables: SVbenchmark minus Truvari discrepancy</div>
  {df_to_html(delta_df)}
</div>

<hr>

<div class="section">
  <div class="section-label">Tables: SURVIVOR-style VCF checks</div>
  <h3>VCF totals</h3>
  {df_to_html(vcf_totals)}
  <h3>SVTYPE counts</h3>
  {df_to_html(svtype_counts)}
  <h3>SVLEN bin counts</h3>
  {df_to_html(svlen_counts)}
  <h3>FILTER counts</h3>
  {df_to_html(filter_counts)}
</div>

<hr>

<div class="section">
  <div class="section-label">SURVIVOR files found by nf-core</div>
  <div class="file-list">
    {files_found_html()}
  </div>
</div>

<hr>

<div class="section">
  <div class="section-label">Files generated by this report</div>
  <ul>
    <li><code>{TABLES / "hg38_truvari_svbenchmark_metrics.csv"}</code></li>
    <li><code>{TABLES / "hg38_mean_metrics_by_caller.csv"}</code></li>
    <li><code>{TABLES / "hg38_method_discrepancy_svbenchmark_minus_truvari.csv"}</code></li>
    <li><code>{TABLES / "hg38_survivor_style_svtype_counts.csv"}</code></li>
    <li><code>{TABLES / "hg38_survivor_style_svlen_bin_counts.csv"}</code></li>
    <li><code>{TABLES / "hg38_filter_counts.csv"}</code></li>
    <li><code>{LOGS / "hg38_survivor_stats_files_found.txt"}</code></li>
  </ul>
</div>

</div>
</body>
</html>
"""

REPORT = HTML / "hg38_complete_sv_benchmarking_analysis.html"
REPORT.write_text(html_text)

print()
print("DONE")
print("HTML report:")
print(REPORT)
print()
print("Plots:")
for p in sorted(PLOTS.glob("*.png")):
    print(p)
print()
print("Tables:")
for p in sorted(TABLES.glob("*.csv")):
    print(p)
