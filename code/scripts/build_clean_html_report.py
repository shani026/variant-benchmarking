from pathlib import Path
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import html
import textwrap

BASE = Path(__file__).resolve().parents[2]
QC = BASE / "reports/chr22" / "qc_tables"
MET = BASE / "reports/chr22" / "metrics"
OUTDIR = BASE / "reports/chr22"
PLOTDIR = OUTDIR / "plots_clean"
OUTHTML = OUTDIR / "final_chr22_sv_benchmark_report_clean.html"

PLOTDIR.mkdir(parents=True, exist_ok=True)
OUTDIR.mkdir(parents=True, exist_ok=True)

# ----------------------------
# Helpers
# ----------------------------
def read_tsv(path):
    path = Path(path)
    if not path.exists():
        return None
    try:
        return pd.read_csv(path, sep="\t")
    except Exception:
        return None

def read_text(path):
    path = Path(path)
    if not path.exists():
        return None
    return path.read_text(errors="replace")

def safe_numeric(series):
    return pd.to_numeric(series, errors="coerce")

def norm_cols(df):
    d = df.copy()
    d.columns = [str(c).strip() for c in d.columns]
    return d

def df_to_html(df, index=False):
    if df is None:
        return "<p><em>Not available.</em></p>"
    return df.to_html(index=index, classes="dataframe", border=0)

def pre_text(txt):
    if txt is None:
        return "<p><em>Not available.</em></p>"
    return f"<pre>{html.escape(txt)}</pre>"

def save_fig(path):
    plt.tight_layout()
    plt.savefig(path, dpi=200, bbox_inches="tight")
    plt.close()

def find_first_existing(paths):
    for p in paths:
        if Path(p).exists():
            return Path(p)
    return None

def extract_build(run):
    run = str(run)
    if run.startswith("hg19"):
        return "hg19"
    if run.startswith("hg38"):
        return "hg38"
    return "unknown"

def extract_caller(run):
    run = str(run).lower()
    if "sniffles" in run:
        return "Sniffles"
    if "dragen" in run:
        return "DRAGEN"
    if "delly" in run:
        return "DELLY"
    return run

# ----------------------------
# Load main tables
# ----------------------------
truvari_main = read_tsv(MET / "truvari_summary_extended.tsv")
if truvari_main is None:
    truvari_main = read_tsv(MET / "truvari_summary.tsv")
if truvari_main is not None:
    truvari_main = norm_cols(truvari_main)

truvari_svtype = read_tsv(MET / "truvari_svtype_tp_fp_fn_counts.tsv")
if truvari_svtype is not None:
    truvari_svtype = norm_cols(truvari_svtype)

truvari_compare_path = find_first_existing([
    MET / "truvari_original_vs_converted_dup_to_ins.tsv",
    MET / "truvari_original_vs_dup_to_ins.tsv"
])
truvari_compare = read_tsv(truvari_compare_path) if truvari_compare_path else None
if truvari_compare is not None:
    truvari_compare = norm_cols(truvari_compare)

bnd_counts = read_tsv(QC / "BND_counts.tsv")
if bnd_counts is not None:
    bnd_counts = norm_cols(bnd_counts)

end_lt_pos = read_tsv(QC / "END_lt_POS_counts.tsv")
if end_lt_pos is not None:
    end_lt_pos = norm_cols(end_lt_pos)

missing_end_svlen = read_tsv(QC / "missing_END_SVLEN_counts.tsv")
if missing_end_svlen is not None:
    missing_end_svlen = norm_cols(missing_end_svlen)

raw_chr22 = read_tsv(QC / "raw_chr22_record_summary.tsv")
if raw_chr22 is not None:
    raw_chr22 = norm_cols(raw_chr22)

benchmark_interval_counts_txt = read_text(QC / "benchmark_ready_interval_counts.txt")
benchmark_dup_to_ins_counts_txt = read_text(QC / "benchmark_ready_dup_to_ins_counts.txt")

# ----------------------------
# Clean / enrich Truvari tables
# ----------------------------
if truvari_main is not None:
    for col in ["precision", "recall", "f1", "TP_base", "TP_comp", "FP", "FN", "base_cnt_evaluated", "comp_cnt_evaluated", "gt_concordance"]:
        if col in truvari_main.columns:
            truvari_main[col] = safe_numeric(truvari_main[col])

    if "run" in truvari_main.columns:
        truvari_main["build"] = truvari_main["run"].apply(extract_build)
        truvari_main["caller"] = truvari_main["run"].apply(extract_caller)

if truvari_compare is not None:
    for col in truvari_compare.columns:
        if col != "run":
            truvari_compare[col] = safe_numeric(truvari_compare[col]) if "run" not in col else truvari_compare[col]
    truvari_compare["build"] = truvari_compare["run"].apply(extract_build)
    truvari_compare["caller"] = truvari_compare["run"].apply(extract_caller)

# ----------------------------
# Plot 1: Precision / Recall / F1 original
# ----------------------------
plot_metric_original = PLOTDIR / "truvari_original_metrics.png"
if truvari_main is not None and {"run", "precision", "recall", "f1"}.issubset(truvari_main.columns):
    df = truvari_main.copy().sort_values(["build", "caller"])
    x = np.arange(len(df))
    w = 0.25

    plt.figure(figsize=(12, 6))
    plt.bar(x - w, df["precision"], width=w, label="Precision")
    plt.bar(x,     df["recall"],    width=w, label="Recall")
    plt.bar(x + w, df["f1"],        width=w, label="F1")
    plt.xticks(x, df["run"], rotation=45, ha="right")
    plt.ylim(0, 1.05)
    plt.ylabel("Metric value")
    plt.title("Truvari original benchmark: Precision, Recall, and F1")
    plt.legend()
    save_fig(plot_metric_original)

# ----------------------------
# Plot 2: Original vs DUP-to-INS precision
# Plot 3: Original vs DUP-to-INS recall
# Plot 4: Original vs DUP-to-INS F1
# ----------------------------
plot_prec_compare = PLOTDIR / "precision_original_vs_dup_to_ins.png"
plot_recall_compare = PLOTDIR / "recall_original_vs_dup_to_ins.png"
plot_f1_compare = PLOTDIR / "f1_original_vs_dup_to_ins.png"

if truvari_compare is not None:
    df = truvari_compare.copy().sort_values(["build", "caller"])
    x = np.arange(len(df))
    w = 0.35

    if {"original_precision", "dup_to_ins_precision"}.issubset(df.columns):
        plt.figure(figsize=(12, 6))
        plt.bar(x - w/2, df["original_precision"], width=w, label="Original")
        plt.bar(x + w/2, df["dup_to_ins_precision"], width=w, label="DUP-to-INS")
        plt.xticks(x, df["run"], rotation=45, ha="right")
        plt.ylim(0, 1.05)
        plt.ylabel("Precision")
        plt.title("Precision: original vs DUP-to-INS sensitivity analysis")
        plt.legend()
        save_fig(plot_prec_compare)

    if {"original_recall", "dup_to_ins_recall"}.issubset(df.columns):
        plt.figure(figsize=(12, 6))
        plt.bar(x - w/2, df["original_recall"], width=w, label="Original")
        plt.bar(x + w/2, df["dup_to_ins_recall"], width=w, label="DUP-to-INS")
        plt.xticks(x, df["run"], rotation=45, ha="right")
        plt.ylim(0, 1.05)
        plt.ylabel("Recall")
        plt.title("Recall: original vs DUP-to-INS sensitivity analysis")
        plt.legend()
        save_fig(plot_recall_compare)

    if {"original_f1", "dup_to_ins_f1"}.issubset(df.columns):
        plt.figure(figsize=(12, 6))
        plt.bar(x - w/2, df["original_f1"], width=w, label="Original")
        plt.bar(x + w/2, df["dup_to_ins_f1"], width=w, label="DUP-to-INS")
        plt.xticks(x, df["run"], rotation=45, ha="right")
        plt.ylim(0, 1.05)
        plt.ylabel("F1")
        plt.title("F1: original vs DUP-to-INS sensitivity analysis")
        plt.legend()
        save_fig(plot_f1_compare)

# ----------------------------
# Plot 5: Precision vs Recall scatter with arrows
# ----------------------------
plot_scatter = PLOTDIR / "precision_recall_scatter_original_vs_dup_to_ins.png"
if truvari_compare is not None and {
    "original_precision", "original_recall",
    "dup_to_ins_precision", "dup_to_ins_recall",
    "run"
}.issubset(truvari_compare.columns):
    df = truvari_compare.copy()

    plt.figure(figsize=(9, 8))
    plt.scatter(df["original_recall"], df["original_precision"], label="Original", s=80)
    plt.scatter(df["dup_to_ins_recall"], df["dup_to_ins_precision"], label="DUP-to-INS", s=80, marker="x")

    for _, row in df.iterrows():
        x1, y1 = row["original_recall"], row["original_precision"]
        x2, y2 = row["dup_to_ins_recall"], row["dup_to_ins_precision"]
        plt.annotate("", xy=(x2, y2), xytext=(x1, y1),
                     arrowprops=dict(arrowstyle="->", lw=1))
        plt.text(x2, y2, str(row["run"]), fontsize=8)

    plt.xlim(0, 1.05)
    plt.ylim(0, 1.05)
    plt.xlabel("Recall")
    plt.ylabel("Precision")
    plt.title("Precision vs Recall: movement from original to DUP-to-INS")
    plt.legend()
    save_fig(plot_scatter)

# ----------------------------
# Plot 6: Delta F1
# ----------------------------
plot_delta_f1 = PLOTDIR / "delta_f1_dup_to_ins_minus_original.png"
if truvari_compare is not None and {"run", "delta_f1"}.issubset(truvari_compare.columns):
    df = truvari_compare.copy().sort_values(["build", "caller"])
    plt.figure(figsize=(12, 6))
    plt.bar(df["run"], df["delta_f1"])
    plt.axhline(0, linestyle="--", linewidth=1)
    plt.xticks(rotation=45, ha="right")
    plt.ylabel("ΔF1 (DUP-to-INS - Original)")
    plt.title("Effect of DUP-to-INS conversion on F1")
    save_fig(plot_delta_f1)

# ----------------------------
# Plot 7: BND counts
# ----------------------------
plot_bnd = PLOTDIR / "bnd_counts.png"
if bnd_counts is not None:
    cols = {c.lower(): c for c in bnd_counts.columns}
    file_col = cols.get("file", bnd_counts.columns[1] if len(bnd_counts.columns) > 1 else bnd_counts.columns[0])
    count_col = None
    for c in bnd_counts.columns:
        if "count" in c.lower():
            count_col = c
            break
    if count_col is None:
        count_col = bnd_counts.columns[-1]

    tmp = bnd_counts.copy()
    tmp[count_col] = safe_numeric(tmp[count_col])

    plt.figure(figsize=(12, 6))
    plt.bar(tmp[file_col], tmp[count_col])
    plt.xticks(rotation=60, ha="right")
    plt.ylabel("BND count")
    plt.title("BND records in raw chr22 files")
    save_fig(plot_bnd)

# ----------------------------
# Plot 8: END<POS
# ----------------------------
plot_endltpos = PLOTDIR / "end_lt_pos_counts.png"
if end_lt_pos is not None:
    count_col = None
    file_col = None
    for c in end_lt_pos.columns:
        if "file" in c.lower():
            file_col = c
        if "count" in c.lower() or "bad" in c.lower():
            count_col = c
    if file_col is None:
        file_col = end_lt_pos.columns[1] if len(end_lt_pos.columns) > 1 else end_lt_pos.columns[0]
    if count_col is None:
        count_col = end_lt_pos.columns[-1]

    tmp = end_lt_pos.copy()
    tmp[count_col] = safe_numeric(tmp[count_col])

    plt.figure(figsize=(12, 6))
    plt.bar(tmp[file_col], tmp[count_col])
    plt.xticks(rotation=60, ha="right")
    plt.ylabel("Count")
    plt.title("Records with END < POS")
    save_fig(plot_endltpos)

# ----------------------------
# Plot 9: Missing END / Missing SVLEN
# ----------------------------
plot_missing = PLOTDIR / "missing_end_svlen_counts.png"
if missing_end_svlen is not None:
    cols = list(missing_end_svlen.columns)
    file_col = None
    end_col = None
    svlen_col = None

    for c in cols:
        lc = c.lower()
        if "file" in lc:
            file_col = c
        if "missing" in lc and "end" in lc:
            end_col = c
        if "missing" in lc and "svlen" in lc:
            svlen_col = c

    if file_col and end_col and svlen_col:
        tmp = missing_end_svlen.copy()
        tmp[end_col] = safe_numeric(tmp[end_col])
        tmp[svlen_col] = safe_numeric(tmp[svlen_col])

        x = np.arange(len(tmp))
        w = 0.35

        plt.figure(figsize=(12, 6))
        plt.bar(x - w/2, tmp[end_col], width=w, label="Missing END")
        plt.bar(x + w/2, tmp[svlen_col], width=w, label="Missing SVLEN")
        plt.xticks(x, tmp[file_col], rotation=60, ha="right")
        plt.ylabel("Count")
        plt.title("Missing END and Missing SVLEN counts")
        plt.legend()
        save_fig(plot_missing)

# ----------------------------
# Summary text
# ----------------------------
summary_items = []

if truvari_main is not None and "f1" in truvari_main.columns:
    best = truvari_main.loc[truvari_main["f1"].idxmax()]
    summary_items.append(
        f"<li><strong>Best original Truvari F1:</strong> {html.escape(str(best['run']))} "
        f"(F1 = {best['f1']:.3f})</li>"
    )

if truvari_compare is not None and "delta_f1" in truvari_compare.columns:
    tmp = truvari_compare.dropna(subset=["delta_f1"]).copy()
    if not tmp.empty:
        best_gain = tmp.loc[tmp["delta_f1"].idxmax()]
        worst_gain = tmp.loc[tmp["delta_f1"].idxmin()]
        summary_items.append(
            f"<li><strong>Largest DUP-to-INS F1 increase:</strong> {html.escape(str(best_gain['run']))} "
            f"(ΔF1 = {best_gain['delta_f1']:.4f})</li>"
        )
        summary_items.append(
            f"<li><strong>Smallest / most negative DUP-to-INS F1 change:</strong> {html.escape(str(worst_gain['run']))} "
            f"(ΔF1 = {worst_gain['delta_f1']:.4f})</li>"
        )

if bnd_counts is not None:
    count_col = None
    for c in bnd_counts.columns:
        if "count" in c.lower():
            count_col = c
            break
    if count_col is None:
        count_col = bnd_counts.columns[-1]
    total_bnd = safe_numeric(bnd_counts[count_col]).fillna(0).sum()
    summary_items.append(
        f"<li><strong>Total BND records counted separately:</strong> {int(total_bnd)}</li>"
    )

summary_html = "<ul>" + "\n".join(summary_items) + "</ul>" if summary_items else "<p>No summary statistics available.</p>"

# ----------------------------
# Plot figure blocks
# ----------------------------
def img_block(path, caption):
    path = Path(path)
    if not path.exists():
        return f"<p><em>{html.escape(caption)}: plot not available.</em></p>"
    return f"""
    <div class="figure">
      <img src="{path.relative_to(OUTDIR).as_posix()}" alt="{html.escape(caption)}">
      <div class="caption">{html.escape(caption)}</div>
    </div>
    """

# ----------------------------
# Build HTML
# ----------------------------
html_parts = []

html_parts.append(f"""
<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="utf-8">
<title>Chromosome 22 SV Benchmarking Report</title>
<style>
body {{
    font-family: Arial, Helvetica, sans-serif;
    margin: 40px;
    line-height: 1.5;
    color: #222;
    background: #fff;
}}
h1, h2, h3 {{
    color: #1f2d3d;
}}
pre {{
    background: #f7f7f7;
    border: 1px solid #ddd;
    padding: 12px;
    overflow-x: auto;
    white-space: pre-wrap;
}}
table.dataframe {{
    border-collapse: collapse;
    width: 100%;
    margin: 10px 0 20px 0;
    font-size: 14px;
}}
table.dataframe th, table.dataframe td {{
    border: 1px solid #ccc;
    padding: 6px 8px;
    text-align: left;
}}
table.dataframe th {{
    background: #efefef;
}}
.figure {{
    margin: 25px 0;
    text-align: center;
}}
.figure img {{
    max-width: 100%;
    border: 1px solid #ddd;
}}
.caption {{
    font-size: 14px;
    color: #444;
    margin-top: 8px;
}}
.note {{
    background: #f8fbff;
    border-left: 5px solid #4a90e2;
    padding: 12px;
    margin: 15px 0;
}}
.small {{
    font-size: 13px;
    color: #555;
}}
</style>
</head>
<body>

<h1>Chromosome 22 Structural Variant Benchmarking Report</h1>

<h2>1. Aim</h2>
<p>
This report summarises the chromosome 22 benchmarking workflow for HG002 structural variant calls from
DRAGEN, Sniffles, and DELLY against the truth set in both hg19 / GRCh37 and hg38 / GRCh38.
The report focuses on diagnosis, preprocessing, benchmark interpretation, and DUP-to-INS sensitivity analysis.
</p>

<h2>2. Key findings</h2>
{summary_html}

<div class="note">
<strong>Main interpretation:</strong>
The core benchmark should be interpreted as an <strong>interval-SV benchmark</strong>. BND records,
END&lt;POS issues, and missing END / SVLEN issues were counted and visualised separately so the reader can see
what was excluded and why.
</div>

<h2>3. Diagnosis</h2>
<p>
The original chr22 subsetting problem came from chromosome naming mismatch:
</p>
<ul>
<li><strong>hg19 / GRCh37:</strong> chromosome label <code>22</code></li>
<li><strong>hg38 / GRCh38:</strong> chromosome label <code>chr22</code></li>
</ul>

<h2>4. Benchmark scope</h2>
<p>The main interval-SV benchmark includes:</p>
<ul>
<li>DEL</li>
<li>DUP</li>
<li>INS</li>
<li>INV</li>
</ul>
<p>The following are handled separately and interpreted as exclusions / QC context rather than direct interval benchmark targets:</p>
<ul>
<li>BND</li>
<li>UNKNOWN</li>
<li>NA</li>
<li>missing SVTYPE</li>
<li>malformed interval records</li>
</ul>

<h2>5. Selected QC tables</h2>

<h3>5.1 Raw chr22 record summary</h3>
{df_to_html(raw_chr22)}

<h3>5.2 Benchmark-ready interval counts</h3>
{pre_text(benchmark_interval_counts_txt)}

<h3>5.3 BND counts</h3>
{df_to_html(bnd_counts)}

<h3>5.4 END &lt; POS counts</h3>
{df_to_html(end_lt_pos)}

<h3>5.5 Missing END / SVLEN counts</h3>
{df_to_html(missing_end_svlen)}

<h2>6. QC interpretation through plots</h2>
<p>
These plots explain what was counted separately and why those records matter for interpretation.
</p>

{img_block(plot_bnd, "BND records counted separately across raw chr22 files")}
{img_block(plot_endltpos, "Records with END < POS")}
{img_block(plot_missing, "Missing END and Missing SVLEN counts")}

<h2>7. Truvari main benchmark</h2>
<p>
Truvari was used as the main structural variant benchmarking tool.
The table below shows the original benchmark summary metrics.
</p>

{df_to_html(truvari_main)}

{img_block(plot_metric_original, "Original Truvari benchmark: precision, recall, and F1")}

<h3>7.1 Per-SVTYPE Truvari counts</h3>
{df_to_html(truvari_svtype)}

<h2>8. DUP-to-INS sensitivity analysis</h2>
<p>
Because tandem duplications may be represented either as <code>SVTYPE=DUP</code> or as inserted duplicated sequence
represented as <code>SVTYPE=INS</code>, a sensitivity analysis was carried out.
The main question is whether changing DUP representation to INS changes precision, recall, or F1.
</p>

<h3>8.1 DUP-to-INS converted input counts</h3>
{pre_text(benchmark_dup_to_ins_counts_txt)}

<h3>8.2 Original vs DUP-to-INS comparison table</h3>
{df_to_html(truvari_compare)}

<h2>9. Metric comparison plots</h2>
<p>
These plots show how the benchmark changes when DUP calls are treated as INS in the sensitivity analysis.
</p>

{img_block(plot_prec_compare, "Precision: original vs DUP-to-INS")}
{img_block(plot_recall_compare, "Recall: original vs DUP-to-INS")}
{img_block(plot_f1_compare, "F1: original vs DUP-to-INS")}
{img_block(plot_scatter, "Precision vs Recall: movement from original to DUP-to-INS")}
{img_block(plot_delta_f1, "ΔF1 (DUP-to-INS minus original)")}

<h2>10. Interpretation</h2>
<p>
The figures above should be used to interpret both performance and representation sensitivity:
</p>
<ul>
<li><strong>Precision / recall / F1 plots</strong> show overall caller performance in the original benchmark.</li>
<li><strong>Original vs DUP-to-INS plots</strong> show whether differences are due to biology / caller performance or due to representation.</li>
<li><strong>BND plot</strong> shows the number of records that were counted separately because they are not directly comparable as interval SVs.</li>
<li><strong>END&lt;POS and missing END / SVLEN plots</strong> show why some records needed QC handling before benchmarking.</li>
</ul>

<p>
In practical terms:
</p>
<ul>
<li>If DUP-to-INS barely changes the metrics, then the original results are stable to representation differences.</li>
<li>If DUP-to-INS clearly improves recall or F1 for a caller, then some apparent misses were likely due to DUP vs INS representation rather than completely wrong calls.</li>
<li>A high BND burden reminds the reader that not all structural variants fit cleanly into interval-style benchmarking.</li>
</ul>

<h2>11. Main conclusion</h2>
<p>
The chromosome 22 workflow was successfully debugged and benchmarked.
The original empty chr22 subsets were caused by chromosome naming mismatch, and the later analyses showed:
</p>
<ul>
<li>which variants entered the interval benchmark,</li>
<li>which categories were handled separately,</li>
<li>how original caller performance compared across runs, and</li>
<li>whether DUP-to-INS conversion materially changed interpretation.</li>
</ul>

<h2>12. Limitations</h2>
<ul>
<li>This is a chromosome 22 subset, not a whole-genome benchmark.</li>
<li>The main benchmark is interval-SV focused.</li>
<li>BND and complex events require separate interpretation.</li>
<li>Representation differences (especially DUP vs INS) can affect apparent performance.</li>
</ul>

<p class="small">
Generated from files already present in <code>reports/chr22/qc_tables/</code> and <code>reports/chr22/metrics/</code>.
</p>

</body>
</html>
""")

OUTHTML.write_text("".join(html_parts), encoding="utf-8")
print(f"Wrote {OUTHTML}")
print(f"Plots saved in {PLOTDIR}")