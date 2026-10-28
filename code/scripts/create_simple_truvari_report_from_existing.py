#!/usr/bin/env python3

import csv
import html
import json
import math
import subprocess
from collections import Counter, defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
TRUVARI_DIR = ROOT / "benchmark-runs/chr22" / "truvari"
OUT_DIR = ROOT / "reports/chr22" / "truvari_simple_report"
PLOTS_DIR = OUT_DIR / "plots"
TABLES_DIR = OUT_DIR / "tables"
HTML_OUT = OUT_DIR / "simple_truvari_report.html"

EXPECTED_RUNS = [
    "hg19_DELLY",
    "hg19_dragen",
    "hg19_sniffles",
    "hg38_delly",
    "hg38_dragen",
    "hg38_sniffles",
]

COLORS = [
    "#4C78A8", "#F58518", "#54A24B", "#E45756",
    "#72B7B2", "#B279A2", "#FF9DA6", "#9D755D",
    "#BAB0AC", "#A0CBE8", "#FFBE7D", "#8CD17D"
]


def ensure_dirs():
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    PLOTS_DIR.mkdir(parents=True, exist_ok=True)
    TABLES_DIR.mkdir(parents=True, exist_ok=True)


def clean_key(key):
    return str(key).strip().lower().replace("-", "_").replace(" ", "_").replace(".", "_")


def load_json(path):
    if not path.exists():
        return {}
    with open(path) as handle:
        data = json.load(handle)
    return {clean_key(k): v for k, v in data.items()}


def get_metric(summary, names):
    for name in names:
        key = clean_key(name)
        if key in summary:
            return summary[key]
    return None


def to_float(x):
    if x is None:
        return None
    try:
        return float(x)
    except Exception:
        return None


def to_int(x):
    if x is None:
        return None
    try:
        return int(float(x))
    except Exception:
        return None


def run_cmd(cmd):
    try:
        res = subprocess.run(
            cmd,
            check=True,
            text=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
        )
        return res.stdout
    except Exception:
        return ""


def bcftools_query(vcf, fmt):
    vcf = Path(vcf)
    if not vcf.exists():
        return []
    out = run_cmd(["bcftools", "query", "-f", fmt, str(vcf)])
    return [x for x in out.splitlines() if x.strip()]


def count_records(vcf):
    vcf = Path(vcf)
    if not vcf.exists():
        return 0
    out = run_cmd(["bcftools", "view", "-H", str(vcf)])
    if not out:
        return 0
    return len(out.splitlines())


def svtype_counts(vcf):
    counts = Counter()
    for line in bcftools_query(vcf, "%INFO/SVTYPE\n"):
        svtype = line.strip()
        if svtype in {"", "."}:
            svtype = "NA"
        counts[svtype] += 1
    return counts


def chrom_counts(vcf):
    counts = Counter()
    for line in bcftools_query(vcf, "%CHROM\n"):
        chrom = line.strip()
        if chrom:
            counts[chrom] += 1
    return counts


def size_values(vcf):
    values = []
    for line in bcftools_query(vcf, "%POS\t%INFO/END\t%INFO/SVLEN\n"):
        parts = line.split("\t")
        if len(parts) != 3:
            continue

        pos, end, svlen = parts
        size = None

        if svlen not in {"", "."}:
            try:
                size = abs(int(float(svlen.split(",")[0])))
            except Exception:
                size = None

        if size is None and end not in {"", "."}:
            try:
                size = abs(int(end) - int(pos))
            except Exception:
                size = None

        if size and size > 0:
            values.append(size)

    return values


def parse_run_name(run):
    build, caller = run.split("_", 1)
    caller_l = caller.lower()

    if caller_l == "delly":
        caller = "DELLY"
    elif caller_l == "dragen":
        caller = "DRAGEN"
    elif caller_l == "sniffles":
        caller = "Sniffles"

    return build, caller


def collect_results():
    metric_rows = []
    svtype_rows = []
    chrom_rows = []
    size_rows = []
    found = []
    missing = []

    for run in EXPECTED_RUNS:
        run_dir = TRUVARI_DIR / run

        if not run_dir.exists():
            missing.append(run)
            continue

        found.append(run)
        build, caller = parse_run_name(run)
        label = f"{build} {caller}"

        summary = load_json(run_dir / "summary.json")

        row = {
            "run": run,
            "build": build,
            "caller": caller,
            "label": label,
            "base_count": to_int(get_metric(summary, ["base cnt", "base_count", "basecnt"])),
            "comparison_count": to_int(get_metric(summary, ["call cnt", "call_count", "comp count", "comp_count", "callcnt"])),
            "tp_base": to_int(get_metric(summary, ["TP-base", "tp_base"])),
            "tp_comp": to_int(get_metric(summary, ["TP-call", "TP-comp", "tp_call", "tp_comp"])),
            "fp": to_int(get_metric(summary, ["FP", "fp"])),
            "fn": to_int(get_metric(summary, ["FN", "fn"])),
            "precision": to_float(get_metric(summary, ["precision"])),
            "recall": to_float(get_metric(summary, ["recall"])),
            "f1": to_float(get_metric(summary, ["f1"])),
            "summary_json": str(run_dir / "summary.json"),
            "truvari_dir": str(run_dir),
        }

        if row["tp_base"] is None:
            row["tp_base"] = count_records(run_dir / "tp-base.vcf.gz")
        if row["tp_comp"] is None:
            row["tp_comp"] = count_records(run_dir / "tp-comp.vcf.gz")
        if row["fp"] is None:
            row["fp"] = count_records(run_dir / "fp.vcf.gz")
        if row["fn"] is None:
            row["fn"] = count_records(run_dir / "fn.vcf.gz")

        metric_rows.append(row)

        role_vcfs = {
            "TP-base": run_dir / "tp-base.vcf.gz",
            "TP-comp": run_dir / "tp-comp.vcf.gz",
            "FP": run_dir / "fp.vcf.gz",
            "FN": run_dir / "fn.vcf.gz",
        }

        for role, vcf in role_vcfs.items():
            for svtype, count in svtype_counts(vcf).items():
                svtype_rows.append({
                    "run": run,
                    "build": build,
                    "caller": caller,
                    "role": role,
                    "svtype": svtype,
                    "count": count,
                })

            if role in {"FP", "FN"}:
                for chrom, count in chrom_counts(vcf).items():
                    chrom_rows.append({
                        "run": run,
                        "build": build,
                        "caller": caller,
                        "role": role,
                        "chrom": chrom,
                        "count": count,
                    })

            for size in size_values(vcf):
                size_rows.append({
                    "run": run,
                    "build": build,
                    "caller": caller,
                    "role": role,
                    "size": size,
                })

    return metric_rows, svtype_rows, chrom_rows, size_rows, found, missing


def write_csv(rows, path):
    if not rows:
        path.write_text("")
        return

    with open(path, "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)


def svg_header(width, height):
    return [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">',
        '<rect width="100%" height="100%" fill="white"/>'
    ]


def svg_footer():
    return ["</svg>"]


def esc(x):
    return html.escape(str(x))


def nice_num(x):
    if x is None:
        return ""
    if isinstance(x, float):
        return f"{x:.4f}"
    return str(x)


def grouped_bar_svg(title, labels, series, ylabel, outpath, ymax=None):
    width = 1100
    height = 620
    left = 90
    right = 30
    top = 70
    bottom = 160
    plot_w = width - left - right
    plot_h = height - top - bottom

    n_groups = len(labels)
    n_series = len(series)

    if n_groups == 0 or n_series == 0:
        return None

    all_vals = []
    for _, vals in series:
        all_vals.extend([v for v in vals if v is not None])

    if not all_vals:
        return None

    max_val = ymax if ymax is not None else max(all_vals)
    if max_val <= 0:
        max_val = 1

    group_w = plot_w / n_groups
    bar_w = group_w * 0.7 / n_series

    svg = svg_header(width, height)
    svg.append(f'<text x="{width/2}" y="35" text-anchor="middle" font-size="22" font-weight="bold">{esc(title)}</text>')
    svg.append(f'<text x="22" y="{top + plot_h/2}" transform="rotate(-90 22 {top + plot_h/2})" text-anchor="middle" font-size="14">{esc(ylabel)}</text>')

    svg.append(f'<line x1="{left}" y1="{top}" x2="{left}" y2="{top+plot_h}" stroke="#333"/>')
    svg.append(f'<line x1="{left}" y1="{top+plot_h}" x2="{left+plot_w}" y2="{top+plot_h}" stroke="#333"/>')

    for i in range(6):
        val = max_val * i / 5
        y = top + plot_h - (val / max_val) * plot_h
        svg.append(f'<line x1="{left}" y1="{y:.1f}" x2="{left+plot_w}" y2="{y:.1f}" stroke="#eee"/>')
        svg.append(f'<text x="{left-8}" y="{y+4:.1f}" text-anchor="end" font-size="11">{val:.2g}</text>')

    for gi, label in enumerate(labels):
        gx = left + gi * group_w + group_w * 0.15
        for si, (sname, vals) in enumerate(series):
            val = vals[gi] if gi < len(vals) and vals[gi] is not None else 0
            h = (val / max_val) * plot_h
            x = gx + si * bar_w
            y = top + plot_h - h
            color = COLORS[si % len(COLORS)]
            svg.append(f'<rect x="{x:.1f}" y="{y:.1f}" width="{bar_w*0.9:.1f}" height="{h:.1f}" fill="{color}"/>')

        lx = left + gi * group_w + group_w / 2
        svg.append(f'<text x="{lx:.1f}" y="{top+plot_h+18}" text-anchor="end" font-size="11" transform="rotate(-45 {lx:.1f} {top+plot_h+18})">{esc(label)}</text>')

    legend_x = left
    legend_y = height - 40
    for si, (sname, _) in enumerate(series):
        x = legend_x + si * 150
        color = COLORS[si % len(COLORS)]
        svg.append(f'<rect x="{x}" y="{legend_y-12}" width="14" height="14" fill="{color}"/>')
        svg.append(f'<text x="{x+20}" y="{legend_y}" font-size="13">{esc(sname)}</text>')

    svg.extend(svg_footer())
    outpath.write_text("\n".join(svg))
    return outpath


def stacked_bar_svg(title, labels, categories, matrix, ylabel, outpath):
    width = 1200
    height = 700
    left = 90
    right = 220
    top = 70
    bottom = 180
    plot_w = width - left - right
    plot_h = height - top - bottom

    if not labels or not categories:
        return None

    totals = [sum(matrix[label].get(cat, 0) for cat in categories) for label in labels]
    max_val = max(totals) if totals else 1
    if max_val <= 0:
        max_val = 1

    bar_w = plot_w / len(labels) * 0.65

    svg = svg_header(width, height)
    svg.append(f'<text x="{width/2}" y="35" text-anchor="middle" font-size="22" font-weight="bold">{esc(title)}</text>')
    svg.append(f'<text x="22" y="{top + plot_h/2}" transform="rotate(-90 22 {top + plot_h/2})" text-anchor="middle" font-size="14">{esc(ylabel)}</text>')

    svg.append(f'<line x1="{left}" y1="{top}" x2="{left}" y2="{top+plot_h}" stroke="#333"/>')
    svg.append(f'<line x1="{left}" y1="{top+plot_h}" x2="{left+plot_w}" y2="{top+plot_h}" stroke="#333"/>')

    for i in range(6):
        val = max_val * i / 5
        y = top + plot_h - (val / max_val) * plot_h
        svg.append(f'<line x1="{left}" y1="{y:.1f}" x2="{left+plot_w}" y2="{y:.1f}" stroke="#eee"/>')
        svg.append(f'<text x="{left-8}" y="{y+4:.1f}" text-anchor="end" font-size="11">{val:.2g}</text>')

    for i, label in enumerate(labels):
        x = left + i * (plot_w / len(labels)) + (plot_w / len(labels) - bar_w) / 2
        y_cursor = top + plot_h

        for ci, cat in enumerate(categories):
            val = matrix[label].get(cat, 0)
            h = (val / max_val) * plot_h
            y_cursor -= h
            color = COLORS[ci % len(COLORS)]
            svg.append(f'<rect x="{x:.1f}" y="{y_cursor:.1f}" width="{bar_w:.1f}" height="{h:.1f}" fill="{color}"/>')

        lx = x + bar_w / 2
        svg.append(f'<text x="{lx:.1f}" y="{top+plot_h+18}" text-anchor="end" font-size="11" transform="rotate(-45 {lx:.1f} {top+plot_h+18})">{esc(label)}</text>')

    legend_x = width - right + 20
    legend_y = top
    for ci, cat in enumerate(categories):
        y = legend_y + ci * 22
        color = COLORS[ci % len(COLORS)]
        svg.append(f'<rect x="{legend_x}" y="{y-12}" width="14" height="14" fill="{color}"/>')
        svg.append(f'<text x="{legend_x+20}" y="{y}" font-size="12">{esc(cat)}</text>')

    svg.extend(svg_footer())
    outpath.write_text("\n".join(svg))
    return outpath


def fp_fn_balance_svg(metric_rows, outpath):
    labels = [r["label"] for r in metric_rows]
    fp = [r["fp"] or 0 for r in metric_rows]
    fn = [r["fn"] or 0 for r in metric_rows]

    max_val = max(fp + fn) if fp + fn else 1
    if max_val <= 0:
        max_val = 1

    width = 1100
    height = 620
    left = 90
    right = 30
    top = 70
    bottom = 160
    plot_w = width - left - right
    plot_h = height - top - bottom
    mid_y = top + plot_h / 2

    group_w = plot_w / len(labels)
    bar_w = group_w * 0.3

    svg = svg_header(width, height)
    svg.append(f'<text x="{width/2}" y="35" text-anchor="middle" font-size="22" font-weight="bold">FP/FN balance</text>')
    svg.append(f'<text x="22" y="{top + plot_h/2}" transform="rotate(-90 22 {top + plot_h/2})" text-anchor="middle" font-size="14">Number of SVs</text>')

    svg.append(f'<line x1="{left}" y1="{top}" x2="{left}" y2="{top+plot_h}" stroke="#333"/>')
    svg.append(f'<line x1="{left}" y1="{mid_y}" x2="{left+plot_w}" y2="{mid_y}" stroke="#333"/>')

    for i, label in enumerate(labels):
        cx = left + i * group_w + group_w / 2

        fp_h = (fp[i] / max_val) * (plot_h / 2)
        fn_h = (fn[i] / max_val) * (plot_h / 2)

        svg.append(f'<rect x="{cx-bar_w:.1f}" y="{mid_y-fp_h:.1f}" width="{bar_w:.1f}" height="{fp_h:.1f}" fill="{COLORS[1]}"/>')
        svg.append(f'<rect x="{cx:.1f}" y="{mid_y:.1f}" width="{bar_w:.1f}" height="{fn_h:.1f}" fill="{COLORS[3]}"/>')

        svg.append(f'<text x="{cx:.1f}" y="{top+plot_h+18}" text-anchor="end" font-size="11" transform="rotate(-45 {cx:.1f} {top+plot_h+18})">{esc(label)}</text>')

    svg.append(f'<rect x="{left}" y="{height-45}" width="14" height="14" fill="{COLORS[1]}"/>')
    svg.append(f'<text x="{left+20}" y="{height-33}" font-size="13">FP above line</text>')
    svg.append(f'<rect x="{left+160}" y="{height-45}" width="14" height="14" fill="{COLORS[3]}"/>')
    svg.append(f'<text x="{left+180}" y="{height-33}" font-size="13">FN below line</text>')

    svg.extend(svg_footer())
    outpath.write_text("\n".join(svg))
    return outpath


def size_hist_svg(size_rows, outpath):
    rows = [r for r in size_rows if r["role"] in {"TP-base", "FP", "FN"} and int(r["size"]) > 0]
    if not rows:
        return None

    vals = [math.log10(int(r["size"])) for r in rows]
    min_v = min(vals)
    max_v = max(vals)

    if min_v == max_v:
        max_v += 1

    bins = 20
    labels = ["TP-base", "FP", "FN"]
    matrix = {role: [0] * bins for role in labels}

    for r in rows:
        role = r["role"]
        v = math.log10(int(r["size"]))
        idx = int((v - min_v) / (max_v - min_v) * bins)
        if idx == bins:
            idx -= 1
        matrix[role][idx] += 1

    x_labels = [f"{min_v + i*(max_v-min_v)/bins:.1f}" for i in range(bins)]
    categories = labels
    m = defaultdict(dict)
    for i, xl in enumerate(x_labels):
        for role in categories:
            m[xl][role] = matrix[role][i]

    return stacked_bar_svg(
        "SV size distribution",
        x_labels,
        categories,
        m,
        "Number of SVs",
        outpath,
    )


def make_plots(metric_rows, svtype_rows, chrom_rows, size_rows):
    plots = []

    labels = [r["label"] for r in metric_rows]

    p = grouped_bar_svg(
        "Truvari precision, recall, and F1",
        labels,
        [
            ("Precision", [r["precision"] or 0 for r in metric_rows]),
            ("Recall", [r["recall"] or 0 for r in metric_rows]),
            ("F1", [r["f1"] or 0 for r in metric_rows]),
        ],
        "Score",
        PLOTS_DIR / "01_precision_recall_f1.svg",
        ymax=1,
    )
    if p:
        plots.append(("Precision, recall, and F1", p))

    p = grouped_bar_svg(
        "Truvari TP, FP, and FN counts",
        labels,
        [
            ("TP-base", [r["tp_base"] or 0 for r in metric_rows]),
            ("FP", [r["fp"] or 0 for r in metric_rows]),
            ("FN", [r["fn"] or 0 for r in metric_rows]),
        ],
        "Number of SVs",
        PLOTS_DIR / "02_tp_fp_fn.svg",
    )
    if p:
        plots.append(("TP, FP, and FN counts", p))

    p = fp_fn_balance_svg(metric_rows, PLOTS_DIR / "03_fp_fn_balance.svg")
    if p:
        plots.append(("FP/FN balance", p))

    fpfn_svtype = [r for r in svtype_rows if r["role"] in {"FP", "FN"}]
    svtypes = sorted(set(r["svtype"] for r in fpfn_svtype))
    sv_labels = sorted(set(f"{r['build']} {r['caller']} {r['role']}" for r in fpfn_svtype))
    sv_matrix = defaultdict(dict)

    for r in fpfn_svtype:
        label = f"{r['build']} {r['caller']} {r['role']}"
        sv_matrix[label][r["svtype"]] = sv_matrix[label].get(r["svtype"], 0) + int(r["count"])

    p = stacked_bar_svg(
        "FP and FN by SVTYPE",
        sv_labels,
        svtypes,
        sv_matrix,
        "Number of SVs",
        PLOTS_DIR / "04_fp_fn_by_svtype.svg",
    )
    if p:
        plots.append(("FP and FN by SVTYPE", p))

    total_chrom = Counter()
    for r in chrom_rows:
        total_chrom[r["chrom"]] += int(r["count"])

    top_chroms = [chrom for chrom, _ in total_chrom.most_common(12)]
    chrom_filtered = [r for r in chrom_rows if r["chrom"] in top_chroms]

    chrom_labels = sorted(set(f"{r['build']} {r['caller']} {r['role']}" for r in chrom_filtered))
    chrom_matrix = defaultdict(dict)

    for r in chrom_filtered:
        label = f"{r['build']} {r['caller']} {r['role']}"
        chrom_matrix[label][r["chrom"]] = chrom_matrix[label].get(r["chrom"], 0) + int(r["count"])

    p = stacked_bar_svg(
        "FP and FN by chromosome",
        chrom_labels,
        top_chroms,
        chrom_matrix,
        "Number of SVs",
        PLOTS_DIR / "05_fp_fn_by_chromosome.svg",
    )
    if p:
        plots.append(("FP and FN by chromosome", p))

    p = size_hist_svg(size_rows, PLOTS_DIR / "06_sv_size_distribution.svg")
    if p:
        plots.append(("SV size distribution", p))

    return plots


def html_table(rows, max_rows=300):
    if not rows:
        return "<p>No data available.</p>"

    rows = rows[:max_rows]
    keys = list(rows[0].keys())

    out = ["<table>", "<thead><tr>"]
    for key in keys:
        out.append(f"<th>{esc(key)}</th>")
    out.append("</tr></thead><tbody>")

    for row in rows:
        out.append("<tr>")
        for key in keys:
            out.append(f"<td>{esc(nice_num(row.get(key, '')))}</td>")
        out.append("</tr>")

    out.append("</tbody></table>")
    return "\n".join(out)


def make_html(metric_rows, svtype_rows, chrom_rows, found, missing, plots):
    display_metrics = []

    for r in metric_rows:
        display_metrics.append({
            "Build": r["build"],
            "Caller": r["caller"],
            "Base count": r["base_count"],
            "Comparison count": r["comparison_count"],
            "TP-base": r["tp_base"],
            "TP-comp": r["tp_comp"],
            "FP": r["fp"],
            "FN": r["fn"],
            "Precision": r["precision"],
            "Recall": r["recall"],
            "F1": r["f1"],
        })

    display_svtype = [
        {
            "Build": r["build"],
            "Caller": r["caller"],
            "Role": r["role"],
            "SVTYPE": r["svtype"],
            "Count": r["count"],
        }
        for r in svtype_rows
        if r["role"] in {"FP", "FN"}
    ]

    display_chrom = [
        {
            "Build": r["build"],
            "Caller": r["caller"],
            "Role": r["role"],
            "Chromosome": r["chrom"],
            "Count": r["count"],
        }
        for r in chrom_rows
    ]

    parts = []

    parts.append("""
<!DOCTYPE html>
<html>
<head>
<meta charset="utf-8">
<title>Simple Truvari SV Benchmark Report</title>
<style>
body {
    font-family: Arial, Helvetica, sans-serif;
    margin: 40px;
    color: #222;
    line-height: 1.5;
}
h1, h2, h3 { color: #111; }
table {
    border-collapse: collapse;
    width: 100%;
    margin: 16px 0 28px 0;
    font-size: 14px;
}
th, td {
    border: 1px solid #ccc;
    padding: 6px 8px;
    text-align: left;
}
th { background: #f2f2f2; }
img {
    max-width: 100%;
    border: 1px solid #ddd;
    margin: 8px 0 28px 0;
}
.note {
    background: #f7f7f7;
    border-left: 5px solid #888;
    padding: 12px 16px;
    margin: 18px 0;
}
code {
    background: #f0f0f0;
    padding: 2px 4px;
}
</style>
</head>
<body>
""")

    parts.append("<h1>Simple Truvari Structural Variant Benchmark Report</h1>")

    parts.append("""
<div class="note">
<b>Scope:</b> This report uses only existing Truvari results from
<code>benchmark-runs/chr22/truvari/</code>. It does not rerun Truvari,
does not use RTGtools, does not use Wittyer, and does not include liftover.
</div>
""")

    parts.append("<h2>1. Aim</h2>")
    parts.append("""
<p>
The aim of this report is to summarize the existing Truvari benchmark results
for HG002 structural variant callsets on hg19 and hg38. The report focuses on
precision, recall, F1, true positives, false positives, false negatives,
SVTYPE patterns, and discordant variant categories.
</p>
""")


    parts.append("<h2>2. Main metrics</h2>")
    parts.append("""
<p>
Precision means how many comparison calls were supported by the truth set.
Recall means how many truth variants were recovered. F1 balances precision
and recall. False positives are extra calls in the comparison set. False
negatives are truth variants missed by the caller.
</p>
""")
    parts.append(html_table(display_metrics))

    parts.append("<h2>3. Plots</h2>")
    for title, path in plots:
        rel = path.relative_to(OUT_DIR)
        parts.append(f"<h3>{esc(title)}</h3>")
        parts.append(f'<img src="{esc(rel)}" alt="{esc(title)}">')

    parts.append("<h2>4. FP and FN by SVTYPE</h2>")
    parts.append(html_table(display_svtype))

    parts.append("<h2>5. FP and FN by chromosome</h2>")
    parts.append(html_table(display_chrom))



    parts.append("</body></html>")

    HTML_OUT.write_text("\n".join(parts))


def main():
    ensure_dirs()

    if not TRUVARI_DIR.exists():
        raise SystemExit(f"ERROR: Missing folder: {TRUVARI_DIR}")

    metric_rows, svtype_rows, chrom_rows, size_rows, found, missing = collect_results()

    if not metric_rows:
        raise SystemExit("ERROR: No Truvari result folders found.")

    write_csv(metric_rows, TABLES_DIR / "truvari_simple_metrics.csv")
    write_csv(svtype_rows, TABLES_DIR / "truvari_svtype_counts.csv")
    write_csv(chrom_rows, TABLES_DIR / "truvari_fp_fn_chromosome_counts.csv")
    write_csv(size_rows, TABLES_DIR / "truvari_size_values.csv")

    plots = make_plots(metric_rows, svtype_rows, chrom_rows, size_rows)

    make_html(metric_rows, svtype_rows, chrom_rows, found, missing, plots)

    print("")
    print("DONE")
    print(f"Report created: {HTML_OUT}")
    print(f"Plots created:  {PLOTS_DIR}")
    print(f"Tables created: {TABLES_DIR}")
    print("")


if __name__ == "__main__":
    main()
