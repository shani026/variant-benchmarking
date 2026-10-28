from pathlib import Path
import csv
import subprocess
import datetime
import html
import statistics

ROOT = Path(__file__).resolve().parents[2]

RTGTOOLS_CSV = ROOT / "benchmark-runs/nfcore/hg38/structural/summary/tables/rtgtools/rtgtools.summary.csv"
SVBENCHMARK_CSV = ROOT / "benchmark-runs/nfcore/hg38/structural/summary/tables/svbenchmark/svbenchmark.summary.csv"
TRUVARI_CSV = ROOT / "benchmark-runs/nfcore/hg38/structural/summary/tables/truvari/truvari.summary.csv"

VCFS = {
    "Truth": ROOT / "Truth_vcfs/hg38/HG002_GRCh38_v5.0q_stvar.vcf.gz",
    "Delly": ROOT / "caller_vcfs/HG002_delly_lr.hg38.renamed.vcf.gz",
    "DRAGEN": ROOT / "caller_vcfs/HG002_dragen.sv.hg38.renamed.vcf.gz",
    "Sniffles": ROOT / "caller_vcfs/HG002_sniffles.hg38.renamed.vcf.gz",
}

OUT_HTML = ROOT / "benchmark-runs/nfcore/hg38/hg38_final_sv_benchmarking_report.html"
OUT_VALUES = ROOT / "analysis/hg38/hg38_values_used.csv"
OUT_SVTYPES = ROOT / "analysis/hg38/hg38_svtype_counts_from_actual_vcfs.csv"

def read_metrics(path, method):
    rows = []
    with open(path, newline="") as f:
        reader = csv.DictReader(f)
        for row in reader:
            rows.append({
                "method": method,
                "tool": row.get("Tool", ""),
                "file": row.get("File", ""),
                "caller": row.get("Caller", ""),
                "threshold": row.get("Threshold", ""),
                "tp_base": int(float(row.get("TP_base", 0) or 0)),
                "tp_comp": int(float(row.get("TP_comp", 0) or 0)),
                "fp": int(float(row.get("FP", 0) or 0)),
                "fn": int(float(row.get("FN", 0) or 0)),
                "precision": float(row.get("Precision", 0) or 0),
                "recall": float(row.get("Recall", 0) or 0),
                "f1": float(row.get("F1", 0) or 0),
            })
    return rows

all_metrics = []
all_metrics += read_metrics(RTGTOOLS_CSV, "RTGtools")
all_metrics += read_metrics(SVBENCHMARK_CSV, "SVbenchmark")
all_metrics += read_metrics(TRUVARI_CSV, "Truvari")

with open(OUT_VALUES, "w", newline="") as f:
    writer = csv.DictWriter(f, fieldnames=list(all_metrics[0].keys()))
    writer.writeheader()
    writer.writerows(all_metrics)

def run(cmd):
    p = subprocess.run(cmd, text=True, capture_output=True)
    if p.returncode != 0:
        return None, p.stderr.strip()
    return p.stdout, None

def svtype_counts(vcf):
    stdout, err = run(["bcftools", "query", "-f", "%INFO/SVTYPE\n", str(vcf)])
    if err:
        return {"ERROR": err}
    counts = {}
    for line in stdout.splitlines():
        svt = line.strip() or "UNKNOWN"
        if svt == ".":
            svt = "UNKNOWN"
        counts[svt] = counts.get(svt, 0) + 1
    return counts

def record_count(vcf):
    stdout, err = run(["bcftools", "view", "-H", str(vcf)])
    if err:
        return "ERROR"
    return len(stdout.splitlines())

def pass_count(vcf):
    stdout, err = run(["bcftools", "view", "-H", "-f", "PASS", str(vcf)])
    if err:
        return "ERROR"
    return len(stdout.splitlines())

svtype_rows = []
for label, vcf in VCFS.items():
    counts = svtype_counts(vcf)
    total = record_count(vcf)
    passed = pass_count(vcf)
    for svtype, count in counts.items():
        svtype_rows.append({
            "source": label,
            "vcf": str(vcf),
            "total_records": total,
            "pass_records": passed,
            "svtype": svtype,
            "count": count,
        })

with open(OUT_SVTYPES, "w", newline="") as f:
    writer = csv.DictWriter(f, fieldnames=["source", "vcf", "total_records", "pass_records", "svtype", "count"])
    writer.writeheader()
    writer.writerows(svtype_rows)

def fmt(x, digits=4):
    try:
        return f"{float(x):.{digits}f}"
    except Exception:
        return str(x)

def fmt_pct(x):
    try:
        return f"{float(x) * 100:.2f}%"
    except Exception:
        return str(x)

def fmt_int(x):
    try:
        return f"{int(x):,}"
    except Exception:
        return str(x)

def caller_title(c):
    return {"delly": "Delly", "dragen": "DRAGEN", "sniffles": "Sniffles"}.get(c, c)

def caller_metrics(caller):
    rows = [r for r in all_metrics if r["caller"] == caller]
    return {
        "precision_mean": statistics.mean([r["precision"] for r in rows]),
        "recall_mean": statistics.mean([r["recall"] for r in rows]),
        "f1_mean": statistics.mean([r["f1"] for r in rows]),
        "tp_base_sum": sum(r["tp_base"] for r in rows),
        "tp_comp_sum": sum(r["tp_comp"] for r in rows),
        "fp_sum": sum(r["fp"] for r in rows),
        "fn_sum": sum(r["fn"] for r in rows),
    }

callers = ["delly", "dragen", "sniffles"]
means = {c: caller_metrics(c) for c in callers}

def metric_table(method):
    rows = [r for r in all_metrics if r["method"] == method]
    body = ""
    for r in rows:
        body += f"""
        <tr>
          <td class="mono">{caller_title(r["caller"])}</td>
          <td class="num">{fmt_int(r["tp_base"])}</td>
          <td class="num">{fmt_int(r["tp_comp"])}</td>
          <td class="num warn">{fmt_int(r["fp"])}</td>
          <td class="num bad">{fmt_int(r["fn"])}</td>
          <td class="num good">{fmt(r["precision"], 6)}</td>
          <td class="num bad">{fmt(r["recall"], 6)}</td>
          <td class="num bad">{fmt(r["f1"], 6)}</td>
        </tr>
        """
    return body

def bar_chart(metric, title, max_value=None, percent=True):
    rows = []
    values = [r[metric] for r in all_metrics]
    if max_value is None:
        max_value = max(values) if values else 1
    if max_value == 0:
        max_value = 1

    for r in all_metrics:
        value = r[metric]
        width = min(100, (value / max_value) * 100)
        label = fmt_pct(value) if percent else fmt_int(value)
        rows.append(f"""
        <div class="bar-row">
          <div class="bar-label">{caller_title(r["caller"])} · {r["method"]}</div>
          <div class="bar-bg"><div class="bar-fill" style="width:{width:.2f}%"></div></div>
          <div class="bar-value">{label}</div>
        </div>
        """)
    return f"""
    <div class="plot-card">
      <h3>{title}</h3>
      {''.join(rows)}
    </div>
    """

def svtype_table():
    grouped = {}
    for r in svtype_rows:
        source = r["source"]
        grouped.setdefault(source, {"total": r["total_records"], "pass": r["pass_records"], "types": {}})
        grouped[source]["types"][r["svtype"]] = r["count"]

    all_types = sorted({r["svtype"] for r in svtype_rows})
    header = "".join(f"<th>{html.escape(t)}</th>" for t in all_types)
    body = ""

    for source in ["Truth", "Delly", "DRAGEN", "Sniffles"]:
        item = grouped.get(source, {"total": "NA", "pass": "NA", "types": {}})
        cells = "".join(f"<td class='num'>{fmt_int(item['types'].get(t, 0))}</td>" for t in all_types)
        body += f"""
        <tr>
          <td class="mono">{source}</td>
          <td class="num">{fmt_int(item["total"])}</td>
          <td class="num">{fmt_int(item["pass"])}</td>
          {cells}
        </tr>
        """

    return f"""
    <table>
      <thead>
        <tr>
          <th>Source</th>
          <th>Total records</th>
          <th>PASS records</th>
          {header}
        </tr>
      </thead>
      <tbody>{body}</tbody>
    </table>
    """

def svtype_bars():
    rows = []
    max_count = max([int(r["count"]) for r in svtype_rows if str(r["count"]).isdigit()] or [1])
    for r in svtype_rows:
        try:
            count = int(r["count"])
        except Exception:
            continue
        width = min(100, (count / max_count) * 100)
        rows.append(f"""
        <div class="bar-row">
          <div class="bar-label">{html.escape(r["source"])} · {html.escape(r["svtype"])}</div>
          <div class="bar-bg"><div class="bar-fill blue" style="width:{width:.2f}%"></div></div>
          <div class="bar-value">{fmt_int(count)}</div>
        </div>
        """)
    return f"""
    <div class="plot-card">
      <h3>SVTYPE counts from actual VCF files</h3>
      {''.join(rows)}
    </div>
    """

cards = ""
for c in callers:
    m = means[c]
    cards += f"""
    <div style="margin-bottom:24px;">
      <div class="caller-title">{caller_title(c)}</div>
      <div class="cards">
        <div class="card"><div class="val good">{fmt(m["precision_mean"], 4)}</div><div class="lbl">Mean precision</div></div>
        <div class="card"><div class="val bad">{fmt(m["recall_mean"], 4)}</div><div class="lbl">Mean recall</div></div>
        <div class="card"><div class="val bad">{fmt(m["f1_mean"], 4)}</div><div class="lbl">Mean F1</div></div>
        <div class="card"><div class="val">{fmt_int(m["tp_base_sum"])}</div><div class="lbl">TP-base total</div></div>
        <div class="card"><div class="val warn">{fmt_int(m["fp_sum"])}</div><div class="lbl">FP total</div></div>
        <div class="card"><div class="val bad">{fmt_int(m["fn_sum"])}</div><div class="lbl">FN total</div></div>
      </div>
    </div>
    """

html_text = f"""<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="UTF-8">
<title>SV Benchmarking Report — hg38</title>
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
  font-family: Arial, sans-serif;
  font-size:13px;
  line-height:1.6;
  padding-bottom:60px;
}}
.header {{
  background:linear-gradient(135deg,#0d1117 0%,#161b22 60%,#1c2a3a 100%);
  border-bottom:1px solid var(--border);
  padding:36px 48px 28px;
}}
.header h1 {{
  margin:0;
  font-size:28px;
  color:var(--accent);
}}
.meta {{
  color:var(--muted);
  margin-top:6px;
  font-size:12px;
}}
.container {{
  max-width:1300px;
  margin:0 auto;
  padding:32px 48px;
}}
.badge {{
  display:inline-block;
  background:rgba(88,166,255,.12);
  color:var(--accent);
  border:1px solid rgba(88,166,255,.3);
  border-radius:4px;
  padding:2px 8px;
  font-size:10px;
  font-weight:700;
  text-transform:uppercase;
  margin-right:6px;
}}
.badge.green {{
  color:var(--green);
  border-color:rgba(63,185,80,.3);
  background:rgba(63,185,80,.1);
}}
.section {{
  margin-bottom:42px;
}}
.section-label {{
  font-size:12px;
  font-weight:700;
  letter-spacing:2px;
  text-transform:uppercase;
  color:var(--muted);
  margin-bottom:14px;
  padding-bottom:6px;
  border-bottom:1px solid var(--border);
}}
.note {{
  background:rgba(88,166,255,.06);
  border-left:3px solid var(--accent);
  padding:12px 16px;
  color:var(--subtext);
  margin-bottom:18px;
}}
.warning {{
  background:rgba(248,81,73,.08);
  border-left:3px solid var(--danger);
  padding:12px 16px;
  color:var(--subtext);
  margin-bottom:18px;
}}
.config-block {{
  background:var(--surface);
  border:1px solid var(--border);
  border-radius:8px;
  padding:16px 20px;
}}
.config-block .row {{
  display:flex;
  gap:12px;
  margin-bottom:5px;
}}
.config-block .key {{
  color:var(--muted);
  min-width:180px;
  font-size:11px;
}}
.config-block .value {{
  color:var(--subtext);
  word-break:break-all;
}}
.cards {{
  display:grid;
  grid-template-columns:repeat(auto-fill,minmax(180px,1fr));
  gap:12px;
}}
.card {{
  background:var(--surface);
  border:1px solid var(--border);
  border-radius:8px;
  padding:16px;
}}
.val {{
  font-size:27px;
  font-weight:800;
  color:var(--accent);
}}
.val.good {{ color:var(--green); }}
.val.warn {{ color:var(--warn); }}
.val.bad {{ color:var(--danger); }}
.lbl {{
  color:var(--muted);
  font-size:10px;
  text-transform:uppercase;
  letter-spacing:1px;
}}
.caller-title {{
  font-size:16px;
  font-weight:800;
  color:var(--subtext);
  margin-bottom:10px;
}}
table {{
  width:100%;
  border-collapse:collapse;
  background:var(--surface);
  border:1px solid var(--border);
  margin-bottom:16px;
}}
th {{
  background:var(--highlight);
  color:var(--accent);
  text-align:left;
  padding:8px 10px;
  font-size:10px;
  text-transform:uppercase;
  letter-spacing:1px;
}}
td {{
  padding:7px 10px;
  border-top:1px solid var(--border);
  color:var(--subtext);
}}
td.num {{
  text-align:right;
  font-variant-numeric:tabular-nums;
}}
td.good {{ color:var(--green); }}
td.warn {{ color:var(--warn); }}
td.bad {{ color:var(--danger); }}
td.mono {{
  font-family:monospace;
}}
.plot-card {{
  background:var(--surface);
  border:1px solid var(--border);
  border-radius:8px;
  padding:18px;
  margin-bottom:20px;
}}
.plot-card h3 {{
  margin-top:0;
  color:var(--subtext);
}}
.bar-row {{
  display:grid;
  grid-template-columns:230px 1fr 90px;
  gap:10px;
  align-items:center;
  margin-bottom:8px;
}}
.bar-label {{
  color:var(--subtext);
  font-size:12px;
}}
.bar-bg {{
  height:12px;
  background:var(--highlight);
  border-radius:10px;
  overflow:hidden;
}}
.bar-fill {{
  height:100%;
  background:var(--danger);
  border-radius:10px;
}}
.bar-fill.blue {{
  background:var(--accent);
}}
.bar-value {{
  text-align:right;
  font-variant-numeric:tabular-nums;
  color:var(--subtext);
  font-size:12px;
}}
.source {{
  color:var(--muted);
  font-size:11px;
  margin-top:8px;
}}
hr {{
  border:0;
  border-top:1px solid var(--border);
  margin:32px 0;
}}
</style>
</head>
<body>
<div class="header">
  <span class="badge">nf-core/variantbenchmarking</span>
  <span class="badge green">v1.5.0</span>
  <span class="badge">hg38</span>
  <h1>SV Benchmarking Report</h1>
  <div class="meta">
    Truth set: HG002 | Reference: hg38 / GRCh38 | Callers: Delly · DRAGEN · Sniffles | Benchmarks: RTGtools · SVbenchmark · Truvari
  </div>
</div>

<div class="container">

<div class="section">
  <div class="section-label">Run summary</div>
  <div class="warning">
    The report was generated by reading the actual local files from <code>benchmark-runs/nfcore/hg38</code>, <code>Truth_vcfs/hg38</code>, and <code>caller_vcfs</code>. Recall is extremely low across all tools; the main driver is the very large FN count.
  </div>
  <div class="config-block">
    <div class="row"><span class="key">Generated</span><span class="value">{datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")}</span></div>
    <div class="row"><span class="key">RTGtools table</span><span class="value">{RTGTOOLS_CSV}</span></div>
    <div class="row"><span class="key">SVbenchmark table</span><span class="value">{SVBENCHMARK_CSV}</span></div>
    <div class="row"><span class="key">Truvari table</span><span class="value">{TRUVARI_CSV}</span></div>
    <div class="row"><span class="key">Truth VCF</span><span class="value">{VCFS["Truth"]}</span></div>
    <div class="row"><span class="key">Output values CSV</span><span class="value">{OUT_VALUES}</span></div>
    <div class="row"><span class="key">SVTYPE counts CSV</span><span class="value">{OUT_SVTYPES}</span></div>
  </div>
</div>

<div class="section">
  <div class="section-label">Per-caller mean summary across benchmarking tools</div>
  {cards}
</div>

<hr>

<div class="section">
  <div class="section-label">Recall-focused plots</div>
  {bar_chart("recall", "Recall by caller and benchmarking tool")}
  {bar_chart("f1", "F1 by caller and benchmarking tool")}
  {bar_chart("precision", "Precision by caller and benchmarking tool", max_value=1)}
  {bar_chart("fn", "False negatives by caller and benchmarking tool", percent=False)}
</div>

<hr>

<div class="section">
  <div class="section-label">Truvari results</div>
  <table>
    <thead>
      <tr>
        <th>Caller</th><th>TP-base</th><th>TP-comp</th><th>FP</th><th>FN</th><th>Precision</th><th>Recall</th><th>F1</th>
      </tr>
    </thead>
    <tbody>{metric_table("Truvari")}</tbody>
  </table>
  <div class="source">Source: {TRUVARI_CSV}</div>
</div>

<div class="section">
  <div class="section-label">SVbenchmark results</div>
  <table>
    <thead>
      <tr>
        <th>Caller</th><th>TP-base</th><th>TP-comp</th><th>FP</th><th>FN</th><th>Precision</th><th>Recall</th><th>F1</th>
      </tr>
    </thead>
    <tbody>{metric_table("SVbenchmark")}</tbody>
  </table>
  <div class="source">Source: {SVBENCHMARK_CSV}</div>
</div>

<div class="section">
  <div class="section-label">RTGtools results</div>
  <table>
    <thead>
      <tr>
        <th>Caller</th><th>TP-base</th><th>TP-comp</th><th>FP</th><th>FN</th><th>Precision</th><th>Recall</th><th>F1</th>
      </tr>
    </thead>
    <tbody>{metric_table("RTGtools")}</tbody>
  </table>
  <div class="source">Source: {RTGTOOLS_CSV}</div>
</div>

<hr>

<div class="section">
  <div class="section-label">SVTYPE counts from actual VCF files</div>
  <div class="note">
    These counts are produced directly from the actual VCF files with:
    <code>bcftools query -f '%INFO/SVTYPE\\n'</code>
  </div>
  {svtype_bars()}
  {svtype_table()}
</div>

<hr>

<div class="section">
  <div class="section-label">Interpretation</div>
  <div class="note">
    The benchmarking outputs show high precision for SVbenchmark and Truvari, but extremely low recall. That means the callers mostly avoid many false positives, but they recover only a tiny fraction of the truth set under the current benchmark settings. RTGtools returned zero TP/FP/FN for all three callers, so RTGtools should be treated separately as a compatibility or representation issue rather than interpreted biologically.
  </div>
</div>

</div>
</body>
</html>
"""

OUT_HTML.write_text(html_text)

print("DONE")
print("HTML report:", OUT_HTML)
print("Values CSV:", OUT_VALUES)
print("SVTYPE counts CSV:", OUT_SVTYPES)
