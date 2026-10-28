#!/usr/bin/env python3
"""
Build a complete HTML report for the chr22 SV benchmarking project.

Run this from the project root:
    python3 make_final_html_report.py

It reads the existing tables/results under reports/chr22/, benchmark-runs/chr22/,
and writes:
    reports/chr22/final_chr22_sv_benchmark_report.html

No heredoc, no manual table copying.
"""

from __future__ import annotations

import csv
import html
import json
import re
from pathlib import Path
from datetime import datetime
from typing import Iterable, Optional


ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "reports/chr22" / "final_chr22_sv_benchmark_report.html"


def rel(path: Path | str) -> str:
    p = Path(path)
    try:
        return str(p.relative_to(ROOT))
    except ValueError:
        return str(p)


def exists(path: str | Path) -> bool:
    return (ROOT / path).exists() if not Path(path).is_absolute() else Path(path).exists()


def read_text(path: str | Path, fallback: str = "Not available.") -> str:
    p = ROOT / path if not Path(path).is_absolute() else Path(path)
    if p.exists():
        return p.read_text(errors="replace")
    return fallback


def esc(x) -> str:
    return html.escape("" if x is None else str(x))


def safe_float(x) -> Optional[float]:
    try:
        if x is None or str(x).strip() == "":
            return None
        return float(x)
    except Exception:
        return None


def fmt_float(x, digits: int = 3) -> str:
    f = safe_float(x)
    return "" if f is None else f"{f:.{digits}f}"


def read_tsv(path: str | Path) -> list[dict[str, str]]:
    p = ROOT / path if not Path(path).is_absolute() else Path(path)
    if not p.exists():
        return []
    with p.open(newline="", errors="replace") as f:
        return list(csv.DictReader(f, delimiter="\t"))


def read_table_any(path: str | Path) -> tuple[list[str], list[list[str]]]:
    """
    Reads TSV/CSV-ish files. Falls back to whitespace table if needed.
    """
    p = ROOT / path if not Path(path).is_absolute() else Path(path)
    if not p.exists():
        return [], []
    text = p.read_text(errors="replace").strip()
    if not text:
        return [], []

    lines = text.splitlines()

    # TSV
    if "\t" in lines[0]:
        reader = csv.reader(lines, delimiter="\t")
        rows = list(reader)
        return rows[0], rows[1:]

    # CSV
    if "," in lines[0] and not lines[0].startswith("##"):
        reader = csv.reader(lines)
        rows = list(reader)
        return rows[0], rows[1:]

    # Fallback: whitespace split if header-like first line
    rows = [re.split(r"\s{2,}|\t", ln.strip()) for ln in lines if ln.strip()]
    if len(rows) >= 2 and len(rows[0]) > 1:
        return rows[0], rows[1:]

    return [], []


def html_table(headers: list[str], rows: Iterable[Iterable[str]], caption: str = "") -> str:
    rows = list(rows)
    if not headers:
        return '<p class="missing">Table not available.</p>'
    cap = f"<caption>{esc(caption)}</caption>" if caption else ""
    thead = "".join(f"<th>{esc(h)}</th>" for h in headers)
    body = []
    for row in rows:
        body.append("<tr>" + "".join(f"<td>{esc(cell)}</td>" for cell in row) + "</tr>")
    return f"""
<table>
{cap}
<thead><tr>{thead}</tr></thead>
<tbody>
{''.join(body)}
</tbody>
</table>
"""


def table_from_file(path: str, caption: str = "") -> str:
    headers, rows = read_table_any(path)
    if headers:
        return html_table(headers, rows, caption=caption)
    txt = read_text(path, "")
    if txt.strip():
        return f"<pre>{esc(txt)}</pre>"
    return f'<p class="missing">Not available: <code>{esc(path)}</code></p>'


def pre_from_file(path: str, max_chars: Optional[int] = None) -> str:
    txt = read_text(path, f"Not available: {path}")
    if max_chars is not None and len(txt) > max_chars:
        txt = txt[:max_chars] + "\n\n[truncated in report preview]"
    return f"<pre>{esc(txt)}</pre>"


def file_list(pattern: str) -> list[str]:
    return sorted(rel(p) for p in ROOT.glob(pattern) if p.is_file())


def list_as_html(items: Iterable[str]) -> str:
    items = list(items)
    if not items:
        return '<p class="missing">No files found.</p>'
    return "<ul>" + "".join(f"<li><code>{esc(x)}</code></li>" for x in items) + "</ul>"


def metric_cards_from_truvari(path: str) -> str:
    rows = read_tsv(path)
    if not rows:
        return '<p class="missing">Truvari summary table not available.</p>'

    rows_sorted = sorted(rows, key=lambda r: safe_float(r.get("f1")) or -1, reverse=True)
    best = rows_sorted[0]
    worst = rows_sorted[-1]

    # Best by build when run names start with hg19/hg38
    by_build = {}
    for row in rows:
        run = row.get("run", "")
        build = "hg19" if run.startswith("hg19") else "hg38" if run.startswith("hg38") else "other"
        if build not in by_build or (safe_float(row.get("f1")) or -1) > (safe_float(by_build[build].get("f1")) or -1):
            by_build[build] = row

    cards = []
    cards.append(card("Best overall F1", f"{best.get('run','')}", f"F1 {fmt_float(best.get('f1'))}, precision {fmt_float(best.get('precision'))}, recall {fmt_float(best.get('recall'))}"))
    cards.append(card("Lowest F1", f"{worst.get('run','')}", f"F1 {fmt_float(worst.get('f1'))}, precision {fmt_float(worst.get('precision'))}, recall {fmt_float(worst.get('recall'))}"))
    for build in ["hg19", "hg38"]:
        if build in by_build:
            row = by_build[build]
            cards.append(card(f"Best {build}", row.get("run", ""), f"F1 {fmt_float(row.get('f1'))}, recall {fmt_float(row.get('recall'))}"))
    return f'<div class="cards">{"".join(cards)}</div>'


def card(title: str, value: str, sub: str = "") -> str:
    return f"""
<div class="card">
  <div class="card-title">{esc(title)}</div>
  <div class="card-value">{esc(value)}</div>
  <div class="card-sub">{esc(sub)}</div>
</div>
"""


def original_vs_dup_summary(path: str) -> str:
    rows = read_tsv(path)
    if not rows:
        return '<p class="missing">Original vs DUP-to-INS comparison table not available.</p>'
    changed = []
    unchanged = []
    for r in rows:
        d = safe_float(r.get("delta_f1"))
        if d is None:
            continue
        if abs(d) > 1e-12:
            changed.append((r.get("run", ""), d))
        else:
            unchanged.append(r.get("run", ""))

    if changed:
        detail = "; ".join(f"{run}: ΔF1={d:+.4f}" for run, d in changed)
        return f"<p>The DUP-to-INS conversion changed at least one F1 score: {esc(detail)}.</p>"
    return "<p>The DUP-to-INS conversion did not change F1 for the runs represented in the comparison table, or the changes were effectively zero.</p>"


def rtg_overview() -> str:
    base = file_list("benchmark-runs/chr22/rtg/vcfeval/*/summary.txt")
    dup = file_list("benchmark-runs/chr22/rtg/vcfeval_dup_to_ins/*/summary.txt")
    return f"""
<div class="cards">
{card("RTG original summaries", str(len(base)), "summary.txt files")}
{card("RTG DUP-to-INS summaries", str(len(dup)), "summary.txt files")}
</div>
"""


def witty_overview() -> str:
    base = file_list("benchmark-runs/chr22/wittyer/**/*")
    dup = file_list("benchmark-runs/chr22/wittyer_dup_to_ins/**/*")
    return f"""
<div class="cards">
{card("Witty.er original outputs", str(len(base)), "files")}
{card("Witty.er DUP-to-INS outputs", str(len(dup)), "files")}
</div>
"""


def source_file_audit() -> str:
    expected = [
        "reports/chr22/qc_tables/raw_chr22_record_summary.tsv",
        "reports/chr22/qc_tables/raw_chr22_vcf_counts.txt",
        "reports/chr22/qc_tables/benchmark_ready_interval_counts.txt",
        "reports/chr22/qc_tables/BND_counts.tsv",
        "reports/chr22/qc_tables/END_lt_POS_counts.tsv",
        "reports/chr22/qc_tables/missing_END_SVLEN_counts.tsv",
        "reports/chr22/qc_tables/benchmark_ready_dup_to_ins_counts.txt",
        "reports/chr22/metrics/truvari_summary.tsv",
        "reports/chr22/metrics/truvari_summary_extended.tsv",
        "reports/chr22/metrics/truvari_svtype_tp_fp_fn_counts.tsv",
        "reports/chr22/metrics/truvari_converted_dup_to_ins_summary.tsv",
        "reports/chr22/metrics/truvari_original_vs_converted_dup_to_ins.tsv",
        "reports/chr22/metrics/rtg_vcfeval_summaries.txt",
        "reports/chr22/metrics/rtg_vcfeval_dup_to_ins_summaries.txt",
        "reports/chr22/metrics/wittyer_metric_like_files.txt",
        "reports/chr22/metrics/wittyer_preview.txt",
        "reports/chr22/metrics/wittyer_dup_to_ins_metric_like_files.txt",
        "reports/chr22/logs/reproducibility_log.txt",
    ]
    rows = []
    for f in expected:
        p = ROOT / f
        rows.append([f, "present" if p.exists() else "missing", str(p.stat().st_size) if p.exists() else ""])
    return html_table(["file", "status", "bytes"], rows)


CSS = """
:root {
  --bg: #f7f7fb;
  --card: #ffffff;
  --text: #222233;
  --muted: #667085;
  --line: #d9dce7;
  --accent: #355caa;
  --accent2: #ecf1ff;
  --warn: #fff4e5;
  --warn-border: #f0b35a;
}
* { box-sizing: border-box; }
body {
  margin: 0;
  font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", Roboto, Helvetica, Arial, sans-serif;
  background: var(--bg);
  color: var(--text);
  line-height: 1.55;
}
header {
  background: linear-gradient(135deg, #233866, #3f65c8);
  color: white;
  padding: 42px 48px;
}
header h1 { margin: 0 0 8px 0; font-size: 34px; }
header p { margin: 6px 0 0 0; opacity: 0.92; max-width: 1100px; }
main { max-width: 1250px; margin: 0 auto; padding: 28px; }
section {
  background: var(--card);
  border: 1px solid var(--line);
  border-radius: 14px;
  padding: 24px;
  margin: 20px 0;
  box-shadow: 0 2px 10px rgba(16, 24, 40, 0.04);
}
h2 { margin: 0 0 14px 0; color: #1d2d54; }
h3 { margin-top: 26px; color: #263b70; }
code {
  background: #eef2ff;
  padding: 2px 5px;
  border-radius: 4px;
}
pre {
  background: #101828;
  color: #f2f4f7;
  padding: 16px;
  border-radius: 10px;
  overflow-x: auto;
  font-size: 13px;
  line-height: 1.4;
}
table {
  border-collapse: collapse;
  width: 100%;
  margin: 12px 0 22px 0;
  font-size: 13px;
}
caption {
  text-align: left;
  font-weight: 700;
  margin-bottom: 8px;
}
th, td {
  border: 1px solid var(--line);
  padding: 8px 10px;
  text-align: left;
  vertical-align: top;
}
th {
  background: var(--accent2);
  color: #1d2d54;
  position: sticky;
  top: 0;
}
tr:nth-child(even) td { background: #fafbff; }
.cards {
  display: grid;
  grid-template-columns: repeat(auto-fit, minmax(210px, 1fr));
  gap: 14px;
  margin: 16px 0;
}
.card {
  background: #ffffff;
  border: 1px solid var(--line);
  border-left: 5px solid var(--accent);
  padding: 14px;
  border-radius: 10px;
}
.card-title { color: var(--muted); font-size: 13px; font-weight: 700; text-transform: uppercase; letter-spacing: .02em; }
.card-value { font-size: 20px; font-weight: 800; margin: 5px 0; }
.card-sub { color: var(--muted); font-size: 13px; }
.note {
  background: var(--warn);
  border-left: 5px solid var(--warn-border);
  padding: 12px 14px;
  border-radius: 8px;
}
.missing {
  color: #9a3412;
  background: #fff7ed;
  border: 1px solid #fed7aa;
  padding: 10px;
  border-radius: 8px;
}
.toc {
  columns: 2;
  line-height: 1.8;
}
.toc a { color: var(--accent); text-decoration: none; }
footer {
  color: var(--muted);
  text-align: center;
  padding: 30px;
}
@media print {
  body { background: white; }
  header { background: #fff; color: #000; border-bottom: 2px solid #000; }
  section { box-shadow: none; break-inside: avoid; }
  pre { white-space: pre-wrap; color: #000; background: #f2f2f2; }
}
"""


def section(id_: str, title: str, body: str) -> str:
    return f'<section id="{esc(id_)}"><h2>{esc(title)}</h2>\n{body}\n</section>'


def main() -> None:
    OUT.parent.mkdir(parents=True, exist_ok=True)

    truvari_summary_path = "reports/chr22/metrics/truvari_summary_extended.tsv"
    if not exists(truvari_summary_path):
        truvari_summary_path = "reports/chr22/metrics/truvari_summary.tsv"

    generated = datetime.now().strftime("%Y-%m-%d %H:%M:%S")

    sections = []

    sections.append(section("toc", "Contents", """
<div class="toc">
<a href="#aim">1. Aim</a><br>
<a href="#diagnosis">2. Diagnosis and preprocessing</a><br>
<a href="#scope">3. Benchmark scope</a><br>
<a href="#qc">4. QC tables</a><br>
<a href="#truvari">5. Truvari primary benchmark</a><br>
<a href="#rtg">6. RTGtools benchmark</a><br>
<a href="#witty">7. Witty.er benchmark</a><br>
<a href="#duptoins">8. DUP-to-INS sensitivity analysis</a><br>
<a href="#interpretation">9. Interpretation</a><br>
<a href="#limitations">10. Limitations</a><br>
<a href="#repro">11. Reproducibility and source file audit</a><br>
</div>
"""))

    sections.append(section("aim", "1. Aim", """
<p>The aim of this analysis was to benchmark HG002 structural variant calls from DRAGEN, Sniffles, and DELLY against the GIAB truth set on chromosome 22 for both hg19/GRCh37 and hg38/GRCh38.</p>
<p>This chromosome 22 benchmark was used as a debugging and validation subset before whole-genome benchmarking.</p>
"""))

    sections.append(section("diagnosis", "2. Diagnosis and preprocessing", """
<h3>Chromosome naming diagnosis</h3>
<p>The first major bug was chromosome naming mismatch. The hg19/GRCh37 files used <code>22</code>, while the hg38/GRCh38 files used <code>chr22</code>. The original zero-record subsets were therefore caused by using the wrong chromosome label, not by missing chromosome 22 data.</p>
<table>
<thead><tr><th>Build</th><th>Correct chromosome label</th></tr></thead>
<tbody>
<tr><td>hg19 / GRCh37</td><td><code>22</code></td></tr>
<tr><td>hg38 / GRCh38</td><td><code>chr22</code></td></tr>
</tbody>
</table>

<h3>Preprocessing workflow</h3>
<ol>
<li>Subset truth and caller VCFs to chromosome 22.</li>
<li>Subset truth BED files to chromosome 22.</li>
<li>Confirm matching chromosome names between VCF, BED, and reference.</li>
<li>Generate raw QC summaries before filtering.</li>
<li>Create benchmark-ready interval-SV VCFs under <code>03_benchmark_inputs/</code>.</li>
<li>Keep raw chr22 VCFs unchanged under <code>02_chr22/</code>.</li>
</ol>
"""))

    sections.append(section("scope", "3. Benchmark scope", """
<p>The main benchmark was restricted to interval-style structural variants:</p>
<ul>
<li>DEL</li>
<li>DUP</li>
<li>INS</li>
<li>INV</li>
</ul>
<p>The following were excluded from the main interval-SV benchmark but counted separately:</p>
<ul>
<li>BND</li>
<li>UNKNOWN</li>
<li>NA</li>
<li>missing SVTYPE</li>
<li>malformed interval records</li>
</ul>
<p>BND records were not silently removed. They were reported separately because breakend/translocation-style records are not directly comparable to interval-style DEL/DUP/INS/INV records in a chromosome 22-only benchmark.</p>
"""))

    qc_body = ""
    qc_body += "<h3>Raw chr22 record summary</h3>" + table_from_file("reports/chr22/qc_tables/raw_chr22_record_summary.tsv")
    qc_body += "<h3>Raw chr22 VCF counts</h3>" + pre_from_file("reports/chr22/qc_tables/raw_chr22_vcf_counts.txt")
    qc_body += "<h3>Benchmark-ready interval counts</h3>" + pre_from_file("reports/chr22/qc_tables/benchmark_ready_interval_counts.txt")
    qc_body += "<h3>BND counts</h3>" + table_from_file("reports/chr22/qc_tables/BND_counts.tsv")
    qc_body += "<h3>END &lt; POS counts</h3>" + table_from_file("reports/chr22/qc_tables/END_lt_POS_counts.tsv")
    qc_body += "<h3>Missing END/SVLEN counts</h3>" + table_from_file("reports/chr22/qc_tables/missing_END_SVLEN_counts.tsv")
    qc_body += "<p>The hg19 DELLY file had records where <code>INFO/END</code> was smaller than <code>POS</code>. Inspection showed these were BND records, not ordinary interval SV records. They were therefore counted separately and excluded from the main interval-SV benchmark.</p>"
    sections.append(section("qc", "4. QC tables", qc_body))

    truvari_body = metric_cards_from_truvari(truvari_summary_path)
    truvari_body += "<h3>Main Truvari summary</h3>" + table_from_file(truvari_summary_path)
    truvari_body += """
<p>Truvari was used as the primary SV benchmarking tool. The main run used the truth VCF as baseline, caller VCF as comparison, matching reference FASTA, chromosome 22 benchmark BED through <code>--includebed</code>, and PASS-only records through <code>--passonly</code>.</p>
<p>The benchmark-ready VCFs contained more records than Truvari finally evaluated because Truvari further restricts variants using BED/high-confidence regions, PASS filtering, size thresholds, and internal matching criteria.</p>
"""
    truvari_body += "<h3>Per-SVTYPE Truvari TP/FP/FN counts</h3>" + table_from_file("reports/chr22/metrics/truvari_svtype_tp_fp_fn_counts.tsv")
    sections.append(section("truvari", "5. Truvari primary benchmark", truvari_body))

    rtg_body = rtg_overview()
    rtg_body += "<h3>RTGtools original summaries</h3>" + pre_from_file("reports/chr22/metrics/rtg_vcfeval_summaries.txt", max_chars=60000)
    rtg_body += "<h3>RTGtools original output file list</h3>" + pre_from_file("reports/chr22/metrics/rtg_vcfeval_output_files.txt", max_chars=30000)
    sections.append(section("rtg", "6. RTGtools benchmark", rtg_body))

    witty_body = witty_overview()
    witty_body += "<h3>Witty.er metric-like output files</h3>" + pre_from_file("reports/chr22/metrics/wittyer_metric_like_files.txt")
    witty_body += "<h3>Witty.er output preview</h3>" + pre_from_file("reports/chr22/metrics/wittyer_preview.txt", max_chars=60000)
    sections.append(section("witty", "7. Witty.er benchmark", witty_body))

    dup_body = """
<p>A DUP-to-INS sensitivity analysis was included because tandem duplications can sometimes be represented either as <code>SVTYPE=DUP</code> or as inserted duplicated sequence with <code>SVTYPE=INS</code>.</p>
<p>The original benchmark preserved the caller representation. A separate sensitivity input set was created where <code>SVTYPE=DUP</code> was rewritten to <code>SVTYPE=INS</code>, and symbolic <code>&lt;DUP&gt;</code> alleles were rewritten to <code>&lt;INS&gt;</code>. This converted set was then used for tools that do not have a native Truvari-style <code>--dup-to-ins</code> option.</p>
"""
    dup_body += "<h3>DUP-to-INS converted input counts</h3>" + pre_from_file("reports/chr22/qc_tables/benchmark_ready_dup_to_ins_counts.txt")
    dup_body += "<h3>Truvari converted DUP-to-INS summary</h3>" + table_from_file("reports/chr22/metrics/truvari_converted_dup_to_ins_summary.tsv")
    dup_body += "<h3>Original vs converted DUP-to-INS Truvari comparison</h3>" + table_from_file("reports/chr22/metrics/truvari_original_vs_converted_dup_to_ins.tsv")
    dup_body += original_vs_dup_summary("reports/chr22/metrics/truvari_original_vs_converted_dup_to_ins.tsv")
    dup_body += "<h3>RTGtools DUP-to-INS sensitivity summaries</h3>" + pre_from_file("reports/chr22/metrics/rtg_vcfeval_dup_to_ins_summaries.txt", max_chars=60000)
    dup_body += "<h3>Witty.er DUP-to-INS output files</h3>" + pre_from_file("reports/chr22/metrics/wittyer_dup_to_ins_metric_like_files.txt")
    sections.append(section("duptoins", "8. DUP-to-INS sensitivity analysis", dup_body))

    sections.append(section("interpretation", "9. Interpretation", """
<p>The primary benchmark result is the Truvari comparison. Sniffles had the highest F1 score in both genome builds. DRAGEN showed high precision but lower recall, especially in hg38. DELLY performed poorly in hg19 but much better in hg38.</p>
<p>The benchmark should be interpreted as an interval-SV benchmark, not as a full assessment of every structural variant representation emitted by the callers. BND records were not silently removed; they were counted and reported separately because they are not directly comparable to interval-style SVs in a chromosome 22-only benchmark.</p>
<p>The DUP-to-INS sensitivity analysis tests whether apparent false negatives or false positives were caused by true biological disagreement or by DUP/INS representation differences.</p>
<p>The observed Truvari ranking by F1 in the original benchmark was hg19 Sniffles, hg38 Sniffles, hg38 DELLY, hg19 DRAGEN, hg38 DRAGEN, and hg19 DELLY.</p>
"""))

    sections.append(section("limitations", "10. Limitations", """
<ul>
<li>This was a chromosome 22 subset, not a whole-genome benchmark.</li>
<li>The main benchmark was restricted to interval-style SVs.</li>
<li>BND and other complex rearrangements require separate evaluation.</li>
<li>Differences between callers may reflect caller performance, but also representation differences, filtering, genotype encoding, symbolic allele handling, missing END/SVLEN fields, and confidence-region overlap.</li>
<li>RTGtools and Witty.er may behave differently from Truvari because they use different matching assumptions and input requirements.</li>
</ul>
"""))

    repro_body = "<h3>Reproducibility log</h3>" + pre_from_file("reports/chr22/logs/reproducibility_log.txt", max_chars=60000)
    repro_body += "<h3>Source file audit</h3>" + source_file_audit()
    repro_body += "<h3>Result directories found</h3>"
    repro_body += "<h4>Truvari original</h4>" + list_as_html(file_list("benchmark-runs/chr22/truvari/*/summary.json"))
    repro_body += "<h4>Truvari DUP-to-INS</h4>" + list_as_html(file_list("benchmark-runs/chr22/truvari_dup_to_ins/*/summary.json"))
    repro_body += "<h4>RTGtools original</h4>" + list_as_html(file_list("benchmark-runs/chr22/rtg/vcfeval/*/summary.txt"))
    repro_body += "<h4>RTGtools DUP-to-INS</h4>" + list_as_html(file_list("benchmark-runs/chr22/rtg/vcfeval_dup_to_ins/*/summary.txt"))
    sections.append(section("repro", "11. Reproducibility and source file audit", repro_body))

    html_doc = f"""<!doctype html>
<html lang="en">
<head>
<meta charset="utf-8">
<title>Chromosome 22 SV Benchmarking Report</title>
<meta name="viewport" content="width=device-width, initial-scale=1">
<style>{CSS}</style>
</head>
<body>
<header>
<h1>Chromosome 22 Structural Variant Benchmarking Report</h1>
<p>HG002 SV benchmark against GIAB truth, covering diagnosis, preprocessing, QC, Truvari, RTGtools, Witty.er, and DUP-to-INS sensitivity analysis.</p>
<p>Generated: {esc(generated)} from project root <code>{esc(str(ROOT))}</code></p>
</header>
<main>
{''.join(sections)}
</main>
<footer>
Generated from local result files under <code>reports/chr22/</code> and <code>benchmark-runs/chr22/</code>.
</footer>
</body>
</html>
"""
    OUT.write_text(html_doc, encoding="utf-8")
    print(f"Wrote {OUT}")
    print(f"Open with: open {OUT}")


if __name__ == "__main__":
    main()
