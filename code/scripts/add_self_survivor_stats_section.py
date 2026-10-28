from pathlib import Path
import subprocess
import shutil
import pandas as pd
import matplotlib.pyplot as plt
import re
import base64
import mimetypes
import html as html_lib

ROOT = Path(__file__).resolve().parents[2]

REPORT_DIR = ROOT / "reports/hg38"
TABLES = REPORT_DIR / "tables"
PLOTS = REPORT_DIR / "plots"
HTML_DIR = REPORT_DIR / "html"
SELF = REPORT_DIR / "self_survivor"
PLAIN = SELF / "plain_vcfs"

for d in [TABLES, PLOTS, HTML_DIR, SELF, PLAIN]:
    d.mkdir(parents=True, exist_ok=True)

STRUCT = ROOT / "benchmark-runs/nfcore/hg38/structural"

PREPARED_VCFS = {
    "Truth": STRUCT / "HG002/preprocess/HG002_GRCh38_v5.0q_stvar.reheader.filtermissing.vcf.gz",
    "Delly": STRUCT / "HG002_delly_hg38/preprocess/HG002_delly_lr.hg38.renamed.reheader.filtermissing.vcf.gz",
    "DRAGEN": STRUCT / "HG002_dragen_hg38/preprocess/HG002_dragen.sv.hg38.renamed.reheader.filtermissing.vcf.gz",
    "Sniffles": STRUCT / "HG002_sniffles_hg38/preprocess/HG002_sniffles.hg38.renamed.reheader.filtermissing.vcf.gz",
}

NFCORE_SURVIVOR_STATS = {
    "Truth": STRUCT / "HG002/stats/survivor/HG002_mqc.stats",
    "Delly": STRUCT / "HG002_delly_hg38/stats/survivor/HG002_delly_hg38.delly_mqc.stats",
    "DRAGEN": STRUCT / "HG002_dragen_hg38/stats/survivor/HG002_dragen_hg38.dragen_mqc.stats",
    "Sniffles": STRUCT / "HG002_sniffles_hg38/stats/survivor/HG002_sniffles_hg38.sniffles_mqc.stats",
}

DOCKER_IMAGE = "informationsea/survivor"

def run_command(cmd):
    p = subprocess.run(cmd, text=True, capture_output=True)
    return p.returncode, p.stdout, p.stderr

def savefig(name):
    plt.tight_layout()
    plt.savefig(PLOTS / name, dpi=220, bbox_inches="tight")
    plt.close()

def slug(x):
    return re.sub(r"[^A-Za-z0-9]+", "_", x).strip("_")

def parse_survivor_stdout(stdout, source):
    rows = []
    lines = [l.strip() for l in stdout.splitlines() if l.strip()]

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
                        "origin": "self_SURVIVOR_stdout",
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

    lines = [l.strip() for l in path.read_text(errors="replace").splitlines() if l.strip()]
    if len(lines) >= 2:
        header = lines[0].split()
        vals = lines[1].split()
        if len(header) == len(vals):
            for h, v in zip(header, vals):
                try:
                    rows.append({
                        "source": source,
                        "metric": h,
                        "value": float(v),
                        "origin": origin,
                        "file": str(path),
                    })
                except Exception:
                    pass

    # fallback parser
    for line in lines:
        parts = line.replace(":", " ").split()
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
                "file": str(path),
            })

    return rows

# clean old self-survivor files except plain vcfs folder
for p in SELF.glob("*"):
    if p.is_file():
        p.unlink()

command_rows = []
self_file_rows = []
self_stat_rows = []

# 1. Decompress prepared VCFs to plain VCF because SURVIVOR stats is more reliable on uncompressed VCF
for source, gzvcf in PREPARED_VCFS.items():
    plain_vcf = PLAIN / f"{slug(source)}.prepared.vcf"

    if not gzvcf.exists():
        command_rows.append({
            "source": source,
            "step": "decompress",
            "command": "",
            "return_code": "",
            "status": "input VCF missing",
            "stdout": "",
            "stderr": str(gzvcf),
        })
        continue

    cmd = ["bcftools", "view", "-o", str(plain_vcf), str(gzvcf)]
    rc, stdout, stderr = run_command(cmd)

    command_rows.append({
        "source": source,
        "step": "decompress",
        "command": " ".join(cmd),
        "return_code": rc,
        "status": "success" if rc == 0 and plain_vcf.exists() else "failed",
        "stdout": stdout[-1000:],
        "stderr": stderr[-1000:],
    })

# 2. Run SURVIVOR stats through Docker with the full correct argument set
docker_ok = shutil.which("docker") is not None

for source in PREPARED_VCFS:
    plain_vcf = PLAIN / f"{slug(source)}.prepared.vcf"
    out_prefix = SELF / f"{slug(source)}_self_survivor"

    if not docker_ok:
        command_rows.append({
            "source": source,
            "step": "SURVIVOR stats",
            "command": "docker",
            "return_code": "",
            "status": "docker not found",
            "stdout": "",
            "stderr": "Docker was not found in PATH.",
        })
        continue

    if not plain_vcf.exists():
        command_rows.append({
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
        "-1",
        "-1",
        "-1",
        container_out,
    ]

    rc, stdout, stderr = run_command(cmd)

    produced = sorted(SELF.glob(f"{out_prefix.name}*"))

    command_rows.append({
        "source": source,
        "step": "SURVIVOR stats",
        "command": " ".join(cmd),
        "return_code": rc,
        "status": "success" if rc == 0 and produced else "failed",
        "stdout": stdout[-2000:],
        "stderr": stderr[-2000:],
    })

    self_stat_rows.extend(parse_survivor_stdout(stdout, source))

    for f in produced:
        if f.is_file():
            self_file_rows.append({
                "source": source,
                "file": str(f),
                "filename": f.name,
                "size_bytes": f.stat().st_size,
            })
            self_stat_rows.extend(parse_survivor_file(f, source, "self_SURVIVOR_file"))

# 3. Parse nf-core SURVIVOR stats too
nfcore_rows = []
for source, f in NFCORE_SURVIVOR_STATS.items():
    nfcore_rows.extend(parse_survivor_file(f, source, "nfcore_SURVIVOR_file"))

commands_df = pd.DataFrame(command_rows)
self_files_df = pd.DataFrame(self_file_rows)
self_df = pd.DataFrame(self_stat_rows)
nfcore_df = pd.DataFrame(nfcore_rows)

commands_df.to_csv(TABLES / "self_survivor_command_log.csv", index=False)
self_files_df.to_csv(TABLES / "self_survivor_output_files.csv", index=False)
self_df.to_csv(TABLES / "self_survivor_parsed_stats.csv", index=False)
nfcore_df.to_csv(TABLES / "nfcore_survivor_parsed_stats.csv", index=False)

# 4. Compare self vs nf-core
diff_df = pd.DataFrame()
if not self_df.empty and not nfcore_df.empty:
    a = self_df.groupby(["source", "metric"], as_index=False)["value"].max()
    b = nfcore_df.groupby(["source", "metric"], as_index=False)["value"].max()

    merged = a.merge(b, on=["source", "metric"], suffixes=("_self", "_nfcore"), how="inner")
    if not merged.empty:
        merged["difference_self_minus_nfcore"] = merged["value_self"] - merged["value_nfcore"]
        diff_df = merged

diff_df.to_csv(TABLES / "self_survivor_vs_nfcore_survivor_differences.csv", index=False)

# 5. Plots
if not commands_df.empty:
    status_counts = commands_df["status"].value_counts().reset_index()
    status_counts.columns = ["status", "count"]

    ax = status_counts.set_index("status")["count"].plot(kind="bar", figsize=(8, 4))
    ax.set_title("Self-run SURVIVOR command status")
    ax.set_xlabel("Status")
    ax.set_ylabel("Count")
    plt.xticks(rotation=25, ha="right")
    savefig("self_survivor_command_status.png")

if not self_df.empty:
    keep = ["Tot", "DEL", "DUP", "INS", "INV", "TRA"]
    plot_df = self_df[self_df["metric"].isin(keep)].copy()

    if not plot_df.empty:
        plot_df = plot_df.groupby(["source", "metric"], as_index=False)["value"].max()
        pivot = plot_df.pivot_table(index="metric", columns="source", values="value", fill_value=0)
        ax = pivot.plot(kind="bar", figsize=(10, 5), logy=True)
        ax.set_title("Self-computed SURVIVOR stats")
        ax.set_xlabel("SURVIVOR metric")
        ax.set_ylabel("Count, log scale")
        ax.legend(title="VCF source", bbox_to_anchor=(1.02, 1), loc="upper left")
        plt.xticks(rotation=0)
        savefig("self_survivor_metrics_log.png")

if not diff_df.empty:
    d = diff_df.copy()
    d["abs_difference"] = d["difference_self_minus_nfcore"].abs()
    d = d.sort_values("abs_difference", ascending=False).head(30)
    d["source_metric"] = d["source"] + " · " + d["metric"]

    ax = d.set_index("source_metric")["difference_self_minus_nfcore"].plot(kind="bar", figsize=(12, 5))
    ax.axhline(0, linewidth=1)
    ax.set_title("Self SURVIVOR minus nf-core SURVIVOR differences")
    ax.set_xlabel("Source · metric")
    ax.set_ylabel("Difference")
    plt.xticks(rotation=35, ha="right")
    savefig("self_survivor_vs_nfcore_differences.png")

# 6. Patch HTML
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
      <h3>{html_lib.escape(title)}</h3>
      <img src="{img_data_uri(path)}" alt="{html_lib.escape(title)}">
    </div>
    """

def df_to_html(df, max_rows=80):
    if df is None or df.empty:
        return "<p>No rows found.</p>"
    d = df.copy()
    if len(d) > max_rows:
        d = d.head(max_rows)
        note = f"<p class='small-note'>Showing first {max_rows} rows. Full table is saved in the tables folder.</p>"
    else:
        note = ""
    return note + d.to_html(index=False, classes="tbl", border=0, escape=True, float_format=lambda x: f"{x:.6g}")

survivor_section = f"""
<section>
  <h2>Self-computed SURVIVOR statistics</h2>
  <div class="note">
    This section decompresses the prepared hg38 VCFs to plain VCF, runs <code>SURVIVOR stats</code> through Docker using <code>-1 -1 -1</code> thresholds, parses the self-computed outputs, and compares them with the nf-core SURVIVOR statistics files.
  </div>

  <h3>Self-run SURVIVOR command log</h3>
  {df_to_html(commands_df, max_rows=50)}

  <h3>Self-run SURVIVOR output files</h3>
  {df_to_html(self_files_df, max_rows=50)}

  <div class="plot-grid">
    {plot_card("self_survivor_command_status.png", "Self-run SURVIVOR command status")}
    {plot_card("self_survivor_metrics_log.png", "Self-computed SURVIVOR statistics")}
    {plot_card("self_survivor_vs_nfcore_differences.png", "Self SURVIVOR vs nf-core SURVIVOR differences")}
  </div>

  <h3>Parsed self-computed SURVIVOR statistics</h3>
  {df_to_html(self_df, max_rows=100)}

  <h3>Parsed nf-core SURVIVOR statistics</h3>
  {df_to_html(nfcore_df, max_rows=100)}

  <h3>Self-computed SURVIVOR vs nf-core SURVIVOR differences</h3>
  {df_to_html(diff_df, max_rows=100)}
</section>
"""

candidate_reports = [
    HTML_DIR / "hg38_sv_benchmarking_full_results_report_EMBEDDED_CLEANED.html",
    HTML_DIR / "hg38_sv_benchmarking_full_results_report_EMBEDDED.html",
    HTML_DIR / "hg38_sv_benchmarking_full_results_report_CLEANED.html",
    HTML_DIR / "hg38_sv_benchmarking_full_results_report.html",
]

report_in = None
for r in candidate_reports:
    if r.exists():
        report_in = r
        break

if report_in is None:
    raise FileNotFoundError("Could not find existing hg38 report HTML.")

report_html = report_in.read_text()

report_html = re.sub(
    r"\n<section>\s*<h2>Self-computed SURVIVOR statistics</h2>.*?</section>\s*",
    "\n",
    report_html,
    flags=re.DOTALL,
)

insert_point = report_html.rfind("</div>\n</body>")
if insert_point == -1:
    insert_point = report_html.rfind("</body>")

report_html = report_html[:insert_point] + survivor_section + "\n" + report_html[insert_point:]

report_out = HTML_DIR / "hg38_sv_benchmarking_full_results_report_WITH_SELF_SURVIVOR.html"
report_out.write_text(report_html)

print()
print("DONE")
print("Report:")
print(report_out)
print()
print("Self SURVIVOR files:")
print(SELF)
print()
print("Main tables:")
print(TABLES / "self_survivor_command_log.csv")
print(TABLES / "self_survivor_output_files.csv")
print(TABLES / "self_survivor_parsed_stats.csv")
print(TABLES / "nfcore_survivor_parsed_stats.csv")
print(TABLES / "self_survivor_vs_nfcore_survivor_differences.csv")
