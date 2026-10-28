from pathlib import Path
import re

REPORTS = [
    Path(__file__).resolve().parents[2] / Path("reports/hg38/html/hg38_sv_benchmarking_full_results_report.html"),
    Path(__file__).resolve().parents[2] / Path("reports/hg38/html/hg38_sv_benchmarking_full_results_report_EMBEDDED.html"),
]

def remove_section(html, h2_text):
    pattern = re.compile(
        r'\n<section>\s*<h2>' + re.escape(h2_text) + r'</h2>.*?</section>\s*',
        re.DOTALL
    )
    return pattern.sub("\n", html)

for report in REPORTS:
    if not report.exists():
        continue

    html = report.read_text()

    # Remove sections you said you do not want
    remove_sections = [
        "3. Methods validation file checklist",
        "4. Nextflow trace status counts",
        "22. Main input paths",
        "23. Generated files",
    ]

    for sec in remove_sections:
        html = remove_section(html, sec)

    # Rename confusing PASS card wording
    html = html.replace("PASS:", "FILTER=PASS records:")

    # Add short explanation under input overview
    html = html.replace(
        "<h2>2. Input record overview</h2>",
        """<h2>2. Input record overview</h2>
  <div class="note">
    These cards show the number of records in the prepared truth and caller VCFs after nf-core preprocessing.
    FILTER=PASS records means records where the VCF FILTER column is exactly PASS.
  </div>"""
    )

    out = report.with_name(report.stem + "_CLEANED.html")
    out.write_text(html)
    print("Wrote:", out)

