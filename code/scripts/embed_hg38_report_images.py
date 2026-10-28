from pathlib import Path
import re
import base64
import mimetypes

REPORT = Path(__file__).resolve().parents[2] / Path("reports/hg38/html/hg38_sv_benchmarking_full_results_report.html")
OUT = Path(__file__).resolve().parents[2] / Path("reports/hg38/html/hg38_sv_benchmarking_full_results_report_EMBEDDED.html")

html = REPORT.read_text()

def embed_img(match):
    src = match.group(1)

    if src.startswith("data:") or src.startswith("http"):
        return match.group(0)

    img_path = (REPORT.parent / src).resolve()

    if not img_path.exists():
        print("Missing image:", img_path)
        return match.group(0)

    mime = mimetypes.guess_type(img_path)[0] or "image/png"
    encoded = base64.b64encode(img_path.read_bytes()).decode("ascii")

    return match.group(0).replace(src, f"data:{mime};base64,{encoded}")

html2 = re.sub(r'<img\s+[^>]*src="([^"]+)"[^>]*>', embed_img, html)

OUT.write_text(html2)

print("Wrote embedded report:")
print(OUT)
