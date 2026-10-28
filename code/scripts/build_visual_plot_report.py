#!/usr/bin/env python3
from pathlib import Path
import csv, html, math, re

ROOT = Path(__file__).resolve().parents[2]
QC = ROOT / 'reports/chr22' / 'qc_tables'
MET = ROOT / 'reports/chr22' / 'metrics'
OUT = ROOT / 'reports/chr22' / 'plot_report'
FIG = OUT / 'figures'
HTML = OUT / 'final_chr22_visual_report.html'
FIG.mkdir(parents=True, exist_ok=True)
OUT.mkdir(parents=True, exist_ok=True)

# ----------------------------
# Basic reading helpers
# ----------------------------
def read_tsv(path):
    path = Path(path)
    if not path.exists():
        return []
    with path.open(newline='', errors='replace') as f:
        return list(csv.DictReader(f, delimiter='\t'))

def to_float(x, default=0.0):
    try:
        if x is None or str(x).strip() == '':
            return default
        return float(x)
    except Exception:
        return default

def to_int(x, default=0):
    try:
        if x is None or str(x).strip() == '':
            return default
        return int(float(x))
    except Exception:
        return default

def esc(x):
    return html.escape(str(x))

def clean_run(run):
    return str(run).replace('hg19_', 'hg19 ').replace('hg38_', 'hg38 ').replace('_', ' ')

def caller_from_run(run):
    r = str(run).lower()
    if 'sniffles' in r:
        return 'Sniffles'
    if 'dragen' in r:
        return 'DRAGEN'
    if 'delly' in r:
        return 'DELLY'
    return str(run)

def build_from_run(run):
    r = str(run).lower()
    if r.startswith('hg19'):
        return 'hg19'
    if r.startswith('hg38'):
        return 'hg38'
    return ''

# ----------------------------
# SVG helpers, no matplotlib needed
# ----------------------------
PALETTE = ['#2563eb', '#16a34a', '#f97316', '#9333ea', '#dc2626', '#0891b2', '#ca8a04', '#4b5563']
BG = '#ffffff'
GRID = '#e5e7eb'
TEXT = '#111827'
SUBT = '#4b5563'

def save_svg(name, svg):
    path = FIG / name
    path.write_text(svg, encoding='utf-8')
    return path

def wrap_label(label, width=18):
    words = str(label).replace('_', ' ').split()
    lines, cur = [], ''
    for w in words:
        if len(cur) + len(w) + 1 <= width:
            cur = (cur + ' ' + w).strip()
        else:
            if cur:
                lines.append(cur)
            cur = w
    if cur:
        lines.append(cur)
    return lines[:3]

def svg_header(w, h, title, subtitle=''):
    return [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{w}" height="{h}" viewBox="0 0 {w} {h}">',
        f'<rect width="{w}" height="{h}" fill="{BG}"/>',
        f'<text x="30" y="34" font-family="Arial" font-size="20" font-weight="700" fill="{TEXT}">{esc(title)}</text>'
    ] + ([f'<text x="30" y="58" font-family="Arial" font-size="13" fill="{SUBT}">{esc(subtitle)}</text>'] if subtitle else [])

def svg_footer(parts, x, y):
    out = []
    xpos = x
    for i, label in enumerate(parts):
        col = PALETTE[i % len(PALETTE)]
        out.append(f'<rect x="{xpos}" y="{y-10}" width="12" height="12" rx="2" fill="{col}"/>')
        out.append(f'<text x="{xpos+17}" y="{y}" font-family="Arial" font-size="12" fill="{TEXT}">{esc(label)}</text>')
        xpos += 17 + len(label)*7 + 25
    return out

def grouped_bar_svg(title, subtitle, rows, category_key, series, filename, y_max=None, value_suffix='', rotate=False):
    # series = [(key, label)]
    w = 1200
    h = 680
    ml, mr, mt, mb = 80, 30, 90, 145
    plot_w = w - ml - mr
    plot_h = h - mt - mb
    if not rows:
        svg = '\n'.join(svg_header(w, 220, title, 'No data available.') + ['</svg>'])
        return save_svg(filename, svg)
    max_val = y_max if y_max is not None else max([to_float(r.get(k)) for r in rows for k, _ in series] + [1])
    if max_val <= 0:
        max_val = 1
    n = len(rows)
    group_w = plot_w / n
    bar_w = min(32, group_w / (len(series) + 1.3))
    out = svg_header(w, h, title, subtitle)
    # grid + y labels
    for tick in range(0, 6):
        val = max_val * tick / 5
        y = mt + plot_h - (val / max_val) * plot_h
        out.append(f'<line x1="{ml}" y1="{y:.1f}" x2="{w-mr}" y2="{y:.1f}" stroke="{GRID}"/>')
        lab = f'{val:.1f}{value_suffix}' if max_val > 1.2 else f'{val:.2f}{value_suffix}'
        if y_max == 1:
            lab = f'{val:.1f}'
        out.append(f'<text x="{ml-10}" y="{y+4:.1f}" font-family="Arial" font-size="11" text-anchor="end" fill="{SUBT}">{esc(lab)}</text>')
    out.append(f'<line x1="{ml}" y1="{mt}" x2="{ml}" y2="{mt+plot_h}" stroke="#9ca3af"/>')
    out.append(f'<line x1="{ml}" y1="{mt+plot_h}" x2="{w-mr}" y2="{mt+plot_h}" stroke="#9ca3af"/>')
    for i, r in enumerate(rows):
        cx = ml + i*group_w + group_w/2
        for j, (key, label) in enumerate(series):
            val = to_float(r.get(key))
            bh = (val / max_val) * plot_h
            x = cx - (len(series)*bar_w)/2 + j*bar_w
            y = mt + plot_h - bh
            out.append(f'<rect x="{x:.1f}" y="{y:.1f}" width="{bar_w-2:.1f}" height="{bh:.1f}" fill="{PALETTE[j%len(PALETTE)]}" rx="3"><title>{esc(r.get(category_key))} {esc(label)}: {val:.4g}</title></rect>')
        lines = wrap_label(r.get(category_key, ''), width=16)
        if rotate:
            out.append(f'<g transform="translate({cx:.1f},{mt+plot_h+18}) rotate(35)">')
            out.append(f'<text font-family="Arial" font-size="11" fill="{TEXT}">')
            for k, line in enumerate(lines):
                out.append(f'<tspan x="0" dy="{13 if k else 0}">{esc(line)}</tspan>')
            out.append('</text></g>')
        else:
            for k, line in enumerate(lines):
                out.append(f'<text x="{cx:.1f}" y="{mt+plot_h+18+k*13}" font-family="Arial" font-size="11" text-anchor="middle" fill="{TEXT}">{esc(line)}</text>')
    out += svg_footer([label for _, label in series], ml, h-25)
    out.append('</svg>')
    return save_svg(filename, '\n'.join(out))

def stacked_percent_svg(title, subtitle, rows, category_key, series, filename):
    w = 1200
    h = 700
    ml, mr, mt, mb = 230, 40, 90, 70
    plot_w = w - ml - mr
    plot_h = h - mt - mb
    if not rows:
        svg = '\n'.join(svg_header(w, 220, title, 'No data available.') + ['</svg>'])
        return save_svg(filename, svg)
    n = len(rows)
    row_h = min(42, plot_h / max(n, 1))
    out = svg_header(w, h, title, subtitle)
    for tick in range(0, 6):
        x = ml + plot_w * tick / 5
        out.append(f'<line x1="{x:.1f}" y1="{mt-10}" x2="{x:.1f}" y2="{mt+n*row_h}" stroke="{GRID}"/>')
        out.append(f'<text x="{x:.1f}" y="{mt-18}" font-family="Arial" font-size="11" text-anchor="middle" fill="{SUBT}">{tick*20}%</text>')
    for i, r in enumerate(rows):
        y = mt + i*row_h
        total = sum(max(0, to_float(r.get(k))) for k, _ in series)
        if total <= 0:
            total = 1
        x = ml
        for j, (key, label) in enumerate(series):
            val = max(0, to_float(r.get(key)))
            ww = plot_w * val / total
            out.append(f'<rect x="{x:.1f}" y="{y:.1f}" width="{ww:.1f}" height="{row_h*0.65:.1f}" fill="{PALETTE[j%len(PALETTE)]}" rx="3"><title>{esc(r.get(category_key))} {esc(label)}: {int(val)} ({100*val/total:.1f}%)</title></rect>')
            x += ww
        out.append(f'<text x="{ml-10}" y="{y+row_h*0.43:.1f}" font-family="Arial" font-size="11" text-anchor="end" fill="{TEXT}">{esc(r.get(category_key))}</text>')
    out += svg_footer([label for _, label in series], ml, h-25)
    out.append('</svg>')
    return save_svg(filename, '\n'.join(out))

def horizontal_bar_svg(title, subtitle, rows, category_key, value_key, filename, log=False):
    w = 1200
    h = max(380, 95 + len(rows)*42)
    ml, mr, mt, mb = 300, 50, 80, 50
    plot_w = w - ml - mr
    if not rows:
        return save_svg(filename, '\n'.join(svg_header(w, 220, title, 'No data available.') + ['</svg>']))
    vals = [to_float(r.get(value_key)) for r in rows]
    plot_vals = [math.log10(v+1) if log else v for v in vals]
    maxv = max(plot_vals + [1])
    out = svg_header(w, h, title, subtitle)
    for i, r in enumerate(rows):
        y = mt + i*42
        v = to_float(r.get(value_key))
        pv = math.log10(v+1) if log else v
        ww = plot_w * pv / maxv if maxv else 0
        out.append(f'<text x="{ml-12}" y="{y+20}" font-family="Arial" font-size="11" text-anchor="end" fill="{TEXT}">{esc(r.get(category_key))}</text>')
        out.append(f'<rect x="{ml}" y="{y}" width="{ww:.1f}" height="26" fill="{PALETTE[i%len(PALETTE)]}" rx="4"><title>{esc(r.get(category_key))}: {int(v)}</title></rect>')
        out.append(f'<text x="{ml+ww+8:.1f}" y="{y+18}" font-family="Arial" font-size="11" fill="{TEXT}">{int(v)}</text>')
    out.append('</svg>')
    return save_svg(filename, '\n'.join(out))

def scatter_arrow_svg(title, subtitle, rows, filename):
    w = 900
    h = 760
    ml, mr, mt, mb = 80, 40, 90, 70
    plot_w = w - ml - mr
    plot_h = h - mt - mb
    out = svg_header(w, h, title, subtitle)
    out.append('<defs><marker id="arrow" markerWidth="10" markerHeight="10" refX="7" refY="3" orient="auto" markerUnits="strokeWidth"><path d="M0,0 L0,6 L9,3 z" fill="#4b5563" /></marker></defs>')
    # axes and grid
    for tick in range(0, 6):
        x = ml + plot_w * tick / 5
        y = mt + plot_h - plot_h * tick / 5
        out.append(f'<line x1="{x:.1f}" y1="{mt}" x2="{x:.1f}" y2="{mt+plot_h}" stroke="{GRID}"/>')
        out.append(f'<line x1="{ml}" y1="{y:.1f}" x2="{ml+plot_w}" y2="{y:.1f}" stroke="{GRID}"/>')
        out.append(f'<text x="{x:.1f}" y="{mt+plot_h+22}" font-family="Arial" font-size="11" text-anchor="middle" fill="{SUBT}">{tick/5:.1f}</text>')
        out.append(f'<text x="{ml-10}" y="{y+4:.1f}" font-family="Arial" font-size="11" text-anchor="end" fill="{SUBT}">{tick/5:.1f}</text>')
    out.append(f'<line x1="{ml}" y1="{mt}" x2="{ml}" y2="{mt+plot_h}" stroke="#9ca3af"/>')
    out.append(f'<line x1="{ml}" y1="{mt+plot_h}" x2="{ml+plot_w}" y2="{mt+plot_h}" stroke="#9ca3af"/>')
    out.append(f'<text x="{ml+plot_w/2}" y="{h-20}" font-family="Arial" font-size="13" text-anchor="middle" fill="{TEXT}">Recall</text>')
    out.append(f'<text transform="translate(22,{mt+plot_h/2}) rotate(-90)" font-family="Arial" font-size="13" text-anchor="middle" fill="{TEXT}">Precision</text>')
    def xy(rec, prec):
        return ml + rec*plot_w, mt + plot_h - prec*plot_h
    for i, r in enumerate(rows):
        ox, oy = xy(to_float(r.get('original_recall')), to_float(r.get('original_precision')))
        dx, dy = xy(to_float(r.get('dup_to_ins_recall')), to_float(r.get('dup_to_ins_precision')))
        out.append(f'<line x1="{ox:.1f}" y1="{oy:.1f}" x2="{dx:.1f}" y2="{dy:.1f}" stroke="#4b5563" stroke-width="1.4" marker-end="url(#arrow)"/>')
        out.append(f'<circle cx="{ox:.1f}" cy="{oy:.1f}" r="5" fill="{PALETTE[0]}"><title>{esc(r.get("run"))} original</title></circle>')
        out.append(f'<rect x="{dx-5:.1f}" y="{dy-5:.1f}" width="10" height="10" fill="{PALETTE[1]}"><title>{esc(r.get("run"))} DUP-to-INS</title></rect>')
        out.append(f'<text x="{dx+8:.1f}" y="{dy-6:.1f}" font-family="Arial" font-size="10" fill="{TEXT}">{esc(r.get("run"))}</text>')
    out += svg_footer(['Original', 'DUP-to-INS'], ml, mt-18)
    out.append('</svg>')
    return save_svg(filename, '\n'.join(out))

# ----------------------------
# Load and prepare data
# ----------------------------
raw = read_tsv(QC / 'raw_chr22_record_summary.tsv')
bnd = read_tsv(QC / 'BND_counts.tsv')
endpos = read_tsv(QC / 'END_lt_POS_counts.tsv')
miss = read_tsv(QC / 'missing_END_SVLEN_counts.tsv')
tru = read_tsv(MET / 'truvari_summary_extended.tsv') or read_tsv(MET / 'truvari_summary.tsv')
comp = read_tsv(MET / 'truvari_original_vs_converted_dup_to_ins.tsv') or read_tsv(MET / 'truvari_original_vs_dup_to_ins.tsv')
svtype = read_tsv(MET / 'truvari_svtype_tp_fp_fn_counts.tsv')

def add_run_order(rows):
    def key(r):
        run = str(r.get('run', ''))
        b = 0 if run.startswith('hg19') else 1
        c = {'sniffles':0, 'dragen':1, 'delly':2, 'DELLY':2}
        cr = 9
        for name, val in c.items():
            if name.lower() in run.lower(): cr = val
        return (b, cr, run)
    return sorted(rows, key=key)

tru = add_run_order(tru)
comp = add_run_order(comp)

# Short labels for visual report
for r in tru:
    r['label'] = clean_run(r.get('run',''))
for r in comp:
    r['label'] = clean_run(r.get('run',''))

# Raw accounting rows: make readable label and numeric cols
raw_plot = []
for r in raw:
    raw_plot.append({
        'label': f"{r.get('build','')} {caller_from_run(r.get('file','')) if 'truth' not in str(r.get('file','')).lower() else 'Truth'}",
        'Interval SV': to_int(r.get('DEL_DUP_INS_INV')),
        'BND': to_int(r.get('BND')),
        'Unknown/missing': to_int(r.get('UNKNOWN_NA_missing')),
    })

# BND labels
bnd_plot = []
for r in bnd:
    filev = r.get('file') or r.get('File') or ''
    count_key = next((k for k in r if 'count' in k.lower()), list(r.keys())[-1])
    bnd_plot.append({'label': f"{r.get('build','')} {caller_from_run(filev) if 'truth' not in str(filev).lower() else 'Truth'}", 'BND': to_int(r.get(count_key))})
bnd_plot = sorted(bnd_plot, key=lambda x: x['BND'], reverse=True)

# QC issue plot from missing/endltpos as log issue count
issue_rows = []
# by file/build from miss table
for r in miss:
    filev = r.get('file','')
    label = f"{r.get('build','')} {caller_from_run(filev) if 'truth' not in str(filev).lower() else 'Truth'}"
    missing_end_key = next((k for k in r if 'missing' in k.lower() and 'end' in k.lower()), None)
    missing_svlen_key = next((k for k in r if 'svlen' in k.lower()), None)
    issue_rows.append({
        'label': label,
        'Missing END': to_int(r.get(missing_end_key)) if missing_end_key else 0,
        'Missing SVLEN': to_int(r.get(missing_svlen_key)) if missing_svlen_key else 0,
        'END < POS': 0
    })
# merge END<POS
for r in endpos:
    filev = r.get('file','')
    label = f"{r.get('build','')} {caller_from_run(filev) if 'truth' not in str(filev).lower() else 'Truth'}"
    count_key = next((k for k in r if 'count' in k.lower() or 'records' in k.lower()), list(r.keys())[-1])
    found = next((x for x in issue_rows if x['label'] == label), None)
    if found:
        found['END < POS'] = to_int(r.get(count_key))
    else:
        issue_rows.append({'label': label, 'Missing END': 0, 'Missing SVLEN': 0, 'END < POS': to_int(r.get(count_key))})
# convert to log scale values for plotting, preserve visual names
issue_log = []
for r in issue_rows:
    issue_log.append({
        'label': r['label'],
        'Missing END': math.log10(r['Missing END'] + 1),
        'Missing SVLEN': math.log10(r['Missing SVLEN'] + 1),
        'END < POS': math.log10(r['END < POS'] + 1),
    })

# FP/FN balance
fpfn_rows = []
for r in tru:
    fpfn_rows.append({'label': r.get('label'), 'FP': to_float(r.get('FP')), 'FN': to_float(r.get('FN'))})

# FN by SVTYPE plot: aggregate only result_file/fn
fn_sv = {}
for r in svtype:
    result = str(r.get('result_file') or r.get('result') or '').lower()
    if result != 'fn':
        continue
    run = r.get('run','')
    sv = r.get('SVTYPE') or r.get('svtype') or 'NA'
    cnt = to_int(r.get('count'))
    fn_sv.setdefault(run, {})[sv] = fn_sv.setdefault(run, {}).get(sv, 0) + cnt
fn_sv_rows = []
all_sv = sorted({sv for d in fn_sv.values() for sv in d})
for run, vals in sorted(fn_sv.items()):
    row = {'label': clean_run(run)}
    for sv in all_sv:
        row[sv] = vals.get(sv, 0)
    fn_sv_rows.append(row)

# ----------------------------
# Build all figures
# ----------------------------
figs = {}
figs['raw_accounting'] = stacked_percent_svg(
    'Raw chr22 record accounting',
    'Each row sums to 100%; shows what is interval-benchmarkable vs BND vs unknown/missing.',
    raw_plot, 'label', [('Interval SV','Interval SV'), ('BND','BND'), ('Unknown/missing','Unknown/missing')],
    '01_raw_chr22_record_accounting.svg'
)
figs['bnd'] = horizontal_bar_svg(
    'BND burden in raw chr22 VCFs',
    'BND was counted separately instead of silently mixed into interval-SV benchmarking.',
    bnd_plot, 'label', 'BND', '02_bnd_burden.svg'
)
figs['qc_issues'] = grouped_bar_svg(
    'Malformed / incomplete interval fields',
    'Log10(count + 1), so large truth-set missing fields and small END<POS signals are visible together.',
    issue_log, 'label', [('Missing END','Missing END'), ('Missing SVLEN','Missing SVLEN'), ('END < POS','END < POS')],
    '03_qc_issue_burden_log.svg', y_max=None, rotate=True
)
figs['main_metrics'] = grouped_bar_svg(
    'Original Truvari performance',
    'Precision, recall, and F1 for the original interval-SV benchmark.',
    tru, 'label', [('precision','Precision'), ('recall','Recall'), ('f1','F1')],
    '04_original_truvari_metrics.svg', y_max=1, rotate=True
)
figs['fpfn'] = grouped_bar_svg(
    'FP/FN balance explains F1 loss',
    'False negatives dominate low-recall runs, especially hg19 DELLY and hg38 DRAGEN.',
    fpfn_rows, 'label', [('FP','FP'), ('FN','FN')],
    '05_fp_fn_balance.svg', y_max=None, rotate=True
)
figs['fn_svtype'] = grouped_bar_svg(
    'False negatives by SVTYPE',
    'Shows which SV types drive missed truth calls in each original Truvari run.',
    fn_sv_rows, 'label', [(sv, sv) for sv in all_sv],
    '06_fn_by_svtype.svg', y_max=None, rotate=True
)
figs['precision_dup'] = grouped_bar_svg(
    'Precision: original vs DUP-to-INS',
    'Representation sensitivity check: precision before and after DUP calls are treated as INS.',
    comp, 'label', [('original_precision','Original'), ('dup_to_ins_precision','DUP-to-INS')],
    '07_precision_original_vs_dup.svg', y_max=1, rotate=True
)
figs['recall_dup'] = grouped_bar_svg(
    'Recall: original vs DUP-to-INS',
    'Representation sensitivity check: recall before and after DUP calls are treated as INS.',
    comp, 'label', [('original_recall','Original'), ('dup_to_ins_recall','DUP-to-INS')],
    '08_recall_original_vs_dup.svg', y_max=1, rotate=True
)
figs['f1_dup'] = grouped_bar_svg(
    'F1: original vs DUP-to-INS',
    'Representation sensitivity check: F1 before and after DUP calls are treated as INS.',
    comp, 'label', [('original_f1','Original'), ('dup_to_ins_f1','DUP-to-INS')],
    '09_f1_original_vs_dup.svg', y_max=1, rotate=True
)
figs['scatter'] = scatter_arrow_svg(
    'Precision-recall movement after DUP-to-INS',
    'Circle = original; square = DUP-to-INS. Arrows show whether representation conversion changes the result.',
    comp, '10_precision_recall_movement.svg'
)
figs['delta_f1'] = grouped_bar_svg(
    'ΔF1 from DUP-to-INS conversion',
    'Positive values mean the DUP-to-INS representation improved F1; zero means stable.',
    comp, 'label', [('delta_f1','ΔF1')],
    '11_delta_f1_dup_minus_original.svg', y_max=None, rotate=True
)

# ----------------------------
# Minimal HTML report: plot-first, no conclusions/references/logs/previews
# ----------------------------
def img(path, caption):
    rel = path.relative_to(OUT)
    return f'''<section class="figure"><img src="{esc(rel.as_posix())}" alt="{esc(caption)}"><p>{esc(caption)}</p></section>'''

# Small metric table only, no raw count tables
metric_rows = ''
for r in tru:
    metric_rows += f"<tr><td>{esc(r.get('label'))}</td><td>{to_float(r.get('precision')):.3f}</td><td>{to_float(r.get('recall')):.3f}</td><td>{to_float(r.get('f1')):.3f}</td><td>{to_int(r.get('TP_base'))}</td><td>{to_int(r.get('FP'))}</td><td>{to_int(r.get('FN'))}</td></tr>\n"
metric_table = f'''
<table>
<thead><tr><th>Run</th><th>Precision</th><th>Recall</th><th>F1</th><th>TP</th><th>FP</th><th>FN</th></tr></thead>
<tbody>{metric_rows}</tbody>
</table>
'''

dup_rows = ''
for r in comp:
    dup_rows += f"<tr><td>{esc(r.get('label'))}</td><td>{to_float(r.get('original_f1')):.3f}</td><td>{to_float(r.get('dup_to_ins_f1')):.3f}</td><td>{to_float(r.get('delta_f1')):.4f}</td><td>{to_float(r.get('original_recall')):.3f}</td><td>{to_float(r.get('dup_to_ins_recall')):.3f}</td></tr>\n"
dup_table = f'''
<table>
<thead><tr><th>Run</th><th>Original F1</th><th>DUP-to-INS F1</th><th>ΔF1</th><th>Original recall</th><th>DUP-to-INS recall</th></tr></thead>
<tbody>{dup_rows}</tbody>
</table>
'''

html_text = f'''<!doctype html>
<html lang="en">
<head>
<meta charset="utf-8">
<title>chr22 SV benchmark visual report</title>
<style>
body {{ font-family: Arial, Helvetica, sans-serif; margin: 36px; color: #111827; background: #fff; line-height: 1.45; }}
h1 {{ font-size: 30px; margin-bottom: 8px; }}
h2 {{ font-size: 22px; margin-top: 34px; border-bottom: 2px solid #e5e7eb; padding-bottom: 5px; }}
.lead {{ color: #374151; max-width: 1100px; }}
.grid {{ display: grid; grid-template-columns: 1fr; gap: 26px; }}
.figure {{ margin: 22px 0 34px 0; }}
.figure img {{ width: 100%; max-width: 1250px; border: 1px solid #e5e7eb; border-radius: 8px; background: white; }}
.figure p {{ font-size: 14px; color: #374151; margin: 8px 0 0 2px; max-width: 1050px; }}
table {{ border-collapse: collapse; width: 100%; max-width: 1100px; margin: 15px 0 30px 0; font-size: 14px; }}
th, td {{ border: 1px solid #d1d5db; padding: 7px 9px; text-align: left; }}
th {{ background: #f3f4f6; }}
.note {{ background: #f9fafb; border-left: 5px solid #2563eb; padding: 12px 16px; margin: 18px 0; max-width: 1100px; }}
</style>
</head>
<body>
<h1>Chromosome 22 SV benchmark — visual report</h1>
<p class="lead">Plot-focused report for the chr22 HG002 structural variant benchmark. This version removes long logs, tool previews, references, and conclusion sections. It keeps only minimal metric tables plus figures that show benchmark behaviour and representation sensitivity.</p>

<h2>Diagnosis and benchmark accounting</h2>
<div class="note"><strong>Scope:</strong> interval-SV benchmark = DEL, DUP, INS, INV. BND, unknown/missing SVTYPE, and malformed interval fields are visualised separately.</div>
{img(figs['raw_accounting'], 'Raw VCF accounting: shows how much of each chr22 file is interval-style SV, BND, or unknown/missing.')}
{img(figs['bnd'], 'BND burden: these were counted separately because breakend-style records are not directly comparable to interval DEL/DUP/INS/INV calls.')}
{img(figs['qc_issues'], 'QC issue burden: log-scaled view of missing END, missing SVLEN, and END < POS records.')}

<h2>Original benchmark performance</h2>
{img(figs['main_metrics'], 'Original Truvari precision, recall, and F1 across builds and callers.')}
{img(figs['fpfn'], 'FP/FN balance: highlights whether low F1 is mainly caused by false positives or false negatives.')}
{img(figs['fn_svtype'], 'FN by SVTYPE: shows which variant type drives missed truth calls.')}
<h3>Minimal original metric table</h3>
{metric_table}

<h2>DUP-to-INS sensitivity</h2>
<div class="note"><strong>Question:</strong> do apparent misses change when DUP records are represented as INS? This tests representation sensitivity, not a new primary benchmark.</div>
{img(figs['precision_dup'], 'Precision before and after DUP-to-INS conversion.')}
{img(figs['recall_dup'], 'Recall before and after DUP-to-INS conversion.')}
{img(figs['f1_dup'], 'F1 before and after DUP-to-INS conversion.')}
{img(figs['scatter'], 'Precision-recall movement: arrows show whether each run moves after DUP-to-INS conversion.')}
{img(figs['delta_f1'], 'ΔF1 from DUP-to-INS conversion: positive means DUP-to-INS improved F1.')}
<h3>Minimal DUP-to-INS metric table</h3>
{dup_table}
</body>
</html>
'''

HTML.write_text(html_text, encoding='utf-8')
print(f'Wrote {HTML}')
print(f'Wrote {len(list(FIG.glob("*.svg")))} SVG plots to {FIG}')
