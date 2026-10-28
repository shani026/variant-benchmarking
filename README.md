# Variant benchmarking

Structural variant benchmarking for HG002 on hg19 and hg38, with retained
final reports and aggregate results. Reference genomes, truth sets, caller
VCFs, variant-level comparison tables, intermediates, execution logs and
Nextflow caches have been removed from the working tree.

## Final reports

- [hg19 complete report](reports/hg19/html/hg19_sv_benchmarking_full_results_report_WITH_SELF_SURVIVOR.html)
- [hg38 complete report](reports/hg38/html/hg38_sv_benchmarking_full_results_report_WITH_SELF_SURVIVOR.html)
- [hg19 MultiQC](reports/hg19/multiqc_report.html)
- [hg38 MultiQC](reports/hg38/multiqc_report.html)
- [chr22 clean report](reports/chr22/final_chr22_sv_benchmark_report_clean.html)
- [chr22 visual report](reports/chr22/plot_report/final_chr22_visual_report.html)
- [chr22 Truvari report](reports/chr22/truvari_simple_report/simple_truvari_report.html)

`reports/` also contains supporting plots and aggregate tables. Historical
paths or commands embedded in reports describe the original runs; the raw
files they describe are no longer present.

## Compact final results

`final-results/hg19/` and `final-results/hg38/` contain the consolidated
benchmark metrics, caller means and tool comparisons copied from the final
analysis tables. `final-results/chr22/` retains the final summary tables and
DUP-to-INS comparison metrics. These preserve the existing calculated values;
no benchmarking has been rerun or results recalculated during cleanup.

`final-results/retained-files.json` lists retained reports and result files
with SHA-256 hashes. Additional native-tool JSON summaries and ROC files
were removed during cleanup following a case-insensitive Results/results
path collision; the aggregate metric tables and final reports remain.

## Source layout

```text
code/
├── run-nfcore.sh                # launch the pinned nf-core workflow
├── config/                     # hg19/hg38 parameter and samplesheet templates
├── nextflow-custom/            # original local Truvari workflow, repaired
├── scripts/                    # report generation and analysis scripts
├── crossbuild-legacy/          # preserved historical liftover implementation
└── requirements.previous-environment.txt
reports/                        # completed reports, figures and tables
final-results/                  # compact final metrics and file inventory
tools/                          # retained RTG and Witty.er code/binaries
```

The legacy crossbuild source is retained because it differs from the main
`../crossbuild-sv` project. Vendor example genomic files have been removed;
vendor test suites that require those fixtures need a fresh upstream checkout.
Installed Nextflow assets and plugins in the home directory are retained.
Git history is unchanged and can still contain historical data.

## Run again with external data

The completed nf-core runs used `nf-core/variantbenchmarking` 1.5.0,
revision `8b21c01749`. The launch script pins that release. Replace the
`/absolute/path/to/...` placeholders in `code/config/params_hg19.yaml`,
`params_hg38.yaml` and their referenced samplesheets with external files.
From this project root:

```bash
code/run-nfcore.sh hg19
code/run-nfcore.sh hg38
```

This requires Nextflow, Docker, indexed reference FASTAs, indexed truth and
caller VCFs and the truth BED. New output goes to `benchmark-runs/nfcore/`
and temporary work to `work/`; both are excluded from Git.

The local Truvari alternative uses the samplesheets whose names contain a
period (`samplesheet.hg19.csv` and `samplesheet.hg38.csv`), with the schema
`sample,caller,vcf`. Update its separate parameter templates before running:

```bash
nextflow run code/nextflow-custom/main.nf \
  -params-file code/nextflow-custom/params.hg19.yaml
```

That workflow stages VCFs, FASTA and their indexes as file inputs, and lets
Truvari create its own output directory. Truvari must be installed in the
execution environment. The former Python environment's package versions are
recorded in `code/requirements.previous-environment.txt`.

Report generators are retained for reuse. Generators that inspect individual
VCFs or native benchmark outputs require those inputs to be supplied again;
the preserved reports can be viewed directly. New raw runs belong under
`benchmark-runs/`, separate from the preserved final deliverables.
