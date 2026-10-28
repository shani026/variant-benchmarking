nextflow.enable.dsl = 2

workflow {
    if (!params.samplesheet) error 'Missing required parameter: samplesheet'
    if (!params.reference) error 'Missing required parameter: reference'
    if (!params.truth_vcf) error 'Missing required parameter: truth_vcf'
    if (!params.truth_bed) error 'Missing required parameter: truth_bed'


    callers = Channel
        .fromPath(params.samplesheet, checkIfExists: true)
        .splitCsv(header: true)
        .map { row ->
            if (!row.sample || !row.caller || !row.vcf) {
                error 'Samplesheet columns must be sample,caller,vcf'
            }
            tuple(row.sample.toString(), row.caller.toString(),
                  file(row.vcf, checkIfExists: true),
                  file("${row.vcf}.tbi", checkIfExists: true))
        }

    TRUVARI_BENCH(
        callers,
        file(params.truth_vcf, checkIfExists: true),
        file("${params.truth_vcf}.tbi", checkIfExists: true),
        file(params.reference, checkIfExists: true),
        file("${params.reference}.fai", checkIfExists: true),
        file(params.truth_bed, checkIfExists: true)
    )
}

process TRUVARI_BENCH {
    tag "${sample}:${caller}"
    publishDir { "${params.outdir}/${caller}" }, mode: 'copy', overwrite: true

    input:
    tuple val(sample), val(caller), path(query_vcf), path(query_index)
    path truth_vcf
    path truth_index
    path reference
    path reference_index
    path truth_bed

    output:
    path 'bench', emit: benchmark

    script:
    """
    set -euo pipefail
    truvari bench \
        -b "${truth_vcf}" \
        -c "${query_vcf}" \
        -f "${reference}" \
        --includebed "${truth_bed}" \
        -o bench
    """
}
