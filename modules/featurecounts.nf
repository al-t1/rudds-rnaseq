/*
 * FEATURECOUNTS: gene counts, all samples in one call
 * In:  all <sample>_RSR.bam, samplesheet (column order), GTF, featurecounts.R
 * Out: <counts_prefix>.rds, <counts_prefix>.summary (MultiQC), versions.txt
 */
process FEATURECOUNTS {
    label 'featurecounts'

    input:
    path bams
    path samplesheet
    path gtf
    path script

    output:
    path "${params.counts_prefix}.rds",     emit: rds
    path "${params.counts_prefix}.summary", emit: summary
    path 'versions.txt',                    emit: versions

    script:
    // Plain GTF, as in the original run
    def gz    = gtf.name.endsWith('.gz')
    def annot = gz ? 'annotation.gtf' : gtf
    def unzip = gz ? "gzip -dc ${gtf} > annotation.gtf" : ''
    """
    ${unzip}
    Rscript ${script} ${samplesheet} ${annot} ${task.cpus} ${params.counts_prefix}
    """

    stub:
    """
    touch ${params.counts_prefix}.rds ${params.counts_prefix}.summary versions.txt
    """
}
