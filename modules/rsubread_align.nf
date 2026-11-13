/*
 * RSUBREAD_ALIGN: paired-end alignment
 * In:  [sample, R1, R2], index files, index prefix, rsubread_align.R
 * Out: [sample, <sample>_RSR.bam], <sample>_RSR.bam.summary, versions.txt
 */
process RSUBREAD_ALIGN {
    tag "${sample}"
    label 'rsubread_align'

    input:
    tuple val(sample), path(r1, stageAs: 'input_1.fq.gz'), path(r2, stageAs: 'input_2.fq.gz')
    path index
    val  index_name
    path script

    output:
    tuple val(sample), path("${sample}_RSR.bam"), emit: bam
    path "${sample}_RSR.bam.summary",             emit: summary, optional: true
    path 'versions.txt',                           emit: versions

    script:
    """
    Rscript ${script} ${index_name} ${sample} ${r1} ${r2} ${task.cpus}
    """

    stub:
    """
    touch ${sample}_RSR.bam ${sample}_RSR.bam.summary versions.txt
    """
}
