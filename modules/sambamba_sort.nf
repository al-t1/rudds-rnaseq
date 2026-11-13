/*
 * SAMBAMBA_SORT: coordinate sort for RNA-SeQC
 * In:  [sample, BAM]
 * Out: [sample, <sample>_RSR.sorted.bam], .bai, versions.txt
 */
process SAMBAMBA_SORT {
    tag "${sample}"
    label 'sambamba'

    input:
    tuple val(sample), path(bam)

    output:
    tuple val(sample), path("${sample}_RSR.sorted.bam"), emit: bam
    path "${sample}_RSR.sorted.bam.bai",                 emit: bai, optional: true
    path 'versions.txt',                                  emit: versions

    script:
    def mem_mb = (task.memory.toMega() * 0.8) as long   // sort buffer; temp files in work dir
    """
    sambamba sort -l 9 -p -t ${task.cpus} -m ${mem_mb}M --tmpdir . -o ${sample}_RSR.sorted.bam ${bam}
    echo "sambamba: \$(sambamba --version 2>&1 | grep -m1 -o '[0-9][0-9.]*')" > versions.txt
    """

    stub:
    """
    touch ${sample}_RSR.sorted.bam versions.txt
    """
}
