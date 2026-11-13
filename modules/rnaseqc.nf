/*
 * RNASEQC: alignment QC, reverse-stranded
 * In:  [sample, sorted BAM], collapsed GTF
 * Out: <sample>.metrics.tsv, <sample>.*.gct (gene/exon reads, TPM), versions.txt
 */
process RNASEQC {
    tag "${sample}"
    label 'rnaseqc'

    input:
    tuple val(sample), path(bam)
    path gtf

    output:
    path "${sample}.metrics.tsv", emit: metrics
    path "${sample}.*.gct",       emit: gct
    path 'versions.txt',          emit: versions

    script:
    """
    rnaseqc ${gtf} ${bam} . -s ${sample} --stranded RF -vv --legacy
    echo "RNA-SeQC: \$(rnaseqc --version 2>&1 | grep -m1 -o '[0-9][0-9.]*')" > versions.txt
    """

    stub:
    """
    touch ${sample}.metrics.tsv ${sample}.gene_reads.gct versions.txt
    """
}
