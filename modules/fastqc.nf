/*
 * FASTQC: read QC per sample
 * In:  [sample, R1, R2]
 * Out: <sample>_{1,2}_fastqc.{zip,html} (names parse_fastqc() expects), versions.txt
 */
process FASTQC {
    tag "${sample}"
    label 'fastqc'

    input:
    tuple val(sample), path(r1, stageAs: 'input_1.fq.gz'), path(r2, stageAs: 'input_2.fq.gz')

    output:
    path '*_fastqc.zip',  emit: zip
    path '*_fastqc.html', emit: html
    path 'versions.txt',  emit: versions

    script:
    """
    ln -s ${r1} ${sample}_1.fq.gz
    ln -s ${r2} ${sample}_2.fq.gz
    fastqc ${sample}_1.fq.gz ${sample}_2.fq.gz --memory ${params.fastqc_memory_mb} -t 1 -o .
    echo "FastQC: \$(fastqc --version | sed 's/^FastQC v*//')" > versions.txt
    """

    stub:
    """
    touch ${sample}_1_fastqc.zip ${sample}_2_fastqc.zip ${sample}_1_fastqc.html ${sample}_2_fastqc.html versions.txt
    """
}
