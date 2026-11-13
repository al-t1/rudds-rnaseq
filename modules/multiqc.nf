/*
 * MULTIQC: one report over FastQC, RNA-SeQC and featureCounts
 * In:  QC files
 * Out: multiqc_report.html, multiqc_report_data/
 */
process MULTIQC {
    label 'multiqc'

    input:
    path qc_files, stageAs: 'qc/*'

    output:
    path 'multiqc_report.html', emit: report
    path '*_data',              emit: data

    script:
    // --filename: else --title renames the outputs
    """
    multiqc --force --title 'RUDDS RNA-seq' --filename multiqc_report.html qc
    """

    stub:
    """
    touch multiqc_report.html
    mkdir multiqc_data
    """
}
