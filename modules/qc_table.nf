/*
 * QC_TABLE: Supplementary Table 1
 * In:  core_functions.R (parsers), ST1 metadata, FastQC zips, RNA-SeQC metrics,
 *      featureCounts RDS, qc_table.R
 * Out: ST1_QC_metrics.tsv
 */
process QC_TABLE {
    label 'r_small'

    input:
    path core_functions
    path metadata
    path fastqc_zips,     stageAs: 'fastqc/*'
    path rnaseqc_metrics, stageAs: 'rnaseqc/*'
    path counts
    path script

    output:
    path 'ST1_QC_metrics.tsv', emit: st1

    script:
    """
    Rscript ${script} --core ${core_functions} --metadata ${metadata} \\
        --counts ${counts} --out ST1_QC_metrics.tsv
    """

    stub:
    """
    touch ST1_QC_metrics.tsv
    """
}
