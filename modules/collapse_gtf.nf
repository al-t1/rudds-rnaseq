/*
 * COLLAPSE_GTF: RNA-SeQC gene model
 * In:  GENCODE GTF, collapse_annotation.py
 * Out: genes_stranded.gtf, versions.txt
 */
process COLLAPSE_GTF {
    tag "${gtf.name}"
    label 'collapse'

    input:
    path gtf
    path script

    output:
    path 'genes_stranded.gtf', emit: gtf
    path 'versions.txt',       emit: versions

    script:
    """
    python3 ${script} --stranded ${gtf} genes_stranded.gtf
    echo "bx-python: \$(python3 -c 'from importlib.metadata import version; print(version("bx-python"))')" > versions.txt
    """

    stub:
    """
    touch genes_stranded.gtf versions.txt
    """
}
