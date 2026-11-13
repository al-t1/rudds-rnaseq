/*
 * RSUBREAD_INDEX: genome index; skipped with --rsubread_index
 * In:  genome FASTA, rsubread_index.R, index prefix
 * Out: <index_name>.*, versions.txt
 */
process RSUBREAD_INDEX {
    tag "${index_name}"
    label 'rsubread_index'

    input:
    path fasta
    path script
    val  index_name

    output:
    path "${index_name}.*", emit: index
    path 'versions.txt',    emit: versions

    script:
    """
    Rscript ${script} ${fasta} ${index_name}
    """

    stub:
    """
    touch ${index_name}.00.b.array ${index_name}.files versions.txt
    """
}
