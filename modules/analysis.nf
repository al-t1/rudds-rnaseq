/*
 * ANALYSIS: scripts/analysis.R, run whole
 * In:  script, functions, counts RDS, ST1, GTF, sample tables, frozen online
 *      resources (CollecTRI, MSigDB, KEGG, ExperimentHub), Open Targets tables,
 *      [BayesPrism reference]; staged under the relative paths the script reads
 * Out: Results/ (figures, tables, sessionInfo), ORA_layout_data.js
 */
process ANALYSIS {
    label 'analysis'

    input:
    path main_script,    stageAs: 'analysis.R'
    path core_functions, stageAs: 'core_functions.R'
    path counts,         stageAs: 'RSR_counts.rds'
    path st1,            stageAs: 'ST1_QC_metrics.tsv'
    path gtf,            stageAs: 'input_*'
    path jin_xlsx,       stageAs: '40779_2025_637_MOESM3_ESM.xlsx'
    path sample_meta,    stageAs: 'ST_safe_Sheets.tsv'
    path deconvolution,  stageAs: 'ST2_deconvolution.tsv'
    path ora_layout,     stageAs: 'ORA_layout_top_5_5.json'
    path collectri,      stageAs: 'CollecTRI_human.tsv'
    path mirdb,          stageAs: 'MSigDB_C3_MIRDB.tsv'
    path kegg,           stageAs: 'KEGG_hsa_pathways.tsv'
    path ehub,           stageAs: 'experimenthub_frozen'
    path ot_files
    path bayesprism_ref, stageAs: 'scrnaseq_reference'

    output:
    path 'Results',            emit: results
    path 'ORA_layout_data.js', emit: ora_editor_data, optional: true   // for ora_layout_editor.html
    path 'omnipathr-log',      emit: logs, optional: true

    script:
    // Script reads the GTF gzipped, under its original name
    def gtf_target = 'gencode.v47.basic.annotation.gtf.gz'
    def place_gtf  = gtf.name.endsWith('.gz') ? "ln -s ${gtf} ${gtf_target}" : "gzip -c ${gtf} > ${gtf_target}"
    """
    ${place_gtf}
    export RUN_BAYESPRISM=${params.run_bayesprism ? 'TRUE' : 'FALSE'}
    # Writable home and caches (container home is read-only)
    export HOME="\$PWD/home" XDG_CACHE_HOME="\$PWD/home/.cache" XDG_CONFIG_HOME="\$PWD/home/.config"
    mkdir -p "\$XDG_CACHE_HOME" "\$XDG_CONFIG_HOME"
    # Offline ExperimentHub: small files copied (BiocFileCache writes), large linked
    mkdir eh_cache
    for f in ${ehub}/*; do
        if [ "\$(stat -Lc %s "\$f")" -gt 100000000 ]; then ln -s "\$(readlink -f "\$f")" eh_cache/
        else cp -L "\$f" eh_cache/; fi
    done
    export EXPERIMENT_HUB_CACHE="\$PWD/eh_cache" EXPERIMENT_HUB_LOCAL=TRUE EXPERIMENT_HUB_ASK=FALSE
    Rscript analysis.R
    """

    stub:
    """
    mkdir Results
    touch Results/stub.txt
    """
}
