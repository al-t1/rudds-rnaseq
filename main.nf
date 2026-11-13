#!/usr/bin/env nextflow
/*
 * RUDDS lumbar disc RNA-seq: FASTQ -> manuscript figures and tables
 *   UPSTREAM  FastQC; Rsubread align -> sambamba sort -> RNA-SeQC;
 *             featureCounts -> ST1 (QC_TABLE); MultiQC
 *   ANALYSIS  scripts/analysis.R: DESeq2, enrichment, F1-F3, SF1-SF5, ST2-ST10
 * --counts <rds> skips UPSTREAM (ST1 from --st1)
 */

include { RSUBREAD_INDEX } from './modules/rsubread_index'
include { COLLAPSE_GTF   } from './modules/collapse_gtf'
include { FASTQC         } from './modules/fastqc'
include { RSUBREAD_ALIGN } from './modules/rsubread_align'
include { SAMBAMBA_SORT  } from './modules/sambamba_sort'
include { RNASEQC        } from './modules/rnaseqc'
include { FEATURECOUNTS  } from './modules/featurecounts'
include { QC_TABLE       } from './modules/qc_table'
include { MULTIQC        } from './modules/multiqc'
include { ANALYSIS       } from './modules/analysis'

// params.<name> as an existing file, else stop
def required(String name) {
    def p = params[name]
    if (!p) error "Missing --${name}"
    def f = file(p)
    if (!f.exists()) error "--${name}: file not found: ${p}"
    return f
}

/*
 * UPSTREAM
 * Take: samplesheet (sample,fastq_1,fastq_2), GENCODE GTF
 * Emit: counts (featureCounts RDS), st1 (ST1_QC_metrics.tsv)
 */
workflow UPSTREAM {
    take:
    samplesheet
    gtf

    main:
    def scripts = "${projectDir}/scripts"
    def versions = channel.empty()
    def index
    def index_name
    def collapsed_gtf

    def reads = channel.fromPath(samplesheet)
        .splitCsv(header: true, strip: true)
        .map { row ->
            if (!row.sample || !row.fastq_1 || !row.fastq_2)
                error "Samplesheet row without sample, fastq_1 or fastq_2: ${row}"
            tuple(row.sample, file(row.fastq_1, checkIfExists: true), file(row.fastq_2, checkIfExists: true))
        }

    // Index: existing prefix or built from --fasta
    if (params.rsubread_index) {
        def prefix = file(params.rsubread_index)
        def index_files = files("${params.rsubread_index}.*")
        if (!index_files) error "--rsubread_index: no files ${params.rsubread_index}.*"
        index      = channel.value(index_files)
        index_name = prefix.name
    } else {
        RSUBREAD_INDEX(required('fasta'), file("${scripts}/rsubread_index.R"), params.index_name)
        index      = RSUBREAD_INDEX.out.index
        index_name = params.index_name
        versions   = versions.mix(RSUBREAD_INDEX.out.versions)
    }

    // RNA-SeQC gene model: existing or built
    if (params.collapsed_gtf) {
        collapsed_gtf = channel.value(required('collapsed_gtf'))
    } else {
        COLLAPSE_GTF(gtf, required('collapse_script'))
        collapsed_gtf = COLLAPSE_GTF.out.gtf
        versions      = versions.mix(COLLAPSE_GTF.out.versions)
    }

    FASTQC(reads)
    RSUBREAD_ALIGN(reads, index, index_name, file("${scripts}/rsubread_align.R"))
    SAMBAMBA_SORT(RSUBREAD_ALIGN.out.bam)
    RNASEQC(SAMBAMBA_SORT.out.bam, collapsed_gtf)

    // Unsorted BAMs, as in the original run
    FEATURECOUNTS(RSUBREAD_ALIGN.out.bam.map { _sample, bam -> bam }.collect(),
                  samplesheet, gtf, file("${scripts}/featurecounts.R"))

    QC_TABLE(required('core_functions'), required('st1_metadata'),
             FASTQC.out.zip.collect(), RNASEQC.out.metrics.collect(),
             FEATURECOUNTS.out.rds,
             file("${scripts}/qc_table.R"))

    if (!params.skip_multiqc)
        MULTIQC(FASTQC.out.zip.collect()
                    .mix(RNASEQC.out.metrics.collect(), FEATURECOUNTS.out.summary)
                    .collect())

    versions
        .mix(FASTQC.out.versions.first(), RSUBREAD_ALIGN.out.versions.first(),
             SAMBAMBA_SORT.out.versions.first(), RNASEQC.out.versions.first())
        .collectFile(name: 'software_versions.txt', storeDir: "${params.outdir}/pipeline_info",
                     sort: true, newLine: false)

    emit:
    counts = FEATURECOUNTS.out.rds
    st1    = QC_TABLE.out.st1
}

workflow {
    def gtf = required('gtf')
    def counts
    def st1

    if (params.counts) {
        counts = channel.value(required('counts'))
        st1    = channel.value(required('st1'))
    } else {
        UPSTREAM(required('input'), gtf)
        counts = UPSTREAM.out.counts
        st1    = UPSTREAM.out.st1
    }

    if (!params.skip_analysis) {
        def ot = files(params.ot_files)
        if (ot.size() != 3) error "--ot_files: expected 3 Open Targets tables, found ${ot.size()}: ${params.ot_files}"
        if (params.run_bayesprism && !params.bayesprism_ref) error "--run_bayesprism needs --bayesprism_ref"

        ANALYSIS(required('analysis_script'), required('core_functions'),
                 counts, st1, gtf, required('jin_xlsx'),
                 required('sample_metadata'), required('deconvolution'),
                 required('ora_layout'), required('collectri'), required('mirdb'),
                 required('kegg'), required('experimenthub'), ot,
                 params.run_bayesprism ? required('bayesprism_ref') : [])
    }
}
