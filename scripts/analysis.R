# RUDDS RNA-seq analysis: DGE, enrichment, figures and tables of the manuscript
# In (relative to the working dir):
#   core_functions.R                 functions
#   RSR_counts.rds                   featureCounts result
#   ST_safe_Sheets.tsv               sample table
#   ST1_QC_metrics.tsv               QC metrics and outlier flags (ST1)
#   ST2_deconvolution.tsv            cell-type fractions (unless RUN_BAYESPRISM)
#   gencode.v47.basic.annotation.gtf.gz (or gencode47_gene_table.tsv)
#   40779_2025_637_MOESM3_ESM.xlsx, *_OT-*-associated-targets-*.tsv
#   ORA_layout_top_5_5.json          Figure 2 layout
#   CollecTRI_human.tsv, MSigDB_C3_MIRDB.tsv, KEGG_hsa_pathways.tsv  frozen
#     online resources (downloaded when absent); ExperimentHub EH3226/EH3228/EH3234
#   scrnaseq_reference/              BayesPrism reference (RUN_BAYESPRISM only)
# Out: Results/  T1, F1-F3, SF1-SF5, ST2-ST10, sessionInfo.txt
# Env: RUN_BAYESPRISM=TRUE recomputes ST2 (hours)
# R 4.5.1, Bioconductor 3.21

# In RStudio: run from this script's folder
if (requireNamespace("rstudioapi", quietly = TRUE) && rstudioapi::isAvailable()) {
  .p <- rstudioapi::getSourceEditorContext()$path
  if (nzchar(.p)) setwd(dirname(.p))
}
RESULTS_DIR <- "Results"
dir.create(RESULTS_DIR, showWarnings = FALSE, recursive = TRUE)
res_path <- function(file) file.path(RESULTS_DIR, file)

source("core_functions.R")

# ── 1. Setup ──────────────────────────────────────────────────────────────────
# Attached
packages <- c("ggplot2", "ggrepel", "patchwork", "dplyr",
              "DESeq2", "sva", "BiocParallel", "org.Hs.eg.db",
              "clusterProfiler", "enrichplot", "decoupleR", "ggraph",
              "grid", "ComplexHeatmap", "circlize",
              "msigdbr", "fgsea", "signatureSearch", "ExperimentHub")
# Used as pkg::fun
packages_ns <- c("AnnotationDbi", "data.table", "tibble", "tidyr", "tidygraph", "igraph",
                 "limma", "matrixStats", "DOSE", "ggpubr", "gridExtra", "ggforce",
                 "scatterpie", "png", "jsonlite", "openxlsx", "rstudioapi",
                 "rtracklayer")

RUN_BAYESPRISM <- as.logical(Sys.getenv("RUN_BAYESPRISM", "FALSE"))
stopifnot(!is.na(RUN_BAYESPRISM))
if (RUN_BAYESPRISM)
  packages_ns <- c(packages_ns, "hdf5r", "Matrix", "Danko-Lab/BayesPrism/BayesPrism")

install_and_load(packages, install_only = packages_ns)

# ── 2. Data ───────────────────────────────────────────────────────────────────
sample_table <- data.table::fread("ST_safe_Sheets.tsv", data.table = FALSE) %>%
  mutate(Sex = as.factor(Sex),
         Condition = relevel(factor(Condition), ref = "control"))
rownames(sample_table) <- sample_table$Name

# Gene table (ID, Transcript.ID, Type, Name, Length): from the GTF, else the saved table
GTF_FILE <- "gencode.v47.basic.annotation.gtf.gz"
if (file.exists(GTF_FILE)) {
  gene_table <- get_gencode_gene_table(GTF_FILE)
  data.table::fwrite(gene_table, res_path("gencode47_gene_table.tsv"), sep = "\t",
                     eol = "\n")
} else {
  gene_table <- data.table::fread("gencode47_gene_table.tsv", data.table = FALSE)
}

# ST1: QC metrics and outlier flags; skip its 3 header lines
st1 <- data.table::fread("ST1_QC_metrics.tsv", skip = 3,
                         data.table = FALSE, check.names = FALSE)
stopifnot(setequal(st1$Name, sample_table$Name),
          all(st1$Outlier %in% c("yes", "no")))

rsr <- readRDS("RSR_counts.rds")
colnames(rsr$counts) <- sub("_RSR[.]bam$", "", colnames(rsr$counts))
stopifnot(setequal(colnames(rsr$counts), sample_table$Name))

# Raw counts: ID + one column per sample
data.table::fwrite(data.frame(ID = rownames(rsr$counts), rsr$counts,
                              check.names = FALSE),
                   res_path("RSR_counts.tsv"), sep = "\t")

# Mean TPM per gene (ST8)
tpm_mat  <- calc_tpm(rsr$counts, rsr$annotation$Length, gene_table)
mean_tpm <- data.frame(Name = rownames(tpm_mat), Mean_TPM = rowMeans(tpm_mat))
rm(tpm_mat)

# ── 3. Deconvolution: BayesPrism (ST2) ────────────────────────────────────────
# Reference: Gan et al. (2021), scrnaseq_reference/
if (RUN_BAYESPRISM) {
  bayesprism <- exec_bayesprism(rsr$counts, gene_table, ref_dir = "scrnaseq_reference")
  saveRDS(bayesprism$bp, res_path("BayesPrism_result.rds"))
  theta_df <- bayesprism$theta
  data.table::fwrite(tibble::rownames_to_column(theta_df, "sample"),
                     res_path("ST2_deconvolution.tsv"), sep = "\t")
} else {
  theta_df <- data.table::fread("ST2_deconvolution.tsv", data.table = FALSE) %>%
    tibble::column_to_rownames("sample")
}

# ── 4. DESeq2: full (72) and clean (57) cohorts ───────────────────────────────
dds <- DESeqDataSetFromMatrix(countData = rsr$counts,
                              colData = sample_table[colnames(rsr$counts), ],
                              design = ~ Sex + Age + Condition)
dds <- estimateSizeFactors(dds)
# Outliers (SF1 PCA cut-offs): ST1 column Outlier
outs <- st1$Name[st1$Outlier == "yes"]
stopifnot(length(outs) == 15)
clean_samples <- setdiff(sample_table$Name, outs)
stopifnot(length(clean_samples) == 57)

dge_fits <- list()
# Full: dds only (PCA)
dge_fits$full <- exec_deseq_downweight(dds = dds,
                                       generate_dds_only = TRUE,
                                       result_name = "full",
                                       keep_samples = NULL)
# Clean: primary DGE
dge_fits$clean <- exec_deseq_downweight(dds = dds,
                                        deseq_formula = ~ SV1 + SV2 + Sex + Age + Condition,
                                        sv_formula_full = ~ Sex + Age + Condition,
                                        sv_formula_null = ~ Sex + Age,
                                        cooks_iqr_threshold = 50,
                                        result_name = "clean",
                                        keep_samples = colnames(dds) %in% clean_samples)

# ── 5. Table 1: cohort description ────────────────────────────────────────────
# Manuscript layout: UTF-8, CRLF, no final newline
cohort_stats <- function(st) {
  mean_sd <- function(x) paste0(round(mean(x, na.rm = TRUE), 2), " ± ", round(sd(x, na.rm = TRUE), 2))
  sapply(split(st, st$Condition)[c("control", "case")], function(d) c(
    "Total Samples" = nrow(d),
    "% Male"        = round(mean(d$Sex == "Male") * 100, 2),
    "Mean Age ± SD" = mean_sd(d$Age_years),
    "Mean BMI ± SD" = mean_sd(d$BMI),
    "Mean RIN ± SD" = mean_sd(d$RIN)))
}
table1 <- cbind(cohort_stats(sample_table),
                cohort_stats(sample_table[clean_samples, ]))
t1_lines <- c("Table 1. RuDDs summary description before quality control.\t\t\t\t",
              "\t\t\t\t",
              "\tBefore outlier deletion\t\tAfter outlier deletion\t",
              paste(c("Condition", colnames(table1)), collapse = "\t"),
              paste(rownames(table1), apply(table1, 1, paste, collapse = "\t"), sep = "\t"))
writeBin(charToRaw(enc2utf8(paste(t1_lines, collapse = "\r\n"))), res_path("T1_cohort.tsv"))

# ── 6. RNA quality control ────────────────────────────────────────────────────

# ── 6.1 Outlier PCA (SF1) ─────────────────────────────────────────────────────
vsd_full  <- vst(dge_fits$full$dds, blind = TRUE)
vsd_clean <- vst(dge_fits$clean$dds, blind = TRUE)
# SV1, SV2 removed; Condition, Sex, Age protected
vsd_clean_rbe <- vsd_clean
assay(vsd_clean_rbe) <- limma::removeBatchEffect(
  assay(vsd_clean),
  covariates = cbind(vsd_clean$SV1, vsd_clean$SV2),
  design = model.matrix(~ Condition + Sex + Age, colData(vsd_clean)))

pca_panel <- function(vsd, title, pcs)
  plot_dynamic_pca(vsd, plot_title = title, title_size = 24, intgroup = "Condition",
                   pcs_to_plot = pcs, ntop = nrow(vsd))
# Dashed outlier cut-off
cut_line <- function(geom, ...)
  geom(..., show.legend = NA, color = ggplot2::alpha("black", 0.5), linetype = "dashed")

# A-C: all 72 with cut-offs; D: clean 57; E: clean 57 without SV1, SV2
sf1 <- list(
  A = pca_panel(vsd_full, "A", c("PC1", "PC2")) + cut_line(geom_abline, slope = -0.75, intercept = -55),
  B = pca_panel(vsd_full, "B", c("PC3", "PC4")) + cut_line(geom_abline, slope = 1, intercept = -50),
  C = pca_panel(vsd_full, "C", c("PC5", "PC6")) + cut_line(geom_vline, xintercept = 50),
  D = pca_panel(vsd_clean, "D", c("PC1", "PC2")),
  E = pca_panel(vsd_clean_rbe, "E", c("PC1", "PC2")))
p_sf1 <- (sf1$A + sf1$B + sf1$C) / (sf1$D + sf1$E) +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom")
save_ggplot_png_pdf(p_sf1, res_path("SF1_PCA"), width = 8000, height = 6000, dpi = 600)

# ── 6.2 Pilot-protocol samples (Supplementary Note 1) ─────────────────────────
# 4 pilot vs. 68 standard: unassigned reads, detected genes; Wilcoxon rank-sum
pilot_samples <- c("45_ShAP", "K20", "N1", "R6")
qc_metrics    <- c("Unassigned_NoFeatures", "Genes Detected")

st1$Protocol <- factor(ifelse(st1$Name %in% pilot_samples, "pilot", "standard"),
                       levels = c("pilot", "standard"))
stopifnot(sum(st1$Protocol == "pilot") == length(pilot_samples))

pilot_stats <- do.call(rbind, lapply(qc_metrics, function(m) {
  x <- st1[[m]][st1$Protocol == "pilot"]
  y <- st1[[m]][st1$Protocol == "standard"]
  data.frame(Metric        = m,
             Mean_pilot    = mean(x, na.rm = TRUE),
             Mean_standard = mean(y, na.rm = TRUE),
             Wilcoxon_p    = wilcox.test(x, y, exact = FALSE)$p.value)
}))
pilot_stats$Mean_pilot    <- round(pilot_stats$Mean_pilot, 4)
pilot_stats$Mean_standard <- round(pilot_stats$Mean_standard, 4)
pilot_stats$Wilcoxon_p    <- signif(pilot_stats$Wilcoxon_p, 3)
write.table(pilot_stats, res_path("QC_pilot_batch_stats.tsv"), sep = "\t",
            quote = FALSE, row.names = FALSE)

# ── 7. Cell-type composition and heterogeneity ────────────────────────────────
# ── 7.1 Deconvolution (SF2) ───────────────────────────────────────────────────
plot_deconvolution(theta_df, outs, file = res_path("SF2_deconvolution"), k = 3)

# ── 7.2 Jin et al. (2025) subtypes (SF3, SF4) ─────────────────────────────────
jin_clusters <- exec_jin_heatmaps(counts = rsr$counts, theta = theta_df,
                                  sample_table = sample_table,
                                  clean_samples = clean_samples,
                                  gene_table = gene_table,
                                  out_dir = RESULTS_DIR,
                                  jin_xlsx = "40779_2025_637_MOESM3_ESM.xlsx",
                                  files = c(signature  = "SF3_Jin_signature",
                                            classifier = "SF4_Jin_classifier"))

# ── 7.3 SVs vs. heterogeneity (SF5) ───────────────────────────────────────────
# A: signature heatmap; B: SV / cluster / Condition correlations; C: SV R^2
sv_correlation   <- exec_sv_correlation(dge_fits$clean$dds, jin_clusters$clean57,
                                        file = res_path("SV_correlation_clean57.tsv"))
sv_heterogeneity <- exec_sv_heterogeneity(dge_fits$clean$dds, jin_clusters$clean57, theta_df,
                                          file = res_path("SV_heterogeneity_clean57.tsv"))
plot_sv_figure(jin_clusters$figures$clean57_signature, sv_correlation$ht, sv_heterogeneity$ht,
               file = res_path("SF5_SV_heterogeneity"))

# ── 8. Differentially expressed genes ─────────────────────────────────────────

# ── 8.1 Annotation, Open Targets overlap (ST3, ST4) ───────────────────────────
# Order breaks ties (get_dge_ot)
OT_FILES <- c(
  lumbar_disc_degeneration         = "lumbar_disc_degeneration_OT-EFO_0004994-associated-targets-6_5_2026-v26_03.tsv",
  lumbar_disc_herniation           = "lumbar_disc_herniation_OT-EFO_1002005-associated-targets-6_5_2026-v26_03.tsv",
  intervertebral_disk_degeneration = "intervertebral_disk_degeneration_OT-HP_0008419-associated-targets-6_5_2026-v26_03.tsv")

dge_ot <- get_dge_ot(res_df = dge_fits$clean$all_results$clean,
                     gene_table = gene_table,
                     gene_types = c("protein_coding", "lncRNA"),
                     ot_files = OT_FILES)
res      <- dge_ot$res
ot_genes <- dge_ot$ot_genes
degs <- res[res$padj < 0.05 & abs(res$logFC) >= log2(1.2), ]

# Fisher's exact test: DEG x Open Targets
n_deg_ot <- sum(degs$Name %in% ot_genes$symbol)
n_deg_not_ot <- sum(!(degs$Name %in% ot_genes$symbol))
n_ot_not_deg <- sum(!(ot_genes$symbol %in% degs$Name))
n_neither <- nrow(res) - n_deg_ot - n_deg_not_ot - n_ot_not_deg
fisher_res <- fisher.test(
  matrix(c(n_deg_ot, n_deg_not_ot, n_ot_not_deg, n_neither), nrow = 2),
  alternative = "two.sided")
cat("OT.DEGs:", n_deg_ot, "\nOR:", fisher_res$estimate,
    "\np.val:", fisher_res$p.value)
fisher_ot <- data.frame(DEG_OT       = n_deg_ot,
                        DEG_notOT    = n_deg_not_ot,
                        notDEG_OT    = n_ot_not_deg,
                        notDEG_notOT = n_neither,
                        Odds_ratio   = unname(fisher_res$estimate),
                        CI95_low     = fisher_res$conf.int[1],
                        CI95_high    = fisher_res$conf.int[2],
                        P_value      = fisher_res$p.value)
write.table(fisher_ot, res_path("Fisher_DEG_OpenTargets.tsv"), sep = "\t",
            quote = FALSE, row.names = FALSE)

write.table(x = res, file = res_path("ST3_DGE_all_genes.tsv"), sep = "\t",
            quote = TRUE, row.names = FALSE)
write.table(x = degs, file = res_path("ST4_DEGs.tsv"), sep = "\t",
            quote = FALSE, row.names = FALSE)

# ── 8.2 Volcano plot (F1) ─────────────────────────────────────────────────────
plot_volcano(degs, ot_genes, file = res_path("F1_volcano"))

# ── 9. Pathway enrichment ─────────────────────────────────────────────────────
# Frozen KEGG (KEGG REST, 2026-11-12) put into clusterProfiler's KEGG cache, so
# enrichKEGG / gseKEGG do not download
KEGG_FILE <- "KEGG_hsa_pathways.tsv"
if (file.exists(KEGG_FILE)) {
  kegg_frozen <- data.table::fread(KEGG_FILE, data.table = FALSE,
                                   colClasses = "character", na.strings = "")
  kegg_env <- clusterProfiler:::get_KEGG_Env()
  assign("organism", "hsa", envir = kegg_env)
  assign("_type_", "KEGG", envir = kegg_env)
  assign("KEGGPATHID2EXTID", with(kegg_frozen[!is.na(kegg_frozen$gene), ],
                                  data.frame(from = pathway, to = gene)), envir = kegg_env)
  assign("KEGGPATHID2NAME", unique(with(kegg_frozen, data.frame(from = pathway, to = name))),
         envir = kegg_env)
}

# ── 9.1 ORA: GO and KEGG (ST5, ST6) ───────────────────────────────────────────
# Up / down DEGs vs. all tested genes; figure: top 5 GO + top 5 KEGG per direction
ora <- exec_ora_go_kegg(degs = degs, background = res, top_n = c(GO = 5, KEGG = 5))

lfc_by_entrez <- ora$lfc$entrez

# ── 9.2 ORA network (F2) ──────────────────────────────────────────────────────
# Writes F2_ORA_network_draft and editor data; with ORA_layout_<fig_id>.json
# (saved from ora_layout_editor.html) also the final F2_ORA_network
ora_net <- ora_network_spec(ora, bold_genes = ot_genes$symbol)   # Open Targets genes in bold
exec_ora_network(ora_net, out_dir = RESULTS_DIR,
                 draft_name = "F2_ORA_network_draft", final_name = "F2_ORA_network")
# ora_unpatch_enrichplot()   # restore enrichplot's pie drawing

write.table(x = ora$tables$up, file = res_path("ST5_ORA_up.tsv"), sep = "\t",
            quote = FALSE, row.names = FALSE)
write.table(x = ora$tables$down, file = res_path("ST6_ORA_down.tsv"), sep = "\t",
            quote = FALSE, row.names = FALSE)

# ── 9.3 GSEA: KEGG (ST7) ──────────────────────────────────────────────────────
# Ranked by Wald statistic
res_sorted        <- res[order(res$stat, decreasing = TRUE), ]
res_sorted_entrez <- res_sorted[res_sorted$Entrez.ID != "0", ]
stat_all_entrez   <- setNames(as.numeric(res_sorted_entrez$stat),
                              res_sorted_entrez$Entrez.ID)
set.seed(326)

gsea_kegg <- gseKEGG(
  geneList      = stat_all_entrez,
  organism      = "hsa",
  minGSSize     = 50,
  by            = "DOSE",
  nPerm         = 10000,
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  seed          = TRUE)
gsea_kegg         <- setReadable(gsea_kegg, "org.Hs.eg.db", "ENTREZID")
gsea_kegg         <- append_kegg_category(gsea_kegg)
gsea_kegg@result  <- subset(gsea_kegg@result,
                            category != c("Human Diseases") &
                              p.adjust < 0.01 & subcategory != c("Digestive system"))
gsea_kegg_results <- gsea_kegg@result
write.table(x = gsea_kegg_results, file = res_path("ST7_GSEA_KEGG.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

# ── 10. TF activity: decoupleR (F3, ST8) ──────────────────────────────────────
# CollecTRI ULM on the Wald statistic. TFs drawn: padj < TF_PADJ, >= TF_MIN_CONSISTENT
# DEG targets consistent with the TF's direction, consistent fraction > TF_MIN_FRAC,
# binomial P < TF_CONSISTENCY_P. ST8_TF_activity: all TFs; *_filtered: TFs drawn
TF_PADJ           <- 0.01
TF_MIN_CONSISTENT <- 5
TF_CONSISTENCY_P  <- 0.05
TF_MIN_FRAC       <- 0.75
TF_NET_TARGETS    <- 10   # targets drawn per TF (largest |Wald stat|)
# Frozen CollecTRI (OmniPath, 2026-06-26)
COLLECTRI_FILE <- "CollecTRI_human.tsv"
collectri <- if (file.exists(COLLECTRI_FILE)) {
  data.table::fread(COLLECTRI_FILE, data.table = FALSE)
} else {
  decoupleR::get_collectri(organism = "human", split_complexes = FALSE)
}
tf <- exec_tf_activity(res = res, degs = degs, mean_tpm = mean_tpm,
                       gene_table = gene_table, net = collectri,
                       padj_cutoff = TF_PADJ, min_consistent = TF_MIN_CONSISTENT,
                       consistency_p_cutoff = TF_CONSISTENCY_P, min_frac = TF_MIN_FRAC)
write.table(x = tf$full, file = res_path("ST8_TF_activity.tsv"), sep = "\t",
            quote = FALSE, row.names = FALSE)
write.table(x = tf$network, file = res_path("ST8.TFEA_decoupleR_filtered.tsv"), sep = "\t",
            quote = FALSE, row.names = FALSE)

# Open Targets genes underlined
p_tfea <- plot_tf_network(tf, degs, underline_genes = ot_genes$symbol, n_targets = TF_NET_TARGETS)
save_ggplot_png_pdf(p_tfea, res_path("F3_TF_network"), width = 6000, height = 4000, dpi = 800)

# ── 11. miRNA targets: GSEA, MSigDB C3 miRDB (ST9) ────────────────────────────
# Frozen MSigDB 2026.1.Hs C3 MIR:MIRDB (msigdbr 26.1.0)
MIRDB_FILE <- "MSigDB_C3_MIRDB.tsv"
mir_raw <- if (file.exists(MIRDB_FILE)) {
  data.table::fread(MIRDB_FILE, data.table = FALSE, colClasses = "character")
} else {
  msigdbr(species = "Homo sapiens", collection = "C3", subcollection = "MIR:MIRDB")
}
mir_sets <- split(as.character(mir_raw$ncbi_gene), mir_raw$gs_name)
# One value per Entrez ID: highest |stat|
res_mir_df <- data.frame(entrez = names(stat_all_entrez), stat = as.numeric(stat_all_entrez))
res_mir_df <- res_mir_df[order(abs(res_mir_df$stat), decreasing = TRUE), ]
res_mir_df <- res_mir_df[!duplicated(res_mir_df$entrez), ]
res_mir_df <- res_mir_df[order(res_mir_df$stat, decreasing = TRUE), ]
res_mir_stats <- setNames(res_mir_df$stat, res_mir_df$entrez)

set.seed(326)
gsea_mir <- fgsea(pathways    = mir_sets,
                  stats       = res_mir_stats,
                  minSize     = 14,
                  maxSize     = 500,
                  nPermSimple = 20000,
                  nproc       = 1)
gsea_mir     <- as.data.frame(gsea_mir)
gsea_mir_sig <- gsea_mir[!is.na(gsea_mir$padj) & gsea_mir$padj < 0.01, ]
gsea_mir_sig <- gsea_mir_sig[order(gsea_mir_sig$NES, decreasing = TRUE), ]
# miRBase names from the description (MIR103A_3P_MIR107 -> miR-103a-3p/miR-107)
mir_desc <- unique(mir_raw[, c("gs_name", "gs_description")])
mir_desc$miRNA_label <- vapply(
  regmatches(mir_desc$gs_description, gregexpr("hsa-(miR|let)-[^ ,;]+", mir_desc$gs_description)),
  function(x) paste(sub("^hsa-", "", x), collapse = "/"), character(1))
gsea_mir_sig$miRNA_label <- mir_desc$miRNA_label[match(gsea_mir_sig$pathway, mir_desc$gs_name)]
gsea_mir_sig$leadingEdge_sym <- vapply(gsea_mir_sig$leadingEdge, function(ids)
  paste(res$Name[match(ids, res$Entrez.ID)], collapse = "/"), character(1))

gsea_mir_out <- gsea_mir_sig[, c("miRNA_label", "pval", "padj",
                                 "log2err", "ES", "NES", "size", "leadingEdge_sym")]
write.table(gsea_mir_out, file = res_path("ST9_GSEA_miRNA.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

# ── 12. Connectivity mapping: LINCS L1000 (ST10) ──────────────────────────────
# ExperimentHub EH3226 (LINCS mod-z scores)
upset_ids   <- names(lfc_by_entrez[lfc_by_entrez > 0])
downset_ids <- names(lfc_by_entrez[lfc_by_entrez < 0])
eh <- ExperimentHub()
lincs_db_path <- eh[["EH3226"]]

qsig_lincs <- qSig(
  query = list(upset = upset_ids, downset = downset_ids),
  gess_method = "LINCS",
  refdb = lincs_db_path
)

lincs_res <- gess_lincs(
  qSig = qsig_lincs,
  sortby = "NCS",
  tau = FALSE,
  workers = 12   # MulticoreParam inside; serial on Windows
)

reversal_candidates <- result(lincs_res) %>%
  arrange(NCS) %>%
  filter(WTCS_FDR < 0.01)
write.table(
  reversal_candidates,
  file = res_path("ST10_CMAP_LINCS.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)

# ── 13. Session info ──────────────────────────────────────────────────────────
writeLines(capture.output(sessionInfo()), res_path("sessionInfo.txt"))
