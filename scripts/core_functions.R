# Functions of analysis.R 

# ── Packages, DESeq2, PCA ─────────────────────────────────────────────────────
# Installs missing packages (CRAN / Bioconductor via BiocManager; 'owner/repo'
# via remotes, named if the package name differs), never updates; attaches pkgs,
# loads install_only. install = FALSE: check only. Stops with a hint on failure
install_and_load <- function(pkgs, install_only = character(), install = TRUE) {
  specs  <- c(pkgs, install_only)
  github <- grepl("/", specs)
  nms    <- names(specs) %||% rep("", length(specs))
  names(specs) <- ifelse(nzchar(nms), nms, ifelse(github, sub("@.*$", "", basename(specs)), specs))
  attach <- names(specs)[seq_along(pkgs)]
  is_installed <- function(p) length(find.package(p, quiet = TRUE)) > 0

  todo <- specs[!vapply(names(specs), is_installed, logical(1))]
  if (length(todo) && !install)
    stop("Not installed: ", paste(names(todo), collapse = ", "), ". Install with\n  ",
         "BiocManager::install(c(", paste0("'", todo, "'", collapse = ", "), "), update = FALSE)",
         call. = FALSE)
  if (length(todo)) {
    message("Installing: ", paste(names(todo), collapse = ", "))
    if (!is_installed("BiocManager"))
      install.packages("BiocManager", repos = "https://cloud.r-project.org")
    from_repos <- unname(todo[!grepl("/", todo)])
    if (length(from_repos))
      tryCatch(BiocManager::install(from_repos, update = FALSE, ask = FALSE),
               error = function(e) message("BiocManager: ", conditionMessage(e)))
    from_github <- todo[grepl("/", todo)]
    if (length(from_github) && !is_installed("remotes"))
      BiocManager::install("remotes", update = FALSE, ask = FALSE)
    for (s in from_github)
      tryCatch(remotes::install_github(s, upgrade = "never"),
               error = function(e) message("GitHub ", s, ": ", conditionMessage(e)))
  }

  hint_install <- function(p) {
    if (grepl("/", specs[[p]]))
      return(paste0("GitHub install of ", specs[[p]], " failed: check owner/repo, the network ",
                    "(GitHub rate limit: set GITHUB_PAT) and build tools"))
    avail <- tryCatch(rownames(available.packages(repos = BiocManager::repositories())),
                      error = function(e) NULL)
    if (!is.null(avail) && !p %in% avail)
      return(paste0("not on CRAN or Bioconductor ", BiocManager::version(), " (R ", getRversion(),
                    "): check the spelling; archived on CRAN -> remotes::install_version(); ",
                    'GitHub only -> list it as "owner/repo"; R-universe -> ',
                    'install.packages("', p, '", repos = "https://<owner>.r-universe.dev"); ',
                    "other Bioconductor release -> other R version"))
    paste0("available but did not install (see the messages above): usually build tools ",
           "(Rtools on Windows, Xcode command line tools on macOS), system libraries (Linux), ",
           "no write access to the library, or a network problem")
  }
  hint_load <- function(p)
    paste0("installed but does not load: often built for another R / Bioconductor version; ",
           'reinstall with BiocManager::install("', p, '", force = TRUE), see BiocManager::valid()')

  problems <- character()
  for (p in names(specs)) {
    if (!is_installed(p)) {
      problems[p] <- hint_install(p)
      next
    }
    ok <- if (p %in% attach) suppressPackageStartupMessages(require(p, character.only = TRUE))
          else requireNamespace(p, quietly = TRUE)
    if (!ok) problems[p] <- hint_load(p)
  }
  if (length(problems))
    stop("Package problems:\n", paste0("  ", names(problems), ": ", problems, collapse = "\n"),
         call. = FALSE)
  invisible(attach)
}

# DESeq2 with 2 SVA surrogate variables; Cook's outliers (> Q3 + cooks_iqr_threshold
# * IQR per gene, or NA) weighted 0 in the final fit
# Out: list(dds, all_results[[result_name]], outlier_matrix)
exec_deseq_downweight <- function(dds,
                                  deseq_formula = ~ SV1 + SV2 + Sex + Age + Condition,
                                  sv_formula_full = ~ Sex + Age + Condition,
                                  sv_formula_null = ~ Sex + Age,
                                  cooks_iqr_threshold = 50,
                                  result_name = paste0("Run_", as.character(Sys.time())),
                                  keep_samples = NULL,
                                  generate_dds_only = FALSE) {

  if (!is.null(keep_samples)) {
    dds <- dds[, keep_samples]
  }

  keep <- rowSums(counts(dds) >= 25) >= 15
  dds <- dds[keep, ]

  if (generate_dds_only) {
    dds <- DESeq(dds)
    return(list(
      dds = dds,
      all_results = list()
    ))
  }

  dat  <- counts(dds, normalized = TRUE)
  mod  <- model.matrix(sv_formula_full, colData(dds))
  mod0 <- model.matrix(sv_formula_null, colData(dds))
  svseq <- svaseq(dat, mod, mod0, n.sv = 2)

  dds$SV1 <- svseq$sv[, 1]
  dds$SV2 <- svseq$sv[, 2]

  design(dds) <- deseq_formula
  dds <- DESeq(dds)   # Cook's distances under the SV design

  cooks_res <- assays(dds)[["cooks"]]
  na_mask   <- is.na(cooks_res)

  n_valid    <- rowSums(!na_mask)
  q1         <- matrixStats::rowQuantiles(cooks_res, probs = 0.25, na.rm = TRUE)
  q3         <- matrixStats::rowQuantiles(cooks_res, probs = 0.75, na.rm = TRUE)
  iqr_vals   <- q3 - q1
  thresholds <- ifelse(iqr_vals == 0, q3 + 1e-10, q3 + cooks_iqr_threshold * iqr_vals)
  thresholds[n_valid < 4] <- Inf   # too few values: only NAs flagged

  iqr_mask             <- sweep(cooks_res, 1, thresholds, `>`)
  iqr_mask[is.na(iqr_mask)] <- FALSE
  outlier_mask         <- na_mask | iqr_mask

  assay(dds, "weights", withDimnames = FALSE) <- matrix(1, nrow(dds), ncol(dds))
  assays(dds)[["weights"]][outlier_mask] <- 0

  dds <- DESeq(dds)

  res <- results(dds,
                 alpha = 0.05,
                 independentFiltering = FALSE,
                 pAdjustMethod = "BH")

  all_results <- list()

  res_df <- res %>%
    as.data.frame() %>%
    na.omit() %>%
    tibble::rownames_to_column(var = "GeneID") %>%
    dplyr::rename(logFC = log2FoldChange)

  result_name_str <- as.character(result_name)

  all_results[[result_name_str]] <- res_df

  list(
    dds = dds,
    all_results = all_results,
    outlier_matrix = outlier_mask
  )
}


# PCA plot of two PCs (e.g. c('PC1', 'PC3')), coloured by intgroup
plot_dynamic_pca <- function(vsd,
                             plot_title = "",
                             title_size = 16,
                             intgroup = "Condition",
                             pcs_to_plot = c("PC1", "PC2"),
                             ntop = 500) {
  if (length(pcs_to_plot) != 2) {
    stop("pcs_to_plot must be two PCs, e.g. c('PC1', 'PC3')", call. = FALSE)
  }
  pcs_numeric <- as.numeric(gsub("PC", "", pcs_to_plot))

  if (any(is.na(pcs_numeric)) || any(pcs_numeric < 1)) {
    stop("pcs_to_plot must be of the form 'PCn', n >= 1, e.g. c('PC1', 'PC3')", call. = FALSE)
  }
  pca_data <- plotPCA(object = vsd,
                      intgroup = intgroup,
                      returnData = TRUE,
                      ntop = ntop,
                      pcsToUse = pcs_numeric)
  percent_var <- round(100 * attr(pca_data, "percentVar"), 1)
  percent_var_x <- percent_var[1]
  percent_var_y <- percent_var[2]

  p <- ggplot(pca_data, aes(x = !!sym(pcs_to_plot[1]),
                            y = !!sym(pcs_to_plot[2]),
                            color = !!sym(intgroup),
                            shape = !!sym(intgroup),
                            label = name)) +
    geom_point(size = 2.5) +
    labs(
      x = paste0(pcs_to_plot[1], ": ", percent_var_x, "% variance"),
      y = paste0(pcs_to_plot[2], ": ", percent_var_y, "% variance"),
      title = plot_title
    ) +
    coord_fixed() +
    # reproducible labels: fixed seed, time for all iterations
    geom_text_repel(show.legend = FALSE, seed = 1, max.time = 15) +
    scale_color_manual(values = c("case" = "#3F84F5", control = "#EC8E12")) +
    theme_minimal() +
    theme(plot.title = element_text(size = title_size, face = "bold"),
          legend.text = element_text(size = 10),
          legend.title = element_blank(),
          legend.background = element_rect(fill = "transparent", color = NA),
          legend.position = "bottom",
          panel.background = element_rect(fill = "white"),
          plot.background = element_rect(fill = "white", color = "white"),
          panel.border = element_rect(colour = "darkgrey", fill = NA)
    )
  p
}

# ── Gene table, TPM, QC parsers ───────────────────────────────────────────────

# GENCODE GTF -> one row per gene: ID (unversioned), Transcript.ID (versioned
# gene_id), Type, Name, Length (end - start, not the exon length)
get_gencode_gene_table <- function(gtf_file) {
  g <- as.data.frame(rtracklayer::import(gtf_file, format = "gtf", feature.type = "gene"))
  data.frame(ID            = sub("[._].*$", "", g$gene_id),
             Transcript.ID = g$gene_id,
             Type          = g$gene_type,
             Name          = g$gene_name,
             Length        = g$end - g$start)
}

# TPM from counts (rows: versioned Ensembl IDs) and merged exon lengths (bp);
# genes not in gene_table dropped, first ID per symbol kept
# Out: genes x samples, rownames = Name
calc_tpm <- function(counts, lengths, gene_table) {
  stopifnot(length(lengths) == nrow(counts))
  ids <- sub("[._].*$", "", rownames(counts))
  o   <- order(ids)
  rpk <- counts[o, , drop = FALSE] / (lengths[o] / 1000)
  tpm <- apply(rpk, 2, function(x) x / sum(x) * 1e6)
  nm   <- gene_table$Name[match(ids[o], gene_table$ID)]
  keep <- !is.na(nm) & !duplicated(nm)
  tpm  <- tpm[keep, , drop = FALSE]
  rownames(tpm) <- nm[keep]
  tpm
}

# featureCounts stat -> per sample: Total of the categories and % of each;
# stops if another category is non-zero
rsubread_stats <- function(stat,
                           categories = c("Assigned", "Unassigned_Unmapped",
                                          "Unassigned_NoFeatures", "Unassigned_Ambiguity")) {
  st <- t(as.matrix(stat[, -1]))
  colnames(st) <- stat$Status
  rownames(st) <- sub("_RSR[.]bam$", "", rownames(st))
  other <- setdiff(colnames(st), categories)
  if (any(st[, other] != 0))
    stop("Non-zero featureCounts categories not in 'categories': ",
         paste(other[colSums(st[, other, drop = FALSE]) != 0], collapse = ", "), call. = FALSE)
  total <- rowSums(st[, categories])
  out <- data.frame(Name = rownames(st), Total = total,
                    round(st[, categories] * 100 / total, 3),
                    check.names = FALSE, row.names = NULL)
  out <- out[order(out$Name), ]
  rownames(out) <- NULL
  out
}

# <sample>_{1,2}_fastqc.zip -> per sample: mean Phred score and %GC of R1, R2
parse_fastqc <- function(fastqc_dir) {
  files <- list.files(fastqc_dir, pattern = "_[12]_fastqc\\.zip$", ignore.case = TRUE, full.names = TRUE)
  if (!length(files)) stop("No *_1_fastqc.zip / *_2_fastqc.zip in ", fastqc_dir, call. = FALSE)
  one <- function(f) {
    qc <- suppressMessages(fastqcr::qc_read(f))
    bs <- qc$basic_statistics
    data.frame(Name    = sub("_[12]_fastqc\\.zip$", "", basename(f), ignore.case = TRUE),
               Read    = paste0("R", sub("^.*_([12])_fastqc\\.zip$", "\\1", basename(f), ignore.case = TRUE)),
               Quality = round(mean(qc$per_base_sequence_quality$Mean), 2),
               GC      = as.numeric(bs$Value[bs$Measure == "%GC"]))
  }
  qc <- do.call(rbind, lapply(files, one))
  r1 <- qc[qc$Read == "R1", ]
  r2 <- qc[qc$Read == "R2", ]
  if (!setequal(r1$Name, r2$Name))
    stop("Samples without both R1 and R2: ",
         paste(setdiff(union(r1$Name, r2$Name), intersect(r1$Name, r2$Name)), collapse = ", "), call. = FALSE)
  out <- merge(r1[, c("Name", "Quality", "GC")], r2[, c("Name", "Quality", "GC")], by = "Name")
  colnames(out) <- c("Name", "Mean Phred Quality Score R1", "GC Content R1",
                     "Mean Phred Quality Score R2", "GC Content R2")
  out
}

# <sample>.metrics.tsv (RNA-SeQC) -> per sample: the chosen metrics
parse_rnaseqc <- function(rnaseqc_dir,
                          metrics = c("High Quality Rate", "Exonic Rate", "Genes Detected")) {
  files <- list.files(rnaseqc_dir, pattern = "\\.metrics\\.tsv$", ignore.case = TRUE, full.names = TRUE)
  if (!length(files)) stop("No *.metrics.tsv in ", rnaseqc_dir, call. = FALSE)
  one <- function(f) {
    m <- data.table::fread(f, sep = "\t", header = TRUE, data.table = FALSE, check.names = FALSE)
    missing <- setdiff(metrics, m[[1]])
    if (length(missing)) stop(basename(f), ": no ", paste(missing, collapse = ", "), call. = FALSE)
    v <- as.numeric(m[[2]][match(metrics, m[[1]])])
    data.frame(Name = colnames(m)[2], t(setNames(v, metrics)), check.names = FALSE)
  }
  out <- do.call(rbind, lapply(files, one))
  out <- out[order(out$Name), ]
  rownames(out) <- NULL
  out
}

# ── DEG annotation: gene table, Entrez IDs, Open Targets ──────────────────────
# In:  res_df (DESeq2 results, GeneID versioned); genes kept if in gene_table and
#      of gene_types; ot_files named by disease (order breaks ties)
# Out: list(res: + Entrez.ID, OT.globalScore, OT.gwasCredibleSets (max over
#      diseases, 0 if absent); ot_genes: top Open Targets row per tested symbol)
get_dge_ot <- function(res_df,
                       gene_table,
                       gene_types = c("protein_coding", "lncRNA"),
                       ot_files,
                       OrgDb = org.Hs.eg.db) {
  res_df <- na.omit(res_df)
  res_df <- dplyr::mutate(res_df, ID = sub("[._].*$", "", GeneID), .before = baseMean)
  if (!is.null(gene_table)) {
    res_df <- na.omit(merge(res_df, gene_table, by = "ID"))
    if (!is.null(gene_types))
      res_df <- res_df[res_df$Type %in% gene_types, ]
  }
  res_df <- res_df[order(res_df$padj), ]
  stopifnot(!anyDuplicated(res_df$ID))

  # Entrez IDs: by Ensembl ID, else by symbol; "0" if neither maps
  res_df$Entrez.ID <- dplyr::coalesce(AnnotationDbi::mapIds(OrgDb, res_df$ID,   "ENTREZID", "ENSEMBL"),
                                      AnnotationDbi::mapIds(OrgDb, res_df$Name, "ENTREZID", "SYMBOL"))
  res_df[is.na(res_df)] <- 0

  ot <- do.call(rbind, lapply(names(ot_files), function(k)
    transform(data.table::fread(ot_files[[k]], data.table = FALSE), key = k)))
  key_prevalence <- dplyr::count(ot, key)   # genes per disease
  # Per symbol: highest globalScore; ties -> disease with the most genes
  ot_genes <- ot %>%
    dplyr::group_by(symbol) %>%
    dplyr::filter(globalScore == max(globalScore)) %>%
    dplyr::left_join(key_prevalence, by = "key") %>%
    dplyr::arrange(symbol, dplyr::desc(n)) %>%
    dplyr::slice(1) %>%
    dplyr::ungroup()
  ot_genes <- ot_genes[ot_genes$symbol %in% res_df$Name, ]
  # Per symbol: highest GWAS credible-sets score; 'No data' = 0
  ot_gwas <- ot %>%
    dplyr::mutate(gwasCredibleSets = suppressWarnings(as.numeric(gwasCredibleSets))) %>%
    dplyr::group_by(symbol) %>%
    dplyr::summarise(gwasCredibleSets = suppressWarnings(max(gwasCredibleSets, na.rm = TRUE))) %>%
    dplyr::mutate(gwasCredibleSets = ifelse(is.finite(gwasCredibleSets), gwasCredibleSets, 0))

  res_df <- merge(res_df, ot_genes[, c("symbol", "globalScore")], by.x = "Name", by.y = "symbol", all.x = TRUE)
  res_df <- merge(res_df, ot_gwas, by.x = "Name", by.y = "symbol", all.x = TRUE)
  res_df[is.na(res_df)] <- 0
  res_df <- dplyr::rename(res_df, OT.globalScore = globalScore, OT.gwasCredibleSets = gwasCredibleSets)

  list(res = res_df[, c("ID", "Transcript.ID", "Entrez.ID", "Name", "Type",
                        "baseMean", "logFC", "lfcSE", "stat", "pvalue", "padj",
                        "OT.globalScore", "OT.gwasCredibleSets")],
       ot_genes = ot_genes)
}


# ── Volcano plot (F1) ─────────────────────────────────────────────────────────
# x: padj (reversed log10), y: log2FC; ot_genes larger, opaque. Labelled: logFC
# outside label_lfc, padj < label_padj or OT score >= label_ot. Points beyond
# the caps: triangles at the cap, true values in the label
# file: without extension (NULL = not saved); size in px
plot_volcano <- function(degs, ot_genes, file = NULL,
                         label_lfc = c(-2, 1.5), label_padj = 5e-05, label_ot = 0.1,
                         padj_cap = 1e-6, lfc_cap = 3,
                         width = 8000, height = 6000, dpi = 800) {
  in_ot <- degs$Name %in% ot_genes$symbol
  goi   <- with(degs, Name[logFC <= label_lfc[1] | logFC > label_lfc[2] |
                           padj < label_padj | OT.globalScore >= label_ot])
  padj_capped <- degs$padj < padj_cap
  lfc_capped  <- abs(degs$logFC) > lfc_cap

  # name + true values of capped coordinates
  cap_note <- paste0(ifelse(padj_capped, paste0("p=", formatC(degs$padj, format = "e", digits = 1)), ""),
                     ifelse(padj_capped & lfc_capped, ", ", ""),
                     ifelse(lfc_capped, paste0("LFC=", formatC(degs$logFC, format = "f", digits = 1)), ""))
  cap_note <- ifelse(padj_capped | lfc_capped, paste0("\n(", cap_note, ")"), "")

  df <- data.frame(plot_padj  = pmax(degs$padj, padj_cap),
                   plot_lfc = pmax(pmin(degs$logFC, lfc_cap), -lfc_cap),
                   plot_label = ifelse(degs$Name %in% goi, paste0(degs$Name, cap_note), ""),
                   pt_color = as.factor(ifelse(degs$logFC >= 0, "#EC8E12", "#3F84F5")),
                   pt_alpha = ifelse(in_ot, 0.75, 0.25),
                   pt_size  = ifelse(in_ot, 2, 1.3),
                   pt_shape = ifelse(padj_capped | lfc_capped, 17, 16))

  p <- ggplot(df, aes(x = plot_padj, y = plot_lfc, label = plot_label)) +
    scale_x_continuous(trans = c("log10", "reverse"),
                       breaks = c(5e-02, 5e-03, 5e-04, 5e-05, 5e-06)) +
    ylab(expression("Log"["2"] ~ "fold-change")) + xlab("Adjusted p-value") +
    geom_point(color = df$pt_color, alpha = df$pt_alpha, size = df$pt_size, shape = df$pt_shape) +
    geom_hline(yintercept = 0, color = "gray50", alpha = 0.5, linewidth = 0.25) +
    # all labels, with leader lines; reproducible
    ggrepel::geom_text_repel(show.legend = FALSE, max.overlaps = Inf, min.segment.length = 0,
                             seed = 1, max.time = 15, size = 4.2, lineheight = 0.9, fontface = "italic",
                             segment.color = "grey50", segment.size = 0.3) +
    theme_classic(base_size = 18) +
    theme(panel.grid = element_blank(),
          panel.background = element_rect(fill = "white"),
          plot.background  = element_rect(fill = "white", color = "white"),
          axis.line  = element_line(colour = "black", linewidth = 0.5),
          axis.ticks = element_line(colour = "black", linewidth = 0.5),
          axis.title = element_text(size = 20),
          axis.text  = element_text(size = 16, colour = "black"),
          legend.position = "none")

  if (!is.null(file))
    save_ggplot_png_pdf(p, file, width = width, height = height, dpi = dpi)
  invisible(p)
}


# ── ORA: GO + KEGG, network figure (F2) ───────────────────────────────────────

# ── Part 1. Enrichment ────────────────────────────────────────────────────────

# GO BP and KEGG ORA of up / down DEGs, merged per direction
# In:  degs, background (ID, Name, Entrez.ID ("0" = none), logFC)
# GO: Count >= go_min_count, simplified at go_simplify_cutoff (go_plot: at
# go_plot_simplify_cutoff); KEGG: excluded (sub)categories dropped, Count >=
# kegg_min_count; top_categories: top_n[[s]] per source s (smallest padj), no repeats
# Out: list(lfc, go, go_plot, kegg, combined, combined_relabel, top_categories,
#      top_n, tables), per direction
exec_ora_go_kegg <- function(degs,
                             background,
                             OrgDb = "org.Hs.eg.db",
                             go_ont = "BP",
                             kegg_organism = "hsa",
                             min_gs_size = 20,
                             p_adjust_method = "BH",
                             pvalue_cutoff = 0.05,
                             go_qvalue_cutoff = 0.05,
                             go_min_count = 16,
                             go_simplify_cutoff = 0.7,
                             go_plot_simplify_cutoff = 0.35,
                             kegg_min_count = 8,
                             kegg_exclude_category = "Human Diseases",
                             kegg_exclude_subcategory = "Digestive system",
                             top_n = c(GO = 5, KEGG = 5)) {
  need <- c("ID", "Name", "Entrez.ID", "logFC")
  for (tab in c("degs", "background")) {
    miss <- setdiff(need, names(get(tab)))
    if (length(miss)) stop(tab, " lacks column(s): ", paste(miss, collapse = ", "))
  }
  orgdb <- if (is.character(OrgDb)) getExportedValue(OrgDb, OrgDb) else OrgDb

  degs_sorted <- degs[order(degs$logFC, decreasing = TRUE), ]
  lfc <- list(ensembl = setNames(as.numeric(degs_sorted$logFC), degs_sorted$ID),
              symbol  = setNames(as.numeric(degs_sorted$logFC), degs_sorted$Name))
  degs_sorted_entrez <- degs_sorted[degs_sorted$Entrez.ID != "0", ]
  lfc$entrez <- setNames(as.numeric(degs_sorted_entrez$logFC), degs_sorted_entrez$Entrez.ID)

  dirs  <- c(up = "up", down = "down")
  genes <- function(v, dir) names(v[if (dir == "up") v > 0 else v < 0])
  # rows kept by `keep` (NA = FALSE, as in subset())
  keep_rows <- function(r, keep) r[keep & !is.na(keep), , drop = FALSE]

  go <- lapply(dirs, function(dir) {
    x <- enrichGO(
      gene = genes(lfc$ensembl, dir),
      universe = background$ID,
      keyType = "ENSEMBL",
      OrgDb = orgdb,
      ont = go_ont,
      minGSSize = min_gs_size,
      pAdjustMethod = p_adjust_method,
      pvalueCutoff  = pvalue_cutoff,
      qvalueCutoff = go_qvalue_cutoff,
      readable = TRUE)
    x@result <- keep_rows(x@result, x@result$Count >= go_min_count)
    x@result$Description <- tools::toTitleCase(x@result$Description)
    x <- pairwise_termsim(x)
    clusterProfiler::simplify(x, cutoff = go_simplify_cutoff, by = "p.adjust")
  })
  go_plot <- lapply(go, function(x)
    clusterProfiler::simplify(x, cutoff = go_plot_simplify_cutoff, by = "p.adjust")[seq_len(top_n[["GO"]]), ])

  kegg <- lapply(dirs, function(dir) {
    x <- enrichKEGG(
      gene = genes(lfc$entrez, dir),
      universe = background$Entrez.ID,
      organism = kegg_organism,
      pvalueCutoff = pvalue_cutoff,
      pAdjustMethod = p_adjust_method,
      minGSSize = min_gs_size)
    keep <- rep(TRUE, nrow(x@result))
    for (v in kegg_exclude_category) keep <- keep & x@result$category != v
    for (v in kegg_exclude_subcategory) keep <- keep & x@result$subcategory != v
    x@result <- keep_rows(x@result, keep)
    x@result <- keep_rows(x@result, x@result$Count >= kegg_min_count)
    x <- pairwise_termsim(x)
    x@result$Description <- tools::toTitleCase(x@result$Description)
    setReadable(x, OrgDb, "ENTREZID")
  })

  combined <- lapply(dirs, function(dir) {
    x <- merge_result(setNames(list(go[[dir]], kegg[[dir]]), paste0(dir, c(".GO", ".KEGG"))))
    x@compareClusterResult$foldChange <- list(lfc$symbol)
    x
  })
  combined_relabel <- lapply(dirs, function(dir) {
    x <- combined[[dir]]
    x@compareClusterResult$Cluster <- factor(
      sub(paste0("^", dir, "\\."), "", as.character(x@compareClusterResult$Cluster)),
      levels = c("GO", "KEGG"))
    x
  })

  top_categories <- lapply(dirs, function(dir) {
    results <- list(GO = go_plot[[dir]], KEGG = kegg[[dir]])
    picked  <- character()
    for (s in names(top_n)) {
      r   <- as.data.frame(results[[s]])
      r   <- r[!r$Description %in% picked, ]
      top <- head(r$Description[order(r$p.adjust)], top_n[[s]])
      if (length(top) < top_n[[s]])
        warning(sprintf("ORA top categories: only %d significant %s categories (asked for %d)",
                        length(top), s, top_n[[s]]), call. = FALSE)
      picked <- c(picked, top)
    }
    picked
  })

  tables <- lapply(combined, function(x)
    data.frame(x) %>% dplyr::select(-foldChange, -category, -subcategory, -Count))

  list(lfc = lfc, go = go, go_plot = go_plot, kegg = kegg,
       combined = combined, combined_relabel = combined_relabel,
       top_categories = top_categories, top_n = top_n, tables = tables)
}

# ── Part 2. Network figure ────────────────────────────────────────────────────
# 1. exec_ora_network(): draft (kk layout, ggrepel labels) + ORA_layout_data.js
# 2. ora_layout_editor.html: arrange, save ORA_layout_<fig_id>.json
# 3. exec_ora_network() again: final figure as arranged

# Pie scaling of the panel being drawn (read by the patched enrichplot code)
.ora_state <- new.env()
ORA_ENRICHPLOT_TESTED <- "1.28.4"

# Figure settings (A = up, B = down; from exec_ora_go_kegg()). max_genes per
# category (largest |log2FC|); pie radius = pie_scale * n_DEGs^pie_exp * width
# Patches enrichplot:::add_node_pie() (translucent pies, opaque ring, no size
# key); repeatable, undone by ora_unpatch_enrichplot(); stops if its source changed
ora_network_spec <- function(ora = NULL,
                             up = ora$combined_relabel$up,
                             down = ora$combined_relabel$down,
                             fold_change = ora$lfc$symbol,
                             top_categories = ora$top_categories,
                             bold_genes = character(),
                             fig_id = paste0("top_", paste(ora$top_n, collapse = "_")),
                             colors = c(GO = "#1B9E77", KEGG = "#E7298A"),
                             fc_names = list(up   = expression("A. Log"["2"] ~ "fold-change"),
                                             down = expression("B. Log"["2"] ~ "fold-change")),
                             layout_dir = ".",
                             editor_data = "ORA_layout_data.js",
                             fig_px = c(width = 8000, height = 5000),
                             dpi = 600,
                             max_genes = 10,
                             pie_scale = 0.012,
                             pie_exp = 0.3) {
  # ── enrichplot pie patch
  if (packageVersion("enrichplot") != ORA_ENRICHPLOT_TESTED)
    warning("ORA pie patch tested with enrichplot ", ORA_ENRICHPLOT_TESTED, ", installed ",
            packageVersion("enrichplot"), ": check the ORA network figure", call. = FALSE)
  current <- enrichplot:::add_node_pie
  orig    <- attr(current, "ora_original") %||% current
  sub_once <- function(pattern, replacement, x) {
    if (!any(grepl(pattern, x, fixed = TRUE)))
      stop('ora_network_spec(): "', pattern, '" not found in enrichplot:::add_node_pie ',
           "(enrichplot ", packageVersion("enrichplot"), ", patch tested with ",
           ORA_ENRICHPLOT_TESTED, ")", call. = FALSE)
    sub(pattern, replacement, x, fixed = TRUE)
  }
  pie_data_txt <- paste0(
    "transform(merge(p$data[, c(\"x\", \"y\", \"label\")], dd, by.x = \"label\", by.y = \"Description\"), ",
    ".pie_r = pathway_size * pie_scale)")
  txt <- sub_once(
    "p <- p %<+% dd + scatterpie::geom_scatterpie(",
    paste0("p <- p %<+% dd + ",
           "ggforce::geom_circle(data = ", pie_data_txt, ", aes(x0 = x, y0 = y, r = .pie_r), ",
           "fill = \"white\", color = NA, inherit.aes = FALSE) + ",
           "scatterpie::geom_scatterpie("),
    deparse(orig))
  txt <- sub_once(
    "legend_name = \"Cluster\", color = NA)",
    paste0("legend_name = \"Cluster\", color = NA, alpha = 0.75) + ",
           "scatterpie::geom_scatterpie(data = ", pie_data_txt, ", aes(x = x, y = y, r = .pie_r), ",
           "cols = as.character(unique(d$Cluster)), legend_name = \"Cluster\", ",
           "color = NA, donut_radius = 0.9)"),
    txt)
  txt <- sub_once("scatterpie::geom_scatterpie_legend(", "(function(...) NULL)(", txt)
  txt <- sub_once("pathway_size[dd$Description]/sum(pathway_size)",
                  "ora_pie_radius(pathway_size[dd$Description])", txt)
  patched <- eval(parse(text = txt))
  env <- new.env(parent = environment(orig))   # enrichplot's namespace, plus ora_pie_radius
  env$ora_pie_radius <- function(n) .ora_state$pie_scale * n^.ora_state$pie_exp * .ora_state$span
  environment(patched) <- env
  attr(patched, "ora_original") <- orig
  assignInNamespace("add_node_pie", patched, ns = "enrichplot")

  list(panels = list(A = list(ora = up,   direction = "up",   tag = "A", fc_name = fc_names$up),
                     B = list(ora = down, direction = "down", tag = "B", fc_name = fc_names$down)),
       fold_change = fold_change, top_categories = top_categories, bold_genes = bold_genes,
       fig_id = fig_id, colors = colors, layout_dir = layout_dir, editor_data = editor_data,
       fig_px = fig_px, dpi = dpi, max_genes = max_genes,
       pie_scale = pie_scale, pie_exp = pie_exp)
}

# Undo the enrichplot patch of ora_network_spec()
ora_unpatch_enrichplot <- function() {
  orig <- attr(enrichplot:::add_node_pie, "ora_original")
  if (!is.null(orig)) assignInNamespace("add_node_pie", orig, ns = "enrichplot")
  invisible(orig)
}

# Draft figure + editor data; final figure if the layout file exists
# categories: list(up, down); NULL = top categories (use a new fig_id)
# relabel = TRUE: new label positions instead of the saved ones
exec_ora_network <- function(spec, out_dir = ".", prefix = "ORA_integrated_",
                             categories = NULL, fig_id = spec$fig_id, relabel = FALSE,
                             draft_name = paste0(prefix, fig_id),
                             final_name = paste0(prefix, fig_id, "_final")) {
  layout_file <- paste0("ORA_layout_", fig_id, ".json")
  if (!identical(spec$layout_dir, ".")) layout_file <- file.path(spec$layout_dir, layout_file)
  panel_names <- setNames(names(spec$panels), names(spec$panels))

  # One panel's saved layout (frame, nodes, label offsets); NULL if none
  read_layout <- function(nm, layout) {
    if (is.list(layout)) return(layout[[nm]])
    if (is.null(layout) || !file.exists(layout)) return(NULL)
    jsonlite::fromJSON(layout)[[nm]]
  }

  # Long category names -> up to 3 balanced lines (~25 characters each)
  wrap_category <- function(x, chars_per_line = 25, max_lines = 3) {
    vapply(x, function(s) {
      words <- strsplit(s, " ", fixed = TRUE)[[1]]
      n_lines <- min(ceiling(nchar(s) / chars_per_line), max_lines, length(words))
      if (n_lines <= 1) return(s)
      breaks <- combn(length(words) - 1, n_lines - 1, simplify = FALSE)
      split_at <- function(b) {
        grp <- rep(seq_len(n_lines), diff(c(0, b, length(words))))
        vapply(split(words, grp), paste, character(1), collapse = " ")
      }
      longest <- vapply(breaks, function(b) max(nchar(split_at(b))), numeric(1))
      paste(split_at(breaks[[which.min(longest)]]), collapse = "\n")
    }, character(1), USE.NAMES = FALSE)
  }

  # Categories per panel; those no longer enriched dropped with a warning
  ora_categories <- function(categories = NULL) {
    cats <- lapply(spec$panels, function(pnl)
      categories[[pnl$direction]] %||% spec$top_categories[[pnl$direction]])
    cats <- lapply(cats, function(x) x[!is.na(x)])   # NA: fewer top categories than asked for
    enriched <- lapply(spec$panels, function(pnl) pnl$ora@compareClusterResult$Description)
    words <- function(x) setdiff(strsplit(tolower(x), "[^a-z0-9]+")[[1]], c("", "of", "and", "by", "to", "the", "in"))
    for (nm in names(cats)) {
      miss <- setdiff(cats[[nm]], enriched[[nm]])
      if (length(miss)) {
        similar <- unique(unlist(lapply(miss, function(m)
          Filter(function(e) length(intersect(words(m), words(e))) >= 2, enriched[[nm]]))))
        warning(sprintf("Panel %s (%s): not among the enriched categories, left out: %s.%s",
                        nm, spec$panels[[nm]]$direction, paste(miss, collapse = ", "),
                        if (length(similar)) paste0(" Similar enriched: ", paste(similar, collapse = ", ")) else ""),
                call. = FALSE)
        cats[[nm]] <- intersect(cats[[nm]], enriched[[nm]])
      }
    }
    cats
  }

  # Both panels. labels: 'repel' (draft), 'fixed' (saved offsets), 'none'
  # (measuring); frames: panel limits; fill_area: stretch to the space (measuring)
  ora_panels <- function(labels, layout, cats, frames = NULL, fill_area = FALSE) {
    labels <- match.arg(labels, c("repel", "fixed", "none"))
    # Padding keeps other labels off category labels
    pad_category <- function(x) {
      x <- gsub("\n", "  \n  ", x, fixed = TRUE)
      paste0(" \n  ", x, "  \n ")
    }
    gene_nodes <- function(d) d[!d$.isCategory, ]
    # Top n genes per category (|log2FC|); pie size still = all DEGs
    top_genes_per_category <- function(ora, n) {
      r <- ora@compareClusterResult
      genes <- strsplit(r$geneID, "/", fixed = TRUE)
      keep <- lapply(split(genes, r$Description), function(g) {
        g <- unique(unlist(g))
        head(g[order(-abs(spec$fold_change[g]))], n)
      })
      r$geneID <- mapply(function(g, d) paste(intersect(g, keep[[d]]), collapse = "/"),
                         genes, r$Description)
      ora@compareClusterResult <- r
      ora
    }
    # kk layout, with saved positions where available
    kk_layout <- function(saved) {
      function(graph, ...) {
        xy <- igraph::layout_with_kk(graph, ...)
        .ora_state$span <- max(apply(xy, 2, function(v) diff(range(v))))
        if (!is.null(saved)) {
          i <- match(igraph::V(graph)$name, saved$nodes$name)
          if (anyNA(i))
            warning("Not in the saved layout, left at the kk position: ",
                    paste(igraph::V(graph)$name[is.na(i)], collapse = ", "), call. = FALSE)
          xy[!is.na(i), ] <- as.matrix(saved$nodes[i[!is.na(i)], c("x", "y")])
        }
        xy
      }
    }

    lapply(panel_names, function(nm) {
      saved <- read_layout(nm, layout)
      if (labels == "fixed" && is.null(saved))
        stop('labels = "fixed" needs a saved layout')
      ora         <- spec$panels[[nm]]$ora
      fold_change <- spec$fold_change
      bold_genes  <- spec$bold_genes
      categories  <- cats[[nm]]
      missing <- setdiff(categories, ora@compareClusterResult$Description)
      if (length(missing))
        stop("Not among the enriched categories: ", paste(missing, collapse = ", "), call. = FALSE)
      .ora_state$pie_scale <- spec$pie_scale
      .ora_state$pie_exp   <- spec$pie_exp
      p <- cnetplot(top_genes_per_category(ora, spec$max_genes),
                    showCategory = categories,
                    foldChange = fold_change,
                    layout = kk_layout(saved),
                    size_category = 1,   # pie radius set by the patched pie code
                    size_item = 0.6,
                    size_edge = 0.1,
                    node_label = "none") +
        # panel letter, top-left
        labs(tag = spec$panels[[nm]]$tag) +
        # gene dots: fold-change fill, up / down border
        geom_point(data = gene_nodes, aes(x = x, y = y, color = foldChange),
                   shape = 16, size = 2.5, alpha = 0.75) +
        geom_point(data = gene_nodes,
                   aes(x = x, y = y, color = I(ifelse(foldChange > 0, "#EC8E12", "#2679ff"))),
                   shape = 21, fill = NA, size = 2.5, stroke = 0.5, alpha = 0.75, show.legend = FALSE) +
        theme(plot.tag = element_text(size = 32, face = "bold", hjust = 0, vjust = 1),
              plot.tag.position = c(0, 1),
              plot.background = element_rect(fill = "white", colour = "white")) +
        scale_color_gradient2(low = "#2679ff", mid = "white", high = "#EC8E12", midpoint = 0,
                              name = spec$panels[[nm]]$fc_name) +
        scale_fill_manual(values = spec$colors, name = "Source")

      # one label layer, so gene and category labels repel each other
      label_aes <- function(label) aes(label    = {{ label }},
                                       size     = I(ifelse(.isCategory, 3.75, 2)),
                                       alpha    = I(ifelse(.isCategory, 0.9, 0.5)),
                                       fontface = ifelse(!.isCategory & label %in% bold_genes, "bold", "plain"))
      if (labels == "repel") {
        p <- p + geom_cnet_label(mapping = label_aes(ifelse(.isCategory, pad_category(wrap_category(label)), label)),
                                 node_label = "all", color = "black", lineheight = 0.9, show.legend = FALSE,
                                 max.overlaps = Inf,               # never drop labels
                                 segment.size = 0.2,               # thin connector lines
                                 max.time = 5, max.iter = 1e5,     # let the layout converge
                                 seed = 1)
      } else if (labels == "fixed") {
        d <- p$data
        i <- match(d$name, saved$nodes$name)
        if (anyNA(i) || is.null(saved$nodes$lx))
          warning("Some label positions missing from the saved layout; those labels sit on their node",
                  call. = FALSE)
        d$lx <- if (is.null(saved$nodes$lx)) 0 else saved$nodes$lx[i]
        d$ly <- if (is.null(saved$nodes$ly)) 0 else saved$nodes$ly[i]
        d$lx[is.na(d$lx)] <- 0
        d$ly[is.na(d$ly)] <- 0
        p <- p + ggrepel::geom_text_repel(data = d, mapping = label_aes(ifelse(.isCategory, wrap_category(label), label)),
                                          position = ggrepel::position_nudge_repel(x = d$lx, y = d$ly),
                                          color = "black", lineheight = 0.9, show.legend = FALSE,
                                          max.overlaps = Inf, segment.size = 0.2,
                                          max.iter = 0, force = 0, force_pull = 0)
      }

      # 1:1 aspect keeps pies round; frames evaluated only here (may draw to measure)
      if (fill_area) {
        p <- suppressMessages(p + coord_cartesian(expand = FALSE))
      } else {
        frame <- frames[[nm]] %||% saved$frame
        if (!is.null(frame))
          p <- suppressMessages(p + coord_fixed(ratio = p$coordinates$ratio, expand = FALSE,
                                                xlim = frame$x, ylim = frame$y))
      }
      p
    })
  }

  # Panels side by side over a one-row legend; panel_theme for measuring
  ora_figure <- function(panels, panel_theme = NULL) {
    divider <- ggplot() + theme_void() +
      theme(plot.background = element_rect(fill = "grey30", colour = NA))
    legend_source  <- ggpubr::get_legend(panels$A + guides(color = "none") + theme(legend.position = "bottom"))
    legend_fc_up   <- ggpubr::get_legend(panels$A + guides(fill = "none") + theme(legend.position = "bottom"))
    legend_fc_down <- ggpubr::get_legend(panels$B + guides(fill = "none") + theme(legend.position = "bottom"))
    legend_row  <- ggpubr::as_ggplot(gridExtra::arrangeGrob(
      legend_fc_up, legend_source, legend_fc_down, ncol = 3))

    p_panels <- (panels$A + theme(legend.position = "none") + panel_theme) +
      divider +
      (panels$B + theme(legend.position = "none") + panel_theme) +
      plot_layout(nrow = 1, widths = c(1, 0.001, 1)) # divider
    p_panels / legend_row + plot_layout(heights = c(8, 0.8))
  }

  # Panel sizes (in) of the saved figure, from magenta panels in a low-res PNG,
  # plus ggrepel's label positions
  measure_ora_figure <- function(panels, dpi = 50) {
    # ggrepel label positions while `expr` draws (hooks ggrepel:::repel_boxes2, tested 0.9.8)
    record_repel <- function(expr) {
      if (!exists("repel_boxes2", envir = asNamespace("ggrepel"), inherits = FALSE))
        stop("record_repel(): ggrepel:::repel_boxes2 not found (ggrepel ", packageVersion("ggrepel"),
             ", tested with 0.9.8)", call. = FALSE)
      orig  <- ggrepel:::repel_boxes2
      calls <- list()
      assignInNamespace("repel_boxes2", function(data_points, ...) {
        res <- orig(data_points, ...)
        calls[[length(calls) + 1]] <<- list(points = data_points, x = res$x, y = res$y)
        res
      }, ns = "ggrepel")
      on.exit(assignInNamespace("repel_boxes2", orig, ns = "ggrepel"))
      force(expr)
      calls
    }

    f <- tempfile(fileext = ".png")
    on.exit(unlink(f))
    fig <- ora_figure(panels, theme(panel.background = element_rect(fill = "#FF00FF", colour = NA)))
    repel <- record_repel(
      ggsave(f, fig, device = "png", dpi = dpi, units = "in",
             width = spec$fig_px[["width"]] / spec$dpi, height = spec$fig_px[["height"]] / spec$dpi))
    img <- png::readPNG(f)
    magenta <- img[, , 1] > 0.9 & img[, , 2] < 0.1 & img[, , 3] > 0.9
    extent <- function(cols) {
      rows <- which(rowSums(magenta[, cols, drop = FALSE]) > 0)
      hit  <- which(colSums(magenta[, cols, drop = FALSE]) > 0)
      c(width = diff(range(hit)) + 1, height = diff(range(rows)) + 1) / dpi
    }
    half <- ncol(magenta) %/% 2
    list(sizes = list(A = extent(1:half), B = extent((half + 1):ncol(magenta))), repel = repel)
  }

  # Panel limits widened (never narrowed) to the shape of each panel's space
  ora_frames <- function(layout, cats) {
    base <- lapply(ora_panels("none", layout, cats), function(p) {
      pp <- ggplot_build(p)$layout$panel_params[[1]]
      list(x = pp$x.range, y = pp$y.range)
    })
    area <- measure_ora_figure(ora_panels("none", layout, cats, fill_area = TRUE))$sizes
    Map(function(f, a) {
      target <- a[["height"]] / a[["width"]]
      w <- diff(f$x)
      h <- diff(f$y)
      if (h / w < target) f$y <- f$y + c(-1, 1) * (target * w - h) / 2
      else                f$x <- f$x + c(-1, 1) * (h / target - w) / 2
      f
    }, base, area[names(base)])
  }

  # Layout as drawn; label offsets saved, or from ggrepel with relabel = TRUE
  ora_current_layout <- function(panels, cats, measured) {
    # Label offsets (data units) from the recorded ggrepel calls, matched to nodes by position
    repel_offsets <- function(repel, d, xr, yr) {
      px <- (d$x - xr[1]) / diff(xr)
      py <- (d$y - yr[1]) / diff(yr)
      fit <- lapply(repel, function(r) {
        if (nrow(r$points) != nrow(d)) return(NULL)
        i <- vapply(seq_len(nrow(r$points)), function(k)
          which.min((px - r$points[k, 1])^2 + (py - r$points[k, 2])^2), integer(1))
        list(i = i, err = sum((px[i] - r$points[, 1])^2 + (py[i] - r$points[, 2])^2), r = r)
      })
      fit <- Filter(Negate(is.null), fit)
      off <- data.frame(lx = rep(0, nrow(d)), ly = rep(0, nrow(d)))
      if (!length(fit)) return(off)
      best <- fit[[which.min(vapply(fit, `[[`, numeric(1), "err"))]]
      off$lx[best$i] <- (best$r$x - best$r$points[, 1]) * diff(xr)
      off$ly[best$i] <- (best$r$y - best$r$points[, 2]) * diff(yr)
      off
    }

    lapply(setNames(names(panels), names(panels)), function(nm) {
      d     <- panels[[nm]]$data
      pp    <- ggplot_build(panels[[nm]])$layout$panel_params[[1]]
      off   <- repel_offsets(measured$repel, d, pp$x.range, pp$y.range)
      saved <- read_layout(nm, layout_file)
      if (!relabel && !is.null(saved$nodes$lx)) {
        i <- match(d$name, saved$nodes$name)
        off$lx <- ifelse(is.na(i), off$lx, saved$nodes$lx[i])
        off$ly <- ifelse(is.na(i), off$ly, saved$nodes$ly[i])
      }
      list(categories = cats[[nm]],
           frame = list(x = pp$x.range, y = pp$y.range),
           nodes = data.frame(name = d$name, x = d$x, y = d$y, lx = off$lx, ly = off$ly))
    })
  }

  # Adds / replaces this figure in the editor data file (one ORA_FIGS["<fig_id>"] line each)
  write_editor_data <- function(panels, layout, sizes, file = spec$editor_data) {
    out <- lapply(setNames(names(panels), names(panels)), function(nm) {
      d <- panels[[nm]]$data
      lay <- layout[[nm]]
      e <- igraph::as_edgelist(get("graph", envir = panels[[nm]]$plot_env))
      list(width_in = sizes[[nm]][["width"]], height_in = sizes[[nm]][["height"]],
           categories = lay$categories,
           frame = lay$frame,
           nodes = data.frame(
             lay$nodes,
             label    = ifelse(d$.isCategory, wrap_category(d$label), d$label),
             category = d$.isCategory,
             r        = ifelse(d$.isCategory, d$pathway_size, NA),   # pie radius, data units
             fc       = ifelse(d$.isCategory, NA, d$foldChange),
             bold     = !d$.isCategory & d$label %in% spec$bold_genes,
             d[, intersect(names(spec$colors), names(d)), drop = FALSE]),
           edges = cbind(match(e[, 1], d$name), match(e[, 2], d$name)) - 1)
    })
    json <- jsonlite::toJSON(list(created = format(Sys.time()), layout_file = basename(layout_file),
                                  colors = as.list(spec$colors), panels = out),
                             auto_unbox = TRUE, digits = NA, na = "null")
    key   <- sprintf('ORA_FIGS["%s"] = ', fig_id)
    lines <- if (file.exists(file)) readLines(file, warn = FALSE) else character()
    lines <- lines[startsWith(lines, 'ORA_FIGS["') & !startsWith(lines, key)]   # other figures
    writeLines(c("window.ORA_FIGS = window.ORA_FIGS || {};", lines, paste0(key, json, ";")), file)
  }

  save_figure <- function(fig, file)
    save_ggplot_png_pdf(fig, sub("\\.png$", "", file), width = spec$fig_px[["width"]],
                        height = spec$fig_px[["height"]], dpi = spec$dpi)

  # ── Draft
  draft    <- file.path(out_dir, paste0(draft_name, ".png"))
  cats     <- ora_categories(categories)
  panels   <- ora_panels("repel", layout_file, cats, ora_frames(layout_file, cats))
  measured <- measure_ora_figure(panels)
  if (length(measured$repel) < length(panels) ||
      !all(vapply(measured$repel, function(r) !is.null(r$x) && !is.null(r$y), logical(1))))
    stop("exec_ora_network(): ggrepel label positions not recorded as expected (ggrepel ",
         packageVersion("ggrepel"), ", tested with 0.9.8)", call. = FALSE)
  layout <- ora_current_layout(panels, cats, measured)
  write_editor_data(panels, layout, measured$sizes)
  save_figure(ora_figure(ora_panels("fixed", layout, cats)), draft)

  # ── Final, from the edited layout
  final <- NULL
  if (file.exists(layout_file)) {
    final <- file.path(out_dir, paste0(final_name, ".png"))
    saved <- jsonlite::fromJSON(layout_file)
    cats  <- ora_categories(list(up = saved$A$categories, down = saved$B$categories))
    save_figure(ora_figure(ora_panels("fixed", layout_file, cats, ora_frames(layout_file, cats))), final)
  } else {
    message("No ", layout_file, " yet: arrange the layout in ora_layout_editor.html, ",
            "save it there, then rerun for the final figure.")
  }
  invisible(list(draft = draft, final = final, layout_file = layout_file,
                 editor_data = spec$editor_data))
}


# ── Saving figures ────────────────────────────────────────────────────────────

# ggplot -> <file>.png / .pdf; width, height in px at dpi
save_ggplot_png_pdf <- function(plot, file, width, height, dpi) {
  ggsave(paste0(file, ".png"), plot, device = "png", dpi = dpi, units = "px",
         width = width, height = height)
  ggsave(paste0(file, ".pdf"), plot, device = cairo_pdf, dpi = dpi, units = "px",
         width = width, height = height)
}

# Grid drawing (draw_fun()) -> <file>.png / .pdf at the final size; retries while
# the PNG is locked (Windows viewer), then skips with a warning
save_grid_png_pdf <- function(file, draw_fun, width, height, dpi = 1200, tries = 5, wait = 2) {
  g <- grid.grabExpr(draw_fun(), width = width / dpi, height = height / dpi, wrap.grobs = TRUE)
  png_file <- paste0(file, ".png")
  for (i in seq_len(tries)) {
    err <- tryCatch({
      ggsave(png_file, g, device = "png", dpi = dpi, units = "px",
             width = width, height = height, bg = "white")
      NULL
    },
    # ragg only warns when the file is locked
    warning = function(w) conditionMessage(w),
    error = function(e) conditionMessage(e))
    if (is.null(err)) break
    if (i < tries) Sys.sleep(wait)
  }
  if (!is.null(err)) {
    warning(sprintf("Could not write %s after %d attempts - figure skipped (%s)", basename(png_file), tries, err),
            call. = FALSE)
    return(invisible(FALSE))
  }
  pdf_file <- paste0(file, ".pdf")
  err <- tryCatch({
    ggsave(pdf_file, g, device = cairo_pdf, dpi = dpi, units = "px",
           width = width, height = height, bg = "white")
    NULL
  },
  error = function(e) conditionMessage(e))
  if (!is.null(err))
    warning(sprintf("Could not write %s (%s)", basename(pdf_file), err), call. = FALSE)
  invisible(TRUE)
}


# ── TF activity: decoupleR ULM, TF network (F3) ───────────────────────────────
# ULM on the Wald statistic. A TF is drawn
# when padj < padj_cutoff and its DEG targets support its direction (sign(score)
# * mor * sign(logFC) > 0): >= min_consistent, fraction > min_frac,
# binomial P < consistency_p_cutoff;
# Out: list(full: all TFs; network: TFs drawn; edges: TF -> DEG; net; acts)
exec_tf_activity <- function(res, degs, mean_tpm, gene_table,
                             net = decoupleR::get_collectri(organism = "human", split_complexes = FALSE),
                             minsize = 15,
                             padj_cutoff = 0.01,
                             min_consistent = 5,
                             consistency_p_cutoff = 0.05,
                             min_frac = 0.75) {
  gtf_names <- gene_table$Name

  res_tf      <- res[order(-res$baseMean), ]
  res_tf      <- res_tf[!duplicated(res_tf$Name), ]
  stat_vector <- as.matrix(data.frame(stat = res_tf$stat, row.names = res_tf$Name))
  sample_acts <- decoupleR::run_ulm(mat = stat_vector, net = net,
                                    .source = "source", .target = "target", .mor = "mor",
                                    minsize = minsize) %>%
    dplyr::mutate(padj = p.adjust(p_value, method = "BH"))

  # Per TF: regulon size, DEG targets, DEG targets consistent with its activity
  tf_targets <- net %>%
    dplyr::filter(target %in% rownames(stat_vector)) %>%
    dplyr::left_join(dplyr::select(degs, Name, logFC), by = c("target" = "Name")) %>%
    dplyr::inner_join(dplyr::select(sample_acts, source, score), by = "source") %>%
    dplyr::mutate(is_DEG     = !is.na(logFC),
                  consistent = is_DEG & sign(score) * mor * sign(logFC) > 0) %>%
    dplyr::group_by(source) %>%
    dplyr::summarize(n_targets                = dplyr::n(),
                     n_DEG_targets            = sum(is_DEG),
                     n_consistent_DEG_targets = sum(consistent),
                     targets                  = dplyr::na_if(paste(target[is_DEG], collapse = "/"), ""),
                     consistent_targets       = dplyr::na_if(paste(target[consistent], collapse = "/"), ""),
                     .groups = "drop") %>%
    dplyr::mutate(frac_consistent = ifelse(n_DEG_targets > 0, n_consistent_DEG_targets / n_DEG_targets, NA),
                  consistency_p   = ifelse(n_DEG_targets > 0,
                                           pbinom(n_consistent_DEG_targets - 1, n_DEG_targets, 0.5,
                                                  lower.tail = FALSE), NA))

  tfea_full <- sample_acts %>%
    dplyr::select(-statistic, -condition) %>%
    dplyr::left_join(tf_targets, by = "source") %>%
    dplyr::arrange(p_value) %>%
    dplyr::rename(Name = source, padj_tfea = padj) %>%
    dplyr::relocate(targets, consistent_targets, .after = dplyr::last_col()) %>%
    dplyr::left_join(dplyr::select(res_tf, Name, ID, Transcript.ID, logFC, padj_dge = padj, OT.globalScore),
                     by = "Name") %>%
    dplyr::left_join(mean_tpm, by = "Name")

  tfea_net <- tfea_full %>%
    dplyr::filter(padj_tfea < padj_cutoff,
                  n_consistent_DEG_targets >= min_consistent,
                  consistency_p < consistency_p_cutoff,
                  frac_consistent > min_frac,
                  Name %in% gtf_names)

  tf_edges <- net %>% dplyr::filter(source %in% tfea_net$Name, target %in% degs$Name)

  list(full = tfea_full, network = tfea_net, edges = tf_edges, net = net, acts = sample_acts)
}

# TF network: each TF with its n_targets DEG targets of largest |stat|; fill =
# log2FC, TFs bold, genes italic, underline_genes underlined
plot_tf_network <- function(tf, degs, underline_genes = character(), n_targets = 10) {
  network_data <- tf$edges %>%
    dplyr::inner_join(degs, by = c("target" = "Name")) %>%
    dplyr::group_by(source) %>%
    dplyr::slice_max(order_by = abs(stat), n = n_targets) %>%
    dplyr::ungroup()

  border_class <- function(logfc) ifelse(is.na(logfc), "None", ifelse(logfc > 0, "Up", "Down"))
  lfc_deg      <- function(g) degs$logFC[match(g, degs$Name)]
  nodes <- dplyr::bind_rows(
    tf$network %>% dplyr::transmute(name = Name, type = "TF",
                                    shape_class = ifelse(score > 0, "TF_Pos", "TF_Neg")),
    network_data %>% dplyr::distinct(name = target) %>% dplyr::mutate(type = "Gene", shape_class = "Gene")
  ) %>%
    dplyr::distinct(name, .keep_all = TRUE) %>%
    dplyr::mutate(value       = lfc_deg(name),
                  border_type = border_class(value),
                  plot_label  = paste0('"', name, '"'),
                  plot_label  = ifelse(type == "TF", paste0("bold(", plot_label, ")"),
                                       paste0("italic(", plot_label, ")")),
                  plot_label  = ifelse(name %in% underline_genes, paste0("underline(", plot_label, ")"), plot_label))

  edges <- network_data %>% dplyr::transmute(from = source, to = target, mor = as.character(mor))

  tidygraph::tbl_graph(nodes = nodes, edges = edges) %>%
    ggraph(layout = "stress") +
    geom_edge_link(aes(color = mor), alpha = 0.4, width = 0.5) +
    geom_node_point(aes(fill = value, color = border_type, size = type, shape = shape_class), stroke = 0.5) +
    geom_node_text(aes(label = plot_label, size = I(ifelse(type == "TF", 3, 2.5))),
                   parse = TRUE, repel = TRUE, max.overlaps = 50, seed = 1, max.time = 15,
                   show.legend = FALSE) +
    scale_edge_color_manual(values = c("-1" = "blue", "1" = "orangered"), guide = "none") +
    scale_fill_gradient2(low = "#2679ff", mid = "white", high = "#EC8E12", midpoint = 0,
                         na.value = "gray90", name = expression("Log"["2"] ~ "fold-change")) +
    scale_color_manual(values = c("Up" = "#EC8E12", "Down" = "#2679ff", "None" = "gray35"),
                       guide = "none") +
    scale_shape_manual(values = c("TF_Pos" = 23, "TF_Neg" = 25, "Gene" = 21),
                       labels = c("TF_Pos" = "Activated TF", "TF_Neg" = "Inhibited TF",
                                  "Gene" = "Target Gene"),
                       name = "Node Type") +
    scale_size_manual(values = c("TF" = 3, "Gene" = 2.5), guide = "none") +
    theme_void() +
    theme(legend.position = "right",
          plot.background = element_rect(fill = "white", colour = "white"))
}


# ── Deconvolution: BayesPrism (ST2) ───────────────────────────────────────────
# Reference: Gan et al. (2021) IVD atlas (GSE160756), 9 cell types
# In:  counts (rows in gene_table order, summed per symbol); ref_dir with
#      GSE160756_RAW/*.loom and Gan_hIVD_cluster_annotations.csv
# Genes: shared protein-coding, without ribosomal, mitochondrial, Hb, sex-chromosome,
# actin, MALAT1. Seeded per sample: independent of n_cores. Runtime: hours
# Out: list(theta: samples x cell types, bp: BayesPrism object)
exec_bayesprism <- function(counts, gene_table, ref_dir,
                            cell_types = c("Chond1", "Chond2", "Chond3", "Noto", "NPPC",
                                           "Stroma", "EC", "Pericyte", "Blood"),
                            n_cores = max(1L, parallel::detectCores() - 1L),
                            seed = 123) {
  # Loom -> sparse cells x genes (keep_genes); hdf5r reads it already transposed
  read_loom <- function(path, keep_genes) {
    f <- hdf5r::H5File$new(path, mode = "r")
    on.exit(f$close_all())
    sp <- Matrix::Matrix(f[["matrix"]][, ], sparse = TRUE)
    rownames(sp) <- f[["col_attrs/cell_names"]][]
    colnames(sp) <- f[["row_attrs/gene_names"]][]
    sp[, intersect(colnames(sp), keep_genes), drop = FALSE]
  }

  stopifnot(identical(sub("[.].*$", "", rownames(counts)), gene_table$ID))
  bulk <- t(rowsum(counts, gene_table$Name))   # samples x symbols
  storage.mode(bulk) <- "double"               # totals can exceed the integer range
  # Protein-coding pre-filter before the dense conversion (memory)
  pc_genes <- intersect(unique(gene_table$Name[gene_table$Type == "protein_coding"]), colnames(bulk))

  annot <- data.table::fread(file.path(ref_dir, "Gan_hIVD_cluster_annotations.csv"), data.table = FALSE)
  names(annot)[1] <- "cell_id"
  looms <- sort(list.files(file.path(ref_dir, "GSE160756_RAW"), pattern = "\\.loom$", full.names = TRUE))
  sc <- do.call(rbind, lapply(looms, read_loom, keep_genes = pc_genes))
  # Cells without a cluster label are dropped
  cluster <- annot$cluster[match(rownames(sc), annot$cell_id)]
  labels  <- factor(cluster[!is.na(cluster)], levels = cell_types)
  sc <- as.matrix(sc[!is.na(cluster), ])
  mode(sc) <- "integer"
  bulk <- bulk[, intersect(colnames(bulk), colnames(sc)), drop = FALSE]

  gene_groups <- c("Rb", "Mrp", "other_Rb", "chrM", "MALAT1", "chrX", "chrY", "hb", "act")
  sc   <- BayesPrism::cleanup.genes(input = sc, input.type = "count.matrix", species = "hs",
                                    gene.group = gene_groups, exp.cells = 5)
  bulk <- BayesPrism::cleanup.genes(input = bulk, input.type = "count.matrix", species = "hs",
                                    gene.group = gene_groups, exp.cells = 5)
  sc   <- BayesPrism::select.gene.type(sc, gene.type = "protein_coding")

  prism <- BayesPrism::new.prism(reference = sc, mixture = bulk, input.type = "count.matrix",
                                 cell.type.labels = labels, cell.state.labels = labels,
                                 key = NULL, outlier.cut = 0.01, outlier.fraction = 0.1)
  message(sprintf("BayesPrism: %d samples, %d cells, %d genes, %d cores",
                  nrow(bulk), nrow(sc), ncol(sc), n_cores))
  rm(sc, bulk)
  bp <- BayesPrism::run.prism(prism = prism, n.cores = n_cores, gibbs.control = list(seed = seed))
  theta <- BayesPrism::get.fraction(bp = bp, which.theta = "final", state.or.type = "type")
  list(theta = as.data.frame(theta), bp = bp)
}


# ── Deconvolution clusters and bar plot (SF2) ─────────────────────────────────

# Cell types in plot order; Paul Tol "muted" palette
DECONV_CELL_COLORS <- c(
  Chond1   = "#882255",
  Chond2   = "#88CCEE",
  Chond3   = "#EE7733",
  Noto     = "#AA4499",
  NPPC     = "#44AA99",
  Stroma   = "#DDCC77",
  EC       = "#117733",
  Pericyte = "#332288",
  Blood    = "#CC6677"
)

# Ward.D2 clusters of the proportions, named by marker: D1 Chond2, D2 Stroma, D3 Blood
DECONV_NAMES  <- c(D1 = "Chondrocyte-rich", D2 = "Stroma-rich", D3 = "Blood-rich")
DECONV_MARKER <- c(D1 = "Chond2", D2 = "Stroma", D3 = "Blood")
DECONV_COLS   <- list(`2` = c(D1 = "#4D5A66", D2 = "#B8C2CC"),
                      `3` = c(D1 = "#40004b", D2 = "#9970ab", D3 = "#e7d4e8"))
deconv_label <- function(d) paste(d, DECONV_NAMES[d])   # e.g. 'D2 Stroma-rich'

deconv_tree <- function(theta, samples, cell_types = names(DECONV_CELL_COLORS))
  hclust(dist(as.matrix(theta[samples, cell_types])), method = "ward.D2")

# sample -> D1 / D2 (/ D3); warns if a name is not on its marker's richest cluster
deconv_clusters <- function(theta, samples, k, tree = deconv_tree(theta, samples)) {
  ct  <- cutree(tree, k = k)
  mu  <- sapply(DECONV_MARKER, function(m) tapply(theta[names(ct), m], ct, mean))   # cluster x name
  lab <- setNames(character(k), rownames(mu))
  left <- rownames(mu)
  for (d in c(if (k == 3) "D3", "D1", "D2")) {
    best <- left[which.max(mu[left, d])]
    lab[best] <- d
    left <- setdiff(left, best)
  }
  for (d in lab)
    if (lab[rownames(mu)[which.max(mu[, d])]] != d)
      warning(sprintf("Deconvolution k=%d: %s is not the cluster richest in %s",
                      k, deconv_label(d), DECONV_MARKER[[d]]), call. = FALSE)
  setNames(unname(lab[as.character(ct)]), names(ct))
}

# Stacked proportions in dendrogram order, outlier and cluster strips, outliers
# bold; file without extension (NULL = not saved)
plot_deconvolution <- function(theta, outs, file = NULL, k = 3,
                               width = 6000, height = 4000, dpi = 600) {
  ratio <- theta %>%
    tibble::rownames_to_column("variable") %>%
    tidyr::pivot_longer(cols = -variable, names_to = "celltype", values_to = "proportion") %>%
    dplyr::mutate(celltype = factor(celltype, levels = names(DECONV_CELL_COLORS)))

  deconv_hc    <- deconv_tree(theta, rownames(theta))
  sample_order <- deconv_hc$labels[deconv_hc$order]
  deconv_k     <- deconv_clusters(theta, rownames(theta), k, deconv_hc)
  ratio$variable <- factor(ratio$variable, levels = sample_order)

  cl_cols    <- DECONV_COLS[[as.character(k)]]
  cluster_df <- data.frame(variable = factor(sample_order, levels = sample_order),
                           cluster  = factor(deconv_k[sample_order], levels = names(cl_cols)))
  outlier_df <- data.frame(variable = factor(sample_order, levels = sample_order),
                           outlier  = factor(ifelse(sample_order %in% outs, "yes", "no"),
                                             levels = c("no", "yes")))
  outlier_col <- c(no = "#f0f0f0", yes = "#252525")
  label_faces <- ifelse(sample_order %in% outs, "bold", "plain")

  p_strip <- ggplot(cluster_df, aes(x = variable, y = 1, fill = cluster)) +
    geom_tile() +
    scale_x_discrete(limits = sample_order, expand = c(0, 0)) +
    scale_fill_manual(values = cl_cols, breaks = names(cl_cols),
                      labels = deconv_label(names(cl_cols)), drop = FALSE,
                      name = paste0("Deconvolution, k=", k)) +
    theme_void() +
    theme(plot.margin = margin(b = 8))

  p_out <- ggplot(outlier_df, aes(x = variable, y = 1, fill = outlier)) +
    geom_tile() +
    scale_x_discrete(limits = sample_order, expand = c(0, 0)) +
    scale_fill_manual(values = outlier_col, drop = FALSE, name = "Outlier") +
    theme_void() +
    theme(plot.margin = margin(b = 2))

  p_bar <- ggplot(ratio, aes(x = variable, y = proportion, fill = celltype)) +
    geom_col(width = 0.8275, colour = NA) +
    scale_fill_manual(values = DECONV_CELL_COLORS, name = "Cell type") +
    scale_y_continuous(expand = c(0, 0), breaks = seq(0, 1, 0.25)) +
    scale_x_discrete(limits = sample_order, expand = c(0, 0)) +
    theme_classic(base_size = 13) +
    labs(y = "Estimated proportion", x = NULL) +
    theme(
      axis.text.x = element_text(face = label_faces, size = 8, angle = 90,
                                 hjust = 1, vjust = 0.5),
      axis.text.y = element_text(size = 11),
      axis.ticks.x = element_blank(),
      axis.line = element_line(linewidth = 0.5),
      plot.margin = margin(t = 0)
    )

  fig <- p_out / p_strip / p_bar +
    plot_layout(heights = c(0.04, 0.04, 1), guides = "collect") &
    theme(legend.title = element_text(face = "bold", size = 13), legend.text = element_text(size = 11))

  if (!is.null(file))
    save_ggplot_png_pdf(fig, file, width = width, height = height, dpi = dpi)
  invisible(fig)
}


# ── Jin et al. (2025) heatmaps (SF3, SF4, SF5 A) ──────────────────────────────
# Samples: Ward.D2 on 1 - Pearson r of the classifier genes (z-scored VST); each
# cluster labelled with its estimated Jin subtype

# Per cohort: k Pearson clusters; dk / dks deconvolution k shown / written; blood:
# Hb genes, proportions, classifier panel; subtypes per cluster (by eye, left to
# right), valid for cluster sizes `sizes` (warns otherwise)
JIN_COHORTS <- list(
  list(tag = "all72",   label = "all 72",   clean = FALSE, k = 4, dk = 3, dks = c(2, 3), blood = TRUE,
       subtypes = c("C4", "C1", "C3", "C2"), sizes = c(26, 25, 5, 16)),
  list(tag = "clean57", label = "clean 57", clean = TRUE,  k = 3, dk = 2, dks = 2,       blood = FALSE,
       subtypes = c("C1", "C4", "C2"),       sizes = c(23, 19, 15)))

JIN_BLOOD_GENES <- c("HBG1", "HBE1", "HBG2", "HBD", "HBQ1", "HBM", "HBA1", "HBA2", "HBB")
JIN_FIG2E <- list( # Jin et al. Fig. 2e signature genes
  C1 = c("ACAN", "COL2A1", "A2M", "FGF2", "PIEZO1", "TRPV4"),
  C2 = c("THBS1", "ADAMTS5", "ASPN", "BMP2", "BMP4", "COL10A1"),
  C3 = c("S100A9", "S100A8", "S100B", "IL13", "HIF1A"),
  C4 = c("LOXL1", "SERPINH1", "ADAMTS4", "COL1A1", "COL1A2", "TNF", "IL1B", "MMP2", "MMP9"))

JIN_SUBTYPES     <- paste0("C", 1:4)
JIN_BLOCK_COLORS <- c(C1 = "#4477AA", C2 = "#228833", C3 = "#CCBB44", C4 = "#EE6677", Blood = "#882255")
JIN_COND_COL     <- c(case = "orange", control = "#3F84F5")
JIN_DIR_COL      <- c(up = "#EC8E12", down = "#2679ff")
JIN_BLOCK_TITLE  <- "Jin's subtypes"
JIN_OUTLIER_COL  <- c(no = "#f0f0f0", yes = "#252525")

# Position (in) of a rotated label right of heatmap ht_name over row slices
# 1..n_row, measured on an off-screen drawing
heatmap_side_label_pos <- function(draw_fun, width, height, n_row, n_col, names_shown, ht_name = "z(VST)") {
  pdf(NULL, width = width, height = height)
  on.exit(dev.off())
  draw_fun()
  box  <- new.env()
  edge <- function(slice, y) {
    code <- quote(assign("loc", lapply(deviceLoc(unit(1, "npc"), unit(y, "npc")),
                                       convertUnit, "in", valueOnly = TRUE), envir = box))
    if (names_shown) decorate_row_names(ht_name, eval(code), slice = slice)
    else decorate_heatmap_body(ht_name, eval(code), row_slice = slice, column_slice = n_col)
    box$loc
  }
  top <- edge(1, 1)
  bot <- edge(n_row, 0)
  list(x = top$x, y = (top$y + bot$y) / 2)
}

draw_side_label <- function(label, pos, sz)
  grid.text(label$text, x = unit(pos$x, "in") + unit(2, "mm"), y = unit(pos$y, "in"), rot = 270,
            just = c("centre", "bottom"), gp = gpar(fontsize = sz$rtitle, fontface = "bold"))

# Jin heatmaps per cohort
# In:  counts (raw, versioned Ensembl IDs x samples), theta, sample_table (Name,
#      Condition), clean_samples, gene_table
# Out: <tag>_cluster_assignments.tsv; SF3 / SF4 (all 72);
#      list(<tag> = assignments, figures$clean57_signature = list(ht, label, sz))
exec_jin_heatmaps <- function(counts, theta, sample_table, clean_samples, gene_table, out_dir,
                              jin_xlsx = "40779_2025_637_MOESM3_ESM.xlsx",
                              files = c(signature = "SF3_Jin_signature", classifier = "SF4_Jin_classifier"),
                              cohorts = JIN_COHORTS,
                              dpi = 1200) {
  # Clusters numbered left to right along the dendrogram
  run_hc <- function(x, k) {
    z  <- t(scale(t(x)))
    hc <- hclust(as.dist(1 - cor(z)), method = "ward.D2")
    ct <- cutree(hc, k = k)
    list(cl = setNames(match(ct, unique(ct[hc$order])), colnames(x)), hc = hc)
  }

  # Cell-type proportions as log(% + 1), drawn under the gene heatmap
  deconv_heatmap <- function(samples, sz) {
    dm <- log1p(t(deconv_mat[samples, , drop = FALSE]) * 100)
    Heatmap(dm, name = "log(% + 1)",
            height = if (is.null(sz$row_mm)) NULL else unit(nrow(dm) * sz$row_mm, "mm"),
            col = colorRamp2(c(0, max(dm) / 2, max(dm)), c("white", "#74c476", "#00441b")),
            cluster_rows = TRUE, clustering_method_rows = "ward.D2", show_row_dend = FALSE,
            row_title = "Cell types", row_title_gp = gpar(fontsize = sz$rtitle, fontface = "bold"),
            row_names_gp = gpar(fontsize = sz$genes),
            show_column_names = TRUE, column_names_gp = gpar(fontsize = sz$samples),
            border = TRUE)
  }

  # z(VST) heatmap of a gene panel split by Pearson cluster; deconv adds cell-type rows
  make_panel <- function(panel, vst_mat, assignments, hc, k, dk, show_rows, sz, deconv = FALSE) {
    samples <- colnames(vst_mat)
    panel <- panel[panel$gene %in% rownames(vst_mat), ]
    mat   <- t(scale(t(vst_mat[panel$gene, samples, drop = FALSE])))
    mat[!is.finite(mat)] <- 0
    blocks <- intersect(c(JIN_SUBTYPES, "Blood"), panel$Jin_block)

    ann  <- data.frame(Condition = factor(assignments[samples, "Condition"]), row.names = samples)
    acol <- list(Condition = JIN_COND_COL[levels(ann$Condition)])
    ann[["Estimated subtype"]]  <- factor(assignments[samples, "Subtype"],
                                          levels = intersect(JIN_SUBTYPES, assignments[samples, "Subtype"]))
    acol[["Estimated subtype"]] <- JIN_BLOCK_COLORS[levels(ann[["Estimated subtype"]])]
    slice_titles <- sapply(seq_len(k), function(i) unique(assignments$Subtype[assignments$HC == i]))
    if (any(assignments$Outlier == "yes")) {
      ann$Outlier  <- factor(assignments[samples, "Outlier"], levels = names(JIN_OUTLIER_COL))
      acol$Outlier <- JIN_OUTLIER_COL
    }
    dn <- paste0("Deconvolution, k=", dk)
    ann[[dn]] <- factor(assignments[samples, paste0("Deconv_k", dk)])
    acol[[dn]] <- DECONV_COLS[[as.character(dk)]][levels(ann[[dn]])]
    rann <- data.frame(factor(panel$Jin_block, levels = blocks),
                       factor(panel$Direction, levels = intersect(names(JIN_DIR_COL), panel$Direction)))
    names(rann) <- c(JIN_BLOCK_TITLE, "Direction")
    rcol <- setNames(list(JIN_BLOCK_COLORS[blocks], JIN_DIR_COL[levels(rann[[2]])]), names(rann))

    ht <- Heatmap(
      mat, name = "z(VST)",
      col = colorRamp2(c(-2, 0, 2), c("#2679ff", "white", "#EC8E12")),
      height = if (is.null(sz$row_mm)) NULL else unit(nrow(mat) * sz$row_mm, "mm"),
      cluster_columns = as.dendrogram(hc), column_split = k,
      column_title = slice_titles, column_title_gp = gpar(fontsize = sz$clust),
      column_gap = unit(sz$gap, "mm"), show_column_dend = TRUE, column_dend_height = unit(sz$dend, "mm"),
      cluster_rows = FALSE,
      row_split = factor(panel$Jin_block, levels = blocks),
      row_title_gp = gpar(fontsize = sz$rtitle, fontface = "bold"),
      row_gap = unit(sz$gap, "mm"),
      row_names_gp = gpar(fontsize = if (nrow(mat) > 60) sz$genes_dense else sz$genes),
      show_row_names = show_rows,
      show_column_names = !deconv,   # with deconv, sample names go under the cell types
      column_names_gp = gpar(fontsize = sz$samples),
      left_annotation = rowAnnotation(
        df = rann, col = rcol,
        annotation_name_gp = gpar(fontsize = sz$ann), simple_anno_size = unit(sz$strip, "mm")),
      top_annotation = HeatmapAnnotation(df = ann, col = acol, annotation_name_side = "left",
                                         annotation_name_gp = gpar(fontsize = sz$ann),
                                         simple_anno_size = unit(sz$strip, "mm"),
                                         annotation_legend_param = setNames(list(list(
                                           at = levels(ann[[dn]]), labels = deconv_label(levels(ann[[dn]])))), dn)),
      border = TRUE)
    if (deconv) ht %v% deconv_heatmap(samples, sz) else ht
  }

  # One heatmap figure with an optional side label (list(text, n_row, n_col, names_shown))
  save_single <- function(ht, file, height, sz, label = NULL) {
    draw_ht <- function() draw(ht, merge_legend = TRUE, ht_gap = unit(sz$gap * 2, "mm"))
    if (!is.null(label)) {
      ht_opt(HEATMAP_LEGEND_PADDING = unit(9, "mm"), ANNOTATION_LEGEND_PADDING = unit(9, "mm"))
      pos <- heatmap_side_label_pos(draw_ht, sz$width / dpi, height / dpi,
                                    label$n_row, label$n_col, label$names_shown)
    }
    save_grid_png_pdf(file, function() {
      draw_ht()
      if (!is.null(label)) draw_side_label(label, pos, sz)
    }, sz$width, height, dpi = dpi)
  }

  # Sizes: text pt, strips / rows mm, width px
  sz_default <- list(width = 11700, clust = 14, rtitle = 10, genes = 9, genes_dense = 4,
                     samples = 7, ann = 9, gap = 1.5, dend = 12, strip = 5, row_mm = NULL, legend = NULL)
  sz_sets <- list(default = sz_default,
                  all72   = modifyList(sz_default, list(width = 12300, clust = 15, rtitle = 11, genes = 9.5,
                                                        samples = 7, ann = 10.5, row_mm = 3.3,
                                                        legend = list(legend_title_gp = gpar(fontsize = 11, fontface = "bold"),
                                                                      legend_labels_gp = gpar(fontsize = 10.5)))))

  st <- sample_table
  rownames(st) <- st$Name

  gtf <- gene_table[!duplicated(gene_table$ID), ]
  rownames(gtf) <- gtf$ID
  cnt <- as.matrix(counts)
  colnames(cnt) <- sub("_RSR[.]bam$", "", colnames(cnt))
  ens <- sub("[.].*$", "", rownames(cnt))
  pc  <- gtf[ens, "Type"] %in% "protein_coding"
  cnt_sym <- rowsum(cnt[pc, st$Name], gtf[ens[pc], "Name"])
  storage.mode(cnt_sym) <- "integer"

  deconv_mat <- as.matrix(theta)

  genes_248 <- openxlsx::read.xlsx(jin_xlsx, sheet = "248 DEGs in the gene classifier")

  panel_a <- data.frame(gene = unlist(JIN_FIG2E), Jin_block = rep(names(JIN_FIG2E), lengths(JIN_FIG2E)),
                        stringsAsFactors = FALSE, row.names = NULL)
  panel_a$Direction <- ifelse(panel_a$gene == "HIF1A", "down", "up")
  panel_a_blood <- rbind(panel_a, data.frame(gene = JIN_BLOOD_GENES, Jin_block = "Blood", Direction = "up"))

  panel_c <- unique(data.frame(gene = genes_248$SYMBOL, Jin_block = genes_248$Cluster,
                               Direction = "up", stringsAsFactors = FALSE))

  # f: file (NULL = not saved); h: height px
  panels_for <- function(co) if (co$blood) list(
    list(p = panel_c, f = files[["classifier"]], rows = FALSE, h = 10080,
         label = "Classifier genes"),
    list(p = panel_a_blood, f = files[["signature"]], rows = TRUE, h = 10320, deconv = TRUE,
         sz = sz_sets$all72, label = "Signature genes")
  ) else list(
    list(p = panel_a, f = NULL, rows = TRUE, h = 7440, label = "Signature genes"))

  out <- list()
  figures <- list()
  for (co in cohorts) {
    samples <- if (co$clean) clean_samples else st$Name

    cnt  <- cnt_sym[, samples]
    cnt  <- cnt[rowSums(cnt) > 0, ]
    vst_mat  <- assay(vst(DESeqDataSetFromMatrix(cnt, data.frame(row.names = samples), ~1), blind = TRUE))
    clf  <- intersect(unique(genes_248$SYMBOL), rownames(vst_mat))

    hc  <- run_hc(vst_mat[clf, ], co$k)
    sizes <- as.vector(table(factor(hc$cl, levels = seq_len(co$k))))
    if (!identical(sizes, as.integer(co$sizes)))
      warning(sprintf("%s: Pearson cluster sizes %s differ from %s - re-check the subtype mapping",
                      co$label, paste(sizes, collapse = "/"), paste(co$sizes, collapse = "/")), call. = FALSE)
    assignments <- data.frame(Name = samples, Condition = st[samples, "Condition"], HC = hc$cl[samples],
                              Subtype = co$subtypes[hc$cl[samples]],
                              Outlier = ifelse(samples %in% clean_samples, "no", "yes"),
                              sapply(co$dks, function(dk) deconv_clusters(theta, samples, dk)[samples]),
                              row.names = samples, stringsAsFactors = FALSE)
    colnames(assignments)[(ncol(assignments) - length(co$dks) + 1):ncol(assignments)] <- paste0("Deconv_k", co$dks)
    data.table::fwrite(assignments, file.path(out_dir, sprintf("%s_cluster_assignments.tsv", co$tag)), sep = "\t")
    out[[co$tag]] <- assignments

    for (pp in panels_for(co)) {
      sz <- if (is.null(pp$sz)) sz_sets$default else pp$sz
      # legend fonts are global ht_opt options: set for this figure only
      if (!is.null(sz$legend)) do.call(ht_opt, sz$legend)
      n_blocks <- length(intersect(c(JIN_SUBTYPES, "Blood"), pp$p$Jin_block[pp$p$gene %in% rownames(vst_mat)]))
      label    <- list(text = pp$label, n_row = n_blocks, n_col = co$k, names_shown = pp$rows)
      ht       <- make_panel(pp$p, vst_mat, assignments, hc$hc, co$k, co$dk, pp$rows, sz,
                             deconv = isTRUE(pp$deconv))
      if (!is.null(pp$f))
        save_single(ht, file.path(out_dir, pp$f), pp$h, sz, label = label)
      ht_opt(RESET = TRUE)
      if (co$tag == "clean57" && pp$label == "Signature genes")
        figures$clean57_signature <- list(ht = ht, label = label, sz = sz)
    }
  }
  out$figures <- figures
  out
}


# ── SVs vs. clusters and heterogeneity (SF5) ──────────────────────────────────
# dds: with SV1, SV2; clusters: exec_jin_heatmaps()$clean57; nominal p;
# file: table (NULL = not written)

stars_p <- function(p) ifelse(is.na(p), "", ifelse(p < 0.001, "***", ifelse(p < 0.01, "**", ifelse(p < 0.05, "*", ""))))

# Panel B: Pearson r of SV1, SV2, D2, C1 / C2 / C4, Condition (SV signs arbitrary)
exec_sv_correlation <- function(dds, clusters, file = NULL) {
  cd  <- as.data.frame(colData(dds))
  smp <- rownames(cd)
  jc  <- clusters
  stopifnot(setequal(smp, jc$Name))

  subtype <- jc[smp, "Subtype"]
  x_vars  <- cbind("SV1"                = cd$SV1,
                   "SV2"                = cd$SV2,
                   "Deconvolution (D2)" = as.numeric(jc[smp, "Deconv_k2"] == "D2"),
                   "C1"                 = as.numeric(subtype == "C1"),
                   "C2"                 = as.numeric(subtype == "C2"),
                   "C4"                 = as.numeric(subtype == "C4"),
                   "Condition (case)"   = as.numeric(cd$Condition == "case"))
  r <- cor(x_vars)
  p <- outer(seq_len(ncol(x_vars)), seq_len(ncol(x_vars)), Vectorize(function(i, j)
    if (i == j) NA_real_ else cor.test(x_vars[, i], x_vars[, j])$p.value))
  dimnames(p) <- dimnames(r)

  if (!is.null(file))
    write.table(data.frame(Variable = rownames(r), setNames(as.data.frame(r), paste0("r_", colnames(r))),
                           setNames(as.data.frame(p), paste0("p_", colnames(p))), check.names = FALSE),
                file, sep = "\t", quote = FALSE, row.names = FALSE)

  lab   <- matrix(paste0(sprintf("%.2f", r), stars_p(p)), nrow(r), dimnames = dimnames(r))
  grp   <- factor(rep(c("SVs", "Clusters", "Condition"), c(2, 4, 1)), levels = c("SVs", "Clusters", "Condition"))
  ht <- Heatmap(
    r, name = "r",
    col = colorRamp2(c(-1, 0, 1), c("#2679ff", "white", "#EC8E12")),
    cluster_rows = FALSE, cluster_columns = FALSE,
    row_split = grp, column_split = grp, row_title = NULL, column_title = NULL,
    row_gap = unit(1.5, "mm"), column_gap = unit(1.5, "mm"),
    row_names_side = "left", row_names_gp = gpar(fontsize = 12),
    column_labels = c("SV1", "SV2", "Deconv", "C1", "C2", "C4", "Condition"),
    column_names_side = "bottom", column_names_rot = 0, column_names_centered = TRUE,
    column_names_gp = gpar(fontsize = c(12, 12, 8, 12, 12, 12, 8)),
    rect_gp = gpar(col = "white", lwd = 2),
    border = TRUE, border_gp = gpar(col = "grey30", lwd = 0.5),
    width = unit(ncol(r) * 14, "mm"), height = unit(nrow(r) * 14, "mm"),
    cell_fun = function(j, i, x, y, w, h, fill)
      grid.text(if (i == j) "" else lab[i, j], x, y,
                gp = gpar(fontsize = 10, col = if (abs(r[i, j]) > 0.6) "white" else "black")),
    heatmap_legend_param = list(at = seq(-1, 1, 0.5), title_gp = gpar(fontsize = 12, fontface = "bold"),
                                labels_gp = gpar(fontsize = 11), legend_height = unit(35, "mm")))

  list(r = r, p = p, ht = ht)
}

# Panel C: R^2 of SV1, SV2 (lm, F test) and the SV plane (canonical correlation,
# Pillai) with Condition (control), deconvolution, subtype, cell-type %
exec_sv_heterogeneity <- function(dds, clusters, theta, file = NULL,
                                  props = c("Chond2", "Stroma", "NPPC", "EC", "Blood")) {
  cd  <- as.data.frame(colData(dds))
  smp <- rownames(cd)
  sv  <- as.matrix(cd[, c("SV1", "SV2")])
  jc  <- clusters
  stopifnot(setequal(smp, jc$Name))
  stopifnot(all(smp %in% rownames(theta)))

  groupings <- list("Condition"         = factor(cd$Condition),
                    "Deconvolution"     = factor(jc[smp, "Deconv_k2"]),
                    "Estimated subtype" = factor(jc[smp, "Subtype"]))
  vars  <- c(groupings, setNames(lapply(props, function(ct) theta[smp, ct] * 100),
                                 paste0(props, " %")))

  cc2  <- function(x, y) cancor(model.matrix(~x)[, -1, drop = FALSE], y)$cor[1]^2
  pval <- function(x, y) {
    if (ncol(y) == 1) anova(lm(y[, 1] ~ x))[["Pr(>F)"]][1]
    else summary(manova(y ~ x), test = "Pillai")$stats[1, "Pr(>F)"]
  }
  sv_sets <- list(SV1 = sv[, "SV1", drop = FALSE], SV2 = sv[, "SV2", drop = FALSE], `SV1 + SV2` = sv)
  r2 <- sapply(sv_sets, function(y) sapply(vars, cc2, y = y))
  p  <- sapply(sv_sets, function(y) sapply(vars, pval, y = y))

  if (!is.null(file))
    write.table(data.frame(Variable = rownames(r2),
                           setNames(as.data.frame(r2), paste0("R2_", colnames(r2))),
                           setNames(as.data.frame(p),  paste0("p_",  colnames(p))),
                           check.names = FALSE),
                file, sep = "\t", quote = FALSE, row.names = FALSE)

  lab   <- matrix(paste0(sprintf("%.2f", r2), stars_p(p)), nrow(r2), dimnames = dimnames(r2))
  ht <- Heatmap(
    r2, name = "R²",
    col = colorRamp2(c(0, 1), c("white", "#08306b")),
    cluster_rows = FALSE, cluster_columns = FALSE,
    row_split = factor(rep(c("Groupings", "Cell-type %"), c(length(groupings), length(props))),
                       levels = c("Groupings", "Cell-type %")),
    column_split = factor(c("Single SV", "Single SV", "Joint"), levels = c("Single SV", "Joint")),
    row_title = NULL, column_title = NULL,
    border = TRUE, border_gp = gpar(col = "grey30", lwd = 0.5),
    row_names_side = "left", row_names_gp = gpar(fontsize = 12),
    column_names_side = "bottom", column_names_rot = 0, column_names_centered = TRUE,
    column_names_gp = gpar(fontsize = 12),
    rect_gp = gpar(col = "white", lwd = 2),
    width = unit(ncol(r2) * 24, "mm"), height = unit(nrow(r2) * 9, "mm"),
    cell_fun = function(j, i, x, y, w, h, fill)
      grid.text(lab[i, j], x, y, gp = gpar(fontsize = 11, col = if (r2[i, j] > 0.5) "white" else "black")),
    heatmap_legend_param = list(at = seq(0, 1, 0.25), title_gp = gpar(fontsize = 12, fontface = "bold"),
                                labels_gp = gpar(fontsize = 11), legend_height = unit(30, "mm")))

  list(r2 = r2, p = p, ht = ht)
}

# SF5: A = signature heatmap (list(ht, label, sz)), B and C below; px at dpi
plot_sv_figure <- function(panel_a, ht_b, ht_c, file, dpi = 1200,
                           widths = c(7560, 6360), heights = c(7440, 5400)) {
  width_px  <- sum(widths)
  height_px <- sum(heights)

  letter <- function(l) grid.text(l, x = unit(1.5, "mm"), y = unit(1, "npc") - unit(1.5, "mm"),
                                  just = c("left", "top"), gp = gpar(fontsize = 30, fontface = "bold"))
  # upViewport (not popViewport) keeps the heatmaps' viewports for decorate_*()
  in_cell <- function(row, col, code) {
    pushViewport(viewport(layout.pos.row = row, layout.pos.col = col))
    code
    upViewport()
  }
  draw_page <- function() {
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(2, 2, heights = unit(heights / dpi, "in"),
                                               widths  = unit(widths / dpi, "in"))))
    in_cell(1, 1:2, {
      draw(panel_a$ht, newpage = FALSE, merge_legend = TRUE, ht_gap = unit(panel_a$sz$gap * 2, "mm"))
      letter("A")
    })
    in_cell(2, 1, {
      draw(ht_b, newpage = FALSE)
      letter("B")
    })
    in_cell(2, 2, {
      draw(ht_c, newpage = FALSE)
      letter("C")
    })
    upViewport()
  }

  # room for panel A's side label, as in its own figure
  ht_opt(HEATMAP_LEGEND_PADDING = unit(9, "mm"), ANNOTATION_LEGEND_PADDING = unit(9, "mm"))
  on.exit(ht_opt(RESET = TRUE), add = TRUE)

  pos <- heatmap_side_label_pos(draw_page, width_px / dpi, height_px / dpi, panel_a$label$n_row, panel_a$label$n_col,
                                names_shown = TRUE)
  save_grid_png_pdf(file, function() {
    draw_page()
    draw_side_label(panel_a$label, pos, panel_a$sz)
  }, width_px, height_px, dpi = dpi)
}
