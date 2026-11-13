# Freezes the analysis' online resources into assets/ (manual step, not in the
# pipeline; needs the analysis R library and internet). Package caches are used
# where present (OmnipathR, msigdbr, ExperimentHub).
# In:  [assets_dir] [--force] [--skip-experimenthub]
# Out: CollecTRI_human.tsv    TF network (source, target, mor)
#      MSigDB_C3_MIRDB.tsv    C3 MIR:MIRDB sets (gs_name, gs_description, ncbi_gene)
#      KEGG_hsa_pathways.tsv  KEGG REST pathways (pathway, name, gene)
#      ExperimentHub/         cache with EH3226 (LINCS), EH3228 (targets), EH3234 (ES null)
#      Existing outputs kept unless --force

args  <- commandArgs(trailingOnly = TRUE)
force <- "--force" %in% args
skip_eh <- "--skip-experimenthub" %in% args
assets <- c(args[!startsWith(args, "--")], "assets")[1]
dir.create(assets, showWarnings = FALSE, recursive = TRUE)

todo <- function(name) {
  path <- file.path(assets, name)
  if (file.exists(path) && !force) {
    message("exists, kept (--force to redo): ", path)
    return(NULL)
  }
  message("writing: ", path)
  path
}

# Write, read back, check identical
write_tsv <- function(df, path) {
  data.table::fwrite(df, path, sep = "\t")
  back <- data.table::fread(path, data.table = FALSE, colClasses = "character", na.strings = "")
  df[] <- lapply(df, function(x) { x <- as.character(x); x[x %in% ""] <- NA; x })
  stopifnot(identical(dim(back), dim(df)),
            all(mapply(identical, back, df, USE.NAMES = FALSE)))
  message(sprintf("  %d rows, verified", nrow(df)))
}

# ── CollecTRI ─────────────────────────────────────────────────────────────────
if (!is.null(path <- todo("CollecTRI_human.tsv"))) {
  net <- decoupleR::get_collectri(organism = "human", split_complexes = FALSE)
  write_tsv(as.data.frame(net)[, c("source", "target", "mor")], path)
}

# ── MSigDB C3 miRDB ───────────────────────────────────────────────────────────
# One row per set and gene
if (!is.null(path <- todo("MSigDB_C3_MIRDB.tsv"))) {
  m <- msigdbr::msigdbr(species = "Homo sapiens", collection = "C3", subcollection = "MIR:MIRDB")
  write_tsv(unique(data.frame(gs_name = m$gs_name, gs_description = m$gs_description,
                              ncbi_gene = as.character(m$ncbi_gene))), path)
}

# ── KEGG ──────────────────────────────────────────────────────────────────────
# As clusterProfiler downloads them for enrichKEGG / gseKEGG
if (!is.null(path <- todo("KEGG_hsa_pathways.tsv"))) {
  k <- get("download.KEGG.Path", asNamespace("clusterProfiler"))("hsa")
  e <- k$KEGGPATHID2EXTID; n <- k$KEGGPATHID2NAME
  out <- data.frame(pathway = e$from, name = n$to[match(e$from, n$from)], gene = e$to)
  nogene <- n[!n$from %in% e$from, ]
  if (nrow(nogene)) out <- rbind(out, data.frame(pathway = nogene$from, name = nogene$to, gene = ""))
  write_tsv(out, path)
}

# ── ExperimentHub ─────────────────────────────────────────────────────────────
# Copied from the default cache if there (keeps its snapshot), else downloaded.
# EH3233 (taurefList) only needed with tau = TRUE
EH_IDS <- c("EH3226", "EH3228", "EH3234")
if (!skip_eh && !is.null(target <- todo("ExperimentHub"))) {
  suppressMessages(library(ExperimentHub))
  if (dir.exists(target)) unlink(target, recursive = TRUE)
  dir.create(target)
  default <- getExperimentHubOption("CACHE")
  eh <- tryCatch(suppressMessages(ExperimentHub(localHub = TRUE)), error = function(e) NULL)
  have <- !is.null(eh) && all(EH_IDS %in% names(eh))
  if (have) {
    message("  copying from the cache ", default, " (snapshot ", snapshotDate(eh), ")")
    info  <- as.data.frame(BiocFileCache::bfcinfo(BiocFileCache::BiocFileCache(default, ask = FALSE)))
    keep  <- grepl("experimenthub\\.(sqlite3|index\\.rds)$", info$rname) |
             sub(" : .*", "", info$rname) %in% EH_IDS
    files <- c(file.path(default, "BiocFileCache.sqlite"), info$rpath[keep])
    stopifnot(all(file.copy(files, target)))
  } else {
    message("  downloading into the new cache (~2.5 GB)")
    setExperimentHubOption("CACHE", normalizePath(target))
    setExperimentHubOption("ASK", FALSE)
    eh <- ExperimentHub()
    for (id in EH_IDS) eh[[id]]
  }
  # Check offline
  setExperimentHubOption("CACHE", normalizePath(target))
  eh <- suppressMessages(ExperimentHub(localHub = TRUE))
  stopifnot(all(EH_IDS %in% names(eh)))
  for (id in EH_IDS) stopifnot(file.exists(suppressMessages(eh[[id]])))
  message("  ", paste(EH_IDS, collapse = ", "), " found offline, snapshot ", snapshotDate(eh))
}

message("done: ", normalizePath(assets))
