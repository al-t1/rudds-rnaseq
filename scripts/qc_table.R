# Supplementary Table 1 
# In:  --core <core_functions.R> --metadata <tsv> --counts <rds>; fastqc/*.zip,
#      rnaseqc/*.metrics.tsv in the working dir
# Out: --out <ST1.tsv>: 3 header lines, one row per metadata sample
#      (metadata | FastQC Phred, %GC | Rsubread read categories | RNA-SeQC genes, HQ rate)
args <- commandArgs(trailingOnly = TRUE)
opt  <- setNames(args[seq(2, length(args), 2)], sub("^--", "", args[seq(1, length(args), 2)]))
for (k in c("core", "metadata", "counts", "out"))
  if (is.na(opt[k])) stop("Missing --", k, call. = FALSE)

source(opt[["core"]])

meta <- data.table::fread(opt[["metadata"]], sep = "\t", data.table = FALSE,
                          check.names = FALSE, colClasses = list(character = "Name"))
fq <- parse_fastqc("fastqc")
rs <- rsubread_stats(readRDS(opt[["counts"]])$stat)
rq <- parse_rnaseqc("rnaseqc", metrics = c("Genes Detected", "High Quality Rate"))

parts <- list(FastQC = fq, Rsubread = rs, `RNA-SeQC` = rq)
for (p in names(parts)) {
  absent <- setdiff(meta$Name, parts[[p]]$Name)
  extra  <- setdiff(parts[[p]]$Name, meta$Name)
  if (length(absent)) stop(p, ": no results for ", paste(absent, collapse = ", "), call. = FALSE)
  if (length(extra))  stop(p, ": samples not in the metadata: ", paste(extra, collapse = ", "), call. = FALSE)
}
st1 <- cbind(meta,
             fq[match(meta$Name, fq$Name), -1],
             rs[match(meta$Name, rs$Name), -1],
             rq[match(meta$Name, rq$Name), -1])
rownames(st1) <- NULL

# Group label above each group's first column
groups <- c(Metadata = ncol(meta), FastQC = ncol(fq) - 1,
            Rsubread = ncol(rs) - 1, `RNA-SeQC` = ncol(rq) - 1)
group_row <- unlist(lapply(names(groups), function(g) c(g, rep("", groups[[g]] - 1))))
pad <- function(x) paste(c(x, rep("", ncol(st1) - 1)), collapse = "\t")
header <- c(pad("Supplementary Table 1. Case-control status, RNA integrity and sequencing quality by sample."),
            pad(paste0("Updated: ", format(Sys.Date(), "%Y%m%d"))),
            paste(group_row, collapse = "\t"))
writeLines(header, opt[["out"]])
data.table::fwrite(st1, opt[["out"]], sep = "\t", na = "NA", quote = FALSE,
                   append = TRUE, col.names = TRUE)

