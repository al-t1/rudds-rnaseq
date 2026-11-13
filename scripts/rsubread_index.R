# Genome index
# In:  <fasta> <index_basename>
# Out: <index_basename>.*, versions.txt
args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 2)
suppressPackageStartupMessages(library(Rsubread))

Rsubread::buildindex(basename = args[2], reference = args[1])
writeLines(paste0("Rsubread: ", packageVersion("Rsubread")), "versions.txt")
