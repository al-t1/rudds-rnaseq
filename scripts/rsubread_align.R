# Paired-end alignment, one sample
# In:  <index_basename> <sample> <R1.fq.gz> <R2.fq.gz> <threads>
# Out: <sample>_RSR.bam (+ .summary, .indel.vcf), versions.txt
args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 5)
suppressPackageStartupMessages(library(Rsubread))

Rsubread::align(index        = args[1],
                readfile1    = args[3],
                readfile2    = args[4],
                type         = "rna",
                input_format = "gzFASTQ",
                output_format = "BAM",
                output_file  = paste0(args[2], "_RSR.bam"),
                nthreads     = as.integer(args[5]))
writeLines(paste0("Rsubread: ", packageVersion("Rsubread")), "versions.txt")
