# featureCounts, all samples 
# In:  <samplesheet.csv> <annotation.gtf> <threads> <out_prefix>; <sample>_RSR.bam
#      in the working dir, columns in samplesheet order
# Out: <out_prefix>.rds (featureCounts list), <out_prefix>.summary (MultiQC), versions.txt
args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 4)
suppressPackageStartupMessages(library(Rsubread))

samples   <- read.csv(args[1], stringsAsFactors = FALSE)$sample
filenames <- paste0(samples, "_RSR.bam")
missing   <- filenames[!file.exists(filenames)]
if (length(missing)) stop("BAM files not found:\n", paste(missing, collapse = "\n"))

RSR <- featureCounts(files = filenames,
                     strandSpecific = 2,
                     isPairedEnd = TRUE,
                     isGTFAnnotationFile = TRUE,
                     GTF.featureType = "exon",
                     GTF.attrType = "gene_id",
                     annot.ext = args[2],
                     nthreads = as.integer(args[3]),
                     verbose = TRUE)

# Save before checks
saveRDS(RSR, file = paste0(args[4], ".rds"))
write.table(RSR$stat, paste0(args[4], ".summary"), sep = "\t", quote = FALSE,
            row.names = FALSE)

# Some Rsubread versions alter column names (make.names); the analysis matches by name
bad <- which(colnames(RSR$counts) != filenames)
if (length(bad))
  warning(length(bad), " featureCounts column names differ from the BAM file names, e.g. '",
          colnames(RSR$counts)[bad[1]], "' for '", filenames[bad[1]],
          "'. The analysis matches samples by name: check the counts.", call. = FALSE)
writeLines(paste0("Rsubread: ", packageVersion("Rsubread")), "versions.txt")
