# Get input data paths
if (interactive()) {
  # If interactive, use data in ./
  infile <- here::here("testdata/HG002.SnpResult.txt")
  kfile <- here::here("assets/kintelligence.alleles.csv")
  outfile <- here::here("tmp.23.txt")
} else {
  # If running in container, use paths in /
  args = commandArgs(trailingOnly=TRUE)
  infile <- args[1]
  kfile <- "/assets/kintelligence.alleles.csv"
  outfile <- "/tmp/tmp.23.txt"
}
if (!file.exists(infile)) stop(paste("File doesn't exist:", infile))
if (!file.exists(kfile)) stop(paste("File doesn't exist:", kfile))

# Read in kintelligence allele data. See build/README.sh
k <- read.csv(kfile, stringsAsFactors = FALSE)
# Read input .SnpResult.txt file
d <- read.delim(infile)

#' @title Translate genotypes
#' @description Takes a genotype, reference, and alternate allele, return letter genotypes. Vectorized over gt, ref, alt.
tgt <- Vectorize(function(gt, ref, alt) switch(as.character(gt), `NA`=NA, "0"=ref, "1"=alt, "0/0"=paste0(ref, ref), "0/1"=paste0(ref, alt), "1/1"=paste0(alt,alt)))
#' @examples
#' tgt(c("0/0", "0/1", "1/1", NA, "0", "1"), "A", "B")
#' transform(data.frame(gt=rep(c("0/0", "0/1", "1/1"),3), ref=LETTERS[1:9], alt=LETTERS[10:18]), tgt=tgt(gt, ref, alt))

# Subset to columns of interest
d <- d[,c("Locus_ID", "Genotype")]
names(d) <- c("rsid", "gt")

# Change . to NA, and 1/0 to 0/1. Check genotypes are expected 0/0, 0/1, 1,1
d$gt[grep("\\.", d$gt)] <- NA
d$gt[d$gt=="1/0"] <- "0/1"
if (!identical(sort(unique(na.omit(d$gt))), c("0/0", "0/1", "1/1"))) stop("Genotypes aren't 0/0, 0/1, 1/1.")

# Join genotypes to kintelligence allele table
d <- merge(d, k)

# Make Y genotypes haploid
d$gt[d$chr=="Y" & d$gt=="0/0"] <- "0"
d$gt[d$chr=="Y" & d$gt=="1/1"] <- "1"

# Translate genotypes see tgt() function
d$tgt <- unlist(tgt(d$gt, d$ref, d$alt))

# Collect relevant columns and write to file
d23 <- d[,c("rsid", "chr", "pos","tgt")]
write.table(d23, file=outfile, quote=FALSE, sep="\t", na="--", row.names=FALSE, col.names=FALSE)
