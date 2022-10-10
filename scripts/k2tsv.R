if (interactive()) {
  infile <- here::here("exampledata/NA24385.SnpResult.txt")
  kfile <- here::here("assets/kintelligence.alleles.csv")
  outfile <- here::here("tmp.23.txt")
} else {
  args = commandArgs(trailingOnly=TRUE)
  infile <- args[1]
  kfile <- "/assets/kintelligence.alleles.csv"
  outfile <- "/tmp/tmp.23.txt"
}
if (!file.exists(infile)) stop(paste("File doesn't exist:", infile))
if (!file.exists(kfile)) stop(paste("File doesn't exist:", kfile))

k <- read.csv(kfile, stringsAsFactors = FALSE)
d <- read.delim(infile)

d <- d[,c("Locus_ID", "Genotype")]
names(d) <- c("rsid", "gt")
d$gt[grep("\\.", d$gt)] <- NA
d$gt[d$gt=="1/0"] <- "0/1"
if (!identical(sort(unique(na.omit(d$gt))), c("0/0", "0/1", "1/1"))) stop("Genotypes aren't 0/0, 0/1, 1/1.")
d <- merge(d, k)

tgt <- Vectorize(function(gt, ref, alt) switch(as.character(gt), `NA`=NA, "0/0"=paste0(ref, ref), "0/1"=paste0(ref, alt), "1/1"=paste0(alt,alt)))
d$tgt <- unlist(tgt(d$gt, d$ref, d$alt))

d23 <- d[,c("rsid", "chr", "pos","tgt")]

write.table(d23, file=outfile, quote=FALSE, sep="\t", na="--", row.names=FALSE, col.names=FALSE)
