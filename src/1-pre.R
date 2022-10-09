if (interactive()) {
  infile <- here::here("exampledata/NA24385.SnpResult.txt")
  k <- read.csv(here::here("assets/kintelligence.alleles.csv"), stringsAsFactors = FALSE)
} else {
  args = commandArgs(trailingOnly=TRUE)
  infile <- args[1]
  k <- read.csv("/assets/kintelligence.alleles.csv", stringsAsFactors = FALSE)
}
if (!file.exists(infile)) stop(paste("File doesn't exist:", infile))

d <- read.delim(infile)
d <- d[,c("Locus_ID", "Genotype")]
names(d) <- c("rsid", "gt")
d$gt[grep("\\.", d$gt)] <- NA
d$gt[d$gt=="1/0"] <- "0/1"
if (!identical(sort(unique(na.omit(d$gt))), c("0/0", "0/1", "1/1"))) stop("Genotypes aren't 0/0, 0/1, 1/1.")
d <- merge(d, k)

tgt <- Vectorize(function(gt, ref, alt) switch(as.character(gt), `NA`="--", "0/0"=paste0(ref, ref), "0/1"=paste0(ref, alt), "1/1"=paste0(alt,alt)))
d$tgt <- unlist(tgt(d$gt, d$ref, d$alt))

d23 <- d[,c("rsid", "chr", "pos","tgt")]
