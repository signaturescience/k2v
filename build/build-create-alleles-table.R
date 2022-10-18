suppressPackageStartupMessages(suppressWarnings(library(tidyverse)))
d <- read_tsv(here::here("build/tmp"), col_names=c("chr", "pos", "rsid", "ref", "alt", "af"), col_types="cicccd")
d <- d %>%
  group_by(chr, pos, rsid) %>%
  slice_max(order_by=af, n=1) %>%
  ungroup()

# Manually add other missing sites
d <- d %>%
  # https://www.ncbi.nlm.nih.gov/snp/rs201326893 aka rs201326893_Y152OCH 16:89986122 C-A
  add_row(chr="16", pos=89986122L, rsid="rs201326893_Y152OCH", ref="C", alt="A", af=0.00145) %>%
  # https://www.ncbi.nlm.nih.gov/snp/rs938283 17:77468498 T-C
  add_row(chr="17", pos=77468498L, rsid="rs938283", ref="T", alt="C", af=0.16352) %>%
  # https://www.ncbi.nlm.nih.gov/snp/rs34228583 X:11139882 G-A
  add_row(chr="X", pos=11139882L, rsid="rs34228583", ref="G", alt="A", af=0.00001)

# You next need to see which SNPs are genotyped in the sample data but NOT yet in the data from Gnomad.
# Get info from dbSNP on these SNPs, and put back into the allele table
fp <- here::here("build/y_snps_to_fill.rds")
if (!file.exists(fp)) {
  sampledata <- read_tsv(here::here("testdata/HG002.SnpResult.txt"), col_types="ciiicciii")
  hg19 <- read_tsv(here::here("build/hg19_mapping.tsv"), col_types="cic")
  hg19 <- hg19 %>% transmute(chr=Chromosome %>% gsub("23", "X", .) %>% gsub("24", "Y", .), pos=Position, rsid=RsId)
  hg19 <- hg19 %>% semi_join(sampledata, by=c("rsid"="Locus_ID"))
  missingysnps <- hg19 %>%
    anti_join(d, by=c("chr", "pos", "rsid")) %>%
    filter(chr=="Y")
  ysnpsfromdbsnp <- rsnps::ncbi_snp_query(unique(missingysnps$rsid))
  y_snps_to_fill <-
    ysnpsfromdbsnp %>%
    select(rsid=query, ref=ancestral_allele, alt=variation_allele) %>%
    filter(nchar(ref)==1) %>%
    filter(nchar(alt)==1) %>%
    distinct() %>%
    inner_join(missingysnps, by="rsid") %>%
    transmute(chr, pos, rsid, ref, alt, af=NA)
  saveRDS(y_snps_to_fill, file=fp)
} else {
  y_snps_to_fill <- readRDS(fp)
}
# Add those SNPs
d <- bind_rows(d, y_snps_to_fill)


# Manually add final missing Y sites
d <- d %>%
  # https://www.ncbi.nlm.nih.gov/snp/rs9785996 Y:17594966 G-C
  add_row(chr="Y", pos=17594966L, rsid="rs9785996", ref="G", alt="C", af=0.47889) %>%
  # https://www.ncbi.nlm.nih.gov/snp/rs9785996 Y:17594966 G-C
  add_row(chr="Y", pos=23655112L, rsid="rs9786064", ref="A", alt="G", af=0.36560)


# Re-sort
d <- d %>%
  mutate(chr=factor(chr, levels=c(as.character(1:22), "X", "Y"))) %>%
  arrange(chr, pos) %>%
  mutate(chr=as.character(chr)) %>%
  distinct()

# Add the N29insertion, but bcftools tsv2vcf will not convert.
# https://www.ncbi.nlm.nih.gov/snp/rs796296176: aka N29insA 16:89985750 C-CA
# d <- d %>% add_row(chr="16", pos=89985750L, rsid="rs796296176", ref="C", alt="CA", af=0.00285)

# Write out
d %>% write_csv(here::here("build/kintelligence.alleles.csv"))

# Overwrite the current rsids
# d$rsid %>% write_lines(here::here("build/kintelligence.rsids.txt"))

