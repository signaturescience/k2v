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

# Re-sort
d <- d %>%
  mutate(chr=factor(chr, levels=c(as.character(1:22), "X", "Y"))) %>%
  arrange(chr, pos) %>%
  mutate(chr=as.character(chr))

# Add the N29insertion, but bcftools tsv2vcf will not convert.
# https://www.ncbi.nlm.nih.gov/snp/rs796296176: aka N29insA 16:89985750 C-CA
# d <- d %>% add_row(chr="16", pos=89985750L, rsid="rs796296176", ref="C", alt="CA", af=0.00285)

# Write out
d %>% write_csv(here::here("build/kintelligence.alleles.csv"))

# Overwrite the current rsids
# d$rsid %>% write_lines(here::here("build/kintelligence.rsids.txt"))

