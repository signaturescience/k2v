suppressPackageStartupMessages(suppressWarnings(library(tidyverse)))
d <- read_tsv("tmp", col_names=c("chr", "pos", "rsid", "ref", "alt", "af"), col_types="cicccd")
d <- d %>%
  group_by(chr, pos, rsid) %>%
  slice_max(order_by=af, n=1) %>%
  ungroup()
d %>% write_csv("kintelligence.alleles.csv")

