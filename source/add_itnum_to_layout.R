#!/usr/bin/env Rscript
options(width=150)
library(tidyverse)

# Load input directory and layout file
args = commandArgs(trailingOnly=TRUE)
idir_name = args[1]
layo_name = args[2]
itnm_name = file.path(idir_name, "sample_itnum.tab")

# Load data
layo = read_tsv(layo_name)
itnm = read_tsv(itnm_name, col_names=c("sample", "itnum"))

# Save layout table with itnum
write_tsv(
  inner_join(itnm, layo) %>% select(itnum, type, name),
  file.path(idir_name, "layout.tab")
)
