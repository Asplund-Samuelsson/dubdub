#!/usr/bin/env Rscript
options(width=150)
library(tidyverse)

# Load input valid samples file name and layout file name
args = commandArgs(trailingOnly=TRUE)
vals_name = args[1]
smit_name = args[2]
layo_name = args[3]

# Determine if valid samples is dedup and name output file accordingly
if (endsWith(vals_name, ".dedup.txt")){
  outl_name = str_replace(layo_name, ".tab$", ".valid.dedup.tab")
} else {
  outl_name = str_replace(layo_name, ".tab$", ".valid.tab")
}

# Load data
vals = scan(vals_name, character())
smit = read_tsv(smit_name, col_names=c("Sample", "itnum"))
layo = read_tsv(layo_name)

# Filter layout and write to outfile
write_tsv(semi_join(layo, filter(smit, Sample %in% vals)), outl_name)
