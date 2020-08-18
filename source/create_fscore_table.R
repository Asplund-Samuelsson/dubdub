#!/usr/bin/env Rscript
options(width=150)
library(tidyverse)

# Load input directory
args = commandArgs(trailingOnly=TRUE)
idir_name = args[1]

# Define infiles
smit_file = file.path(idir_name, "sample_itnum.tab")
frag_file = file.path(idir_name, "fscore", "fscore_base.tsv")
fdir = file.path(idir_name, "fscore")
fscr_files = grep("IT.*fscore.tsv", list.files(fdir), value=T)

# Load data
smit = read_tsv(smit_file, col_names = c("sample", "itnum"))
frag = read_tsv(frag_file)
fscr = bind_rows(lapply(
  fscr_files,
  # Load each gene score file
  function (infile) {
    # Extract the itnum from the filename
    itnum = str_split(infile, "\\.")[[1]][1]
    # Load the infile and add the itnum
    read_tsv(file.path(fdir, infile)) %>% mutate(itnum = itnum)
  }
))

# Combine gene annotations table with gene score table
fscr = fscr %>%
  inner_join(frag) %>%
  # Add sample name
  inner_join(smit)

# Save combined table
write_tsv(fscr, file.path(idir_name, "fragment_scores.tab"))
