#!/usr/bin/env Rscript
options(width=150)
library(tidyverse)

# Load input directory
args = commandArgs(trailingOnly=TRUE)
idir_name = args[1]

# Define infiles
smit_file = file.path(idir_name, "sample_itnum.tab")
gene_file = file.path(idir_name, "gscore", "gscore_base.tsv")
gdir = file.path(idir_name, "gscore")
gscr_files = grep("IT.*gscore.tsv", list.files(gdir), value=T)

# Load data
smit = read_tsv(smit_file, col_names = c("sample", "itnum"))
gene = read_tsv(gene_file)
gscr = bind_rows(lapply(
  gscr_files,
  # Load each gene score file
  function (infile) {
    # Extract the itnum from the filename
    itnum = str_split(infile, "\\.")[[1]][1]
    # Load the infile and add the itnum
    read_tsv(file.path(gdir, infile)) %>% mutate(itnum = itnum)
  }
))

# Combine gene annotations table with gene score table
gscr = gscr %>%
  inner_join(gene) %>%
  # Remove superfluous columns
  select(-index, -gene_index, -gene_name) %>%
  # Add sample name
  inner_join(smit) %>%
  # Reorder columns
  select(sample, itnum, contig_id, locus_tag, everything())

# Save combined table
write_tsv(gscr, file.path(idir_name, "gene_scores.tab"))
