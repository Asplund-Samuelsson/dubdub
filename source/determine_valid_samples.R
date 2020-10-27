#!/usr/bin/env Rscript
options(width=150)
library(tidyverse)

# Load input directory name, bpag file name, and layout file name
args = commandArgs(trailingOnly=TRUE)
idir_name = args[1]
bpag_name = args[2]
layo_name = args[3]

# Determine if bpag is dedup and name output file accordingly
if (endsWith(bpag_name, ".dedup.bpag.tab")){
  vals_name = file.path(idir_name, "valid_samples.dedup.txt")
} else {
  vals_name = file.path(idir_name, "valid_samples.txt")
}

# Load bpag table
bpag = read_tsv(bpag_name)

# Load layout
layo = read_tsv(layo_name)

# Determine samples and file names
bdir = file.path(idir_name, "barseq")
bars_files = list.files(bdir, pattern="*bstat.tsv", recursive=T)
samples = unlist(lapply(str_split(bars_files, "/"), "[[", 1))

# Load all sample barcode counts
bars = bind_rows(lapply(
  1:length(samples),
  function(x){
    bars_file = file.path(bdir, bars_files[x])
    sample = samples[x]
    read_tsv(bars_file) %>% mutate(Sample = sample)
  }
))

# Filter barcodes to those recommended in bpag and at least 10 reads in time0
vals = bars %>%
  # Remove barcodes that are not recommended in the bpag table
  filter(barcode %in% filter(bpag, recommended == "+")$barcode_up) %>%
  # Remove barcodes that fail to reach 10 reads in time0
  filter(
    barcode %in% (
      bars %>%
        filter(Sample %in% filter(layo, type == "Time0")$sample) %>%
        filter(reads_count >= 10) %>%
        pull(barcode) %>%
        unique()
    )
  ) %>%
  # Select unique barcodes
  pull(barcode) %>%
  unique()

# Filter samples to those with at least one valid barcode
bars %>%
  filter(barcode %in% vals) %>%
  pull(Sample) %>%
  unique() %>%
  # Write valid samples to file
  writeLines(vals_name)
