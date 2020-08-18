#!/usr/bin/env Rscript
options(width=150)
library(tidyverse)

# Load file arguments
args = commandArgs(trailingOnly=TRUE)
fscr_file = args[1]
mult_file = args[2]
gene_file = args[3]

# Load data
fscr = read_tsv(fscr_file)
mult = read_tsv(mult_file)
gene = read_tsv(
  gene_file,
  col_names = c(
    "seqname", "source", "feature", "start", "end",
    "score", "strand", "frame", "attribute"
  )
)

# Add fragment counts
fscr = mult %>%
  select(barcode_up) %>%
  group_by(barcode_up) %>%
  summarise(fragment_count = length(barcode_up)) %>%
  rename(barcode = barcode_up) %>%
  right_join(fscr) %>%
  mutate(fragment_count = replace_na(fragment_count, 1))

# Add missing fragments
fscr = mult %>%
  rename(barcode = barcode_up, pos_to = pos_end, contig_id = up_contig_id) %>%
  select(barcode, pos_from, pos_to, contig_id) %>%
  inner_join(select(fscr, -contig_id, -pos_from, -pos_to)) %>%
  bind_rows(filter(fscr, !(barcode %in% mult$barcode_up)))

# Add locus_tag to gene
gene = gene %>%
  mutate(
    locus_tag = unlist(lapply(
      attribute,
      function(x){
        str_remove(
          grep("locus_tag=", unlist(str_split(x, ";")), value=T), "locus_tag="
        )
      }
    ))
  )

# Extract fragments
frag = fscr %>%
  select(barcode, contig_id, pos_from, pos_to) %>%
  distinct()

# Add genes to fragment scores
fscr = gene %>%
  select(seqname, start, end, locus_tag) %>%
  rename(contig_id = seqname) %>%
  inner_join(frag) %>%
  filter(start >= pos_from, end <= pos_to) %>%
  rename(gene_start = start, gene_end = end) %>%
  right_join(fscr)

# Write new fragment score table
write_tsv(fscr, fscr_file)
