#!/usr/bin/env Rscript
options(width=150)
library(tidyverse)

# Load arguments
args = commandArgs(trailingOnly=TRUE)
opts = commandArgs(trailingOnly=FALSE)

# Name of new library
olib = args[1]

# Name of libraries to combine
libs = args[2:length(args)]

# Directory that contains the DubDub pipeline
dubd = dirname(dirname(sub("--file=", "", opts[grep("--file=", opts)])))

# Load the libraries to combine
bpag = bind_rows(lapply(
  libs,
  function (libr) {
    read_tsv(
      file.path(dubd, "libraries", libr, paste(libr, "bpag.tab", sep="."))
    ) %>% mutate(
      up_contig_id = paste(libr, up_contig_id, sep="_"),
      dn_contig_id = paste(libr, dn_contig_id, sep="_")
    )
  }
))

gene = bind_rows(lapply(
  libs,
  function (libr) {
    read_tsv(
      file.path(dubd, "libraries", libr, paste(libr, "gff", sep=".")),
      col_names = c(
        "seqname", "source", "feature", "start", "end",
        "score", "strand", "frame", "attribute"
      ),
      col_types = cols(
        seqname = col_character(),
        source = col_character(),
        feature = col_character(),
        start = col_integer(),
        end = col_integer(),
        score = col_character(),
        strand = col_character(),
        frame = col_integer(),
        attribute = col_character()
      )
    ) %>% mutate(
      lib = libr,
      seqname = paste(libr, seqname, sep="_")
    )
  }
))

# If locus_tag or ID is non-unique, add library name to them
gdup = gene %>%
  mutate(
    locus_tag = str_remove(
      unlist(lapply(str_split(attribute, ";"), "[[", 2)),
      "^locus_tag="
    ),
    ID = str_remove(
      unlist(lapply(str_split(attribute, ";"), "[[", 3)),
      "^ID="
    ),
    product = str_remove(
      unlist(lapply(str_split(attribute, ";"), "[[", 1)),
      "^product=")
    )

gdup = gdup %>%
  group_by(locus_tag) %>%
  mutate(Count = length(lib)) %>%
  ungroup() %>%
  mutate(
    locus_tag = ifelse(Count > 1, paste(locus_tag, lib, sep="_"), locus_tag)
  ) %>%
  group_by(ID) %>%
  mutate(Count = length(lib)) %>%
  ungroup() %>%
  mutate(
    ID = ifelse(Count > 1, paste(ID, lib, sep="_"), ID)
  )

# Reconstruct attribute
gene = gdup %>%
  mutate(
    attribute = paste(
      paste("product=", product, sep=""),
      paste("locus_tag=", locus_tag, sep=""),
      paste("ID=", ID, sep=""),
      sep=";"
    )
  ) %>%
  select(-lib, -locus_tag, -ID, -product, -Count)

# Determine duplicated UP barcodes
dupl = bpag %>%
  group_by(barcode_up) %>%
  summarise(Count = length(barcode_up)) %>%
  filter(Count > 1) %>%
  pull(barcode_up)

# Create new directory for combined library
odir = file.path(dubd, "libraries", olib)
dir.create(odir)

# If there are duplicates, save additional BPAG tables
if (length(dupl) > 0) {

  # Save duplicated barcodes to new library in separate file
  write_tsv(
    filter(bpag, barcode_up %in% dupl),
    file.path(odir, paste(olib, "multi", "bpag.tab", sep="."))
  )

  # Save deduplicated barcodes to new library in separate file
  write_tsv(
    bpag %>%
      group_by(barcode_up) %>%
      slice_head(n=1),
    file.path(odir, paste(olib, "dedup", "bpag.tab", sep="."))
  )

}

# Save unique barcodes to new library
write_tsv(
  filter(bpag, !(barcode_up %in% dupl)),
  file.path(odir, paste(olib, "bpag.tab", sep="."))
)

# Save combined GFF to new library
write_tsv(gene, file.path(odir, paste(olib, "gff", sep=".")), col_names=F)
