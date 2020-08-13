#!/usr/bin/env Rscript

# Load libraries
library(tidyverse)
library(ggrepel)
library(scales)

# Load input directory
args = commandArgs(trailingOnly=TRUE)
idir_name = args[1]

# Load data
gscr_file = file.path(idir_name, "gene_scores.tab")
gscr = read_tsv(gscr_file) %>%
  select(locus_tag, sample, score_cnnls)

# Center the data
gscr_wide = gscr %>%
  spread(sample, score_cnnls) %>%
  as.data.frame()

rownames(gscr_wide) = gscr_wide$locus_tag

gscr_matrix = select(gscr_wide, -locus_tag) %>%
  as.matrix() %>%
  t() %>%
  asinh() %>%
  scale(center=T, scale=F) %>%
  t()

# Perform PCA
gscr_pca = prcomp(gscr_matrix)

# Create plotting dataframes
plot_gscr = as.data.frame(gscr_pca$rotation)
plot_gscr$Sample = rownames(plot_gscr)

# Calculate fraction of variance per PC
var_pc = percent(gscr_pca$sdev^2 / sum(gscr_pca$sdev^2))[1:3]

gp = ggplot(plot_gscr, aes(x=PC1, y=PC2, colour=PC3, label=Sample))
gp = gp + geom_point(aes(size=PC3))
gp = gp + scale_size(range=c(1,3))
gp = gp + scale_colour_gradient2(low="#d0d1e6", mid="#3690c0", high="#014636")
gp = gp + geom_text_repel(force=3, size=4)
gp = gp + labs(
            x=paste("PC1 (", var_pc[1], ")", sep=""),
            y=paste("PC2 (", var_pc[2], ")", sep=""),
            colour=paste("PC3 (", var_pc[3], ")", sep="")
          )
gp = gp + theme_bw()

outfile = file.path(idir_name, "gene_score_PCA.pdf")

ggsave(outfile, gp, height=15/2.54, width=18/2.54)
