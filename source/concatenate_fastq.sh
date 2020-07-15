#!/usr/bin/env bash

# Define input and output directories
INDIR=$1 # Directory with fastq.gz files
OUTDIR=$2 # Directory with concatenated output fastq.gz files for each sample

# Create associative array to store lists of files
declare -A INFILES

# For each input file
while read INFILE; do
  # Cut out the sample ID
  SAMPLE=`echo $INFILE | rev | cut -f 1 -d / | cut -f 5- -d _ | rev`;
  # Define the sample outfile
  OUTFILE="${OUTDIR}/${SAMPLE}.fastq";
  # Make list of individual infiles in a string in the associative array
  INFILES[$OUTFILE]+="$INFILE ";
# List input files (has to be done this end of the loop...)
done < <(find $INDIR -name *.fastq.gz)

# List the outfiles, which are keys in the associative array
for OUTFILE in "${!INFILES[@]}"; do
  # Create command to concatenate infiles to the correct outfiles and then gzip
  echo "zcat ${INFILES[$OUTFILE]}> $OUTFILE; gzip $OUTFILE";
# Pipe the commands to GNU parallel for faster processing
done | parallel --no-notice
