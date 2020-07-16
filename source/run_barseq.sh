#!/usr/bin/env bash

# Define input and output directories
INDIR=$1 # Directory with fastq.gz files
OUTDIR=$2 # Directory with barseq output

# Define pre- and post-sequences flanking the barcodes left and right
LSEQ=$3 # "ACGGATCCT"
LPOS=$4 # "16"
RSEQ=$5 # "ACTAAACAT"

# Change to dubseq directory in order to run scripts there
cd `dirname $0`/../DubSeq/dubseq

# List all input files
ls ${INDIR}/* | while read INFILE; do
  # Determine sample
  SAMPLE=`echo "$INFILE" | rev | cut -f 1 -d / | cut -f 3 -d . | rev`
  # Create individual output directories
  mkdir ${OUTDIR}/${SAMPLE}
  # Create barseq.py command for each infile, targeting one output directory
  echo "./barseq.py -l $LSEQ -p $LPOS -r $RSEQ -i $INFILE -o ${OUTDIR}/${SAMPLE}"
# Run barseq.py in parallel for each input file
done | parallel --no-notice
