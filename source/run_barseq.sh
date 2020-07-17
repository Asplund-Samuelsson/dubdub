#!/usr/bin/env bash

# Define input and output directories
INDIR=$(realpath $1) # Directory with fastq.gz files
OUTDIR=$(realpath $2) # Directory with barseq output

# Define pre- and post-sequences flanking the barcodes left and right
LSEQ=$3 # 9 nt perfectly matching the sequence upstream (left) of the barcode
LPOS=$4 # Position of LSEQ start
RSEQ=$5 # 9 nt perfectly matching the sequence downstream (right) of the barcode

# Create output directory
mkdir $OUTDIR

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
