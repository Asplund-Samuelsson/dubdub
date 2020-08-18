#!/usr/bin/env bash

# Define input and output directories
INDIR=$(realpath $1) # Directory with bstat files
OUTDIR=$(realpath $2) # Directory with fscore output

# Define layout and library files
LIBRARY=$3 # The library to use
LAYOUT=$(realpath $4) # Barseq layout file

# Create fscore output directory
mkdir $OUTDIR

# Define bpag infile; "file with barcode pairs mapped to a genome"
BPAG="`dirname $(realpath $0)`/../libraries/${LIBRARY}/${LIBRARY}.dedup.bpag.tab"

# Change to dubseq directory in order to run scripts there
cd `dirname $0`/../DubSeq/dubseq

# Run fscore.py
./fscore.py -i $INDIR -l $LAYOUT -p $BPAG -o $OUTDIR
