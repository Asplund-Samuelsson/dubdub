#!/usr/bin/env bash

# Define input directory
INDIR=$(realpath $1) # DubDub directory to clean up

# Remove fastq.gz files
echo "Removing fastq.gz files..."
rm -rf ${INDIR}/fastq
echo "Done."

# Compress intermediate output
echo "Compressing intermediate output..."
cd ${INDIR}
tar cf - barseq bstat gscore | gzip > intermediate_data.tar.gz
rm -rf barseq bstat gscore
echo "Done."