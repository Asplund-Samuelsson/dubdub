#!/usr/bin/env bash

# Define input directory
INDIR=$1 # DubDub analysis results directory

# Create bstat folder
mkdir ${INDIR}/bstat

# Reset counter for itnum
N=0

# For each sample with a directory in the barseq directory
ls ${INDIR}/barseq/ | while read SAMPLE; do
  # Increment counter for itnum
  ((N++))
  # Create itnum string
  ITNUM=`printf "IT%03d" $N`
  # Determine source and target for link
  SOURCE=$(realpath ${INDIR}/barseq/${SAMPLE}/${SAMPLE}.fastq.gz.bstat.tsv)
  TARGET=$(realpath ${INDIR}/bstat/${ITNUM}.bstat.tsv)
  # Create link from source to target
  ln -s $SOURCE $TARGET
  # Save the sample-itnum association
  echo -e "${SAMPLE}\t${ITNUM}" >> ${INDIR}/sample_itnum.tab
done
