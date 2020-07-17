#!/usr/bin/env bash

# Find dubdub directory
DUBDIR=`dirname $0`

# Load configuration
source ${DUBDIR}/config.sh

# Read input flags
while getopts d:f:l:L:U:P:D: flag
	do
	    case "${flag}" in
	        d) WORKDIR=${OPTARG};;  # Working directory in which results appear
	        f) FASTQDIR=${OPTARG};; # Directory with input fastq files
	        l) LAYOUT=${OPTARG};;   # Barseq experiment layout file
          L) LIBRARY=${OPTARG};;  # DubSeq library to use
          U) LSEQ=${OPTARG};;     # If necessary, supply alternative LSEQ
          P) LPOS=${OPTARG};;     # If necessary, supply alternative LPOS
          D) RSEQ=${OPTARG};;     # If necessary, supply alternative RSEQ
	    esac
	done

# Remove trailing slash from inputs
WORKDIR=`echo $WORKDIR | sed -e 's|/$||'`
FASTQDIR=`echo $FASTQDIR | sed -e 's|/$||'`

# Create working directory
mkdir $WORKDIR

# Set step counter to zero
S=0

# Concatenate FASTQ.gz files
((S++))
echo -e "\n\e[94mStep $S: Concatenating fastq.gz files...\e[0m\n"
${DUBDIR}/source/concatenate_fastq.sh $FASTQDIR ${WORKDIR}/fastq
echo -e "\n\e[92mStep $S: Done.\e[0m\n"

# Count barcodes
((S++))
echo -e "\n\e[94mStep $S: Counting barcodes...\e[0m\n"
${DUBDIR}/source/run_barseq.sh ${WORKDIR}/fastq ${WORKDIR}/barseq \
$LSEQ $LPOS $RSEQ
echo -e "\n\e[92mStep $S: Done.\e[0m\n"

# Create links to bstat files in new directory
((S++))
echo -e "\n\e[94mStep $S: Linking barseq statistics files...\e[0m\n"
${DUBDIR}/source/link_bstat_files_with_itnums.sh $WORKDIR
echo -e "\n\e[92mStep $S: Done.\e[0m\n"

# Add IT numbering (itnum) to layout
((S++))
echo -e "\n\e[94mStep $S: Adding IT numbers to layout...\e[0m\n"
${DUBDIR}/source/add_itnum_to_layout.R $WORKDIR $LAYOUT
echo -e "\n\e[92mStep $S: Done.\e[0m\n"

# Calculate fragment fitness values
((S++))
echo -e "\n\e[94mStep $S: Calculating fragment fitness values...\e[0m\n"
${DUBDIR}/source/run_fscore.sh ${WORKDIR}/bstat ${WORKDIR}/fscore \
$LIBRARY ${WORKDIR}/layout.tab
echo -e "\n\e[92mStep $S: Done.\e[0m\n"

# Calculate gene fitness values
# ...