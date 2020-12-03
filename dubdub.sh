#!/usr/bin/env bash

# Find dubdub directory
DUBDIR=`dirname $0`

# Load configuration
source ${DUBDIR}/config.sh

# Read input flags
while getopts d:f:l:L:U:P:D:t: flag
	do
	    case "${flag}" in
	        d) WORKDIR=${OPTARG};;  # Working directory in which results appear
	        f) FASTQDIR=${OPTARG};; # Directory with input fastq files
	        l) LAYOUT=${OPTARG};;   # Barseq experiment layout file
          L) LIBRARY=${OPTARG};;  # DubSeq library to use
          U) LSEQ=${OPTARG};;     # If necessary, supply alternative LSEQ
          P) LPOS=${OPTARG};;     # If necessary, supply alternative LPOS
          D) RSEQ=${OPTARG};;     # If necessary, supply alternative RSEQ
					t) ZEROCUT=${OPTARG};;  # If necessary, supply time0 read count cutoff
	    esac
	done

# Remove trailing slash from inputs
WORKDIR=`echo $WORKDIR | sed -e 's|/$||'`
FASTQDIR=`echo $FASTQDIR | sed -e 's|/$||'`

# Create working directory if it does not exist
[ ! -d "$WORKDIR" ] && mkdir -p $WORKDIR

# Exit if any of the required inputs do not exist
if [ ! -d "$WORKDIR" ]; then
	echo -e "\n\e[31mError: Output directory does not exist.\e[0m\n"; exit 1
fi

if [ ! -d "$FASTQDIR" ]; then
	echo -e "\n\e[31mError: FASTQ directory does not exist.\e[0m\n"; exit 1
fi

if [ ! -f "$LAYOUT" ]; then
	echo -e "\n\e[31mError: Barseq layout file does not exist.\e[0m\n"; exit 1
fi

BPAG="`dirname $(realpath $0)`/libraries/${LIBRARY}/${LIBRARY}.bpag.tab"
if [ ! -f "$BPAG" ]; then
	echo -e "\n\e[31mError: DubSeq library (BPAG) does not exist.\e[0m\n"; exit 1
fi

# Exit if working directory is not empty
if [ "$(ls -A $WORKDIR)" ]; then
	echo -e "\n\e[31mError: Output directory is not empty.\e[0m\n"; exit 1
fi

# Set step counter to zero
S=0

# Count total number of steps
T=`grep -cP "^\(\(S\+\+\)\)$" $0`

# Concatenate FASTQ.gz files
((S++))

echo -e "\n\e[94mStep $S of $T: Concatenating fastq.gz files...\e[0m\n"

${DUBDIR}/source/concatenate_fastq.sh $FASTQDIR ${WORKDIR}/fastq

FILECOUNT=`find ${WORKDIR}/fastq -name *fastq.gz | wc -l`
if [ "$FILECOUNT" -lt "1" ]; then
	echo -e "\n\e[31mError: Step $S failed.\e[0m\n"; exit 1
fi
FILECOUNT="0"

echo -e "\n\e[92mStep $S of $T: Done.\e[0m\n"

# Count barcodes
((S++))

echo -e "\n\e[94mStep $S of $T: Counting barcodes...\e[0m\n"

${DUBDIR}/source/run_barseq.sh ${WORKDIR}/fastq ${WORKDIR}/barseq \
$LSEQ $LPOS $RSEQ

FILECOUNT=`find ${WORKDIR}/barseq -name *bstat.tsv | wc -l`
if [ "$FILECOUNT" -lt "1" ]; then
	echo -e "\n\e[31mError: Step $S failed.\e[0m\n"; exit 1
fi
FILECOUNT="0"

echo -e "\n\e[92mStep $S of $T: Done.\e[0m\n"

# Create links to bstat files in new directory
((S++))

echo -e "\n\e[94mStep $S of $T: Linking barseq statistics files...\e[0m\n"
${DUBDIR}/source/link_bstat_files_with_itnums.sh $WORKDIR

FILECOUNT=`find ${WORKDIR}/bstat -name *bstat.tsv | wc -l`
if [ "$FILECOUNT" -lt "1" ]; then
	echo -e "\n\e[31mError: Step $S failed.\e[0m\n"; exit 1
fi
FILECOUNT="0"

echo -e "\n\e[92mStep $S of $T: Done.\e[0m\n"

# Add IT numbering (itnum) to layout
((S++))

echo -e "\n\e[94mStep $S of $T: Adding IT numbers to layout...\e[0m\n"

${DUBDIR}/source/add_itnum_to_layout.R $WORKDIR $LAYOUT

FILECOUNT=`find ${WORKDIR} -name layout.tab | wc -l`
if [ "$FILECOUNT" -lt "1" ]; then
	echo -e "\n\e[31mError: Step $S failed.\e[0m\n"; exit 1
fi
FILECOUNT="0"

echo -e "\n\e[92mStep $S of $T: Done.\e[0m\n"

# Filter layout file to valid samples (has >= 10 t0 reads recommended barcodes)
((S++))

echo -e "\n\e[94mStep $S of $T: Filtering layout to valid samples...\e[0m\n"

${DUBDIR}/source/determine_valid_samples.R $WORKDIR $BPAG $LAYOUT

${DUBDIR}/source/filter_layout_table_to_valid_samples.R \
${WORKDIR}/valid_samples.txt ${WORKDIR}/sample_itnum.tab ${WORKDIR}/layout.tab

FILECOUNT=`find ${WORKDIR} -name layout.valid.tab | wc -l`
if [ "$FILECOUNT" -lt "1" ]; then
	echo -e "\n\e[31mError: Step $S failed.\e[0m\n"; exit 1
fi
FILECOUNT="0"

echo -e "\n\e[92mStep $S of $T: Done.\e[0m\n"

# Perform fscore calculation if deduplicated barcodes are present in library
DEDUP="`dirname $(realpath $0)`/libraries/${LIBRARY}/${LIBRARY}.dedup.bpag.tab"
if [ -f "$DEDUP" ]; then
	echo -e "\n\e[95mDedup detour: Filtering layout to valid samples...\e[0m\n"
	${DUBDIR}/source/determine_valid_samples.R $WORKDIR $DEDUP $LAYOUT
	${DUBDIR}/source/filter_layout_table_to_valid_samples.R \
	${WORKDIR}/valid_samples.dedup.txt ${WORKDIR}/sample_itnum.tab \
	${WORKDIR}/layout.tab
	echo -e "\n\e[95mDedup detour: Calculating fragment fitness values...\e[0m\n"
	${DUBDIR}/source/run_fscore.sh ${WORKDIR}/bstat ${WORKDIR}/fscore \
	$LIBRARY ${WORKDIR}/layout.valid.dedup.tab $ZEROCUT
	echo -e "\n\e[95mDedup detour: Creating combined fragment score table...\e[0m\n"
	${DUBDIR}/source/create_fscore_table.R $WORKDIR
	echo -e "\n\e[95mDedup detour: Adding genes and multi-barcodes...\e[0m\n"
	${DUBDIR}/source/finalize_fscore_table.R ${WORKDIR}/fragment_scores.tab \
	${DUBDIR}/libraries/${LIBRARY}/${LIBRARY}.multi.bpag.tab \
	${DUBDIR}/libraries/${LIBRARY}/${LIBRARY}.gff
	echo -e "\n\e[93mDedup detour: Done.\e[0m\n"
fi

# Calculate gene fitness values
((S++))

echo -e "\n\e[94mStep $S of $T: Calculating gene fitness values...\e[0m\n"

${DUBDIR}/source/run_gscore.sh ${WORKDIR}/bstat ${WORKDIR}/gscore \
$LIBRARY ${WORKDIR}/layout.valid.tab $ZEROCUT

FILECOUNT=`find ${WORKDIR}/gscore -name *gscore.tsv | wc -l`
if [ "$FILECOUNT" -lt "1" ]; then
	echo -e "\n\e[31mError: Step $S failed.\e[0m\n"; exit 1
fi
FILECOUNT="0"

echo -e "\n\e[92mStep $S of $T: Done.\e[0m\n"

# Create a combined gene score table
((S++))

echo -e "\n\e[94mStep $S of $T: Creating combined gene score table...\e[0m\n"

${DUBDIR}/source/create_gscore_table.R $WORKDIR

FILECOUNT=`find ${WORKDIR} -name gene_scores.tab | wc -l`
if [ "$FILECOUNT" -lt "1" ]; then
	echo -e "\n\e[31mError: Step $S failed.\e[0m\n"; exit 1
fi
FILECOUNT="0"

echo -e "\n\e[92mStep $S of $T: Done.\e[0m\n"

# Perform principal component analysis on gene scores
((S++))

echo -e "\n\e[94mStep $S of $T: Performing PCA on gene scores...\e[0m\n"

${DUBDIR}/source/gscore_pca.R $WORKDIR

FILECOUNT=`find ${WORKDIR} -name gene_score_PCA.pdf | wc -l`
if [ "$FILECOUNT" -lt "1" ]; then
	echo -e "\n\e[31mError: Step $S failed.\e[0m\n"; exit 1
fi
FILECOUNT="0"

echo -e "\n\e[92mStep $S of $T: Done.\e[0m\n"

# Clean up the working directory
((S++))

echo -e "\n\e[94mStep $S of $T: Cleaning up output directory...\e[0m\n"

${DUBDIR}/source/cleanup.sh $WORKDIR

echo -e "\n\e[92mStep $S of $T: Done.\e[0m\n"
