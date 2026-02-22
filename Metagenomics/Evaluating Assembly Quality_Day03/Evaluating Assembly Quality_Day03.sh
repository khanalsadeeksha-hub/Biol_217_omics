#!/bin/bash
#SBATCH --nodes=1
#SBATCH --cpus-per-task=12
#SBATCH --mem=60G
#SBATCH --time=48:00:00  # Increased time as mapping and profiling take a while
#SBATCH --job-name=day3_binning
#SBATCH --output=day3_binning.out
#SBATCH --error=day3_binning.err
#SBATCH --partition=base
#SBATCH --reservation=biol217

# ------------------------------------------------------------------------------
# 0. Environment Setup
# ------------------------------------------------------------------------------
module load gcc12-env/12.1.0
module load micromamba 2> /dev/null
eval "$(micromamba shell hook --shell=bash)"
export MAMBA_ROOT_PREFIX=$WORK/.micromamba

cd $WORK
micromamba activate $WORK/.micromamba/envs/00_anvio/

# Define directories
WORK_DIR="$WORK/metagenomics"
CLEAN_DIR="$WORK_DIR/2_clean_reads"
ASSEMBLY="$WORK_DIR/3_coassembly/final.contigs.fa"
QUAST_DIR="$WORK_DIR/4_quast_eval"
MAPPING_DIR="$WORK_DIR/5_mapping"
ANVIO_DIR="$WORK_DIR/6_anvio_profiles"
BINNING_DIR="$WORK_DIR/7_binning"

mkdir -p $QUAST_DIR $MAPPING_DIR $ANVIO_DIR $BINNING_DIR

echo "Starting Automated Day 3 Workflow..."

# ------------------------------------------------------------------------------
# 1. Evaluating Assembly Quality
# ------------------------------------------------------------------------------
echo "Running MetaQUAST..."
# metaquast assesses assembly quality, ignoring contigs < 1000bp.
metaquast $ASSEMBLY -o $QUAST_DIR/ -t 12 -m 1000

# ------------------------------------------------------------------------------
# 2. Preparing Contigs for Mapping
# ------------------------------------------------------------------------------
echo "Reformatting and Indexing Contigs..."
CONTIGS_ANVIO="$MAPPING_DIR/contigs.anvio.fa"

# Simplify sequence IDs and filter out short contigs (< 1000 nts).
anvi-script-reformat-fasta $ASSEMBLY -o $CONTIGS_ANVIO --min-len 1000 --simplify-names --report-file $MAPPING_DIR/name_conversion.txt

# Index the reformatted contigs to drastically speed up mapping.
bowtie2-build $CONTIGS_ANVIO $MAPPING_DIR/contigs.anvio.fa.index

# ------------------------------------------------------------------------------
# 3. Anvi'o Contigs Database Generation
# ------------------------------------------------------------------------------
echo "Generating Anvi'o Contigs Database..."
CONTIGS_DB="$ANVIO_DIR/contigs.db"

# Convert the fasta into an anvi'o database, calculating k-mer frequencies and identifying ORFs with Prodigal.
anvi-gen-contigs-database -f $CONTIGS_ANVIO -o $CONTIGS_DB -n biol217

# Perform HMM search against gene collections to annotate potential biological functions.
anvi-run-hmms -c $CONTIGS_DB --num-threads 12

# ------------------------------------------------------------------------------
# 4. Read Mapping and Profiling (Automated Loop)
# ------------------------------------------------------------------------------
# We will collect the profile paths in an array to easily merge them later.
PROFILE_LIST=()

for CLEAN_R1 in $CLEAN_DIR/*_R1_clean.fastq.gz; do
    
    # Extract the base sample name dynamically
    BASENAME=$(basename $CLEAN_R1 _R1_clean.fastq.gz)
    CLEAN_R2="$CLEAN_DIR/${BASENAME}_R2_clean.fastq.gz"
    
    echo "Mapping and Profiling sample: $BASENAME"
    
    SAM_FILE="$MAPPING_DIR/${BASENAME}.sam"
    BAM_FILE="$MAPPING_DIR/${BASENAME}.bam"
    SORTED_BAM="$MAPPING_DIR/${BASENAME}_sorted.bam"
    PROFILE_OUT="$ANVIO_DIR/${BASENAME}_profile"
    
    # Map reads to the indexed contigs using the 'very fast' mode.
    bowtie2 -1 $CLEAN_R1 -2 $CLEAN_R2 -x $MAPPING_DIR/contigs.anvio.fa.index -S $SAM_FILE --very-fast
    
    # Convert SAM to binary format (BAM) for faster machine processing.
    samtools view -Sb $SAM_FILE > $BAM_FILE
    
    # Sort and index the BAM file.
    anvi-init-bam $BAM_FILE -o $SORTED_BAM
    
    # Delete the large intermediate SAM and unsorted BAM files to save cluster storage space
    rm $SAM_FILE $BAM_FILE
    
    # Create the anvi'o profile for this sample.
    anvi-profile -i $SORTED_BAM -c $CONTIGS_DB --output-dir $PROFILE_OUT
    
    # Add the newly created PROFILE.db to our list for merging
    PROFILE_LIST+=("$PROFILE_OUT/PROFILE.db")
    
done


# ------------------------------------------------------------------------------
# 5. Merging Profiles
# ------------------------------------------------------------------------------
echo "Merging Anvi'o Profiles..."
MERGED_DIR="$ANVIO_DIR/merged_profiles"

# Merge profiles from all samples to analyze and compare them together.
anvi-merge ${PROFILE_LIST[@]} -o $MERGED_DIR -c $CONTIGS_DB --enforce-hierarchical-clustering

# ------------------------------------------------------------------------------
# 6. Binning Contigs into MAGs
# ------------------------------------------------------------------------------
echo "Running MetaBAT2 and MaxBin2..."
# Remember: We do not create the summary output folders beforehand, Anvi'o will do it.

# MetaBAT2 workflow
anvi-cluster-contigs -p $MERGED_DIR/PROFILE.db -c $CONTIGS_DB -C METABAT2 --driver metabat2 --just-do-it --log-file $BINNING_DIR/metabat2.log
anvi-summarize -p $MERGED_DIR/PROFILE.db -c $CONTIGS_DB -o $BINNING_DIR/SUMMARY_METABAT2 -C METABAT2

# MaxBin2 workflow
anvi-cluster-contigs -p $MERGED_DIR/PROFILE.db -c $CONTIGS_DB -C MAXBIN2 --driver maxbin2 --just-do-it --log-file $BINNING_DIR/maxbin2.log
anvi-summarize -p $MERGED_DIR/PROFILE.db -c $CONTIGS_DB -o $BINNING_DIR/SUMMARY_MAXBIN2 -C MAXBIN2



echo "Day 3 Pipeline Complete!"
echo "Check $QUAST_DIR for your N50 and assembly stats."
echo "Check $BINNING_DIR/SUMMARY_METABAT2/index.html to see your Archaea bins!"