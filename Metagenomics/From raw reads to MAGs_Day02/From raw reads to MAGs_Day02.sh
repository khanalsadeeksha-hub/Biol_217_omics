#!/bin/bash
#SBATCH --nodes=1
#SBATCH --cpus-per-task=12
#SBATCH --mem=60G
#SBATCH --time=24:00:00
#SBATCH --job-name=day2_assembly
#SBATCH --output=day2_assembly.out
#SBATCH --error=day2_assembly.err
#SBATCH --partition=base
#SBATCH --reservation=biol217

# ------------------------------------------------------------------------------
# 1. Environment Setup
# ------------------------------------------------------------------------------
# Load required modules and activate the Anvi'o environment for the CAU cluster.
module load gcc12-env/12.1.0
module load micromamba 2> /dev/null
eval "$(micromamba shell hook --shell=bash)"
export MAMBA_ROOT_PREFIX=$WORK/.micromamba

cd $WORK
micromamba activate $WORK/.micromamba/envs/00_anvio/

# Define our working directories dynamically based on the cluster setup
RAW_DIR="$WORK/metagenomics/0_raw_reads"
QC_RAW_DIR="$WORK/metagenomics/1_qc_raw"
CLEAN_DIR="$WORK/metagenomics/2_clean_reads"
QC_CLEAN_DIR="$WORK/metagenomics/1_qc_clean"
ASSEMBLY_DIR="$WORK/metagenomics/3_coassembly"

# Build the output directories
mkdir -p $QC_RAW_DIR $CLEAN_DIR $QC_CLEAN_DIR

echo "Starting Automated Day 2 Workflow..."

# Initialize empty arrays to hold the list of clean files for the MEGAHIT co-assembly
R1_CLEAN_LIST=()
R2_CLEAN_LIST=()

# ------------------------------------------------------------------------------
# 2. Quality Control & Trimming (Automated Loop)
# ------------------------------------------------------------------------------
# This loop finds every forward read in the raw directory and processes it alongside its matching reverse read.
for R1 in $RAW_DIR/*_R1.fastq.gz; do
    
    # Extract the base sample name (e.g., if file is "sampleA_R1.fastq.gz", this isolates "sampleA")
    BASENAME=$(basename $R1 _R1.fastq.gz)
    R2="$RAW_DIR/${BASENAME}_R2.fastq.gz"
    
    echo "Processing sample: $BASENAME"
    
    # FastQC on raw reads
    fastqc $R1 $R2 -o $QC_RAW_DIR/
    
    # Define paths for the clean outputs
    CLEAN_R1="$CLEAN_DIR/${BASENAME}_R1_clean.fastq.gz"
    CLEAN_R2="$CLEAN_DIR/${BASENAME}_R2_clean.fastq.gz"
    
    # Run fastp to trim and filter reads
    fastp -i $R1 \
          -I $R2 \
          -o $CLEAN_R1 \
          -O $CLEAN_R2 \
          -t 6 -q 20 \
          -h $QC_RAW_DIR/${BASENAME}_fastp.html \
          -R "${BASENAME} Fastp Report"
          
    # FastQC on the newly cleaned reads
    fastqc $CLEAN_R1 $CLEAN_R2 -o $QC_CLEAN_DIR/
    
    # Add the clean files to our arrays for the next step
    R1_CLEAN_LIST+=("$CLEAN_R1")
    R2_CLEAN_LIST+=("$CLEAN_R2")
    
done

# ------------------------------------------------------------------------------
# 3. Metagenomic Co-Assembly
# ------------------------------------------------------------------------------
echo "Formatting files for MEGAHIT Co-assembly..."

# MEGAHIT accepts comma-separated lists for multiple inputs. We join our bash arrays here.
R1_COMMA_SEPARATED=$(IFS=,; echo "${R1_CLEAN_LIST[*]}")
R2_COMMA_SEPARATED=$(IFS=,; echo "${R2_CLEAN_LIST[*]}")

# We intentionally do not create the assembly output folder beforehand, as MEGAHIT requires it to not exist.
megahit -1 $R1_COMMA_SEPARATED \
        -2 $R2_COMMA_SEPARATED \
        -o $ASSEMBLY_DIR \
        --min-contig-len 1000 \
        --presets meta-large \
        -m 0.85 \
        -t 12

# ------------------------------------------------------------------------------
# 4. Graph Conversion for Bandage
# ------------------------------------------------------------------------------
echo "Converting FASTA to FASTG for Bandage visualization..."

# Convert the final assembly (using k-mer size 99) to a graph format.
megahit_toolkit contig2fastg 99 $ASSEMBLY_DIR/final.contigs.fa > $ASSEMBLY_DIR/final.contigs.fastg

echo "Pipeline complete. Ready for visualization."