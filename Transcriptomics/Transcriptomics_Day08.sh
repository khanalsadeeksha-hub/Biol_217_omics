#!/bin/bash
#SBATCH --nodes=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=64G
#SBATCH --time=04:00:00
#SBATCH --job-name=day8_rnaseq
#SBATCH --output=day8_rnaseq.out
#SBATCH --error=day8_rnaseq.err
#SBATCH --partition=base
#SBATCH --reservation=biol217

# ------------------------------------------------------------------------------
# 1. Environment & Proxy Setup
# ------------------------------------------------------------------------------
module load gcc12-env/12.1.0
module load micromamba/1.4.2
eval "$(micromamba shell hook --shell=bash)"

# Set proxy environment to allow downloading data from NCBI/SRA on the backend nodes
export http_proxy=http://relay:3128
export https_proxy=http://relay:3128
export ftp_proxy=http://relay:3128

# Define base working directory
WORK_DIR="$WORK/rnaseq"
mkdir -p $WORK_DIR
cd $WORK_DIR

echo "--- Day 8 Transcriptomics Pipeline Started ---"

# ------------------------------------------------------------------------------
# 2. Download Raw RNA-Seq Reads
# ------------------------------------------------------------------------------
echo "Activating grabseqs environment and downloading SRA data..."
micromamba activate $WORK/.micromamba/envs/10_grabseqs

mkdir -p $WORK_DIR/fastq_raw
cd $WORK_DIR/fastq_raw

# Download the sequence data using SRR numbers from Prasse et al. 2017
# We use an empty metadata.csv to satisfy the -m flag requirement.
touch metadata.csv
grabseqs sra -t 4 -m ./metadata.csv SRR4018514
grabseqs sra -t 4 -m ./metadata.csv SRR4018515
grabseqs sra -t 4 -m ./metadata.csv SRR4018516
grabseqs sra -t 4 -m ./metadata.csv SRR4018517

# Rename each SRR file according to the sample name and biological replicate.
# Using * catches any suffixes grabseqs might add.
mv SRR4018514*.fastq.gz wt_R1.fastq.gz
mv SRR4018515*.fastq.gz wt_R2.fastq.gz
mv SRR4018516*.fastq.gz mut_R1.fastq.gz
mv SRR4018517*.fastq.gz mut_R2.fastq.gz

micromamba deactivate

# ------------------------------------------------------------------------------
# 3. Setup READemption Project & Reference Genome
# ------------------------------------------------------------------------------
echo "Activating READemption environment and setting up project..."
micromamba activate $WORK/.micromamba/envs/08_reademption
cd $WORK_DIR

# Create the strictly required folder structure for READemption.
reademption create --project_path READemption_analysis_2 --species methanosarcina="Methanosarcina mazei Gö1"

echo "Downloading reference genome and annotations..."
# Download reference genome sequence and annotation directly into READemption's expected input folders.
wget -O READemption_analysis_2/input/methanosarcina_reference_sequences/GCF_000007065.1_ASM706v1_genomic.fna.gz https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/007/065/GCF_000007065.1_ASM706v1/GCF_000007065.1_ASM706v1_genomic.fna.gz
wget -O READemption_analysis_2/input/methanosarcina_annotations/GCF_000007065.1_ASM706v1_genomic.gff.gz https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/007/065/GCF_000007065.1_ASM706v1/GCF_000007065.1_ASM706v1_genomic.gff.gz

# Unzip them as READemption requires plain text fasta/gff files.
gunzip READemption_analysis_2/input/methanosarcina_reference_sequences/GCF_000007065.1_ASM706v1_genomic.fna.gz
gunzip READemption_analysis_2/input/methanosarcina_annotations/GCF_000007065.1_ASM706v1_genomic.gff.gz

# ------------------------------------------------------------------------------
# 4. Move Reads & Run the READemption Pipeline
# ------------------------------------------------------------------------------
echo "Moving reads into READemption input folder..."
# The reads do not need to be unzipped here, but they must be in the input/reads folder.
mv $WORK_DIR/fastq_raw/*.fastq.gz READemption_analysis_2/input/reads/

echo "Step 4.1: Aligning reads to reference..."
# The --fastq flag is required because our inputs are .fastq.gz, not standard .fa files.
reademption align -p 32 --poly_a_clipping --project_path READemption_analysis_2 --fastq

echo "Step 4.2: Calculating read coverage..."
reademption coverage -p 32 --project_path READemption_analysis_2

echo "Step 4.3: Quantifying gene expression..."
reademption gene_quanti -p 32 --features CDS,tRNA,rRNA --project_path READemption_analysis_2

echo "Step 4.4: Calculating differential expression (DESeq2)..."
# Setting up DESeq2 contrasts. Note: -r 1,2,1,2 indicates replicates 1 and 2 for mut, and 1 and 2 for wt.
reademption deseq \
    -l mut_R1.fastq.gz,mut_R2.fastq.gz,wt_R1.fastq.gz,wt_R2.fastq.gz \
    -c mut,mut,wt,wt \
    -r 1,2,1,2 \
    --libs_by_species methanosarcina=mut_R1,mut_R2,wt_R1,wt_R2 \
    --project_path READemption_analysis_2

echo "Step 4.5: Generating visual reports..."
reademption viz_align --project_path READemption_analysis_2
reademption viz_gene_quanti --project_path READemption_analysis_2
reademption viz_deseq --project_path READemption_analysis_2

# ------------------------------------------------------------------------------
# 5. Cleanup
# ------------------------------------------------------------------------------
micromamba deactivate
module purge
jobinfo

echo "Day 8 Pipeline Complete! Check the READemption_analysis_2 folder for your results."