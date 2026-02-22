#!/bin/bash
#SBATCH --nodes=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=12:00:00
#SBATCH --job-name=day6_genomics
#SBATCH --output=day6_genomics.out
#SBATCH --error=day6_genomics.err
#SBATCH --partition=base
#SBATCH --reservation=biol217

# ------------------------------------------------------------------------------
# 0. Environment Setup
# ------------------------------------------------------------------------------
# Loading the base environment modules required by the cluster.
module load gcc12-env/12.1.0
module load miniconda3/4.12.0
module load micromamba/1.4.2
eval "$(micromamba shell hook --shell=bash)"

# Setting up my main working directories
GENOMICS_DIR="$WORK/genomics"
RAW_DIR="$GENOMICS_DIR/0_raw_reads"

echo "Starting Day 6 Genomics Pipeline..."

# ------------------------------------------------------------------------------
# 1. Short Read QC & Trimming
# ------------------------------------------------------------------------------
echo "--- Starting Short Read QC ---"
micromamba activate .micromamba/envs/01_short_reads_qc

# Making directories for raw FastQC, cleaned reads, and clean FastQC.
mkdir -p $GENOMICS_DIR/1_short_reads_qc/1_fastqc_raw
mkdir -p $GENOMICS_DIR/1_short_reads_qc/2_cleaned_reads
mkdir -p $GENOMICS_DIR/1_short_reads_qc/3_fastqc_cleaned

# Finding the raw short reads. Assuming they are in the base raw_reads folder.
# We will just grab the first R1 file we find to process the single isolate.
SHORT_R1=$(ls $RAW_DIR/*_R1.fastq.gz | head -n 1)
BASENAME=$(basename $SHORT_R1 _R1.fastq.gz)
SHORT_R2="$RAW_DIR/${BASENAME}_R2.fastq.gz"

# 1.1 FastQC on raw short reads
fastqc $SHORT_R1 $SHORT_R2 -o $GENOMICS_DIR/1_short_reads_qc/1_fastqc_raw -t 16

# 1.2 fastp for cleaning and adapter trimming
CLEAN_R1="$GENOMICS_DIR/1_short_reads_qc/2_cleaned_reads/${BASENAME}_R1_clean.fastq.gz"
CLEAN_R2="$GENOMICS_DIR/1_short_reads_qc/2_cleaned_reads/${BASENAME}_R2_clean.fastq.gz"

fastp -i $SHORT_R1 -I $SHORT_R2 \
      -o $CLEAN_R1 -O $CLEAN_R2 \
      -R "$GENOMICS_DIR/1_short_reads_qc/fastp_report" \
      -h "$GENOMICS_DIR/1_short_reads_qc/report.html" \
      -t 6 -q 25

# 1.3 FastQC on cleaned short reads
fastqc $CLEAN_R1 $CLEAN_R2 -o $GENOMICS_DIR/1_short_reads_qc/3_fastqc_cleaned -t 16


micromamba deactivate

# ------------------------------------------------------------------------------
# 2. Long Read QC & Filtering
# ------------------------------------------------------------------------------
echo "--- Starting Long Read QC ---"
micromamba activate .micromamba/envs/02_long_reads_qc

LONG_RAW_DIR="$RAW_DIR/long_reads"
mkdir -p $GENOMICS_DIR/2_long_reads_qc/1_nanoplot_raw
mkdir -p $GENOMICS_DIR/2_long_reads_qc/2_cleaned_reads
mkdir -p $GENOMICS_DIR/2_long_reads_qc/3_nanoplot_cleaned

LONG_RAW=$(ls $LONG_RAW_DIR/*.gz | head -n 1)
LONG_BASE=$(basename $LONG_RAW .fastq.gz)
CLEAN_LONG="$GENOMICS_DIR/2_long_reads_qc/2_cleaned_reads/${LONG_BASE}_cleaned_filtlong.fastq.gz"

# 2.1 NanoPlot on raw long reads
NanoPlot --fastq $LONG_RAW -o $GENOMICS_DIR/2_long_reads_qc/1_nanoplot_raw -t 16 \
         --maxlength 40000 --minlength 1000 --plots kde --format png \
         --N50 --dpi 300 --store --raw --tsv_stats --info_in_report

# 2.2 Filtlong to keep the best 90% of reads over 1000bp
filtlong --min_length 1000 --keep_percent 90 $LONG_RAW | gzip > $CLEAN_LONG

# 2.3 NanoPlot on cleaned long reads
NanoPlot --fastq $CLEAN_LONG -o $GENOMICS_DIR/2_long_reads_qc/3_nanoplot_cleaned -t 16 \
         --maxlength 40000 --minlength 1000 --plots kde --format png \
         --N50 --dpi 300 --store --raw --tsv_stats --info_in_report


micromamba deactivate

# ------------------------------------------------------------------------------
# 3. Hybrid Assembly
# ------------------------------------------------------------------------------
echo "--- Starting Unicycler Hybrid Assembly ---"
micromamba activate .micromamba/envs/03_unicycler

ASSEMBLY_DIR="$GENOMICS_DIR/3_hybrid_assembly"
mkdir -p $ASSEMBLY_DIR

# Using Unicycler to combine the accuracy of short reads with the structure of long reads.
unicycler -1 $CLEAN_R1 -2 $CLEAN_R2 -l $CLEAN_LONG -o $ASSEMBLY_DIR -t 16

micromamba deactivate

# ------------------------------------------------------------------------------
# 4. Assembly Quality Checks
# ------------------------------------------------------------------------------
echo "--- Starting Assembly Quality Checks ---"

ASSEMBLY_FASTA="$ASSEMBLY_DIR/assembly.fasta"

# 4.1 QUAST
micromamba activate .micromamba/envs/04_quast
mkdir -p $ASSEMBLY_DIR/quast
quast.py $ASSEMBLY_FASTA --circos -L --conserved-genes-finding --rna-finding \
         --glimmer --use-all-alignments --report-all-metrics -o $ASSEMBLY_DIR/quast -t 16
micromamba deactivate

# 4.2 CheckM
micromamba activate .micromamba/envs/04_checkm
mkdir -p $ASSEMBLY_DIR/checkm
checkm lineage_wf $ASSEMBLY_DIR $ASSEMBLY_DIR/checkm -x fasta --tab_table --file $ASSEMBLY_DIR/checkm/checkm_results -r -t 16
checkm tree_qa $ASSEMBLY_DIR/checkm
checkm qa $ASSEMBLY_DIR/checkm/lineage.ms $ASSEMBLY_DIR/checkm/ -o 1 > $ASSEMBLY_DIR/checkm/Final_table_01.csv
checkm qa $ASSEMBLY_DIR/checkm/lineage.ms $ASSEMBLY_DIR/checkm/ -o 2 > $ASSEMBLY_DIR/checkm/final_table_checkm.csv
micromamba deactivate

# 4.3 CheckM2
micromamba activate .micromamba/envs/04_checkm2
mkdir -p $ASSEMBLY_DIR/checkm2
# Note: CheckM2 predict uses a single thread parameter explicitly in the tutorial
checkm2 predict --threads 16 --input $ASSEMBLY_FASTA --output-directory $ASSEMBLY_DIR/checkm2 
micromamba deactivate

# ------------------------------------------------------------------------------
# 5. Genome Annotation
# ------------------------------------------------------------------------------
echo "--- Starting Prokka Annotation ---"
micromamba activate .micromamba/envs/05_prokka

# Prokka creates its own output directory.
prokka $ASSEMBLY_FASTA --outdir $GENOMICS_DIR/4_annotated_genome --kingdom Bacteria --addgenes --cpus 16

micromamba deactivate

# ------------------------------------------------------------------------------
# 6. Taxonomic Classification
# ------------------------------------------------------------------------------
echo "--- Starting GTDB-Tk Classification ---"
micromamba activate .micromamba/envs/06_gtdbtk

# Setting the database path required for GTDB-Tk
conda env config vars set GTDBTK_DATA_PATH="$WORK/databases/gtdbtk/release220"
mkdir -p $GENOMICS_DIR/5_gtdb_classification

# Using cpus=1 as suggested by the tutorial to manage memory, but the SLURM job has 64GB.
gtdbtk classify_wf --cpus 1 \
                   --genome_dir $GENOMICS_DIR/4_annotated_genome/ \
                   --out_dir $GENOMICS_DIR/5_gtdb_classification \
                   --extension .fna \
                   --skip_ani_screen


micromamba deactivate

# ------------------------------------------------------------------------------
# 7. MultiQC Summary
# ------------------------------------------------------------------------------
echo "--- Compiling Reports with MultiQC ---"
micromamba activate .micromamba/envs/01_short_reads_qc

# Combining all the QC reports into a single dashboard.
multiqc -d $GENOMICS_DIR/ -o $GENOMICS_DIR/6_multiqc

micromamba deactivate

echo "Day 6 Pipeline Completed Successfully!"