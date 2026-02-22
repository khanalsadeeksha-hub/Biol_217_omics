#!/bin/bash
#SBATCH --nodes=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH --time=03:00:00 # Upped to 3 hours just to be safe for 52 genomes!
#SBATCH --job-name=anvio_pangenomics
#SBATCH --output=anvio_pangenomics.out
#SBATCH --error=anvio_pangenomics.err
#SBATCH --partition=base
#SBATCH --reservation=biol217

# ------------------------------------------------------------------------------
# 1. Environment Setup
# ------------------------------------------------------------------------------
# Loading modules and activating the Anvi'o environment.
module load gcc12-env/12.1.0
module load micromamba/1.4.2
eval "$(micromamba shell hook --shell=bash)"
cd $WORK
micromamba activate .micromamba/envs/00_anvio/

# Create a dedicated working folder.
PANGENOME_DIR="$WORK/pangenomics/01_anvio_pangenomics"
mkdir -p $PANGENOME_DIR
cd $PANGENOME_DIR

echo "--- Day 7: Vibrio jascida Pangenomics Pipeline Started ---"

# ------------------------------------------------------------------------------
# 2. Download Data
# ------------------------------------------------------------------------------
echo "Downloading the 52 Vibrio jascida genomes..."
# Fetch the dataset and unzip it.
curl -L https://ndownloader.figshare.com/files/28965090 -o V_jascida_genomes.tar.gz
tar -zxvf V_jascida_genomes.tar.gz

# Move into the extracted directory to do our work
cd V_jascida_genomes
PROJECT_NAME="V_jascida"

# ------------------------------------------------------------------------------
# 3. Prep fasta files and build contigs databases
# ------------------------------------------------------------------------------
echo "Formatting fasta files..."
# Rename any .fna files to .fasta just in case.
for file in *.fna; do mv "$file" "${file%.fna}.fasta" 2>/dev/null || true; done

# Make a list of all genome names by stripping everything after the first period.
ls *fasta | awk 'BEGIN{FS="."}{print $1}' > genomes.txt

# Reformat the fasta files to remove tiny contigs (<2500nt) and simplify names.
for g in `cat genomes.txt`; do
    echo "Reformatting $g..."
    anvi-script-reformat-fasta ${g}.fasta \
                               --min-len 2500 \
                               --simplify-names \
                               -o ${g}_2.5K.fasta
done

echo "Building contigs databases..."
# Convert the cleaned fasta files into Anvi'o contigs databases.
for g in `cat genomes.txt`; do
    echo "Creating DB for $g..." 
    anvi-gen-contigs-database -f ${g}_2.5K.fasta \
                              -o ${PROJECT_NAME}_${g}.db \
                              --num-threads 16 \
                              -n ${PROJECT_NAME}_${g}
done

# ------------------------------------------------------------------------------
# 4. Annotate the Contigs Databases
# ------------------------------------------------------------------------------
echo "Annotating all databases (This will take a while!)..."
# Run HMMs, COGs, tRNAs, and single-copy core gene taxonomy on every database.
for g in *.db; do
    echo "Annotating $g..."
    anvi-run-hmms -c $g --num-threads 16
    anvi-run-ncbi-cogs -c $g --num-threads 16
    anvi-scan-trnas -c $g --num-threads 16
    anvi-run-scg-taxonomy -c $g --num-threads 16
done

# ------------------------------------------------------------------------------
# 5. Build External Genomes File & Estimate Contamination
# ------------------------------------------------------------------------------
echo "Generating external genomes file..."
# Link all the individual databases into one master file.
anvi-script-gen-genomes-file --input-dir . -o external-genomes.txt

echo "Estimating genome completeness and contamination..."
# Assess the quality of the downloaded genomes.
anvi-estimate-genome-completeness -e external-genomes.txt -o genome_completeness.txt

# NOTE: The tutorial suggests manual refinement here if bins are bad. 
# Since this is a batch script, we are assuming the reference genomes are already high-quality and proceeding.

# ------------------------------------------------------------------------------
# 6. Compute Pangenome and ANI
# ------------------------------------------------------------------------------
echo "Computing the Pangenome..."
# Generate a combined storage database.
anvi-gen-genomes-storage -e external-genomes.txt -o ${PROJECT_NAME}-GENOMES.db

# Compute the actual pangenome using 16 threads.
anvi-pan-genome -g ${PROJECT_NAME}-GENOMES.db \
                --project-name ${PROJECT_NAME} \
                --num-threads 16

echo "Calculating Average Nucleotide Identity (ANI)..."
# Calculate genome similarities using pyANI.
anvi-compute-genome-similarity --external-genomes external-genomes.txt \
                               --program pyANI \
                               --output-dir ANI \
                               --num-threads 16 \
                               --pan-db ${PROJECT_NAME}/${PROJECT_NAME}-PAN.db 



# ------------------------------------------------------------------------------
# 7. Phylogenomics (Optional Task Automated)
# ------------------------------------------------------------------------------
echo "Computing Phylogenomics Tree..."
# Extract single-copy core genes (SCGs) that occur in at least 5 genomes, max 1 per genome, and concatenate them.
anvi-get-sequences-for-gene-clusters -p ${PROJECT_NAME}/${PROJECT_NAME}-PAN.db \
                                     -g ${PROJECT_NAME}-GENOMES.db \
                                     --min-num-genomes-gene-cluster-occurs 5 \
                                     --max-num-genes-from-each-genome 1 \
                                     --concatenate-gene-clusters \
                                     --output-file ${PROJECT_NAME}/${PROJECT_NAME}-SCGs.fa

# Trim poor alignments.
trimal -in ${PROJECT_NAME}/${PROJECT_NAME}-SCGs.fa \
       -out ${PROJECT_NAME}/${PROJECT_NAME}-SCGs-clean.fa \
       -gt 0.5

# Build the phylogenetic tree using IQTREE.
iqtree -s ${PROJECT_NAME}/${PROJECT_NAME}-SCGs-clean.fa \
       -m WAG \
       -bb 1000 \
       -nt 16

# Add the tree as an organizational layer to the pangenome database.
echo -e "item_name\tdata_type\tdata_value" > ${PROJECT_NAME}/${PROJECT_NAME}-phylogenomic-layer-order.txt
echo -e "SCGs_Bayesian_Tree\tnewick\t`cat ${PROJECT_NAME}/${PROJECT_NAME}-SCGs-clean.fa.treefile`" >> ${PROJECT_NAME}/${PROJECT_NAME}-phylogenomic-layer-order.txt

anvi-import-misc-data -p ${PROJECT_NAME}/${PROJECT_NAME}-PAN.db \
                      -t layer_orders ${PROJECT_NAME}/${PROJECT_NAME}-phylogenomic-layer-order.txt



echo "Day 7 Pipeline Complete! You are now ready to visualize the pangenome."