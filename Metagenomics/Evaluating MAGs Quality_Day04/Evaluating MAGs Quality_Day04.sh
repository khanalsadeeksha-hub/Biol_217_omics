#!/bin/bash
#SBATCH --nodes=1
#SBATCH --cpus-per-task=12
#SBATCH --mem=40G
#SBATCH --time=06:00:00  # Should be faster than mapping
#SBATCH --job-name=day4_gunc
#SBATCH --output=day4_gunc.out
#SBATCH --error=day4_gunc.err
#SBATCH --partition=base
#SBATCH --reservation=biol217

# ------------------------------------------------------------------------------
# 1. Environment Setup
# ------------------------------------------------------------------------------
module load gcc12-env/12.1.0
module load micromamba 2> /dev/null
eval "$(micromamba shell hook --shell=bash)"
export MAMBA_ROOT_PREFIX=$WORK/.micromamba

# IMPORTANT: GUNC uses a different environment than Anvi'o
cd $WORK
micromamba activate $WORK/.micromamba/envs/00_gunc/

# Define directories based on the suggested tutorial structure
WORK_DIR="$WORK/metagenomics"
REFINE_DIR="$WORK_DIR/day4_refine"
GUNC_DB="$WORK/databases/gunc/gunc_db_progenomes2.1.dmnd"

echo "Starting Automated Day 4 GUNC Workflow..."

# ------------------------------------------------------------------------------
# 2. Chimera Detection and Plotting (Automated Loop)
# ------------------------------------------------------------------------------
# This loop dynamically finds any fasta file inside the subdirectories of REFINE_DIR
for MAG in $REFINE_DIR/*/*.fa; do
    
    # Get the directory containing the specific MAG fasta file
    BIN_DIR=$(dirname "$MAG")
    BASENAME=$(basename "$MAG" .fa)
    
    echo "Running GUNC on: $BASENAME"
    
    # Run GUNC to check for chimeras and potential contamination
    gunc run -i "$MAG" \
             -r "$GUNC_DB" \
             --out_dir "$BIN_DIR/gunc_out" \
             --detailed_output \
             --threads 12
             
    echo "Plotting GUNC results for: $BASENAME"
    
    # Create interactive plots of the chimeras
    gunc plot -d $BIN_DIR/gunc_out/diamond_output/*.diamond.progenomes_2.1.out \
              -g $BIN_DIR/gunc_out/gene_calls/gene_counts.json \
              --out_dir "$BIN_DIR/gunc_out"
              
done

echo "GUNC analysis complete! Check the CSS scores and 'PASS GUNC' columns in the output tables."