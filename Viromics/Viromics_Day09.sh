#!/bin/bash
# ==============================================================================
# Day 9: Viromics Pipeline and Data Extraction Script
# Note: The MVP and iPHoP pipeline commands (Section 1) were executed by 
# instructors due to runtime/database size, but are included here for reproducibility.
# ==============================================================================

# ------------------------------------------------------------------------------
# SECTION 1: MVP and iPHoP PIPELINE (Instructor Executed)
# ------------------------------------------------------------------------------
# module load gcc12-env/12.1.0
# module load micromamba/1.3.1
# micromamba activate MVP
# cd $WORK/MVP_test

# MVP Pipeline Steps
# mvip MVP_00_set_up_MVP -i ./WORKING_DIRECTORY/ -m input_file_timeseries_final.csv  --genomad_db_path ./WORKING_DIRECTORY/00_DATABASES/genomad_db/ --checkv_db_path ./WORKING_DIRECTORY/00_DATABASES/checkv-db-v1.5/
# mvip MVP_00_set_up_MVP -i ./WORKING_DIRECTORY/ -m input_file_timeseries_final.csv  --skip_install_databases
# mvip MVP_01_run_genomad_checkv -i WORKING_DIRECTORY/ -m input_file_timeseries_final.csv
# mvip MVP_02_filter_genomad_checkv -i WORKING_DIRECTORY/ -m input_file_timeseries_final.csv
# mvip MVP_03_do_clustering -i WORKING_DIRECTORY/ -m input_file_timeseries_final.csv
# mvip MVP_04_do_read_mapping -i WORKING_DIRECTORY/ -m input_file_timeseries_final.csv --delete_files
# mvip MVP_05_create_vOTU_table -i WORKING_DIRECTORY/ -m input_file_timeseries_final.csv
# mvip MVP_06_do_functional_annotation -i WORKING_DIRECTORY/ -m input_file_timeseries_final.csv
# mvip MVP_07_do_binning -i WORKING_DIRECTORY/ -m input_file_timeseries_final.csv --force_outputs
# mvip MVP_100_summarize_outputs -i WORKING_DIRECTORY/ -m input_file_timeseries_final.csv

# iPHoP Pipeline Steps
# conda activate GTDBTk
# export GTDBTK_DATA_PATH=./GTDB_db/GTDB_db
# gtdbtk de_novo_wf --genome_dir fa_all/ --bacteria --outgroup_taxon p__Patescibacteria --out_dir output/ --cpus 12 --force --extension fa
# gtdbtk de_novo_wf --genome_dir fa_all/ --archaea --outgroup_taxon p__Altarchaeota --out_dir output/ --cpus 12 --force --extension fa
# micromamba activate iphop_env
# iphop add_to_db --fna_dir fa_all/ --gtdb_dir ./output/ --out_dir ./MAGs_iPHoP_db --db_dir iPHoP_db/
# iphop predict --fa_file ./MVP_07_Filtered_conservative_Prokaryote_Host_Only_best_vBins_Representative_Unbinned_vOTUs_Sequences_iPHoP_Input.fasta --db_dir ./MAGs_iPHoP_db --out_dir ./iphop_output -t 12

# ------------------------------------------------------------------------------
# SECTION 2: DATA EXTRACTION & ANALYSIS QUERIES
# ------------------------------------------------------------------------------
echo "--- Virus Identification and Quality Control ---"

echo "1) Number of free viruses in BGR_140717:"
grep ">" -c 01_GENOMAD/BGR_140717/BGR_140717_Viruses_Genomad_Output/BGR_140717_modified_summary/BGR_140717_modified_virus.fna

echo "2) Number of proviruses in BGR_140717:"
grep ">" -c 01_GENOMAD/BGR_140717/BGR_140717_Proviruses_Genomad_Output/proviruses_summary/proviruses_virus.fna

echo "3) Number of Caudoviricetes viruses across all samples:"
grep -c "Caudoviricetes" 02_CHECK_V/BGR_*/MVP_02_BGR_*_Filtered_Relaxed_Merged_Genomad_CheckV_Virus_Proviruses_Quality_Summary.tsv

echo "4) Number of Unclassified viruses across all samples:"
grep -c "Unclassified" 02_CHECK_V/BGR_*/MVP_02_BGR_*_Filtered_Relaxed_Merged_Genomad_CheckV_Virus_Proviruses_Quality_Summary.tsv

echo "5) Other taxonomies excluding Caudoviricetes and Unclassified:"
grep -v "Caudoviricetes" 02_CHECK_V/BGR_*/MVP_02_BGR_*_Filtered_Relaxed_Merged_Genomad_CheckV_Virus_Proviruses_Quality_Summary.tsv | grep -v "Unclassified" | grep -v "Sample" | cut -f 3,5,6 | sort | uniq

echo "6) Number of High-quality and Complete viruses combined:"
cut -f 8 02_CHECK_V/BGR_*/MVP_02_BGR_*_Filtered_Relaxed_Merged_Genomad_CheckV_Virus_Proviruses_Quality_Summary.tsv | grep -e "High-quality" -e "Complete" | grep -c ""

echo "7) Quality category breakdown per sample:"
echo "| Sample | Low-quality | Medium-quality | High-quality | Complete |"
echo "|--------|-------------|----------------|--------------|----------|"
for x in BGR_130305 BGR_130527 BGR_130708 BGR_130829 BGR_130925 BGR_131021 BGR_131118 BGR_140106 BGR_140121 BGR_140221 BGR_140320 BGR_140423 BGR_140605 BGR_140717 BGR_140821 BGR_140919 BGR_141022 BGR_150108; do 
    LOW=$(cut -f 8 02_CHECK_V/"$x"/MVP_02_"$x"_Filtered_Relaxed_Merged_Genomad_CheckV_Virus_Proviruses_Quality_Summary.tsv | grep -c "Low-quality")
    MED=$(cut -f 8 02_CHECK_V/"$x"/MVP_02_"$x"_Filtered_Relaxed_Merged_Genomad_CheckV_Virus_Proviruses_Quality_Summary.tsv | grep -c "Medium-quality")
    HIGH=$(cut -f 8 02_CHECK_V/"$x"/MVP_02_"$x"_Filtered_Relaxed_Merged_Genomad_CheckV_Virus_Proviruses_Quality_Summary.tsv | grep -c "High-quality")
    COMP=$(cut -f 8 02_CHECK_V/"$x"/MVP_02_"$x"_Filtered_Relaxed_Merged_Genomad_CheckV_Virus_Proviruses_Quality_Summary.tsv | grep -c "Complete")
    echo "| $x | $LOW | $MED | $HIGH | $COMP |"
done

echo "--- Clustering and Abundance ---"

echo "11) Number of cluster representatives:"
grep -c "" 03_CLUSTERING/MVP_03_All_Sample_Filtered_Relaxed_Merged_Genomad_CheckV_Representative_Virus_Proviruses_Quality_Summary.tsv | awk '{print $1 - 1}'

echo "12) Number of cluster representatives that are proviruses:"
cut -f 5 03_CLUSTERING/MVP_03_All_Sample_Filtered_Relaxed_Merged_Genomad_CheckV_Representative_Virus_Proviruses_Quality_Summary.tsv | grep -c "Yes"

echo "13) Size of clusters for complete viruses:"
for x in BGR_131021_NODE_96_length_46113_cov_32.412567 BGR_140121_NODE_54_length_34619_cov_66.823718 BGR_140717_NODE_168_length_31258_cov_37.020094; do 
    echo -n "$x cluster size: "
    MEMBERS=$(cut -f 2 03_CLUSTERING/MVP_03_All_Sample_Filtered_Relaxed_Merged_Genomad_CheckV_Representative_Virus_Proviruses_Quality_Summary.tsv | grep "$x" | grep "," -o | grep "" -c)
    echo "$(($MEMBERS + 1))"
done

echo "16) Merging CoverM output files..."
grep --no-filename -v "Sample" 04_READ_MAPPING/BGR_*/*_CoverM.tsv > 04_READ_MAPPING/tmp_CoverM_output.tsv
grep "Sample" 04_READ_MAPPING/BGR_150108/BGR_150108_CoverM.tsv > 04_READ_MAPPING/tmp.tsv
cat 04_READ_MAPPING/tmp.tsv 04_READ_MAPPING/tmp_CoverM_output.tsv > 04_READ_MAPPING/merged_CoverM_output.tsv
rm 04_READ_MAPPING/tmp*

echo "17) Abundance (RPKM) for complete viruses:"
for x in "BGR_131021_NODE_96_length_46113_cov_32.412567" "BGR_140121_NODE_54_length_34619_cov_66.823718"  "BGR_140717_NODE_168_length_31258_cov_37.020094" ; do 
    grep "$x" 04_READ_MAPPING/merged_CoverM_output.tsv | cut -f1,2,11
done > 04_READ_MAPPING/merged_viruses_RPKM.tsv

echo "--- Annotation ---"

echo "19) Gene types (PHROGS Category) for complete viruses:"
for x in "BGR_131021_NODE_96_length_46113_cov_32.412567" "BGR_140121_NODE_54_length_34619_cov_66.823718"  "BGR_140717_NODE_168_length_31258_cov_37.020094" ; do 
    echo "Virus: $x"
    grep "$x" 06_FUNCTIONAL_ANNOTATION/MVP_06_All_Sample_Filtered_Conservative_Merged_Genomad_CheckV_Representative_Virus_Proviruses_Gene_Annotation_GENOMAD_PHROGS_PFAM_Filtered.tsv | cut -f23 | sort | uniq
done

echo "--- Host Prediction ---"

echo "27) Host predictions for complete viruses:"
for x in "BGR_131021_NODE_96_length_46113_cov_32.412567" "BGR_140121_NODE_54_length_34619_cov_66.823718"  "BGR_140717_NODE_168_length_31258_cov_37.020094" ; do 
    echo "Searching for: $x"
    grep "$x" 08_iPHoP/Host_prediction_to_genome_m90.csv
done

echo "Script complete."