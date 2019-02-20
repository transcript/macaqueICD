#!/bin/bash

####################################################################
#
# Adapted from Functions_by_organism.sh
# Created June 17, 2017 by Sam Westreich, github.com/transcript
# This version by Danielle Lemay 10/23/18
#
####################################################################
#
# Purpose: identify functions associated with transcripts mapped to Haemophilus influenzae
#
####################################################################
#
# VARIABLES

python_programs=/share/lemaylab-backedup/milklab/macaque_revisit/programs/samsa2/python_scripts

RefSeq_db="/share/lemaylab-backedup/milklab/sam/databases/bct"

RefSeq_results_location=/share/lemaylab-backedup/milklab/sam/final_macaque_files/4_DIAMOND_annotations/bacteria_RefSeq/

raw_counts_file=/share/lemaylab-backedup/milklab/sam/macaque_MTs/unix_cleaned_file_read_counts.txt

R_programs=/share/lemaylab-backedup/milklab/macaque_revisit/programs/samsa2/R_scripts

#
####################################################################
#
echo "STEP 1 - run Python to subset the RefSeq results for organism of interest"

for file in $RefSeq_results_location/*RefSeq_results.m8
do
    python $python_programs/DIAMOND_specific_organism_retriever.py -I $file -SO "Haemophilus influenzae" -D $RefSeq_db
done

#
####################################################################
#
echo "STEP 2 - run Python to aggregate these results"
for file in $RefSeq_results_location/*Haemophilus\ influenzae.tsv
do
    python $python_programs/DIAMOND_analysis_counter.py -I $file -D $RefSeq_db -F
done

#
####################################################################
#
echo "STEP 2.9 - move results to their own directories"

mkdir $RefSeq_results_location/Haemophilus_influenzae
mkdir $RefSeq_results_location/Haemophilus_influenzae/aggregated

mv $RefSeq_results_location/*Haemophilus\ influenzae.tsv $RefSeq_results_location/Haemophilus_influenzae/
mv $RefSeq_results_location/*Haemophilus\ influenzae*function.tsv $RefSeq_results_location/Haemophilus_influenzae/aggregated

echo "rename as control and experimental files so that the R scripts work"
agg_DIR=$RefSeq_results_location/Haemophilus_influenzae/aggregated
mv $agg_DIR/1_fullFile.RefSeq_result_Haemophilus_influenzae_function.tsv $agg_DIR/experimental_1_fullFile.RefSeq_result_Haemophilus_influenzae_function.tsv
mv $agg_DIR/10_fullFile.RefSeq_result_Haemophilus_influenzae_function.tsv $agg_DIR/experimental_10_fullFile.RefSeq_result_Haemophilus_influenzae_function.tsv
mv $agg_DIR/11_fullFile.RefSeq_result_Haemophilus_influenzae_function.tsv $agg_DIR/experimental_11_fullFile.RefSeq_result_Haemophilus_influenzae_function.tsv
mv $agg_DIR/12_fullFile.RefSeq_result_Haemophilus_influenzae_function.tsv $agg_DIR/experimental_12_fullFile.RefSeq_result_Haemophilus_influenzae_function.tsv
mv $agg_DIR/2_fullFile.RefSeq_result_Haemophilus_influenzae_function.tsv $agg_DIR/experimental_2_fullFile.RefSeq_result_Haemophilus_influenzae_function.tsv
mv $agg_DIR/3_fullFile.RefSeq_result_Haemophilus_influenzae_function.tsv $agg_DIR/experimental_3_fullFile.RefSeq_result_Haemophilus_influenzae_function.tsv
mv $agg_DIR/4_fullFile.RefSeq_result_Haemophilus_influenzae_function.tsv $agg_DIR/experimental_4_fullFile.RefSeq_result_Haemophilus_influenzae_function.tsv
mv $agg_DIR/5_fullFile.RefSeq_result_Haemophilus_influenzae_function.tsv $agg_DIR/experimental_5_fullFile.RefSeq_result_Haemophilus_influenzae_function.tsv
mv $agg_DIR/6_fullFile.RefSeq_result_Haemophilus_influenzae_function.tsv $agg_DIR/experimental_6_fullFile.RefSeq_result_Haemophilus_influenzae_function.tsv
mv $agg_DIR/7_fullFile.RefSeq_result_Haemophilus_influenzae_function.tsv $agg_DIR/experimental_7_fullFile.RefSeq_result_Haemophilus_influenzae_function.tsv
mv $agg_DIR/8_fullFile.RefSeq_result_Haemophilus_influenzae_function.tsv $agg_DIR/experimental_8_fullFile.RefSeq_result_Haemophilus_influenzae_function.tsv
mv $agg_DIR/9_fullFile.RefSeq_result_Haemophilus_influenzae_function.tsv $agg_DIR/experimental_9_fullFile.RefSeq_result_Haemophilus_influenzae_function.tsv
mv $agg_DIR/13_fullFile.RefSeq_result_Haemophilus_influenzae_function.tsv $agg_DIR/control_13_fullFile.RefSeq_result_Haemophilus_influenzae_function.tsv
mv $agg_DIR/14_fullFile.RefSeq_result_Haemophilus_influenzae_function.tsv $agg_DIR/control_14_fullFile.RefSeq_result_Haemophilus_influenzae_function.tsv
mv $agg_DIR/15_fullFile.RefSeq_result_Haemophilus_influenzae_function.tsv $agg_DIR/control_15_fullFile.RefSeq_result_Haemophilus_influenzae_function.tsv
mv $agg_DIR/16_fullFile.RefSeq_result_Haemophilus_influenzae_function.tsv $agg_DIR/control_16_fullFile.RefSeq_result_Haemophilus_influenzae_function.tsv
mv $agg_DIR/17_fullFile.RefSeq_result_Haemophilus_influenzae_function.tsv $agg_DIR/control_17_fullFile.RefSeq_result_Haemophilus_influenzae_function.tsv
mv $agg_DIR/18_fullFile.RefSeq_result_Haemophilus_influenzae_function.tsv $agg_DIR/control_18_fullFile.RefSeq_result_Haemophilus_influenzae_function.tsv
mv $agg_DIR/19_fullFile.RefSeq_result_Haemophilus_influenzae_function.tsv $agg_DIR/control_19_fullFile.RefSeq_result_Haemophilus_influenzae_function.tsv
mv $agg_DIR/20_fullFile.RefSeq_result_Haemophilus_influenzae_function.tsv $agg_DIR/control_20_fullFile.RefSeq_result_Haemophilus_influenzae_function.tsv
mv $agg_DIR/21_fullFile.RefSeq_result_Haemophilus_influenzae_function.tsv $agg_DIR/control_21_fullFile.RefSeq_result_Haemophilus_influenzae_function.tsv
mv $agg_DIR/22_fullFile.RefSeq_result_Haemophilus_influenzae_function.tsv $agg_DIR/control_22_fullFile.RefSeq_result_Haemophilus_influenzae_function.tsv
mv $agg_DIR/23_fullFile.RefSeq_result_Haemophilus_influenzae_function.tsv $agg_DIR/control_23_fullFile.RefSeq_result_Haemophilus_influenzae_function.tsv
mv $agg_DIR/24_fullFile.RefSeq_result_Haemophilus_influenzae_function.tsv $agg_DIR/control_24_fullFile.RefSeq_result_Haemophilus_influenzae_function.tsv

#
####################################################################
#
# STEP 3 - run R to get DESeq results for the new files
#
# NOTE: Make sure that files have the appropriate prefixes ('control_' or 'experimental_') so that R can recognize them!
echo "installing R packages"
Rscript $R_programs/install_packages_local.R
echo "running DESeq2"
Rscript $R_programs/run_DESeq_stats_local.R -I $RefSeq_results_location/Haemophilus_influenzae/aggregated/ -O RefSeq.Haemophilus_influenzae.func_DESeq_results.tab -R $raw_counts_file


