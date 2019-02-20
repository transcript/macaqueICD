#!/bin/bash

####################################################################
# campy_Subsys_functions.sh
# Danielle Lemay & Sam Westreich 10/23/18
####################################################################
# Purpose: Determine campylobacter-specific Subsystems functional annotations
#
# Step 1: Annotate all reads against RefSeq and Subsystems using DIAMOND.
# (previously done)
# Step 2: Use DIAMOND_specific_organism_retriever.py to retrieve all RefSeq annotations matching Campylobacter.
# Step 3: Use db_results_swapper.py on the Campylobacter-specific list of RefSeq results and the full Subsystems results.
#  This finds all the read IDs in the Subsystems results that are also present in the Campylobacter-specific RefSeq results.
# Step 4: Run DIAMOND_Subsystems_analysis_counter on the filtered Subsystems results (filtered to only Campylobacter-specific reads).
# Setp 5: Run Subsystems_DESeq_stats.R
#
####################################################################
#
# Subsystems Functional search by Campylobacter 
#
####################################################################
#
# VARIABLES

python_programs=/share/lemaylab-backedup/milklab/macaque_revisit/programs/samsa2/python_scripts

specific_organism="Campylobacter"

RefSeq_db="/share/lemaylab-backedup/milklab/sam/databases/bct"
Subsys_db="/share/lemaylab-backedup/milklab/sam/databases/full_subsys_db.fa"

RefSeq_results_location=/share/lemaylab-backedup/milklab/sam/final_macaque_files/4_DIAMOND_annotations/bacteria_RefSeq/

Subsys_results_location=/share/lemaylab-backedup/milklab/sam/final_macaque_files/5_aggregated_files/Subsystems/receipts/

raw_counts_file=/share/lemaylab-backedup/milklab/sam/macaque_MTs/unix_cleaned_file_read_counts.txt

export R_LIBS="/share/lemaylab-backedup/milklab/macaque_revisit/programs/samsa2/R_scripts/packages"
R_programs=/share/lemaylab-backedup/milklab/macaque_revisit/programs/samsa2/R_scripts


####################################################################
#
# STEP 1 - previously done
#
# STEP 2 - subset the RefSeq results for organism of interest

for file in $RefSeq_results_location/*RefSeq_results.m8
do
    python $python_programs/DIAMOND_specific_organism_retriever.py -I $file -SO "Campylobacter" -D $RefSeq_db
done

#########################################################################

# STEP 3 - Use db_results_swapper.py on the Campylobacter-specific list of RefSeq results and the full Subsystems results.
#  This finds all the read IDs in the Subsystems results that are also present in the Campylobacter-specific RefSeq results.
echo "Running Step 3\n"
for index in `seq 1 24`
do
    python $python_programs/db_results_swapper.py -I $RefSeq_results_location/Campylobacter/${index}_fullFile.RefSeq_result_Campylobacter.tsv -A $Subsys_results_location/${index}_#fullFile.subsys_results.m8.receipt -O ../data/campylobacter/${index}_results_swapper_output.txt
done

####################################################################
#
# STEP 4 - run Python to aggregate these results
echo "Running Step 4\n"
for index in `seq 1 24`
do
    python $python_programs/DIAMOND_subsystems_analysis_counter.py -I ../data/campylobacter/${index}_results_swapper_output.txt -D $Subsys_db -O ../data/campylobacter/${index}_subsys_counter_output.txt
done

#
####################################################################
#
# STEP 4.9 - move results to their own directories

mkdir ../data/campylobacter/aggregated

echo "rename as control and experimental files so that the R scripts work"
agg_DIR=../data/campylobacter/

for index in `seq 1 12`
do
    mv $agg_DIR/${index}_subsys_counter_output.txt $agg_DIR/aggregated/experimental_${index}_subsys_counter_output.txt
done
for index in `seq 13 24`
do
    mv $agg_DIR/${index}_subsys_counter_output.txt $agg_DIR/aggregated/control_${index}_subsys_counter_output.txt
done


####################################################################
#
# STEP 5 - run R to get DESeq results for the new files
#
# NOTE: Make sure that files have the appropriate prefixes ('control_'
#or 'experimental_') so that R can recognize them!
echo "installing R packages"
Rscript $R_programs/install_packages_local.R
echo "running DESeq2"
echo "$agg_DIR/aggregated/"
Rscript --no-save --no-restore --vanilla $R_programs/Subsystems_DESeq_stats_local.R -I $agg_DIR/aggregated/ -O Subsys.L1.Campylobacter.DESeq_results.tab -R $raw_counts_file -L 1



