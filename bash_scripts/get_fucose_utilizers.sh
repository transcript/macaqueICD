#!/bin/sh

# D.Lemay using SAMSA2 10/18/2018

# which microbes are using fucose?  Subsystems annotation = L_fucose_utilization
# NOTE L-fucose_utlization *include* fucosidase

# SEED results
SEED_DIR=/share/lemaylab-backedup/milklab/sam/macaque_MTs/5_aggregated_files/Subsystems/receipts

# Refseq results
REFSEQ_results_DIR=/share/lemaylab-backedup/milklab/sam/final_macaque_files/4_DIAMOND_annotations/bacteria_RefSeq/

# python scripts location
scripts_DIR=/share/lemaylab-backedup/milklab/macaque_revisit/programs/samsa2/python_scripts

# R scripts directory
R_scripts_DIR=/share/lemaylab-backedup/milklab/macaque_revisit/programs/samsa2/R_scripts

# output directory
OUT_DIR=/share/lemaylab-backedup/milklab/macaque_revisit/results/

#Refseq database
RefSeq_db="/share/lemaylab-backedup/milklab/sam/databases/bct"

# STEP 1 - Python, mapping the Subsystems results matching the
#"L-fucose_utlization" to hits to get the read IDs.

    echo "Begin Step 1"

for file in $SEED_DIR/*.m8.receipt
do
    # gets file name without directory name
    file_name=`echo $file | awk -F "/" '{print $NF}'`

    shortname=`echo $file_name | awk -F ".receipt" '{print $1}'`
    shortername=`echo $shortname | awk -F ".subsys_results.m8" '{print $1}'`

    RefSeq_results_name=$shortername.RefSeq_results.m8
    output_name=$shortername.fucose_utilization.txt

    python $scripts_DIR/Subsys_to_RefSeq_mapper.py -S $file -I $REFSEQ_results_DIR/$RefSeq_results_name -O $OUT_DIR/$output_name -T L-fucose_utilization
done

    # STEP 2 - Python, aggregate these results

    echo "Begin Step 2"

    for file in $OUT_DIR/*.fucose_utilization.txt
    do
	python $scripts_DIR/DIAMOND_analysis_counter.py -I $file -D $RefSeq_db -O
    done

    # move files
    mkdir $OUT_DIR/fucose_utilization
    mkdir $OUT_DIR/fucose_utilization/aggregated

    mv $OUT_DIR/*.fucose_utilization.txt $OUT_DIR/fucose_utilization/.
    mv $OUT_DIR/*.fucose_utilization_organism.tsv $OUT_DIR/fucose_utilization/aggregated/.

    #rename as control and experimental files so that the R scripts work
    agg_DIR=$OUT_DIR/fucose_utilization/aggregated/
    mv $agg_DIR/1_fullFile.fucose_utilization_organism.tsv $agg_DIR/experimental_1_fullFile.fucose_utilization_organism.tsv
    mv $agg_DIR/10_fullFile.fucose_utilization_organism.tsv $agg_DIR/experimental_10_fullFile.fucose_utilization_organism.tsv
    mv $agg_DIR/11_fullFile.fucose_utilization_organism.tsv $agg_DIR/experimental_11_fullFile.fucose_utilization_organism.tsv
    mv $agg_DIR/12_fullFile.fucose_utilization_organism.tsv $agg_DIR/experimental_12_fullFile.fucose_utilization_organism.tsv
    mv $agg_DIR/2_fullFile.fucose_utilization_organism.tsv $agg_DIR/experimental_2_fullFile.fucose_utilization_organism.tsv
    mv $agg_DIR/3_fullFile.fucose_utilization_organism.tsv $agg_DIR/experimental_3_fullFile.fucose_utilization_organism.tsv
    mv $agg_DIR/4_fullFile.fucose_utilization_organism.tsv $agg_DIR/experimental_4_fullFile.fucose_utilization_organism.tsv
    mv $agg_DIR/5_fullFile.fucose_utilization_organism.tsv $agg_DIR/experimental_5_fullFile.fucose_utilization_organism.tsv
    mv $agg_DIR/6_fullFile.fucose_utilization_organism.tsv $agg_DIR/experimental_6_fullFile.fucose_utilization_organism.tsv
    mv $agg_DIR/7_fullFile.fucose_utilization_organism.tsv $agg_DIR/experimental_7_fullFile.fucose_utilization_organism.tsv
    mv $agg_DIR/8_fullFile.fucose_utilization_organism.tsv $agg_DIR/experimental_8_fullFile.fucose_utilization_organism.tsv
    mv $agg_DIR/9_fullFile.fucose_utilization_organism.tsv $agg_DIR/experimental_9_fullFile.fucose_utilization_organism.tsv
    mv $agg_DIR/13_fullFile.fucose_utilization_organism.tsv $agg_DIR/control_13_fullFile.fucose_utilization_organism.tsv
    mv $agg_DIR/14_fullFile.fucose_utilization_organism.tsv $agg_DIR/control_14_fullFile.fucose_utilization_organism.tsv
    mv $agg_DIR/15_fullFile.fucose_utilization_organism.tsv $agg_DIR/control_15_fullFile.fucose_utilization_organism.tsv
    mv $agg_DIR/16_fullFile.fucose_utilization_organism.tsv $agg_DIR/control_16_fullFile.fucose_utilization_organism.tsv
    mv $agg_DIR/17_fullFile.fucose_utilization_organism.tsv $agg_DIR/control_17_fullFile.fucose_utilization_organism.tsv
    mv $agg_DIR/18_fullFile.fucose_utilization_organism.tsv $agg_DIR/control_18_fullFile.fucose_utilization_organism.tsv
    mv $agg_DIR/19_fullFile.fucose_utilization_organism.tsv $agg_DIR/control_19_fullFile.fucose_utilization_organism.tsv
    mv $agg_DIR/20_fullFile.fucose_utilization_organism.tsv $agg_DIR/control_20_fullFile.fucose_utilization_organism.tsv
    mv $agg_DIR/21_fullFile.fucose_utilization_organism.tsv $agg_DIR/control_21_fullFile.fucose_utilization_organism.tsv
    mv $agg_DIR/22_fullFile.fucose_utilization_organism.tsv $agg_DIR/control_22_fullFile.fucose_utilization_organism.tsv
    mv $agg_DIR/23_fullFile.fucose_utilization_organism.tsv $agg_DIR/control_23_fullFile.fucose_utilization_organism.tsv
    mv $agg_DIR/24_fullFile.fucose_utilization_organism.tsv $agg_DIR/control_24_fullFile.fucose_utilization_organism.tsv
    
    raw_counts_file=/share/lemaylab-backedup/milklab/sam/macaque_MTs/unix_cleaned_file_read_counts.txt

    echo "Begin Step 3"

    Rscript $R_scripts_DIR/run_DESeq_stats_local.R -I $OUT_DIR/fucose_utilization/aggregated/ -O RefSeq.fucose_utilization.org_DESeq_results.tab -R $raw_counts_file

    echo "Files have been analyzed, pulling out all organism results to fucose_utilization hits in the Subsystems results."
    echo "The raw results have been saved in $OUT_DIR/fucose_utilization ."
    echo "The aggregated results have been saved in $OUT_DIR/fucose_utilization/aggregated ."
    echo "The DESeq analysis has been saved in the above directory as RefSeq.fucose_utilization.org_DESeq_results.tab ."


    

