#!/bin/bash
#SBATCH --mem=90000                             # sets total memory at 128 Gb
#SBATCH --time=7-0:0:0
#SBATCH --mail-user=swestreich@gmail.com
#SBATCH --mail-type=END

# Macaque host reads pipeline
# created 12 Oct 2017

base_loc="/share/milklab/sam/macaque_MTs"

# Running DIAMOND

diamond_database="/share/milklab/sam/databases/macaque_proteins"
diamond_location="/home/swestreich/programs"

for filename in $base_loc/3_ribodepleted_files/*.gz
do
	echo "Performing annotation search on " $filename
    date
    $diamond_location/diamond blastx --db $diamond_database -q $filename -a $filename.host_matches.daa -t ./ -k 1

    echo "Converting file " $filename.host_matches.daa " to readable format"
    date
    $diamond_location/diamond view --daa $filename.host_matches.daa -o $filename.host_results.m8 -f tab
done

mv $base_loc/3_ribodepleted_files/*host_results.m8 $base_loc/4_DIAMOND_annotations/macaque_host_RefSeq/.
mkdir $base_loc/4_DIAMOND_annotations/macaque_host_RefSeq/daa_files/
mv $base_loc/3_ribodepleted_files/*.daa $base_loc/4_DIAMOND_annotations/macaque_host_RefSeq/daa_files/.

# Aggregation

for file in $base_loc/4_DIAMOND_annotations/macaque_host_RefSeq/*.m8
do
	python /share/milklab/sam/SAMSA_pipeline_v2/aggregation/DIAMOND_nonbacteria_analysis_counter.py -I $file -F -R -D /share/milklab/sam/databases/macaque_proteins.fa
done

mv $base_loc/4_DIAMOND_annotations/macaque_host_RefSeq/*function.tsv $base_loc/5_aggregated_files/macaque_host_RefSeq/.

# Analysis

Rscript /share/milklab/sam/SAMSA_pipeline_v2/R_analysis/run_DESeq_stats.R -I $base_loc/5_aggregated_files/macaque_host_RefSeq/ -O $base_loc/5_aggregated_files/macaque_host_RefSeq/R_analysis/macaque_host_functions_DESeq_11_Oct.tab