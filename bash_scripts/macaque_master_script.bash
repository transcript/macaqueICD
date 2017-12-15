#!/bin/bash
#SBATCH --mem=128000
#SBATCH --time=7-0:0:0
#SBATCH --mail-user=swestreich@gmail.com
#SBATCH --mail-type=END

####################################################################
#
# macaque_master_script.bash
# Created March 2017 by Sam Westreich, github.com/transcript
#
####################################################################
#
# This script sets up and runs through ALL steps in the SAMSA pipeline
# before the analysis (which is done in R, likely in RStudio).  Each
# step is set up below.
#
# The steps are:
#	1. Merging with PEAR, if applicable
#	2. Read cleaning with Trimmomatic
#	3. rRNA removal with SortMeRNA
#	4. Annotation using DIAMOND (by default against the RefSeq database)
#	5. Aggregation using analysis_counter.py
#
####################################################################
#
# VARIABLES - set these paths for each step.
#
#	0. Starting files location
starting_location=/share/milklab/sam/macaque_MTs

#	1. PEAR
pear_location=~/programs/pear-0.9.6

# 	2. Trimmomatic
trimmomatic_location=~/programs/Trimmomatic-0.33

#	3. SortMeRNA
sortmerna_location=~/programs/sortmerna-2.1-linux-64

#	4. DIAMOND
diamond_database="/share/milklab/sam/databases/bct"
diamond_subsys_db="/share/milklab/sam/databases/full_subsys_db"
diamond_location="/software/diamond/0.7.9/x86_64-linux-ubuntu14.04/bin"

#	5. Aggregation
analysis_counter_location=/share/milklab/sam/python_scripts

####################################################################
#
# STEP 1: MERGING OF PAIRED-END FILES USING PEAR
# Note: paired-end files are usually named using R1 and R2 in the name.
# Note: if using single-end sequencing, skip this step (comment out).

#for file in $starting_location/*R1*
#do
#	file1=$file
#	file2=`echo $file1 | awk -F"R1" '{print $1 "R2" $2}'`
#	out_path=`echo $file | awk -F"R1" '{print $1 "merged"}'`
#	out_name=`echo ${out_path##*/}`

#	$pear_location/pear-0.9.6 -f $file1 -r $file2 -o $out_name
#done

#mkdir $starting_location/step_1_output/
#mv $starting_location/*merged* $starting_location/step_1_output/

####################################################################
#
# STEP 2: CLEANING FILES WITH TRIMMOMATIC
# Note: if skipping FLASH, make sure that all starting files are in the
# $starting_location/step_1_output/ folder!

#for file in $starting_location/step_1_output/*.assembled*
#do
#	shortname=`echo $file | awk -F"assembled" '{print $1 "cleaned.fastq"}'`
#	
#	java -jar $trimmomatic_location/trimmomatic-0.33.jar SE -phred33 $file $shortname SLIDINGWINDOW:4:15 MINLEN:99
#done

#mkdir $starting_location/step_2_output/
#mv $starting_location/step_1_output/*cleaned.fastq $starting_location/step_2_output/

####################################################################
#
# STEP 3: REMOVING RIBOSOMAL READS WITH SORTMERNA
# Note: this step assumes that the SortMeRNA databases are indexed.  If not,
# do that first (see the SortMeRNA user manual for details).

for file in $starting_location/merged_files/*.cleaned.fastq.gz
do
	unzip_name=`echo $file | awk -F".gz" '{print $1}'`
	shortname=`echo $file | awk -F"cleaned" '{print $1 "ribodepleted"}'`
	
	gunzip $file
	
	$sortmerna_location/sortmerna --ref $sortmerna_location/rRNA_databases/silva-bac-16s-id90.fasta,$sortmerna_location/index/silva-bac-16s-db --reads $unzip_name --aligned $unzip_name.ribosomes --other $shortname --fastx --num_alignments 0 --log -v

done

####################################################################
#
# STEP 4: ANNOTATING WITH DIAMOND AGAINST REFSEQ
# Note: this step assumes that the DIAMOND database is already built.  If not,
# do that first before running this step.

echo "Now starting on DIAMOND org annotations at: "; date

for file in $starting_location/merged_files/*.ribodepleted.fastq
do
	shortname=`echo $file | awk -F"ribodepleted" '{print $1 "org_annotated"}'`
	echo "Now starting on " $file
	
	$diamond_location/diamond blastx --db $diamond_database -q $file -a $file.dmd -t ./ -k 1 --sensitive
	$diamond_location/diamond view --daa $file.dmd -o $shortname -f tab
done

mv $starting_location/merged_files/*annotated* $starting_location/DIAMOND_annotations/
mv $starting_location/merged_files/*.daa $starting_location/merged_files/DIAMOND_annotations/

echo "DIAMOND annotations completed at: "; date

####################################################################
#
# STEP 5: AGGREGATING WITH ANALYSIS_COUNTER

for file in $starting_location/DIAMOND_annotations/*annotated*
do
	python $analysis_counter_location/DIAMOND_analysis_counter.py -I $file -D $diamond_database -O
	python $analysis_counter_location/DIAMOND_analysis_counter.py -I $file -D $diamond_database -F
done

mv $starting_location/DIAMOND_annotations/*.tsv $starting_location/DIAMOND_annotations/RefSeq_results/

####################################################################
#
# STEP 4.1: ANNOTATING WITH DIAMOND AGAINST SUBSYSTEMS

echo "Now starting on DIAMOND Subsystems annotations at: "; date

for file in $starting_location/merged_files/*.ribodepleted.fastq
do
	shortname=`echo $file | awk -F"ribodepleted" '{print $1 "subsys_annotated"}'`
	echo "Now starting on Subsystems annotations for file: " $file

	$diamond_location/diamond blastx --db $diamond_subsys_db -q $file -a $file.dmd -t ./ -k 1 --sensitive
	$diamond_location/diamond view --daa $file.dmd -o $shortname -f tab
done

mv $starting_location/step_3_output/*annotated* $starting_location/DIAMOND_annotations/
mv $starting_location/step_3_output/*.daa $starting_location/merged_files/DIAMOND_annotations/

echo "DIAMOND Subsystems annotations completed at: "; date

for file in $starting_location/DIAMOND_annotations/*subsys_annotated*
do
	python $analysis_counter_location/DIAMOND_subsystems_analysis_counter.py -I $file -D $diamond_subsys_db.fa -P $file.receipt
done

mkdir $starting_location/DIAMOND_annotations/Subsystems_results/
mv $starting_location/DIAMOND_annotations/*.tsv $starting_location/DIAMOND_annotations/Subsystems_results/

echo "Master bash script finished running at: "; date
exit
####################################################################
