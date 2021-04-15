#! /bin/bash

###################
# Filtering metagenomic samples with fastp
# For use on WEI GLBRC servers running HT Condor
# Elizabeth McDaniel
##################

# set path where fastp is installed in local home directory bin
FASTPATH=/home/GLBRCORG/emcdaniel/bin

# queueing r1 r2 metagenomic reads and output folder/file names
file1=$1
file2=$2

sample=$(basename $file1 _1.fastq)

# move to cleaned directory

cd /home/GLBRCORG/emcdaniel/EBPR/AcDiv/metagenomes/Danish-WWTPs/cleaned_fastqs

$FASTPATH/fastp --in1 $file1 --in2 $file2 --out1 $sample-R1-cleaned.fastq --out2 $sample-R2-cleaned.fastq --detect_adapter_for_pe  --cut_tail --cut_tail_window_size 10 --cut_tail_mean_quality 20 -h $sample.html
