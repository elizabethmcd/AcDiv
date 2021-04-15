#! /bin/bash

file1=$1
file2=$2
name=$3
output=$4
reformat=$5

# fastp QC job

~/bin/fastp -i $file1 -I $file2 -o "$output".qced.R1.fastq -O "$output".qced.R2.fastq --cut_tail -h $name.html

# reformat PE to interleaved

/opt/bifxapps/bbmap-38.32/reformat.sh in="$output".qced.R1.fastq in2="$output".qced.R2.fastq out="$reformat".fastq
