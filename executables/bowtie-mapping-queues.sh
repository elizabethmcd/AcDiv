#! /bin/bash

# arguments 
fastq=$1 
ref=$2
output=$3 

samplename=$(basename $fastq .qced.fastq)
refname=$(basename $ref .fasta)

cd $output

# mapping command
bowtie2 -p 4 -x $ref -q $fastq > $refname-vs-$samplename.sam 

