#! /bin/bash

# arguments 
fastq=$1 
ref=$2
output=$3 

samplename=basename($file $fastq .qced.fastq)
refname=basename($file $ref .fasta)

# mapping command
bowtie2 -p 4 -x $ref -q $fastq > $output/$refname-vs-$samplename.sam 

