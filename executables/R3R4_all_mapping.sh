#! /bin/bash

# arguments
fastq=$1
ref=$2

samplename=$(basename $fastq .cleaned.fastq)
refname=$(basename $ref .fasta)

cd /home/GLBRCORG/emcdaniel/EBPR/AcDiv/metagenomes/R3R4/all_combined

# mapping command
/opt/bifxapps/bowtie2-2.3.5.1/bowtie2 --threads 4 -x $ref --interleaved $fastq > $refname-vs-$samplename.sam

