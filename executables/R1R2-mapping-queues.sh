#! /bin/bash

# arguments 
fastq=$1 
ref=$2
output=$3 

samplename=basename($file $fastq .qced.fastq)
refname=basename($file $ref .fasta)
outname=$refname-vs-$samplename

# mapping command
bowtie2 -p 4 -x $ref -q $fastq > $output/$outname.sam 

# BAM, sort, index
/opt/bifxapps/samtools-1.9/bin/samtools view -S -b $outname.sam > $outname.bam
/opt/bifxapps/samtools-1.9/bin/samtools sort $outname.bam -o $outname.sorted.bam
/opt/bifxapps/samtools-1.9/bin/samtools index $outname.sorted.bam $outname.sorted.bam.bai

