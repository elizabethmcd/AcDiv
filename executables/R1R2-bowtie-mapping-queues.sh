#! /bin/bash

# arguments
fastq=$1
ref=$2

samplename=$(basename $fastq .qced.fastq)
refname=$(basename $ref .fasta)

samplename=basename($file $fastq .qced.fastq)
refname=basename($file $ref .fasta)
outname=$refname-vs-$samplename

cd /home/GLBRCORG/emcdaniel/EBPR/AcDiv/mappingResults/R1R2

# mapping command
/opt/bifxapps/bowtie2-2.3.5.1/bowtie2 --threads 4 -x $ref --interleaved $fastq > $outname.sam


# BAM, sort, index
/opt/bifxapps/samtools-1.9/bin/samtools view -S -b $outname.sam > $outname.bam
/opt/bifxapps/samtools-1.9/bin/samtools sort $outname.bam -o $outname.sorted.bam
/opt/bifxapps/samtools-1.9/bin/samtools index $outname.sorted.bam $outname.sorted.bam.bai