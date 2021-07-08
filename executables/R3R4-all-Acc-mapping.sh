#! /bin/bash

# Environment with samtools
export PATH=/home/GLBRCORG/emcdaniel/anaconda3/bin:$PATH
unset PYTHONPATH
source activate coverM
PYTHONPATH=''

cd /home/GLBRCORG/emcdaniel/EBPR/AcDiv/mappingResults/R3R4

# mapping command
/opt/bifxapps/bowtie2-2.3.5.1/bowtie2 --threads 4 -x /home/GLBRCORG/emcdaniel/EBPR/AcDiv/mappingResults/bowtie-refs/R3R4_Acc_bt2/R3R4_Acc_fasta --interleaved /home/GLBRCORG/emcdaniel/EBPR/AcDiv/metagenomes/R3R4/cleaned_fastqs/EBPRWIReactor34.cleaned.fastq > R3R4-All-Acc.sam

# BAM, sort, index
samtools view -S -b R3R4-All-Acc.sam > R3R4-All-Acc.bam
samtools sort R3R4-All-Acc.bam -o R3R4-All-Acc.sorted.bam
samtools index R3R4-All-Acc.sorted.bam R3R4-All-Acc.sorted.bam.bai