#! /bin/bash

# arguments
r1=$1
r2=$2
samplename=$3

cd /home/GLBRCORG/emcdaniel/EBPR/AcDiv/mappingResults/WWTPs

# mapping command
/opt/bifxapps/bowtie2-2.3.5.1/bowtie2 --threads 4 -x /home/GLBRCORG/emcdaniel/EBPR/AcDiv/ref_genomes/relative_abundance/bt2/Danish-sp-rep-genomes.fasta -1 $r1 -2 $r2 > $samplename-spRep.sam


# BAM, sort, index
/opt/bifxapps/samtools-1.9/bin/samtools view -S -b  $samplename-spRep.sam >  $samplename-spRep.bam
/opt/bifxapps/samtools-1.9/bin/samtools sort  $samplename-spRep.bam -o  $samplename-spRep.sorted.bam
/opt/bifxapps/samtools-1.9/bin/samtools index $samplename-spRep.sorted.bam $samplename-spRep.sorted.bam.bai
