#! /bin/bash

# arguments
sam=$1

outname=$(basename $sam .sam)

cd /home/GLBRCORG/emcdaniel/EBPR/AcDiv/mappingResults/WWTPs

# samtools BAM, sort, index
/opt/bifxapps/samtools-1.9/bin/samtools view -S -b $sam > $outname.bam
/opt/bifxapps/samtools-1.9/bin/samtools sort $outname.bam -o $outname.sorted.bam
/opt/bifxapps/samtools-1.9/bin/samtools index $outname.sorted.bam $outname.sorted.bam.bai
