#! /bin/bash

export PATH=/home/GLBRCORG/emcdaniel/anaconda3/bin:$PATH
unset PYTHONPATH
source activate coverM
PYTHONPATH=''

# arguments
sam=$1

outname=$(basename $sam .sam)

cd /home/GLBRCORG/emcdaniel/EBPR/AcDiv/metagenomes/R1R2/R1R2_UW1_UW3_mappingResults/

# samtools BAM, sort, index
/opt/bifxapps/samtools-1.9/bin/samtools view -S -b $sam > $outname.bam
/opt/bifxapps/samtools-1.9/bin/samtools sort $outname.bam -o $outname.sorted.bam
/opt/bifxapps/samtools-1.9/bin/samtools index $outname.sorted.bam $outname.sorted.bam.bai
