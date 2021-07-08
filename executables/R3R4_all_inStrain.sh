#! /bin/bash 

# load inStrain environment
export PATH=/home/GLBRCORG/emcdaniel/anaconda3/bin:$PATH
unset PYTHONPATH
source activate inStrain
PYTHONPATH=''# arguments

genome=$1
genomeName=$(basename $genome .fasta)

genes=/home/GLBRCORG/emcdaniel/EBPR/AcDiv/metagenomes/R3R4/R3R4_all_inStrain/ref_genomes/$genomeName.genes.fna

# cd to mapping folder

cd /home/GLBRCORG/emcdaniel/EBPR/AcDiv/mappingResults/R3R4

# profile command

inStrain profile R3R4-All-Acc.sorted.bam $genome -o /home/GLBRCORG/emcdaniel/EBPR/AcDiv/metagenomes/R3R4/R3R4_all_inStrain/$genomeName.IS -p 8 -g $genes -s /home/GLBRCORG/emcdaniel/EBPR/AcDiv/metagenomes/R3R4/R3R4_all_inStrain/R3R4_genomes.stb