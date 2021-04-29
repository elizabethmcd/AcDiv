#! /bin/bash

# load inStrain environment
export PATH=/home/GLBRCORG/emcdaniel/anaconda3/bin:$PATH
source activate inStrain
PYTHONPATH=''
unset PERL5LIB

# arguments
mapping=$1
genome=$2
code=$3

sample=$(basename $mapping -mapping.sorted.bam)
genomePath=/home/GLBRCORG/emcdaniel/EBPR/AcDiv/ref_genomes/POB_genomes

# cd to mapping directory
cd /home/GLBRCORG/emcdaniel/EBPR/AcDiv/metagenomes/POB/mappingResults

inStrain profile $mapping $genomePath/$genome.fa -o ../inStrain/$sample/$code.IS -p 8 -g $genomePath/$genome.genes.fna -s ../inStrain/POB-scaffolds-to-bins.stb