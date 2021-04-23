#! /bin/bash 

# load inStrain environment
export PATH=/home/GLBRCORG/emcdaniel/anaconda3/bin:$PATH
unset PYTHONPATH
source activate pipenv
PYTHONPATH=''

# arguments

mapping=$1
genome=$2

plant=$(basename $mapping .sorted.bam)
genomeName=$(basename $genome .fa)

fasta=/home/GLBRCORG/emcdaniel/EBPR/AcDiv/ref_genomes/Singleton_SP_REP_genomes/$genomeName.fa
genes=/home/GLBRCORG/emcdaniel/EBPR/AcDiv/ref_genomes/Singleton_SP_REP_genomes/$genomeName.genes.fna

# cd to mapping folder

cd /home/GLBRCORG/emcdaniel/EBPR/AcDiv/metagenomes/Danish-WWTPs/mappingResults

# profile command

inStrain profile $mapping $fasta -o ../inStrain/$plant/$genomeName.IS -p 8 -g genes -s ../inStrain/Singleton-SPREP.stb
