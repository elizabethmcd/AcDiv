#! /bin/bash 

# load inStrain environment
export PATH=/home/GLBRCORG/emcdaniel/anaconda3/bin:$PATH
unset PYTHONPATH
source activate inStrain
PYTHONPATH=''

# arguments 

mapping=$1
genome=$2
reactor=$3
scaffolds=$4

fasta=/home/GLBRCORG/emcdaniel/EBPR/AcDiv/ref_genomes/annotations/$genome/$genome.fna
genes=/home/GLBRCORG/emcdaniel/EBPR/AcDiv/ref_genomes/annotations/$genome/$genome.genes.fna

# cd to location of inStrain results

cd /home/GLBRCORG/emcdaniel/EBPR/AcDiv/metagenomes/R3R4/inStrain/$reactor

# clear preexisting runs
rm -rf *

# profile command

inStrain profile $mapping $fasta -o $genome.IS -p 8 -g $genes -s /home/GLBRCORG/emcdaniel/EBPR/AcDiv/metagenomes/R3R4/inStrain/$scaffolds.stb