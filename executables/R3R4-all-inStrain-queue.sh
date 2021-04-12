#! /bin/bash 

# load inStrain environment
export PATH=/home/GLBRCORG/emcdaniel/anaconda3/bin:$PATH
unset PYTHONPATH
source activate inStrain
PYTHONPATH=''

# arguments 

mapping=$1
genome=$2
scaffolds=$3

fasta=/home/GLBRCORG/emcdaniel/EBPR/AcDiv/mappingResults/bowtie-refs/ref_genomes/$genome.fasta
genes=/home/GLBRCORG/emcdaniel/EBPR/AcDiv/ref_genomes/annotations/$genome/$genome.genes.fna

# cd to location of inStrain results

cd /home/GLBRCORG/emcdaniel/EBPR/AcDiv/metagenomes/R3R4/all_combined

# profile command

inStrain profile $mapping $fasta -o $genome.IS -p 8 -g $genes -s $scaffolds.stb