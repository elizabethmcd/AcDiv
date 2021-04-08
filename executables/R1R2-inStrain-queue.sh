#! /bin/bash 

#! /bin/bash 

# load inStrain environment
export PATH=/home/GLBRCORG/emcdaniel/anaconda3/bin:$PATH
unset PYTHONPATH
source activate inStrain
PYTHONPATH=''

# arguments 

mapping=$1
genome=$2
sampleDate=$3
scaffolds=$4

fasta=/home/GLBRCORG/emcdaniel/EBPR/AcDiv/mappingResults/bowtie-refs/ref_genomes/$genome.fasta
genes=/home/GLBRCORG/emcdaniel/EBPR/AcDiv/ref_genomes/annotations/$genome/$genome.genes.fna

# cd to location of inStrain results

cd $sampleDate

# clear preexisting runs
rm -rf *

# profile command

inStrain profile $mapping $fasta -o $genome.IS -p 8 -g $genes -s /home/GLBRCORG/emcdaniel/EBPR/AcDiv/metagenomes/R1R2/inStrain/$scaffolds.stb