#! /bin/bash 

# load inStrain environment
export PATH=/home/GLBRCORG/emcdaniel/anaconda3/bin:$PATH
unset PYTHONPATH
source activate inStrain
PYTHONPATH=''


# arguments
sam=$1
outfolder=$(basename $sam .sorted.bam)-quick-profile

# cd to mapping results folder

cd /home/GLBRCORG/emcdaniel/EBPR/AcDiv/mappingResults/WWTPs

# inStrain quick profile command


inStrain quick_profile -p 2 -s ../../metagenomes/Danish-WWTPs/inStrain/Singleton-SPREP.stb -o ../../metagenomes/Danish-WWTPs/quick_profiles/$outfolder $sam ../../ref_genomes/relative_abundance/all-Danish-sp-rep-genomes.fasta