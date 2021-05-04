#! /bin/bash 

# load environment where BLAST installed
export PATH=/home/GLBRCORG/emcdaniel/anaconda3/bin:$PATH
unset PYTHONPATH
source activate pipenv
PYTHONPATH=''

# group fasta file 
fasta=$1
group=$(basename $fasta .faa)

# cd to results directory for all subdirectories
cd /home/GLBRCORG/emcdaniel/EBPR/AcDiv/ref_genomes/pangenomics/archived/test_PID

# setup subdirectory
mkdir $group
cd $group
cp $fasta .

# BLAST commands and parsing 
makeblastdb -in $group.faa -out $group.db -dbtype prot
blastp -seg yes -soft_masking true -use_sw_tback -evalue 1e-3 -outfmt 11 -db $group.db -query $group.faa -out $group.blast 
blast_formatter -archive $group.blast -outfmt "6 std qcovs" -out $group.outfmt6
awk -F "\t" '{print $1"_vs_"$2"\t"$3}' $group.outfmt6 > $group.PID.txt