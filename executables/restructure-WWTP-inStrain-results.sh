#! /bin/bash 

# arguments
genome=$1
sample=$2
code=$3
instrain=$(basename $genome .fa).IS
sample_name=$(basename $sample )-spRep
outfile=$sample_name-vs-$code

# move to filtered results directory 
cd /home/GLBRCORG/emcdaniel/EBPR/AcDiv/metagenomes/Danish-WWTPs/inStrain/filtered_results

# get gene and genome files for each of the filtered genomes
cp ../$sample_name/$instrain/output/$instrain*_gene_info.tsv $code/$outfile-gene-info.tsv

cp ../$sample_name/$instrain/output/$instrain*_genome_info.tsv $code/$outfile-genome-info.tsv

cp ../$sample_name/$instrain/output/$instrain*_SNVs.tsv $code/$outfile-SNVs.tsv