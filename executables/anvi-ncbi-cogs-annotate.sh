#! /bin/bash

source /opt/bifxapps/miniconda3/etc/profile.d/conda.sh
unset $PYTHONPATH

# activate environment
conda activate anvio-6.2

# cd where data is
cd /home/GLBRCORG/emcdaniel/EBPR/AcDiv/anvio/contig_dbs

# arguments
contigs=$1

# run ncbi-cogs step
anvi-run-ncbi-cogs -c $contigs --cog-data-dir /home/GLBRCORG/emcdaniel/bin/ncbi_cogs/ -T 10 --search-with diamond --sensitive 