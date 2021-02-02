#! /bin/bash

source /opt/bifxapps/miniconda3/etc/profile.d/conda.sh
unset $PYTHONPATH

# activate environment
conda activate anvio-6.2

# cd where data is
cd /home/GLBRCORG/emcdaniel/EBPR/AcDiv/anvio/high_qual_refs

# arguments
fasta=$1

# run make contigs DB script
anvi-script-FASTA-to-contigs-db $fasta


