#! /bin/bash

source /opt/bifxapps/miniconda3/etc/profile.d/conda.sh
unset $PYTHONPATH

# activate environment
conda activate anvio-6.2

# cd where data is
cd /home/GLBRCORG/emcdaniel/EBPR/AcDiv/anvio/contig_dbs

# run the construction of pangenome
anvi-pan-genome -g ACC-OUTS-GENOMES.db --project-name "Accumulibacter_Pangenome" --output-dir ACCUMPANG --num-threads 10 --minbit 0.5 --mcl-inflation 5 --sensitive 