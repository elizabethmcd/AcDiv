#! /bin/bash 

in=$1
repair=$2
clean=$3


# repair PE interleaved file and output interleaved file to repair, and then QC with fastp

# repair
/opt/bifxapps/bbmap-38.32/repair.sh in=$in out=$repair 

# filter

~/bin/fastp -i $repair -o $clean --cut_tail
