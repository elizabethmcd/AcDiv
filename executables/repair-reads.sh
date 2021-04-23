#! /bin/bash 

# arguments

in1=$1
in2=$2

in1name=$(basename $in1 .fastq).rp.qced.fastq
in2name=$(basename $in2 .fastq).rp.qced.fastq

cd /home/GLBRCORG/emcdaniel/EBPR/AcDiv/metagenomes/Danish-WWTPs/2018_repaired_fastqs

/opt/bifxapps/bbmap-38.32/repair.sh in=$in1 in2=$in2 out=$in1name out2=$in2name 