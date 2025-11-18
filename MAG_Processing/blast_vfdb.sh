# BLAST MAGS against VFDB
# Author: Dr. Haley Sanderson haley.sanderson@agr.gc.ca
# Copyright: Government of Canada
# License: MIT
# Version 0.1


#!/bin/bash -l


line=$1
#activate conda
source ~/miniconda3/etc/profile.d/conda.sh
#activate conda environment with BLAST
conda activate BLAST
#make custom database (only needs to be done once)
#makeblastdb -in VFDB_setB_nt.fas.gz -dbtype nucl -out vfdb_full
#run BLAST
blastn -query "$line".fa -db vfdb_full -out "$line"_vfdb_output -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qseq sseq"
