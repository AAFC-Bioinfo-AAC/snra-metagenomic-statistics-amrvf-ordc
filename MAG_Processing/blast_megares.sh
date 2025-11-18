#blast for MEGARES and VFDB
# Author: Dr. Haley Sanderson haley.sanderson@agr.gc.ca
# Copyright: Government of Canada
# License: MIT
# Version 0.1


#!/bin/bash -l

line=$1

#activate conda
source ~/miniconda3/etc/profile.d/conda.sh
#activate conda environment with BLAST installed
conda activate BLAST
#make blast custom database (only need to do once)
makeblastdb -in megares_database_v3.00.fasta -dbtype nucl -out megares
#run blast
blastn -query "$line".fa -db megares -out "$line"_megares_output -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qseq sseq"
