#nt_db database search
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
#database should have already been downloaded
#run BLAST
blastn -query "$line".fa -db nt_db -out "$line"_nt_db_output -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qseq sseq"
