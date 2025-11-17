# REad-based detection of AMR genes and VFs: read mapping with KMA and MEGARes
# Author: Dr. Haley Sanderson haley.sanderson@agr.gc.ca
# Copyright: Government of Canada
# License: MIT
# Version 0.1

#!/bin/bash -l


#activate conda
source ~/miniconda3/etc/profile.d/conda.sh
#activate conda environment with kma installed
conda activate metareadcounts
#input variable in assisting script corresponds to line variable in this script
line=$1
#make index of MEGARes (only needs to be done once)
#kma index -i megares_database_v3.00.fasta -o megaresv3_kma
#Run KMA
kma -ipe "$line"_R1.atria.fq.gz "$line"_R2.atria.fq.gz -o output/"$line" -t_db megaresv3_kma -mem_mode -ID 80
