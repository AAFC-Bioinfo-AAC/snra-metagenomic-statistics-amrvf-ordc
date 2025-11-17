# REad-based detection of AMR genes and VFs: read mapping with VFDB and KMA
# Author: Dr. Haley Sanderson haley.sanderson@agr.gc.ca
# Copyright: Government of Canada
# License: MIT
# Version 0.1

#mapping_VFDB_KMA.sh

#!/bin/bash -l

#activate conda
source ~/miniconda3/etc/profile.d/conda.sh
#activate conda environment with KMA installed
conda activate metareadcounts
#input variable in assisting script corresponds to line variable in this script
line=$1
#make index for VFDB database (only needs to be done once)
#kma index -i VFDB_setA_nt.fas.gz -o vfdb_exp_kma
#run KMA
kma -ipe "$line"_R1.atria.fq.gz "$line"_R2.atria.fq.gz -o KMA_output/"$line"_VFDB -t_db vfdb_exp_kma -mem_mode -ID 80
