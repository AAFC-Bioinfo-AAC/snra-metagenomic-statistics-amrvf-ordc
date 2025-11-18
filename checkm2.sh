# Run checkm2 on the MAGs
# Author: Dr. Haley Sanderson haley.sanderson@agr.gc.ca
# Copyright: Government of Canada
# License: MIT
# Version 0.1

#!/bin/bash -l



#get completeness and contamination for MAGs from checkm2
#activate conda
source ~/miniconda3/etc/profile.d/conda.sh
#activate conda environment with checkm2
conda activate checkm2
#run checkm2 on all the MAGs in the input folder
checkm2 predict --threads 30 --input mags-fasta --output-directory checkm2_results -x .fa --force
