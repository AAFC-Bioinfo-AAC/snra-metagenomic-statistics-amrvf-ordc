# Running fastqc on metagenomic samples (fastq)
# Author: Dr. Haley Sanderson haley.sanderson@agr.gc.ca
# Copyright: Government of Canada
# License: MIT
# Version 0.1


#!/bin/bash -l


#activate conda
source ~/miniconda3/etc/profile.d/conda.sh
#activate conda environment with fastqc installed
conda activate multifastqc
#input variable from assisting script corresponds to line variable in this script
line=$1
#Run fastqc
fastqc "$line"_R1.atria.fq.gz
fastqc "$line"_R2.atria.fq.gz
