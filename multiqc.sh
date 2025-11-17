# Running multiqc for metagenomic samples (fastq)
# Author: Dr. Haley Sanderson haley.sanderson@agr.gc.ca
# Copyright: Government of Canada
# License: MIT
# Version 0.1


#!/bin/bash -l

#activate conda
source ~/miniconda3/etc/profile.d/conda.sh
#activate conda environment with multiqc installed
conda activate multifastqc
#Run multiqc
multqc .
