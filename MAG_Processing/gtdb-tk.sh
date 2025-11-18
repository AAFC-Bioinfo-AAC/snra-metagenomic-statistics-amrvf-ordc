#taxonomic Classification of quality MAGs with GTDB-tk
# Author: Dr. Haley Sanderson haley.sanderson@agr.gc.ca
# Copyright: Government of Canada
# License: MIT
# Version 0.1


#!/bin/bash -l

#activate conda
source ~/miniconda3/etc/profile.d/conda.sh
#activate conda environment with gtdb-tk installed
conda activate GTDB-tk

#run GTDB-tk on all the mags in the input folder
gtdbtk classify_wf --genome_dir mags-fasta --out_dir GTDB_class_output --skip ani_screen
