# Getting phylogentic tree for identified MAGS
# Author: Dr. Haley Sanderson haley.sanderson@agr.gc.ca
# Copyright: Government of Canada
# License: MIT
# Version 0.1

#!/bin/bash -l


#redirect the temporary files to a specific tmp subdirectory in your current directory
export TMPDIR=tmp
#activate conda
source ~/miniconda3/etc/profile.d/conda.sh
#activate conda environment with gtotree
conda activate gtotree
#Run gototree on the selected genomes listed in the txt file
GToTree -f bacterialbins_class.txt -H Bacteria -D -j 4 -L Domain,Phylum,Class,Order,Family,Genus,Species -F
