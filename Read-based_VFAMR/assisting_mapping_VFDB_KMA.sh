# REad-based detection of AMR genes and VFs: assisting script for mapping_VFDB_KMA.sh
# Author: Dr. Haley Sanderson haley.sanderson@agr.gc.ca
# Copyright: Government of Canada
# License: MIT
# Version 0.1



#!/bin/bash -l
cat "list_MG_samples.txt" | while IFS= read -r line; do
    echo "Processing line: $line"
    sbatch mapping_VFDB_KMA.sh "$line"
done 
