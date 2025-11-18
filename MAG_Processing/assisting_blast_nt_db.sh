# Assisting script to blast_nt_db.sh for BLASTing MAGS
# Author: Dr. Haley Sanderson haley.sanderson@agr.gc.ca
# Copyright: Government of Canada
# License: MIT
# Version 0.1



#!/bin/bash -l
cat "qualitybins_noid.txt" | while IFS= read -r line; do
    echo "Processing line: $line"
    sbatch blast_nt_db.sh "$line"
done
