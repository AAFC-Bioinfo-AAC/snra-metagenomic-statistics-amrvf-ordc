# assisting script for blast_megares.sh
# Author: Dr. Haley Sanderson haley.sanderson@agr.gc.ca
# Copyright: Government of Canada
# License: MIT
# Version 0.1


#!/bin/bash -l
cat "bacterialbins.txt" | while IFS= read -r line; do
    echo "Processing line: $line"
    sbatch blast_megares.sh "$line"
done
