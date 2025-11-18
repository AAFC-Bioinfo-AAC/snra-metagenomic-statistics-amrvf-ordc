# Assisting script read mapping all the MAGS to get relative abundances where we set up combinations of MAG and sample
# Author: Dr. Haley Sanderson haley.sanderson@agr.gc.ca
# Copyright: Government of Canada
# License: MIT
# Version 0.1

#!/bin/bash -l
cat "list_MG_samples.txt" | while IFS= read -r line; do
    echo "Processing line: $line"
#put name of the MAG at <MAG>
    sbatch MAG_read_mapping.sh  <MAG> "$line"
done
