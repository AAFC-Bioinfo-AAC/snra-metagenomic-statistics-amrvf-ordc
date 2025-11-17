# Running fastqc on metagenomic samples (fastq): assisting script
# Author: Dr. Haley Sanderson haley.sanderson@agr.gc.ca
# Copyright: Government of Canada
# License: MIT
# Version 0.1

#!/bin/bash -l
cat "list_MG_samples.txt" | while IFS= read -r line; do
    echo "Processing line: $line"
    sbatch fastqc.sh "$line"
done
