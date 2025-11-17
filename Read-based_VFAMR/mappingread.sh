# REad-based detection of AMR genes and VFs using bowtie2 and samtools
# Author: Dr. Haley Sanderson haley.sanderson@agr.gc.ca
# Copyright: Government of Canada
# License: MIT
# Version 0.1



mappingread.sh

#!/bin/bash -l


#activate conda
source ~/miniconda3/etc/profile.d/conda.sh
#activate conda environment with bowtie2 and samtools installed
conda activate metareadcounts
#input variable in assisting script corresponds to line variable in this script
line=$1
#build index of MEGARes database (only needs to be done once)
#bowtie2-build megares_database_v3.00.fasta megaresv3

bowtie2 -x megaresv3 -1 "$line"_R1.atria.fq.gz -2 "$line"_R2.atria.fq.gz -S "$line"_AMR_output.sam
samtools view -bS "$line"_AMR_output.sam > "$line"_AMR_output.bam
samtools sort "$line"_AMR_output.bam -o "$line"_AMR_output_sorted.bam
samtools index "$line"_AMR_output_sorted.bam
samtools idxstats "$line"_AMR_output_sorted.bam > "$line"_AMR_read_counts.txt
