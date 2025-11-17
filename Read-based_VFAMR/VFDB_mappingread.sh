# REad-based detection of AMR genes and VFs: mapping reads to VFDB
# Author: Dr. Haley Sanderson haley.sanderson@agr.gc.ca
# Copyright: Government of Canada
# License: MIT
# Version 0.1


#VFDB_mappingread.sh

#!/bin/bash -l


#activate conda
source ~/miniconda3/etc/profile.d/conda.sh
#activate conda environment with bowtie2 and samtools installed
conda activate metareadcounts
#input variable from assiting script correspond to variable line in this script
line=$1
#build index of VFDB database (needs to be done only once)
#bowtie2-build VFDB_setB_nt.fas.gz vfdb

bowtie2 -x vfdb -1 "$line"_R1.atria.fq.gz -2 "$line"_R2.atria.fq.gz -S "$line"_vfdb_output.sam
samtools view -bS "$line"_vfdb_output.sam > "$line"_vfdb_output.bam
samtools sort "$line"_vfdb_output.bam -o "$line"_vfdb_output_sorted.bam
samtools index "$line"_vfdb_output_sorted.bam
samtools idxstats "$line"_vfdb_output_sorted.bam > "$line"_vfdb_read_counts.txt
