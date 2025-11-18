# Read mapping all the MAGS to get relative abundances
# Author: Dr. Haley Sanderson haley.sanderson@agr.gc.ca
# Copyright: Government of Canada
# License: MIT
# Version 0.1


#MAG_read_mapping.sh

#!/bin/bash -l



# the input variables from the assisting script where the first one corresponds to the variable mag and the second one corresponds to the variable sample in this script
mag=$1
sample=$2
#redirect temporary files to a subdirectory in the current directory
export TMPDIR=tmp
#activate conda
source ~/miniconda3/etc/profile.d/conda.sh
#activate conda environement with bowtie2, bwa, and samtools installed
conda activate mapping_reads

#changed from bowtie2 to bwa

#bowtie2-build bacterial_mags/"$mag".fa "$mag"_index
#bowtie2 -x "$mag"_index -1 "$sample"_R1.atria.fq.gz -2 "$sample"_R2.atria.fq.gz -S "$mag"_"$sample".sam>
bwa index bacterial_mags/"$mag".fa
bwa mem bacterial_mags/"$mag".fa "$sample"_R1.atria.fq.gz "$sample"_R2.atria.fq.gz > "$mag"_"$sample".sam
samtools view -b -o "$mag"_"$sample".bam "$mag"_"$sample".sam
samtools sort -o sorted_"$mag"_"$sample".bam "$mag"_"$sample".bam
samtools index sorted_"$mag"_"$sample".bam
samtools flagstat sorted_"$mag"_"$sample".bam > Read_mapping_counts_"$mag"_"$sample".txt
