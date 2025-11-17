# REad-based detection of AMR genes and VFs
# Author: Dr. Haley Sanderson haley.sanderson@agr.gc.ca
# Copyright: Government of Canada
# License: MIT
# Version 0.1



conda activate metareadcounts
(contains bowtie2, samtools, bedtools, kma)

Download VFDB (http://www.mgc.ac.cn/VFs/download.htm; lastupdate: Fri Dec 1 19:30:02 2023)
 and MegaRESv3 (https://www.meglab.org/megares/download/)

mappingread.sh

#!/bin/bash -l
#SBATCH --export=USER,LOGNAME,HOME,MAIL,PATH=/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin
#SBATCH --job-name=map_amr
#SBATCH --account=grdi_genarcc
#SBATCH --partition=standard
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=4000M
#SBATCH --time=24:00:00

source ~/miniconda3/etc/profile.d/conda.sh
conda activate metareadcounts
line=$1
#bowtie2-build megares_database_v3.00.fasta megaresv3
bowtie2 -x megaresv3 -1 "$line"_R1.atria.fq.gz -2 "$line"_R2.atria.fq.gz -S "$line"_AMR_output.sam
samtools view -bS "$line"_AMR_output.sam > "$line"_AMR_output.bam
samtools sort "$line"_AMR_output.bam -o "$line"_AMR_output_sorted.bam
samtools index "$line"_AMR_output_sorted.bam
samtools idxstats "$line"_AMR_output_sorted.bam > "$line"_AMR_read_counts.txt

launch_mappingreads.sh

#!/bin/bash -l
cat "list_MG_samples.txt" | while IFS= read -r line; do
    echo "Processing line: $line"
    sbatch mappingreads.sh "$line"
done

VFDB_mappingread.sh

#!/bin/bash -l
#SBATCH --export=USER,LOGNAME,HOME,MAIL,PATH=/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin
#SBATCH --job-name=map_VFDB
#SBATCH --account=grdi_genarcc
#SBATCH --partition=standard
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=4000M
#SBATCH --time=24:00:00

source ~/miniconda3/etc/profile.d/conda.sh
conda activate metareadcounts
line=$1
#bowtie2-build VFDB_setB_nt.fas.gz vfdb
bowtie2 -x vfdb -1 "$line"_R1.atria.fq.gz -2 "$line"_R2.atria.fq.gz -S "$line"_vfdb_output.sam
samtools view -bS "$line"_vfdb_output.sam > "$line"_vfdb_output.bam
samtools sort "$line"_vfdb_output.bam -o "$line"_vfdb_output_sorted.bam
samtools index "$line"_vfdb_output_sorted.bam
samtools idxstats "$line"_vfdb_output_sorted.bam > "$line"_vfdb_read_counts.txt

launch_VF_DB_mappingreads.sh

#!/bin/bash -l
cat "list_MG_samples.txt" | while IFS= read -r line; do
    echo "Processing line: $line"
    sbatch VFDB_mappingreads.sh "$line"
done

#Use experimentally verified VFDB database version
VFDB_exp_mappingreads.sh

#!/bin/bash -l
#SBATCH --export=USER,LOGNAME,HOME,MAIL,PATH=/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin
#SBATCH --job-name=map_VFDB
#SBATCH --account=grdi_genarcc
#SBATCH --partition=standard
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=4000M
#SBATCH --time=24:00:00

source ~/miniconda3/etc/profile.d/conda.sh
conda activate metareadcounts
line=$1conda install -c bioconda kma
#bowtie2-build VFDB_setA_nt.fas.gz vfdbA
bowtie2 -x vfdbA -1 "$line"_R1.atria.fq.gz -2 "$line"_R2.atria.fq.gz -S "$line"_vfdb_output.sam
samtools view -bS "$line"_vfdb_output.sam > "$line"_vfdb_output.bam
samtools sort "$line"_vfdb_output.bam -o "$line"_vfdb_output_sorted.bam
samtools index "$line"_vfdb_output_sorted.bam
samtools idxstats "$line"_vfdb_output_sorted.bam > "$line"_vfdb_exp_read_counts.txt

launch_VFDB_exp_mappingreads.sh

#!/bin/bash -l
cat "list_MG_samples.txt" | while IFS= read -r line; do
    echo "Processing line: $line"
    sbatch VFDB_exp_mappingreads.sh "$line"
done

mapping_AMR_KMA.sh
#!/bin/bash -l
#SBATCH --export=USER,LOGNAME,HOME,MAIL,PATH=/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin
#SBATCH --job-name=map_VFDB
#SBATCH --account=grdi_genarcc
#SBATCH --partition=standard
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=4000M
#SBATCH --time=24:00:00

source ~/miniconda3/etc/profile.d/conda.sh
conda activate metareadcounts
line=$1

#kma index -i megares_database_v3.00.fasta -o megaresv3_kma
kma -ipe "$line"_R1.atria.fq.gz "$line"_R2.atria.fq.gz -o output/"$line" -t_db megaresv3_kma -mem_mode -ID 80


launch_AMR_KMA.sh

#!/bin/bash -l
cat "list_MG_samples.txt" | while IFS= read -r line; do
    echo "Processing line: $line"
    sbatch mapping_AMR_KMA.sh "$line"
done

mapping_VFDB_KMA.sh
#!/bin/bash -l
#SBATCH --export=USER,LOGNAME,HOME,MAIL,PATH=/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin
#SBATCH --job-name=map_VFDB
#SBATCH --account=grdi_genarcc
#SBATCH --partition=standard
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=4000M
#SBATCH --time=24:00:00

source ~/miniconda3/etc/profile.d/conda.sh
conda activate metareadcounts
line=$1

#kma index -i VFDB_setA_nt.fas.gz -o vfdb_exp_kma
kma -ipe "$line"_R1.atria.fq.gz "$line"_R2.atria.fq.gz -o KMA_output/"$line"_VFDB -t_db vfdb_exp_kma -mem_mode -ID 80

launch_VFDB_KMA.sh

#!/bin/bash -l
cat "list_MG_samples.txt" | while IFS= read -r line; do
    echo "Processing line: $line"
    sbatch mapping_VFDB_KMA.sh "$line"
done 
