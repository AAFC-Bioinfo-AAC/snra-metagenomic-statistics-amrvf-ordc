#!/usr/bin/env bash


##########################
# USER EDIT THESE
##########################
RAWDATA_DIR="/path/to/Raw_data"           # raw fastq location
WKDIR="/path/to/wdir"                     # working directory
REF_DIR="/path/to/Reference/greengene2"   # greengene2 reference dir (qza & taxonomy)
MANIFEST="${WKDIR}/manifest.csv"          # full path to manifest.csv
CONDA_ENV="qiime2-2023.7"                 # qiime2 conda env name
CPUS=16
##########################

LOGDIR="${WKDIR}/logs"
mkdir -p "${LOGDIR}"
mkdir -p "${WKDIR}/qiime2_output"
cd "${WKDIR}"

echo "Starting pipeline at $(date)" | tee "${LOGDIR}/pipeline_start.log"

# load lab config if available (your environment)

source /path/to/.config.sh
#or source ~/.bashrc


# Activate qiime2 environment

 conda activate "${CONDA_ENV}"



##############################
# Step 1: detect & trim adapters with atria
##############################
echo "Step1: Atria adapter detect & trim" | tee -a "${LOGDIR}/pipeline.log"
# run detect (non-fatal)
cd "${RAWDATA_DIR}"
atria --detect-adapter -r *_R1.fastq.gz -R *_R2.fastq.gz > "${LOGDIR}/atria_detect.txt" 2>&1 || true

# remove adapters to output directory
OUT_ATRIA="${WKDIR}/16S_atria_trimmed_sequence_file"
mkdir -p "${OUT_ATRIA}"
atria -a AGATCGGAAGAGCACA -A AGATCGGAAGAGCGTC -r *_R1.fastq.gz -R *_R2.fastq.gz -o "${OUT_ATRIA}" \
  > "${LOGDIR}/atria_trim.log" 2>&1

##############################
# Step 2: Import to QIIME2 (manifest)
# Manifest example lines (header + sample rows):
# sample-id,absolute-filepath,direction
# S1,/full/path/to/Raw_data/S1_R1.fastq.gz,forward
# S1,/full/path/to/Raw_data/S1_R2.fastq.gz,reverse
##############################
echo "Step2: Import into QIIME2" | tee -a "${LOGDIR}/pipeline.log"
mkdir -p qiime2_output
qiime tools import \
  --type 'SampleData[PairedEndSequencesWithQuality]' \
  --input-path "${MANIFEST}" \
  --output-path qiime2_output/16S_PE.qza \
  --input-format PairedEndFastqManifestPhred33 &> "${LOGDIR}/qiime_import.log"

##############################
# Step 3: Trim primers (two-step using cutadapt)
##############################
echo "Step3: Primer trimming (cutadapt)" | tee -a "${LOGDIR}/pipeline.log"
cd "${WKDIR}/qiime2_output"

# primers
FWD="GTGYCAGCMGCCGCGGTAA"
REV="CGYCAATTYMTTTRAGTTT"
REVCOMP_FWD="TTACCGCGGCKGCTGRCAC"
REVCOMP_REV="AAACTAAAKRAATTRCG"

qiime cutadapt trim-paired \
  --i-demultiplexed-sequences 16S_PE.qza \
  --p-cores ${CPUS} \
  --p-front-f "${FWD}" \
  --p-front-r "${REV}" \
  --o-trimmed-sequences 16S_PE.primer.trimmed.qza \
  --verbose &> "${LOGDIR}/primer_trimming_step1.log"

qiime cutadapt trim-paired \
  --i-demultiplexed-sequences 16S_PE.primer.trimmed.qza \
  --p-cores ${CPUS} \
  --p-front-f "${REVCOMP_FWD}" \
  --p-front-r "${REVCOMP_REV}" \
  --o-trimmed-sequences 16S_PE.primer.trimmed2.qza \
  --verbose &> "${LOGDIR}/primer_trimming_step2.log"

##############################
# Step 4: Inspect demux quality and decide truncation
##############################
echo "Step4: Summarize demux for QC" | tee -a "${LOGDIR}/pipeline.log"
qiime demux summarize \
  --i-data 16S_PE.primer.trimmed2.qza \
  --o-visualization 16S_PE.primer.trimmed2.qzv &> "${LOGDIR}/demux_summary.log"

echo "Demux summary created: ${WKDIR}/qiime2_output/16S_PE.primer.trimmed2.qzv"
echo "Open the .qzv with qiime tools view or upload to view.qiime2.org to choose truncation lengths." | tee -a "${LOGDIR}/pipeline.log"

##############################
# Step 5: DADA2 denoise (defaults here; edit truncation lengths below)
##############################
echo "Step5: DADA2 denoising" | tee -a "${LOGDIR}/pipeline.log"
# Default truncation values (edit as needed after inspecting qzv)
TRIM_LEFT_F=0
TRIM_LEFT_R=0
TRUNC_F=229
TRUNC_R=187

mkdir -p DADA2_denoising_output
qiime dada2 denoise-paired \
  --p-n-threads ${CPUS} \
  --i-demultiplexed-seqs 16S_PE.primer.trimmed2.qza \
  --p-trim-left-f ${TRIM_LEFT_F} \
  --p-trunc-len-f ${TRUNC_F} \
  --p-trim-left-r ${TRIM_LEFT_R} \
  --p-trunc-len-r ${TRUNC_R} \
  --output-dir DADA2_denoising_output \
  --verbose &> "${LOGDIR}/DADA2_denoising.log"

##############################
# Step 6: Export DADA2 outputs
##############################
echo "Step6: Exporting outputs" | tee -a "${LOGDIR}/pipeline.log"
qiime metadata tabulate --m-input-file DADA2_denoising_output/denoising_stats.qza --o-visualization DADA2_denoising_output/denoising_stats.qzv

# export representative seqs and table
qiime tools export --input-path DADA2_denoising_output/representative_sequences.qza --output-path representative_sequences
qiime tools export --input-path DADA2_denoising_output/table.qza --output-path feature-table

# convert biom to tsv if biom exists
if [ -f feature-table/feature-table.biom ]; then
  biom convert -i feature-table/feature-table.biom -o feature-table/feature-table.tsv --to-tsv
fi

##############################
# Step 7: Taxonomy assignment (train classifier on V4-V5 or use pre-trained)
##############################
echo "Step7: Taxonomy assignment (greengene2 V4-V5)" | tee -a "${LOGDIR}/pipeline.log"
# Names under REF_DIR: ref fasta variable names may differ on your system
REF_FASTA="${REF_DIR}/ref_seqs.fasta"        # edit to actual file name
REF_TAX="${REF_DIR}/ref_taxonomy.tsv"       # edit to actual file name
REF_PREFIX="greengenes2"         # prefix for outputs

# Extract reads from the full ref (region-specific)
qiime feature-classifier extract-reads \
  --i-sequences "${REF_FASTA}" \
  --p-f-primer "${FWD}" \
  --p-r-primer "${REV}" \
  --p-min-length 300 \
  --p-max-length 600 \
  --o-reads "${REF_PREFIX}.V4-V5.qza" &> "${LOGDIR}/extract_ref_reads.log"

qiime tools export --input-path "${REF_PREFIX}.V4-V5.qza" --output-path "${REF_PREFIX}.V4-V5"
qiime tools export --input-path "${REF_TAX}" --output-path "${REF_PREFIX}.taxonomy_export" &> "${LOGDIR}/export_ref_tax.log"

# Train classifier (can take hours)
qiime feature-classifier fit-classifier-naive-bayes \
  --i-reference-reads "${REF_PREFIX}.V4-V5.qza" \
  --i-reference-taxonomy "${REF_TAX}" \
  --o-classifier "${REF_PREFIX}.classifier_V4-V5.qza" &> "${LOGDIR}/classifier_training.log"

# import representative sequences qza (we exported earlier)
qiime tools import --type 'FeatureData[Sequence]' --input-path representative_sequences/dna-sequences.fasta --output-path representative_sequences.qza &> "${LOGDIR}/import_repseqs.log"

qiime feature-classifier classify-sklearn \
  --i-classifier "${REF_PREFIX}.classifier_V4-V5.qza" \
  --i-reads representative_sequences.qza \
  --o-classification rep_fasta.classified.gg2_V4-V5.qza &> "${LOGDIR}/classify.log"

qiime tools export --input-path rep_fasta.classified.gg2_V4-V5.qza --output-path rep_fasta.classified.gg2_V4-V5 &> "${LOGDIR}/export_classify.log"

##############################
# Step 8: Phylogenetic tree
##############################
echo "Step8: Phylogeny" | tee -a "${LOGDIR}/pipeline.log"
qiime phylogeny align-to-tree-mafft-fasttree \
  --i-sequences DADA2_denoising_output/representative_sequences.qza \
  --output-dir phylogenetic_tree \
  --p-n-threads ${CPUS} &> "${LOGDIR}/phylogenetic_tree_generation.log"

echo "Pipeline finished at $(date)" | tee -a "${LOGDIR}/pipeline_end.log"
