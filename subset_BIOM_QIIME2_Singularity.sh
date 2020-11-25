#!/bin/bash
#SBATCH --account=def-cottenie
#SBATCH --time=01:30:00
#SBATCH --mem=64000 # requested memory (in MB)
#SBATCH --ntasks-per-node=16 # number of threads

module load singularity
# feature-table requires a .qza file so doing that conversion here
singularity exec -B /project/def-cottenie/ttockovs/outputs:/outputs \
  -B /project/def-cottenie/ttockovs/QIIME2_installation:/inputs qiime2-2019.10.sif \
  qiime tools import \
  --input-path /inputs/emp_deblur_90bp.qc_filtered.biom \
  --type 'FeatureTable[Frequency]' \
  --input-format BIOMV210Format \
  --output-path /outputs/feature-table.qza

# now subset the .qza file
singularity exec -B /project/def-cottenie/ttockovs/outputs:/outputs \
  -B /project/def-cottenie/ttockovs/QIIME2_installation:/inputs qiime2-2019.10.sif \
  qiime feature-table filter-samples \
  --i-table /outputs/feature-table.qza \
  --m-metadata-file /inputs/SampleIDs_Excreta_Birds.tsv \
  --o-filtered-table /outputs/excreta.qza
  
# Drop zeros
singularity exec -B /project/def-cottenie/ttockovs/outputs:/outputs \
  -B /project/def-cottenie/ttockovs/QIIME2_installation:/inputs qiime2-2019.10.sif \
  qiime feature-table filter-features \
  --i-table /outputs/excreta.qza \
  --p-min-frequency 1 \
  --o-filtered-table /outputs/excreta_2.qza
