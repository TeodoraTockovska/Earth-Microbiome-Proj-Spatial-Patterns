#!/bin/bash
#SBATCH --account=def-cottenie
#SBATCH --time=01:30:00
#SBATCH --mem=64000 # requested memory (in MB)
#SBATCH --ntasks-per-node=16 # number of threads

module load singularity
tsvfile=$1
basename=`echo $tsvfile | sed 's/egg_samples_s__/egg_samples_/' | sed 's/.tsv//'`

# Subset the .qza file
singularity exec -B /project/def-cottenie/ttockovs/outputs:/outputs \
  -B /project/def-cottenie/ttockovs/QIIME2_installation:/inputs qiime2-2019.10.sif \
  qiime feature-table filter-samples \
  --i-table /outputs/feature-table.qza \
  --m-metadata-file /inputs/$tsvfile \
  --o-filtered-table /outputs/$basename.qza
  
# Drop zeros
singularity exec -B /project/def-cottenie/ttockovs/outputs:/outputs \
  -B /project/def-cottenie/ttockovs/QIIME2_installation:/inputs qiime2-2019.10.sif \
  qiime feature-table filter-features \
  --i-table /outputs/$basename.qza \
  --p-min-frequency 1 \
  --o-filtered-table /outputs/${basename}_filtered.qza