# Earth-Microbiome-Proj-Spatial-Patterns
This repository stores my scripts and data files for my Master project. I analyzed the Earth Microbiome Project metadata data to determine whether microbiome diversity varied based on distance. In other words, my work involved investigating the spatial patterns in the metadata.

## EDA_EMP_dataset.R
#### I performed the exploratory data analysis using the entire EMP metadata using this script.

## Functions.R
#### I developed many functions that were used in different scripts. I created this main script that stores all of my functions. I maintained organization and script management using this script.

## hosts_EMP_spatial_pattern_analysis.R
#### Host-associated microbiome data and metadata were analyzed within this script, as well as filtering, quality control, and the calculation of the PERMANOVA analyses.

## updated_env_PERMANOVA_R.R
#### Free-living microbiome data and metadata were analyzed within this script, as well as filtering, quality control, and the calculation of the PERMANOVA analyses.

## Statistical_Models_and_Plotting.R
#### The script contains code to perform the Fisher exact test, logistic regression models, and the generation of figures part of my thesis report.

## subset_BIOM_TSV.sh
#### This script is used to create BIOM file subsets using sample TSV files.

## subset_BIOM_QIIME2_Singularity.sh
#### I retrieve the feature BIOM table in the QZA artifact format, and subset the feature table using TSV files.
