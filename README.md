# Earth-Microbiome-Proj-Spatial-Patterns
This repository stores my scripts and data files for my Master project. I analyzed the Earth Microbiome Project metadata data to determine whether microbiome diversity varied based on distance. In other words, my work involved investigating the spatial patterns in the metadata. The file descrptions are displayed below. The first 5 scripts are displayed in the order in which they were used. 

Sample IDs were extracted for host-associated and free-living datasets from scripts #3 and #4, which were saved in TSV files. Script #6 was used to extract the microbiome files for each TSV. Then, the PERMANOVA analyses were calculated in scripts #3 and #4. Script #7 was created first but script #6 was used the most.

### 1. EDA_EMP_dataset.R
Contains the exploratory data analysis on the entire EMP metadata. The metadata is split into host-associated data and free-living data. <br />

### 2. Functions.R
This main script that stores all of my functions. It is sourced in the next 3 scripts. <br />

### 3. hosts_EMP_spatial_pattern_analysis.R
Host-associated microbiome data and metadata were analyzed within this script, as well as filtering, quality control, and the calculation of the PERMANOVA analyses. <br />

### 4. updated_env_PERMANOVA_R.R
Free-living microbiome data and metadata were analyzed within this script, as well as filtering, quality control, and the calculation of the PERMANOVA analyses. <br />

### 5. Statistical_Models_and_Plotting.R
The script contains code to perform the Fisher exact test, logistic regression models, and the generation of figures part of my thesis report. <br />

### 6. subset_BIOM_TSV.sh
This script is used to create BIOM file subsets using sample TSV files. <br />

### 7. subset_BIOM_QIIME2_Singularity.sh
Retrieve the feature BIOM table in the QZA artifact format, and subset the feature table using TSV files. <br />
