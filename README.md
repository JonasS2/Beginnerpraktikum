## Analysis
learning R Analysis

### Load.R

read in the raw data (csv) and save as rds

### SampleSheet.R

creates a useful sample sheet with columns for all possible conditions (saved as rds)

### hcluster.R

script for hierachical clustering of the samples on log scale

### QC.R

creates a heatmap with hierachical clustering

### PCA.R

creates a bunch of PCA plots with DESeq2 transformation and dataset
excluded part cpmputes the deseqdataset which is saved as rld (takes up to 40 mins)

### ggplot.R

creates some PCA plots as well with ggplot2 as tool


### DEAnalysis.R

does the differential Analysis of a subset (Sim1/2) matches it to gene names from ensembl database via biomaRT and saves the result as an xlsx file.
In this script interactions where used to differ between Genotype, Treatment and Genotype::Treatment.


### DEAnalysis_Wu4.R

since this subset contains 5 different treatments I am not sure if the interaction terms are correct. This needs to be checked. 



The makefile contains all scripts with the necessary dependencies.
