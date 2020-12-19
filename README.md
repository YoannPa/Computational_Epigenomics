# Computational Epigenomics
This git repository contains functions and data useful for computational epigenomics data analysis.  

## HM450K Analysis
[Click here for HM450K content]("https://github.com/YoannPa/Computational_Epigenomics/tree/master/HM450K_Analysis").  

## WGBS Data Analysis
* `extract_genomic_regions.R` - Reads a file and extracts 'chromosome', 'start', and 'end' columns.It supports almost any format. The file can contain any amount of columns. This function does not limit to extraction of genomic regions in a bed format.  
* `get_dataset_metrics.R` - Computes metrics from a coverage and a beta-values matrices.  
* `get_meth_data.R` - Retrieves methylation calls from a given directory, automatically remove SNPs from the data, format and save them in bed files and compute various statistics on it.  
* `handle_directories.R` - dependency of `get_meth_data.R`. Retrieve samples directories into a matrix.  
* `load_Meth_Data.R` - Loads Methylation Data from a specific bisulfite sequencing dataset.  
* `make_beta_matrix.R` - Creates a matrix of beta values from a bisulfite sequencing dataset.  
* `make_coverage_matrix.R` - Creates a matrix of coverage values from a bisulfite sequencing dataset.  
* `methylome_tiling.R` - Calculates the average methylation value of genome tiles.  
* `OTP_Quality_Control_Reports` - Shiny App giving quality control metrics on the bisulfite sequencing dataset.  
