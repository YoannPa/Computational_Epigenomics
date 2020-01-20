# Computational Epigenomics
Contains various functions related to computational epigenomics data analysis.

## HM450K Analysis
* `probeseq_match.R` - Counts HM450K probe sequence matches on hg19 genome assembly scaffolds.  

## WGBS Data Analysis
* `extract_genomic_regions.R` - Reads a file and extracts 'chromosome', 'start', and 'end' columns.  
* `get_dataset_metrics.R` - Computes metrics from a coverage and a beta-values matrices.  
* `get_meth_data.R` - Retrieves methylation calls from a given directory, automatically remove SNPs from the data, format and save them in bed files and compute various statistics on it.  
* `handle_directories.R` - dependency of `get_meth_data.R`.  
* `load_Meth_Data.R` - Loads Methylation Data from a specific bisulfite sequencing dataset.  
* `make_beta_matrix.R` - Creates a matrix of beta values from a bisulfite sequencing dataset.  
* `make_coverage_matrix.R` - Creates a matrix of coverage values from a bisulfite sequencing dataset.  
* `methylome_tiling.R` - Calculates the average methylation value of genome tiles.  
* `OTP_Quality_Control_Reports` - Shiny App giving quality control metrics on the bisulfite sequencing dataset.  
