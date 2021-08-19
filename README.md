# Computational Genomics & Epigenomics
![GitHub last commit](https://img.shields.io/github/last-commit/YoannPa/Computational_Epigenomics)
![GitHub repo size](https://img.shields.io/github/repo-size/YoannPa/Computational_Epigenomics)
![GitHub issues](https://img.shields.io/github/issues-raw/YoannPa/Computational_Epigenomics)
![GitHub closed issues](https://img.shields.io/github/issues-closed-raw/YoannPa/Computational_Epigenomics)
![GitHub](https://img.shields.io/github/license/YoannPa/Computational_Epigenomics)  

This git repository contains functions and data useful for computational genomics and epigenomics data analysis.  

## [HM450K data analysis](HM450K_data_analysis/)

## WGBS data analysis
* `extract_genomic_regions.R` - Reads a file and extracts 'chromosome', 'start', and 'end' columns.It supports almost any format. The file can contain any amount of columns. This function does not limit to extraction of genomic regions in a bed format.  
* `get_dataset_metrics.R` - Computes metrics from a coverage and a beta-values matrices.  
* `get_meth_data.R` - Retrieves methylation calls from a given directory, automatically remove SNPs from the data, format and save them in bed files and compute various statistics on it.  
* `handle_directories.R` - dependency of `get_meth_data.R`. Retrieve samples directories into a matrix.  
* `load_Meth_Data.R` - Loads Methylation Data from a specific bisulfite sequencing dataset.  
* `make_beta_matrix.R` - Creates a matrix of beta values from a bisulfite sequencing dataset.  
* `make_coverage_matrix.R` - Creates a matrix of coverage values from a bisulfite sequencing dataset.  
* `methylome_tiling.R` - Calculates the average methylation value of genome tiles.  
* `OTP_Quality_Control_Reports` - Shiny App giving quality control metrics on the bisulfite sequencing dataset.  

## NCBI BLAST submission

You can submit DNA sequences to NCBI BLAST servers directly from the console in R using the package I developped [**NCBI.BLAST2DT**](https://github.com/YoannPa/NCBI.BLAST2DT). For more details about the package click on the previous link.  

## [Subtelomeric regions annotation](Subtelomeric_regions_annotation/)

