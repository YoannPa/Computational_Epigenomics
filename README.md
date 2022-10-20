# Computational Genomics & Epigenomics
![GitHub last commit](https://img.shields.io/github/last-commit/YoannPa/Computational_Epigenomics)
![GitHub repo size](https://img.shields.io/github/repo-size/YoannPa/Computational_Epigenomics)
![GitHub issues](https://img.shields.io/github/issues-raw/YoannPa/Computational_Epigenomics)
![GitHub closed issues](https://img.shields.io/github/issues-closed-raw/YoannPa/Computational_Epigenomics)
![GitHub](https://img.shields.io/github/license/YoannPa/Computational_Epigenomics)  

Hi there!  
My name is Yoann Pageaud. I am a researcher and a developper in cancer bioinformatics.  
Welcome to my repository dedicated to computational genomics and epigenomics, where I share all productions that I can make publicly available.  
Her you will find:  
* Some custom annotation tracks,  
* Some tables containing data related to topics I cover,  
* Some scripts useful for genomics and epigenomics data analysis,  
* And of course, direct links to all the tools and packages I have developed.  

## [Annotations BED files](hg19_annotations/)

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

## My visualization R package
_[**BiocompR**](https://github.com/YoannPa/BiocompR) is an R package built upon ggplot2, and using data.table. It improves some visualisations commonly used in biology and genomics for data comparison and dataset exploration, introduces new kind of plots, provides a toolbox of functions to work with ggplot2 and grid objects, and ultimately, allows users to customize plots produced into publication ready figures._  

## My NCBI BLAST submission tool
_[**NCBI.BLAST2DT**](https://github.com/YoannPa/NCBI.BLAST2DT) is an R package allowing you to submit DNA sequences to NCBI BLAST servers directly from the console, to retrieve potential hits on a genome or sequence database, and to collect all results within an R data.table._  
_It makes use of the R package hoardeR to submit sequences to the NCBI BLAST API, and then parses the XML BLAST results returned to load them as an R data.table to make it more easy to query, sort, order and subset the resulting hits._  

## My methylation array quality control tool
_[**methview.qc**](https://github.com/YoannPa/methview.qc) allows you to run quality control analysis on your methylation array dataset, and to collect all results in neat ready-to-publish plots._  

## My clinical data collector tool
_[**biotab.manager**](https://github.com/YoannPa/biotab.manager) is an R package allowing you to download, manage, subset, and aggregate TCGA patients clinical data (biotabs) from the GDC portal. The package is built upon TCGAbiolinks to query TCGA databases, and makes use of R data.table handle queries results._  

## My formatting tool for large data tables
_[**DTrsiv**](https://github.com/YoannPa/DTrsiv) is an R package containing a collection of R data.table functions available to quickly and easily clean your data._  

## [Pericentromeric regions annotation](Pericentromeric_regions_annotation/)

## [Subtelomeric regions annotation](Subtelomeric_regions_annotation/)
