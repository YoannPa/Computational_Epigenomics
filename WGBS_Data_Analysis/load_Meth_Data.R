
##IMPORTS
library(data.table)

##FUNCTIONS

# read.chrs.meth ###############################################################

#' Parallel read chromosomes methylation from a specific sample.
#' 
#' @param sample      A \code{character} matching a sample name folder
#'                    name.
#' @param autosomes   An \code{integer} vector listing autosomes for which the
#'                    methylation data should be read.
#' @value A \code{list} of data.tables of CpGs, 1 dataframe by chromosome. Each
#'                    data.table in the list contains 4 columns:
#'                    1 - the start position (0-based) of the CpG on the
#'                    chromosome, 2 - the SNP probability of the CpG,
#'                    3 - the depth of coverage of methylated reads,
#'                    4 - the depth of coverage of unmethylated reads.
#' @author Yoann Pageaud.

read.chrs.meth<-function(sample, autosomes){
  lapply(autosomes, function(chr){
    fread(paste0(sample,"/",chr,".bed.gz"))
  })
}

# load.meth.data ###############################################################

#' Loads Methylation Data from a specific bisulfite sequencing dataset.
#' 
#' @param dataset.dir Full path to a folder containing the samples of interest
#'                    as folders.
#' @param samples     A \code{character} vector of sample identifiers matching
#'                    subfolders names.
#' @param autosomes   An \code{integer} vector listing autosomes for which the
#'                    methylation data should be loaded.
#' @param ncores      An \code{integer} specifying the number of cores or
#'                    threads to speed up the loading. By default equals 1.
#' @value A \code{list} of data.tables of CpGs, 1 data.table by sample. Each
#'                    data.table contains 4 columns:
#'                    chr - the chromosome, pos - the positionof the CpG on the
#'                    chromosome  (0-based), N - the total depth of coverage
#'                    (methylated + unmethylated reads), X - the depth of
#'                    coverage of methylated reads.
#' @author Yoann Pageaud.

load.meth.data<-function(dataset.dir, samples, autosomes){
  #Get current working directory
  origin_wd<-getwd()
  #Set the dataset directory as the working directory
  setwd(dataset.dir)
  cat("Loading methylation data for...\n")
  List.sample.df<-lapply(samples, function(sample){
    cat(paste0("\t",sample,"\n"))
    #Read meth data chromosome wize
    List.chr.df<-read.chrs.meth(sample = sample, autosomes = autosomes)
    #Rbind chromsome dataframes into one
    raw.df<-rbindlist(l = List.chr.df,idcol = "chr")
    #Format final data.table
    raw.df[,c("chr","pos","N","X"):=list(chr,V1+1,V3+V4,V3)][
      ,c("chr","pos","N","X"),with=FALSE]
  })
  #Set back the original working directory
  setwd(origin_wd)
  #return list of dataframes
  return(List.sample.df)
}
