
##IMPORTS
library(doParallel)

##FUNCTIONS

# Read.Chroms.Meth #############################################################

#' Read Chromosomes Methylation from a specific sample.
#' 
#' @param sample      A \code{character} matching a sample identifier folder
#'                    name.
#' @param autosomes   An \code{integer} vector listing autosomes for which the
#'                    methylation data should be loaded.
#' @param ncores      An \code{integer} specifying the number of cores or
#'                    threads to speed up the loading. By default equals 1.
#' @return a \code{list} of dataframes of CpGs, 1 dataframe by chromosome. Each
#'                    dataframe in the list contains 4 columns:
#'                    1 - the position of the CpG on the chromosome,
#'                    2 - the probability of SNP of the CpG, 3 - the depth of
#'                    coverage of methylated reads, 4 - the depth of coverage of
#'                    unmethylated reads.
#' @author Yoann Pageaud.

Read.Chroms.Meth<-function(sample,autosomes,ncores=1){
  cat(paste0(sample,"\n"))
  
  #Set Number of cores to be used
  registerDoParallel(cores=ncores)
  
  List.chr.df<-foreach(chr=autosomes) %dopar%
    do.call("rbind",lapply(lapply(lapply(list(readLines(paste0(sample,"/",chr,
                                                               ".bed.gz")))[[1]],
                                         strsplit," "),unlist),as.numeric))
  return(List.chr.df)
}

# Load.Meth.Data ###############################################################

#' Load Methylation Data from a specific dataset.
#' 
#' @param dataset.dir Full path to a folder containing the samples of interest
#'                    as folders.
#' @param samples     A \code{character} vector of sample identifiers matching
#'                    subfolders names.
#' @param autosomes   An \code{integer} vector listing autosomes for which the
#'                    methylation data should be loaded.
#' @param ncores      An \code{integer} specifying the number of cores or
#'                    threads to speed up the loading. By default equals 1.
#' @return a \code{list} of dataframes of CpGs, 1 dataframe by sample. Each
#'                    dataframe in the list contains 4 columns:
#'                    chr - the chromosome, pos - the position of the CpG on the
#'                    chromosome, N - the total depth of coverage (methylated +
#'                    unmethylated reads), X - the depth of coverage of
#'                    methylated reads.
#' @author Yoann Pageaud.

Load.Meth.Data<-function(dataset.dir, samples, autosomes, ncores=1){
  #Get current working directory
  origin_wd<-getwd()
  
  #Set the dataset directory as the working directory
  setwd(dataset.dir)

  cat("Loading methylation data for...\n")
  
  List.sample.df<-list()
  #For each sample go through files by chromosomes and create list of chromosome
  # dataframes
  for(sample in samples){

    List.chr.df<-Read.Chroms.Meth(sample = sample, autosomes = autosomes,
                                  ncores = ncores)
    
    #Add Chromosome column
    List.chr.df<-Map(cbind,chr = autosomes,List.chr.df)
    
    #Rbind chromsome dataframes into one
    raw.df<-do.call("rbind",List.chr.df)
    
    #Add it to list of sample dataframes
    List.sample.df[[sample]]<-data.frame(chr = raw.df[,1], pos = raw.df[,2]+1,
                                         N = raw.df[,4] + raw.df[,5],
                                         X = raw.df[,4])
  }
  #Set back the original working directory
  setwd(origin_wd)
  
  #return list of dataframes
  return(List.sample.df)
}
