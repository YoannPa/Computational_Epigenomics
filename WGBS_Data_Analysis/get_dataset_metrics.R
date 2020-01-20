
##IMPORTS
Imports = c("doParallel","matrixStats")
lapply(Imports, library, character.only = T)

##FUNCTIONS

# get.dset.metrics #############################################################

#' @description Computes metrics from a coverage and a beta-values matrices.
#'
#' @param cov.mat   A \code{matrix} containing the coverage values, with
#'                  samples by columns and positions by rows.
#' @param beta.mat  A \code{matrix} containing the beta values, with samples by
#'                  columns and positions by rows.
#' @param autosomes An \code{integer} vector listing autosomes for which the
#'                  metrics should be computed.
#' @param ncores    An \code{integer} specifying the number of cores or threads
#'                  to speed up the loading (Default: ncores = 1).
#' @value a \code{data.frame} contaning the metrics by columns and the positions
#'        by rows. Metrics are stored in this order from left to right:
#'        1) "av.cov": the average coverage of the position across samples.
#'        2) "av.pos.cov": the average coverage of the position across samples
#'        for which coverage > 0.
#'        3) "pos.cov": The number of samples for which coverage > 0 for the
#'        position.
#'        4) "min.cov": the minimum coverage value of the position found across
#'        the samples.
#'        5) "perc05.cov": the 5th percentile value on the distribution of
#'        coverage values across samples for the position.
#'        6) "perc50.cov": the median of coverage values across samples for the
#'        position.
#'        7) "perc95.cov": the 95th percentile value on the distribution of
#'        coverage values across samples for the position.
#'        8) "pos.cov10": the number of samples for which coverage >= 10 for the
#'        position.
#'        9) "av.meth": the average methylation of the position across samples.
#'        10) "sd.meth": the standard deviation of methylation values of the
#'        positon across samples.
#' @author Yoann Pageaud.

get.dset.metrics<-function(cov.mat,beta.mat,autosomes,ncores=1){
  #Set Number of cores to be used
  registerDoParallel(cores=ncores)
  cat("Calculating coverage metrics...")
  #Create list of coverage matrices by chromosomes
  mat.cov_list<-foreach(chr=autosomes) %dopar%
    cov.mat[grepl(paste0("^",chr, ":[0-9]+"),rownames(cov.mat)),]
  #Replace all 0 in the coverage matrices by NAs
  mat.cov_list_0toNA<-foreach(chr=seq(length(autosomes))) %dopar%
    replace(mat.cov_list[[chr]], mat.cov_list[[chr]] == 0, NA)
  #Remove Loaded Matrix
  rm(cov.mat)
  #Create coverage metrics table
  cov.metrics<-foreach(chr=seq(length(autosomes))) %dopar%
    data.frame(av.cov = rowMeans(mat.cov_list[[chr]]),
               av.pos.cov = rowMeans(mat.cov_list_0toNA[[chr]],na.rm = T),
               pos.cov = rowSums(mat.cov_list[[chr]]>0),
               min.cov = rowMins(mat.cov_list[[chr]]),
               perc05.cov = apply(mat.cov_list[[chr]],1,quantile, probs=.05),
               perc50.cov = apply(mat.cov_list[[chr]],1,quantile, probs=.50),
               perc95.cov = apply(mat.cov_list[[chr]],1,quantile, probs=.95),
               pos.cov10 = rowSums(mat.cov_list[[chr]]>=10))
  #Remove lists
  rm(mat.cov_list, mat.cov_list_0toNA)
  cat("Done.\n")
  cat("Calculating methylation metrics...")
  #Create list of beta matrices by chromosomes
  mat.beta_list<-foreach(chr=autosomes) %dopar%
    beta.mat[grepl(paste0("^",chr, ":[0-9]+"),rownames(beta.mat)),]
  #Create methylation metrics table
  beta.metrics<-foreach(chr=seq(length(autosomes))) %dopar%
    data.frame(av.meth = rowMeans(mat.beta_list[[chr]],na.rm = T),
               sd.meth = rowSds(mat.beta_list[[chr]],na.rm = T))
  #Remove Beta matrix
  rm(beta.mat)
  cat("Done.\n")
  cat("Merging results...")
  #Merge lists of metrics as dataframes
  cov.metrics<-do.call(rbind,cov.metrics)
  beta.metrics<-do.call(rbind,beta.metrics)
  #Merge the two metrics dataframes into one.
  dataset.metrics<-cbind(cov.metrics,beta.metrics)
  #Remove the merged metrics
  rm(cov.metrics,beta.metrics)
  cat("Done.\n")
  #Return final metrics dataframe
  return(dataset.metrics)
}
