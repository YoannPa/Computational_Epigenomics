
##IMPORTS
source("Load_Meth_Data.R")

##FUNCTIONS

# Make.cov.matrix ################################################################

#' Creates a matrix of coverage values from a bisulfite sequencing dataset.
#'
#' @param dataset.dir Full path to a folder containing the samples of interest
#'                    as folders.
#' @param samples     A \code{character} vector of sample identifiers matching
#'                    subfolders names.
#' @param autosomes   An \code{integer} vector listing autosomes for which the
#'                    methylation data should be loaded.
#' @param ncores      An \code{integer} specifying the number of cores or
#'                    threads to speed up the loading. By default equals 1.
#' @return a \code{matrix} of coverage values with samples by columns and CpG
#'         positions by rows.
#' @author Yoann Pageaud.

Make.cov.matrix<-function(dataset.dir,samples,autosomes,ncores=1){
  #Move to MC dir
  setwd(dataset.dir)
  #Create Matrix
  list_smpl_cov<-list()
  for (smpl in samples) {
    #Read Chromosomes Methylation for each sample
    list_mc_tbl<-Read.Chroms.Meth(sample = smpl,autosomes = autosomes,
                                  ncores = ncores)

    #If 1st sample add chromosome column
    if(smpl == samples[1]) {
      list_mc_tbl<-Map(cbind, list_mc_tbl, chromosome = as.integer(autosomes))
    }
    #Concatenate chromosome matrices
    cov_mat<-do.call("rbind",list_mc_tbl)
    #If 1st sample add rownames
    if(smpl == samples[1]) {
      rownames(cov_mat)<-paste0(as.character(cov_mat[,5]),":",cov_mat[,1])
    }
    #Keep coverage values
    cov_mat<-cov_mat[,c(3,4)]
    #Sum M and U coverage and add result to list of matrices
    list_smpl_cov[[smpl]]<-cov_mat[,1]+cov_mat[,2]
  }
  #Concatenate by columns
  mat_coverage<-do.call("cbind",list_smpl_cov)
}
