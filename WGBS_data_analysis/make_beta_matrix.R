
##IMPORTS
source("Load_Meth_Data.R")

##FUNCTIONS

# Make.b.matrix ################################################################

#' Creates a matrix of beta values from a bisulfite sequencing dataset.
#' 
#' @param dataset.dir Full path to a folder containing the samples of interest
#'                    as folders.
#' @param samples     A \code{character} vector of sample identifiers matching
#'                    subfolders names.
#' @param autosomes   An \code{integer} vector listing autosomes for which the
#'                    methylation data should be loaded.
#' @param ncores      An \code{integer} specifying the number of cores or
#'                    threads to speed up the loading. By default equals 1.
#' @return a \code{matrix} of beta values with samples by columns and CpG
#'         positions by rows.
#' @author Yoann Pageaud.

Make.b.matrix<-function(dataset.dir,samples,autosomes,ncores=1){
  #Move to MC dir
  setwd(dataset.dir)
  
  #Create Matrix
  list_smpl_beta<-list()
  for (smpl in samples) {
    #Read Chromosomes Methylation for each sample
    list_mc_tbl<-Read.Chroms.Meth(sample = smpl,autosomes = autosomes,
                                  ncores = ncores)
    
    #If 1st sample add chromosome column
    if(smpl == samples[1]) {
      list_mc_tbl<-Map(cbind, list_mc_tbl, chromosome = as.integer(autosomes))
    }
    
    #Concatenate chromosome matrices
    beta_mat<-do.call("rbind",list_mc_tbl)
    
    #If 1st sample add rownames
    if(smpl == samples[1]) {
      rownames(beta_mat)<-paste0(as.character(beta_mat[,5]),":",beta_mat[,1])
    }
    
    #Keep coverage values
    beta_mat<-beta_mat[,c(3,4)]
    
    #Calculate beta value and add result to list of matrices
    list_smpl_beta[[smpl]]<-beta_mat[,1]/(beta_mat[,1]+beta_mat[,2])
  }
  
  #Concatenate samples data by columns
  mat_beta<-do.call("cbind",list_smpl_beta)
}
