
##IMPORTS

##FUNCTIONS

# format.path ##################################################################

#' Format a path by removing given parts
#' 
#' @param path     A \code{character} corresponding to the path you want to
#'                 format.
#' @param rm.parts A \code{character} vector containing all parts you want to
#'                 remove from the path. 
#' @return a \code{character} corresponding to the formated path wanted.
#' @author Yoann Pageaud.

format.path<-function(path,rm.parts){
  
  lapply(rm.parts, function(part){
    path<<-gsub(pattern = part,replacement = "",x = path)
  })
  return(path)  
}

# get.sample.dirs ##############################################################

#' Retrieve samples directories into a matrix
#' 
#' @param PID.dir       A \code{character} matching the path to the PID folder
#'                      where samples data are stored.
#' @param sample.subdir A \code{character} matching the name of the subdirectory
#'                      below the sample directory
#' @param paired        A \code{logical} specifying whether the WGBS data are
#'                      containing paired-end reads or not.
#'                      (Supported: "blood"; "tumor00"; many others possible..).
#' @param data.type     A \code{character} matching the type of data that the
#'                      function should look for
#'                      (Default: data.type = "merged").
#' @param skip.amplicon A \code{logical] to specify whether or not to skip the
#'                      AmpliconSeq samples (Default: skip.amplicon = TRUE).
#' @param rm.parts      A \code{character} vector containing all parts you want
#'                      to remove from the path.
#' @param verbose       A \code{logical} specifying whether or not to display
#'                      sample names (Default: verbose = TRUE).
#' @return a \code{matrix} containing all the formated paths to samples data.
#' @author Yoann Pageaud.
 
get.sample.dirs<-function(PID.dir, sample.subdir, paired=T, data.type="merged",
                          skip.amplicon=TRUE, rm.parts,verbose=T){
  #Get current WD
  main.wd<-getwd()
  #Get Sample IDs
  Sample.IDs <- list.files(PID.dir)
    #Set Working Directory
  setwd(PID.dir)
  
  if(skip.amplicon == T) { #Remove Amplicon-Seq 
    Sample.IDs<-Sample.IDs[!(grepl(".+-1$",Sample.IDs))]
  }
  if(paired){
    reads.type<-"paired"
  } else {stop("reads type not supported yet.")}
  if(data.type == "merged"){data.type <-"merged-alignment"} else {
    stop("sample type not supported yet.")}
  
  list.valid.paths<-lapply(Sample.IDs, function(smpl_id){
    ##Define sample paths
    smplpath<-file.path(PID.dir, smpl_id, file.path(sample.subdir,reads.type),
                        data.type)
    #Main subdir
    smpl.main.path <- file.path(PID.dir, smpl_id, file.path(sample.subdir,
                                                            reads.type))
    if(dir.exists(smplpath)){
      setwd(smplpath)
      #Check if directory contain any file or folder
      if(length(list.files()) != 0){
        if(verbose){
          cat(paste0(smpl_id,"\n")) #Print Sample ID  
        }
        #Add the sample main path to the list of valid paths
        matrix(format.path(path = smpl.main.path,rm.parts = rm.parts),nrow = 1)
      } else {
        warning(paste("No file found for the sample",smpl_id))
        NULL
      } #Close list.file() condition on the input directory
    }
  })
  #re-set the working directory
  setwd(main.wd)
  #Make matrix
  mat.sample.dir<-do.call(rbind,list.valid.paths)
  colnames(mat.sample.dir)<- "Sample.Directories"
  mat.sample.dir
}
