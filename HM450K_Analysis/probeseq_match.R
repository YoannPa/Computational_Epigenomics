##IMPORTS
Imports = c("parallel", "stringr")
lapply(Imports, library, character.only = T)

##FUNCTIONS

#' Counts HM450K probe sequence matches on hg19 genome assembly scaffolds.
#'
#' @param data         A \code{data.frame} containing at least a column 'chr'
#'                     with chromosome names and a column
#'                     'forward.genomic.sequence' with the sequence the probe is
#'                     expected to match.
#' @param ls.scaffold  A \code{character} list of scaffold sequences. Each
#'                     scaffold is usually a chromosome. The list must be named.
#'                     Each scaffold sequence must be given in a single string.
#' @param unique.scaff A \code{logical} to specify whether matches should be
#'                     looked for on the expected scaffold (unique.match = TRUE)
#'                     , or if all scaffolds should be checked
#'                     (unique.match = FALSE).
#' @param ncores       An \code{integer} specifying the number of cores/threads
#'                     to be used for parallel computing (Default: ncores = 1).
#' @return A \code{list} of the same length than the number of HM450K probes
#'         provided. The list contain the number of match found on each scaffold
#'         for each probe.
#' @author Yoann Pageaud.
#' @export
#' @examples
#' @references

probeseq.match<-function(data, ls.scaffold, unique.scaff, ncores = 1){
  mclapply(X = seq(nrow(data)), mc.cores = ncores, FUN = function(i){
    matches<-lapply(X = seq_along(ls.scaffold), FUN = function(j){
      if(unique.scaff){
        if(names(ls.scaffold)[j] == data[["chr"]][i]){
          str_count(string = ls.scaffold[[j]], pattern =
                      data[["forward.genomic.sequence"]][i])
        }
      } else {
        str_count(string = ls.scaffold[[j]], pattern =
                    data[["forward.genomic.sequence"]][i])
      }
    })
  })
}
