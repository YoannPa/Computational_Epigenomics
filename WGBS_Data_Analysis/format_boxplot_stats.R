##IMPORTS
library(parallel)

##FUNCTIONS

# list.boxplot.stats ###########################################################
#
#' @description Apply the boxplot.stats() function on groups of columns of a
#'              matrix.
#'
#' @param mat        A \code{matrix} with column names.
#' @param sample.tbl A \code{data.frame} containing sample annotations. The
#'                   first column must match the matrix colnames.
#' @param grp.column A \code{character} matching the name of the column in
#'                   sample.tbl to use for grouping matrix columns before
#'                   applying the boxplot.stats() function.
#' @param ncores     An \code{integer} specifying the number of cores or threads
#'                   to speed up the process (Default: ncores = 1).
#' @param fun        A \code{character} specifying a function to use on the
#'                   matrix columns, after grouping, before applying the
#'                   boxplot.stats() function
#'                   (Default: fun = NULL ; Supported: fun = "rowSums").
#' @value A \code{list} of vectors containing the boxplot.stats() quantile
#'        values.
#' @author Yoann Pageaud.

list.boxplot.stats<-function(mat, samples.tbl, grp.column, ncores, fun=NULL){
  mclapply(unique(samples.tbl[[grp.column]]),mc.cores = ncores,function(i){
    selection<-mat[,samples.tbl[samples.tbl[[grp.column]] == i,1]]
    if(is.null(fun)){
      boxplot.stats(as.vector(selection), do.conf = F, do.out = F)$stats
    } else {
      if(fun == "rowSums"){
        if(is.null(ncol(selection))){
          boxplot.stats(selection, do.conf = F, do.out = F)$stats
        } else {
          boxplot.stats(rowSums(selection, na.rm = T),
                        do.conf = F,do.out = F)$stats
        }
      } else { stop("Unsupported function on matrix rows.") }
    }
  })
}


# df.boxplot.stats #############################################################

#' @description Apply the boxplot.stats() function on groups of columns of a
#'              matrix and return the quantiles for each group in a melted
#'              dataframe with annotations from a samples table.
#'
#' @param mat        A \code{matrix} with column names.
#' @param sample.tbl A \code{data.frame} containing sample annotations. The
#'                   first column must match the matrix colnames.
#' @param grp.column A \code{character} matching the name of the column in
#'                   sample.tbl to use for grouping matrix columns before
#'                   applying the boxplot.stats() function.
#' @param ncores     An \code{integer} specifying the number of cores or threads
#'                   to speed up the process (Default: ncores = 1).
#' @param fun        A \code{character} specifying a function to use on the
#'                   matrix columns, after grouping, before applying the
#'                   boxplot.stats() function
#'                   (Default: fun = NULL ; Supported: fun = "rowSums").
#' @param id.vars    A \code{character} vector to specify the annotations to add
#'                   to the final dataframe by their column names in the samples
#'                   table.
#' @value A \code{data.frame} with the selected annotations and the
#'         boxplot.stats() quantiles, all by columns.
#' @author Yoann Pageaud.

df.boxplot.stats<-function(mat, samples.tbl, grp.column, ncores, fun=NULL,
                           id.vars){
  list_boxplot_stats<-list.boxplot.stats(
    mat = mat, samples.tbl = samples.tbl,grp.column = grp.column,
    ncores = ncores, fun=fun)
  ids<-rev(id.vars)
  lapply(seq_along(ids), function(i){
    if(i == 1){
      list_boxplot_stats<<-Map(data.frame,
                               samples.tbl[
                                 !duplicated(samples.tbl[[grp.column]]),
                                 ids[i]],list_boxplot_stats)
    } else {
      list_boxplot_stats<<-Map(cbind,
                               samples.tbl[
                                 !duplicated(samples.tbl[[grp.column]]),
                                 ids[i]],list_boxplot_stats)
    }
  })
  list_boxplot_stats<-lapply(list_boxplot_stats , setNames ,
                             nm = c(id.vars,"Stats"))
  df_box_stats<-do.call("rbind",list_boxplot_stats)
  return(df_box_stats)
}
