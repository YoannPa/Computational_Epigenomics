
##IMPORTS
Imports = c("feather","data.table","simsalapar")
lapply(Imports, library, character.only = T)

##FUNCTIONS

# extract_bed ##################################################################

#' @description Read a file and extracts 'chromosome', 'start', and 'end'
#'              columns.
#'
#' @param file          A \code{character} specifying the path to the file
#'                      containing the genomic regions of interest.
#' @param sep           A \code{character} to be used for splitting file rows.
#'                      Defaults to the character in the set [,\t |;:] that
#'                      separates the sample of rows into the most number of
#'                      lines with the same number of fields. Use NULL or "" to
#'                      specify no separator; i.e. each line a single character
#'                      column like base::readLines does.
#' @param header        A \code{logical} value indicating whether the file
#'                      contains the names of the variables as its first line.
#'                      If missing, the value is determined from the file
#'                      format: header is set to TRUE if and only if the first
#'                      row contains one fewer field than the number of columns.
#' @param regex.bed.col A \code{character} to be used as a regular expression
#'                      for identifying columns to be kept following their
#'                      names (Default: regex.bed.col =
#'                      "^(Chromosome|chr|Start|start|End|end)$").
#' @param expected.ncol An \code{integer} specifying the number of columns
#'                      expected to be extracted by the function
#'                      (Default: expected.ncol = 3).
#' @param new.colnames  A \code{character} vector to be used for renaming the
#'                      final table
#'                      (Default: new.colnames = c("chr","start","end")).
#' @value a \code{tibble} or a \code{data.table} object.
#' @author Yoann Pageaud.

extract_bed<-function(file,sep = "auto", header = "auto",
                      regex.bed.col = "^(Chromosome|chr|Start|start|End|end)$",
                      expected.ncol = 3, new.colnames = c("chr","start","end")){
  file.reg<-tryCatch.W.E(read_feather(file))
  if(!is.data.frame(file.reg$value)){
    file.reg<-tryCatch.W.E(fread(file, sep = sep, header = header))
    if(!is.data.frame(file.reg$value)){
      stop("The regions file is neither a feather file nor a csv file.")
    } else {
      regions<-file.reg$value[,colnames(file.reg$value)[
        grepl(pattern = regex.bed.col, x = colnames(file.reg$value))],
        with = FALSE]
    }
  } else {
    regions<-file.reg$value[,colnames(file.reg$value)[
      grepl(pattern = regex.bed.col, x = colnames(file.reg$value))]]
  }
  if(ncol(regions) < expected.ncol){
    stop("missing columns! Check if the regular expression is matching all
       'chromosome', 'start' and 'end' types of columns.")
  } else if(ncol(regions) > expected.ncol){
    stop("too many columns! Check if the regular expression is matching all
       'chromosome', 'start' and 'end' types of columns.")
  } else { colnames(regions)<-new.colnames }
  return(regions)
}
