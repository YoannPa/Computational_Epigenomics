#!/usr/bin/env Rscript
library("optparse")

opt_parser = OptionParser(description = "
##_Make_subtelomeres_annotation_on_hg19_########################################
Version = '0.0.1'
Date = '2021-05-14'
Author = 'Yoann PAGEAUD'
Maintainer = 'Yoann PAGEAUD (yoann.pageaud@gmail.com)'
Dependencies = c('R version 4.0.2 (2020-06-22),
'RStudio Version 1.2.5019 – © 2009-2019','!!Add dependencies here!!')
Description = 'script description here'
################################################################################
");
if (is.null(parse_args(opt_parser)$file) != TRUE){print_help(opt_parser);quit()}

##IMPORTS
setwd("~/Computational_Epigenomics/Subtelomeric_regions_annotation/")
Imports = c("data.table", "seqinr", "hoardeR", "NCBI.BLAST2DT")
invisible(lapply(Imports, library, character.only = T))

##PARAMETERS
#Load hg19 genomic gap annotation
gaps <- fread("gap.txt.gz")
colnames(gaps)[2:4] <- c("chromosomes", "start", "end")
#Keep telomeres and centromeres
telomeres <- gaps[V8 == "telomere", c(2:4), ]
#No Telomeres annotated in hg19 on chr 17:
# https://www.biostars.org/p/76193/#76324
# https://www.biostars.org/p/72730/
# https://www.ncbi.nlm.nih.gov/grc/human/issues/HG-1049
# https://www.nature.com/articles/ng0198-76

#Save table of BLAST results
subtelomere_blasts <- fread("subtelomere_selected_BLAST_hits.csv")

#Get chromosome name
subtelomere_blasts[, chromosome := gsub(
  pattern = "^Homo sapiens chromosome (\\d*[0-9XY]), GRCh37.p13 Primary Assembly",
  replacement = "chr\\1", x = hits_def)]
#Get chromosome arm information
subtelomere_blasts[, arm := gsub(
  pattern = "^[0-9XY]\\d*(pY|*)([pq])tel(_short|*)_1(-|_)500[K]*(_new|shortallele|*)_[14]_(12|3)_12(v2|*)_(hg19orientation_|*)(proximal|distal)\\/1.xml$",
  replacement = "\\2", x = file_name)]
#Get subtelomere sequence extremity
subtelomere_blasts[, extremity := gsub(
  pattern = "^[0-9XY]\\d*(pY|*)[pq]tel(_short|*)_1(-|_)500[K]*(_new|shortallele|*)_[14]_(12|3)_12(v2|*)_(hg19orientation_|*)(proximal|distal)\\/1.xml$",
  replacement = "\\8", x = file_name)]

#Set chromosomes as factor
chrs <- c(unique(telomeres$chromosomes),"chr17")
chrs <- sort(factor(x = chrs, levels = paste0("chr", c(1:22,"X","Y"))))

#Create subtelomere bed annotation
ls.subtel <- lapply(X = chrs, function(i){
  #Get subtelomere from p arm
  p.dt <- subtelomere_blasts[
    chromosome == i & arm == "p" & subject_frame == 1,
    c("chromosome", "subject_start", "subject_end", "extremity")]
  #For distal part of the sequence
  if(nrow(p.dt[extremity == "distal"]) == 1){
    #Check if the distal position is before all proximal positions on the p arm
    if(p.dt[, .SD[which.min(subject_end)]]$extremity == "distal"){
      p.start <- p.dt[extremity == "distal"]$subject_start - 1
    } else { #Unexploitable alignment, get p arm telomere end position
      p.start <- telomeres[chromosomes == i][1]$end
    }
  } else if(nrow(p.dt[extremity == "distal"]) > 1){
    p.start <- min(p.dt[extremity == "distal"]$subject_start) - 1
  } else {
    if(nrow(telomeres[chromosomes == i]) == 0){
      #Chr17 in hg19 does not have any telomere defined
      p.start <- 0
    } else { p.start <- telomeres[chromosomes == i][1]$end }  
  }
  #For proximal part of the sequence
  if(nrow(p.dt[extremity == "proximal"]) == 1){
    p.end <- p.dt[extremity == "proximal"]$subject_end
  } else if(nrow(p.dt[extremity == "proximal"]) > 1){
    p.dt[extremity == "proximal", dist := abs(subject_end - p.start - 500000)]
    p.end <- p.dt[extremity == "proximal", .SD[which.min(dist)]]$subject_end
    # p.end <- max(p.dt[extremity == "proximal"]$subject_end)
  } else { p.end <- p.start + 500000 }
  #Get subtelomere from q arm
  q.dt <- subtelomere_blasts[
    chromosome == i & arm == "q" & subject_frame == -1,
    c("chromosome", "subject_start", "subject_end", "extremity")]
  #For distal part of the sequence
  if(nrow(q.dt[extremity == "distal"]) == 1){
    q.end <- q.dt[extremity == "distal"]$subject_start
  } else if(nrow(q.dt[extremity == "distal"]) > 1){
    q.end <- max(q.dt[extremity == "distal"]$subject_start)
  } else { 
    if(nrow(telomeres[chromosomes == i]) == 0){
      #Chr17 in hg19 does not have any telomere defined, set chr17 end position
      q.end <- 81195210
    } else { q.end <- telomeres[chromosomes == i][2]$start }
  }
  #For proximal part of the sequence
  if(nrow(q.dt[extremity == "proximal"]) == 1){
    q.start <- q.dt[extremity == "proximal"]$subject_end - 1
  } else if(nrow(q.dt[extremity == "proximal"]) > 1){
    q.dt[extremity == "proximal",
         dist := abs(q.end - (subject_end - 1) - 500000)]
    q.start <- q.dt[extremity == "proximal",
                    .SD[which.min(dist)]]$subject_end - 1
  } else { q.start <- q.end - 500000 }
  #Create chromosome subtelomeres annotation
  sub.dt <- data.table("chromosome" = i, "start" = c(p.start, q.start),
                       "end" = c(p.end, q.end))
  return(sub.dt)
})
#Get widths of subtelomeric regions
ls.subtel <- lapply(ls.subtel, function(i){i[, width := end - start]})
hg19.subtelomeric.bed <- rbindlist(ls.subtel)
fwrite(
  x = hg19.subtelomeric.bed, file = "hg19_subtelomeres.bed", sep = "\t",
  scipen = 2)

#Remove gaps between telomeres and subtelomeres
invisible(lapply(X = chrs, function(i){
  if(nrow(telomeres[chromosomes == i]) != 0){
    #Remove gaps between p telomere end and subtelomere start
    if(telomeres[chromosomes == i]$end[1] !=
       hg19.subtelomeric.bed[chromosome == i]$start[1]){
      hg19.subtelomeric.bed[chromosome == i][["start"]][1] <<- telomeres[
        chromosomes == i]$end[1]
    }
    #Remove gaps between q subtelomere end and telomere start
    if(telomeres[chromosomes == i]$start[2] !=
       hg19.subtelomeric.bed[chromosome == i]$end[2]){
      hg19.subtelomeric.bed[chromosome == i][["end"]][2] <<- telomeres[
        chromosomes == i]$start[2]
    }
  }
}))
#Recompute new width of extended subtelomeric regions
hg19.subtelomeric.bed[, width := end - start]
#Save extended version of subtelomeric region bed file without gaps between telomeres
fwrite(
  x = hg19.subtelomeric.bed, file = "hg19_extended_subtelomeres.bed",
  sep = "\t", scipen = 2)
