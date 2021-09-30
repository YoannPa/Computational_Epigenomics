#!/usr/bin/env Rscript
library("optparse")

opt_parser = OptionParser(description = "
##_Create pericentromeric annotation from satellites and PICS_##################
Version = '0.0.1'
Date = '2021-09-06'
Author = 'Yoann PAGEAUD'
Maintainer = 'Yoann PAGEAUD (yoann.pageaud@gmail.com)'
Dependencies = c('R version 4.0.5 (2021-03-31)',
'RStudio Version 1.4.1106 - © 2009-2021','!!Add dependencies here!!')
Description = 'script description here'
################################################################################
");
if (is.null(parse_args(opt_parser)$file) != TRUE){print_help(opt_parser);quit()}

##IMPORTS
setwd("~/")
Imports = c("data.table")
invisible(lapply(Imports, library, character.only = T))

##PARAMETERS
#Load DT.selected.satellites
DT.selected.satellites <- readRDS("~/DT.selected.satellites_without_beta.RDS")
#Load selected.PIR4
selected.PIR4 <- fread(file = "~/hg19_PICS.bed")
#Load centromeres
centromeres <- fread("~/hg19_centromeres.bed")


##ANALYSIS
#Bring Satellites and Pericentromeric Interspersed Composite Sequences (PICS)
selected.Sat.and.PICS <- rbind(DT.selected.satellites[, c(
  "chromosome", "start", "end", "name"), ], selected.PIR4[, c(
    "chromosome", "start", "end", "name"), ])
fwrite(
  x = selected.Sat.and.PICS, file =
    "~/hg19_selected_Satellites_and_PICS_to_define_pericentromeric_regions.bed",
  sep = "\t")

#Compute coordinates of pericentromeric regions
ls.pericentromeric <- lapply(X = centromeres$chromosomes, FUN = function(chr){
  #Get start on p arm
  p.start <- selected.Sat.and.PICS[chromosome == chr][
    selected.Sat.and.PICS[chromosome == chr, .I[which.min(start)]]]$start
  #Get end on p arm
  p.end <- centromeres[chromosomes == chr]$start
  #Get start on q arm
  q.start <- centromeres[chromosomes == chr]$end
  #Get end on q arm
  q.end <- selected.Sat.and.PICS[chromosome == chr][
    selected.Sat.and.PICS[chromosome == chr, .I[which.max(end)]]]$end
  #Make chromosome data.table
  chr.dt <- data.table(
    "Chromosome" = chr, "Start" = c(p.start, q.start), "End" = c(p.end, q.end))
})
#Remove all coordinates where start is higher than end
DT.pericentromeric <- rbindlist(ls.pericentromeric)
DT.pericentromeric <- DT.pericentromeric[Start < End]
#Save pericentromeric annotation
fwrite(x = DT.pericentromeric, file = "~/hg19_pericentromeric_regions.bed",
       sep = "\t")