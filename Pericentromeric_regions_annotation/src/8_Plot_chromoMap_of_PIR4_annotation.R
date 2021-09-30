#!/usr/bin/env Rscript
library("optparse")

opt_parser = OptionParser(description = "
##_Plot chromoMaps of PIR4 annotation_##########################################
Version = '0.0.1'
Date = '2021-08-27'
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
Imports = c("data.table", "chromoMap")
invisible(lapply(Imports, library, character.only = T))

##PARAMETERS
#Load chromosomes lengths
chr.length <- fread("~/hg19_chromosomes.tsv")
#Load chromosomes centers
chr.centers <- fread(
  "~/hg19_chromosomes_centers.bed")
#Load PIR4 annotation
PIR4 <- fread("~/hg19_PIR4.bed")
#Load PICS annotation after PICS-specific cut-off
PICS <- fread("~/hg19_PICS.bed")


##ANALYSIS
#Make chromoMap data.frames
chr.length[, start := 1]
chr.chromoMap <- merge(
  x = chr.length, y = chr.centers, by = "chromosomes", all.y = TRUE)[
    order(match(chromosomes, chr.length$chromosomes)), c(1,3,2,4), ]
#Save chromoMap for hg19
fwrite(x = chr.chromoMap, file = "~/hg19_chromoMap.tsv", sep = "\t")

#Plot chromoMap PIR4 sequences
PIR4.chromoMap <- PIR4[, c(6,1:3, 6), ]
chromoMap(
  ch.files = list(chr.chromoMap), data.files = list(PIR4.chromoMap),
  canvas_height = 700, canvas_width = 1100, chr_color = "lightgrey",
  chr_length = 10, legend = TRUE, lg_x = 50, lg_y = 300,
  data_based_color_map = T, text_font_size = 14,
  data_type = "categorical", data_colors = list(
    c("red", "royalblue", "orange", "darkgreen")),
  title = "PICS annotated along hg19 chromosomes",
  title_font_size = 16)

#Plot chromoMap heatmap alignment score
PIR4.chromoMap <- PIR4[, c(6,1:3, 7), ]
chromoMap(
  ch.files = list(chr.chromoMap), data.files = list(PIR4.chromoMap),
  canvas_height = 700, canvas_width = 1100, chr_color = "lightgrey",
  chr_length = 10, legend = TRUE, lg_x = 50, lg_y = 300,
  data_based_color_map = TRUE, text_font_size = 14,
  data_type = "numeric", data_colors = list(c("orange", "red", "purple")),
  heat_map = TRUE, chr_text = TRUE,
  title = "PICS alignments scores along hg19 chromosomes",
  title_font_size = 16)

#Plot chromoMap PIR4 sequences aligment frames
PIR4.chromoMap <- PIR4[, c(6,1:3, 5), ]
PIR4.chromoMap[frame == "-", frame := "Reverse"]
PIR4.chromoMap[frame == "+", frame := "Forward"]
chromoMap(
  ch.files = list(chr.chromoMap), data.files = list(PIR4.chromoMap),
  canvas_height = 700, canvas_width = 1100, chr_color = "lightgrey",
  chr_length = 10, legend = TRUE, lg_x = 50, lg_y = 300,
  data_based_color_map = T, text_font_size = 14,
  data_type = "categorical", data_colors = list(
    c("blue", "red")),
  title = "PICS alignments frame along hg19 chromosomes",
  title_font_size = 16)

#Plot chromoMap PICS
PICS.chromoMap <- PICS[, c(6,1:3, 6), ]
chromoMap(
  ch.files = list(chr.chromoMap), data.files = list(PICS.chromoMap),
  canvas_height = 700, canvas_width = 1100, chr_color = "lightgrey",
  chr_length = 10, legend = TRUE, lg_x = 50, lg_y = 300,
  data_based_color_map = T, text_font_size = 14,
  data_type = "categorical", data_colors = list(
    c("red", "royalblue", "orange", "darkgreen")),
  title = "PICS annotated along hg19 chromosomes after PICS-specific score cut-off",
  title_font_size = 16)
