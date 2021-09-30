#!/usr/bin/env Rscript
library("optparse")

opt_parser = OptionParser(description = "
##_Plot chromoMaps of pericentromeric repeats & regions_########################
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
Imports = c("data.table", "chromoMap")
invisible(lapply(Imports, library, character.only = T))

##PARAMETERS
#Load hg19 chromoMap
chr.chromoMap <- fread("~/hg19_chromoMap.tsv")
#Load selected PICS and Satellites
PICS_satellites <- fread(
  "~/hg19_selected_Satellites_and_PICS_to_define_pericentromeric_regions.bed")
#Load pericentromeric annotation
pericentromeric <- fread("~/hg19_pericentromeric_regions.bed")

##ANALYSIS
#Plot chromoMap of selected satellites after the cut-off
chromoMap(
  ch.files = list(chr.chromoMap),
  data.files = list(PICS_satellites[!name %like% "PIR"][, c(4, 1:4), ]),
  canvas_height = 700, canvas_width = 1100, chr_color = "lightgrey",
  chr_length = 10, legend = TRUE, lg_x = 50, lg_y = 300,
  data_based_color_map = TRUE, text_font_size = 14, data_type = "categorical",
  data_colors = list(c("red", "darkgreen")),
  title = "Selected satellites annotated along hg19 chromosomes",
  title_font_size = 16)

#Plot chromoMap of PICS and Satellites
chromoMap(
  ch.files = list(chr.chromoMap),
  data.files = list(PICS_satellites[, c(4, 1:4), ]), canvas_height = 700,
  canvas_width = 1100, chr_color = "lightgrey", chr_length = 10, legend = TRUE,
  lg_x = 50, lg_y = 300, data_based_color_map = TRUE, text_font_size = 14,
  data_type = "categorical", data_colors = list(
    c("red", "darkgreen", "steelblue", "royalblue", "blue", "purple")
  ), title = "Selected satellites and PICS annotated along hg19 chromosomes",
  title_font_size = 16, segment_annotation = TRUE)

#Plot chromoMap of pericentromeric regions
pericentromeric[, ID := .I]
chromoMap(
  ch.files = list(chr.chromoMap),
  data.files = list(pericentromeric[, c(4,1:3), ]), canvas_height = 700,
  canvas_width = 1100, chr_color = "lightgrey", chr_length = 10, legend = TRUE,
  lg_x = 50, lg_y = 300, text_font_size = 14,
  title = "Pericentromeric regions annotated along hg19 chromosomes",
  title_font_size = 16, segment_annotation = TRUE)
