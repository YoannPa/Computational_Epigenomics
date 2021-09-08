#!/usr/bin/env Rscript
library("optparse")

opt_parser = OptionParser(description = "
##_Plot chromoMaps of alpha, beta, CAGGG, HSATII and GGGCAAAAGCCG_##############
Version = '0.0.1'
Date = '2021-09-02'
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
#Load GGGCAAAAGCCG BLAST results
dt.GGGCAAAAGCCG <- fread("~/GGGCAAAAGCCG_BLAST_results.csv")
#Load alpha, beta, CAGGG, HSATII coordinates
ls.repeat.interest <- readRDS("~/hg19_alpha_beta_CAGGG_HSATII.RDS")

##ANALYSIS
#Only keep hits on GRCh37.p13 Primary Assembly
dt.GGGCAAAAGCCG <- dt.GGGCAAAAGCCG[grepl(
  pattern = "[0-9XY], GRCh37.p13 Primary Assembly$", x = hits_def)]
dt.GGGCAAAAGCCG <- dt.GGGCAAAAGCCG[, c("file_name","hits_def", "subject_start",
                                       "subject_end", "subject_frame"), ]
#Convert hits_def as factor
dt.GGGCAAAAGCCG[, hits_def := as.factor(hits_def)]
setattr(x = dt.GGGCAAAAGCCG$hits_def, name = "levels", value = paste0(
  "chr", gsub(
    pattern = "^Homo\\ssapiens\\schromosome\\s(\\d\\d*||X||Y),\\sGRCh37\\.p13\\sPrimary\\sAssembly$",
    replacement = "\\1", x = levels(dt.GGGCAAAAGCCG$hits_def))))
#Reorder levels
dt.GGGCAAAAGCCG[, hits_def := factor(hits_def, levels = levels(hits_def)[
  order(match(levels(hits_def), chr.chromoMap$chromosomes))])]

#Change frame integer to signs
dt.GGGCAAAAGCCG[, subject_frame := as.character(subject_frame)]
dt.GGGCAAAAGCCG[subject_frame == "-1", subject_frame := "Reverse"]
dt.GGGCAAAAGCCG[subject_frame == "1", subject_frame := "Forward"]
#Extract names from files
dt.GGGCAAAAGCCG[, file_name := gsub(pattern = "\\/1.xml", replacement = "",
                                    x = file_name)]
#Plot chromoMap
chromoMap(
  ch.files = list(chr.chromoMap), data.files = list(dt.GGGCAAAAGCCG),
  canvas_height = 700, canvas_width = 1100, chr_color = "lightgrey",
  chr_length = 10, legend = TRUE, lg_x = 50, lg_y = 300,
  data_based_color_map = T, text_font_size = 14,
  data_type = "categorical", data_colors = list(
    c("blue", "red")),
  title = "GGGCAAAAGCCG sequences alignments frame along hg19 chromosomes",
  title_font_size = 16)
#GGGCAAAAGCCG is not specific to centromeric ou pericentromeric regions

#Merge other repeats of interest
DT.repeat.interest <- rbindlist(l = ls.repeat.interest)
#Change strand data
DT.repeat.interest[strand == "-", strand := "Reverse"]
DT.repeat.interest[strand == "+", strand := "Forward"]
#Rename chromosomes
DT.repeat.interest[, genoName := paste0("chr", genoName)]
DT.repeat.interest <- DT.repeat.interest[, c(5, 1:3, 5), ]

#Plot all selected repeats chromoMap
chromoMap(
  ch.files = list(chr.chromoMap), data.files = list(DT.repeat.interest),
  canvas_height = 700, canvas_width = 1100, chr_color = "lightgrey",
  chr_length = 10, legend = TRUE, lg_x = 50, lg_y = 300,
  data_based_color_map = TRUE, text_font_size = 14,
  data_type = "categorical", data_colors = list(
    c("red", "royalblue", "orange", "darkgreen")),
  title = "Selected repeats sequences annotated along hg19 chromosomes",
  title_font_size = 16)
#CAGGG doesn't seem to be specific to pericentromeric regions
#alpha, beta and HSATII satellites seems to be more densely located around
# centromeres. HSATII seems very specific to pericentromeric regions

#Remove CAGGG repeats from the analysis
DT.repeat.interest <- DT.repeat.interest[repName != "(CAGGG)n"][, c(2:5), ]

#Plot chromoMap for Satellites of interest
chromoMap(
  ch.files = list(chr.chromoMap),
  data.files = list(DT.repeat.interest[, c(4, 1:4), ]), canvas_height = 700,
  canvas_width = 1100, chr_color = "lightgrey", chr_length = 10, legend = TRUE,
  lg_x = 50, lg_y = 300, data_based_color_map = TRUE, text_font_size = 14,
  data_type = "categorical", data_colors = list(
    c("red", "royalblue", "darkgreen")),
  title = "Selected satellite sequences annotated along hg19 chromosomes",
  title_font_size = 16)

#Save repeats of interest
fwrite(x = DT.repeat.interest, file = "~/hg19_alpha_beta_HSATII.bed",
       sep = "\t")
