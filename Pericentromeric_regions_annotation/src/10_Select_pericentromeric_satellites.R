#!/usr/bin/env Rscript
library("optparse")

opt_parser = OptionParser(description = "
##_Select pericentromeric satellites_###########################################
Version = '0.0.1'
Date = '2021-08-17'
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
Imports = c("data.table", "NCBI.BLAST2DT", "parallel")
invisible(lapply(Imports, library, character.only = T))

##PARAMETERS
#BLAST "GGGCAAAAGCCG" on hg19
# (get.NCBI.BLAST2DTdoesn't seem to work on too short sequences)
dt.GGGCAAAAGCCG <- aggregate_NCBI_BLAST_XMLs2DT(
  dir.to.xmls = "~/GGGCAAAAGCCG_BLAST_results/", ncores = 1)
#Load RepeatMasker annotation
hg19_repeats <- fread("~/rmsk_2020-03-22.txt.gz", nThread = 7)
#Load selected.PIR4
selected.PIR4 <- fread(file = "~/hg19_PIR4.bed")
#Load annotation of alpha, beta, and HSATII satellites
DT.alpha_beta_HSATII <- fread("~/hg19_alpha_beta_HSATII.bed")
#Load chr.length
chr.length <- fread(file = "~/hg19_chromosomes.tsv")
#Load chr.centers
chr.centers <- fread("~/hg19_chromosomes_centers.bed")
#Load chr.proximal.half
chr.proximal.half <- fread("~/hg19_chromosomes_proximal_half.bed")

##ANALYSIS
colnames(hg19_repeats) <- c(
  "bin","swScore","milliDiv","milliDel","milliIns","genoName","genoStart",
  "genoEnd","genoLeft","strand","repName","repClass","repFamily","repStart",
  "repEnd","repLeft","id")
#Keep genoName, genoStart, geneEnd, strand, repName, repClass and repFamily.
hg19_repeats <- hg19_repeats[, c(6:8,10:13), with = FALSE]
#Replace SINE tRNA Family by tRNA_SINE
hg19_repeats[repClass == "SINE" & repFamily == "tRNA", repFamily := "SINE_tRNA"]

#Filter out unmapped scaffolds and uncertain repeats families
hg19_repeats <- hg19_repeats[like(genoName, "^chr[0-9XY][0-9]*$")]
hg19_repeats <- hg19_repeats[like(repFamily, "^.*[^?]$")]
#5212846 repeats after cleaning

#Replace chromosome names
hg19_repeats <- hg19_repeats[, genoName := as.factor(genoName)]
hg19_repeats[, genoName := factor(genoName, levels = levels(genoName)[
  order(match(levels(genoName), paste0("chr", c(1:22, "X", "Y"))))])]
setattr(x = hg19_repeats$genoName, name = "levels",
        value = c(seq(1:22), "X", "Y"))
#Order repeats by chromosome, start and end
hg19_repeats <- hg19_repeats[order(genoName, genoStart, genoEnd)]

#Save BLAST results
fwrite(x = dt.GGGCAAAAGCCG, file = "~/GGGCAAAAGCCG_BLAST_results.csv")

#Subset sequences of interest from hg19_repeats
#Satellite centr ALR/Alpha
dt.alpha <- hg19_repeats[repClass == "Satellite" & repFamily == "centr" &
                           repName == "ALR/Alpha"]
# #Satellite Satellite BSR/Beta
# dt.beta <- hg19_repeats[repClass == "Satellite" & repFamily == "Satellite" &
#                           repName == "BSR/Beta"]
#Simple_repeat Simple_repeat (CAGGG)n
dt.CAGGG <- hg19_repeats[repClass == "Simple_repeat" &
                           repFamily == "Simple_repeat" & repName == "(CAGGG)n"]
#Satellite Satellite HSATII
dt.HSATII <- hg19_repeats[repClass == "Satellite" & repFamily == "Satellite" &
                            repName == "HSATII"]

# ls.repeat.interest <- list("Alpha" = dt.alpha,"Beta" = dt.beta,
#                            "CAGGG" = dt.CAGGG, "HSATII" = dt.HSATII)
ls.repeat.interest <- list(
  "Alpha" = dt.alpha, "CAGGG" = dt.CAGGG, "HSATII" = dt.HSATII)
saveRDS(object = ls.repeat.interest,
        file = 
          "~/hg19_alpha_CAGGG_HSATII.RDS")

#Merge annotations with PIR sub-sequences 
DT.pericentr.repeat <- rbind(
  selected.PIR4[, c(1:3, 6), ], DT.alpha_beta_HSATII, use.names = FALSE)

#Define chromosomes and names as factors
DT.pericentr.repeat[, c("chromosome", "name"):= .(
  as.factor(chromosome), as.factor(name))]
#Add chromosomes centers
DT.pericentr.repeat <- merge(x = DT.pericentr.repeat, y = chr.centers,
                             by.x = "chromosome", by.y = "chromosomes",
                             all.x = TRUE)
setnames(x = DT.pericentr.repeat, old = "centers", new = "chrom.centers")
#Compute centers of repeats
DT.pericentr.repeat[, rep.centers := (start + end)/2]
#Add chromsomes length
DT.pericentr.repeat <- merge(
  x = DT.pericentr.repeat, y = chr.length, by.x = "chromosome",
  by.y = "chromosomes", all.x = TRUE)
DT.pericentr.repeat[, chrom.start := 1]
setnames(x = DT.pericentr.repeat, old = "length", new = "chrom.end")
#Compute distances relative to chrom.centers
DT.pericentr.repeat[, c("distance", "chrom.start", "chrom.end") := .(
  rep.centers - chrom.centers, chrom.start - chrom.centers,
  chrom.end - chrom.centers)]

#Calculate distance from centers of chromosomes proximal half intervals
chr.proximal.half[, c("p.prox.half", "q.prox.half") := .(
  p.prox.half - centers, q.prox.half - centers)]

#Remove the PIR sequences
DT.pericentr.repeat <- DT.pericentr.repeat[!name %like% "PIR4"]
#Remove beta satellites
DT.pericentr.repeat <- DT.pericentr.repeat[name != "BSR/Beta"]

#For intervals decreasing from 10Kb to 20Mb keep intervals that contain strictly more than
# 4 repeats and check if all intervals kept fall within the proximal half of each
# chromosomes, if not decrease the interval size.
intervals <- rev(
  c(10, 50, 100, 200, 500, 1000, 1500, 3000, 10000, 20000)*1000)
DT.chr.uniq <- unique(DT.pericentr.repeat, by = "chromosome")

ls.all.intervals <- lapply(X = intervals, FUN = function(i){
  cat(paste("interval size:", i, "\n"))
  ls.interval.res <- mclapply(
    X = DT.chr.uniq$chromosome, mc.cores = 4, FUN = function(chr){
      #Define number of bins n
      n = round(chr.length[chromosomes == chr]$length/i)
      #Get n+1 breakpoints
      breakpts <- seq(DT.chr.uniq[chromosome == chr]$chrom.start,
                      DT.chr.uniq[chromosome == chr]$chrom.end,
                      length.out = n + 1)
      #Create intervals
      breakpts <- breakpts[2:(length(breakpts)-1)]
      #Create intervals based on breaks
      dt <- data.table::data.table(
        c(DT.chr.uniq[chromosome == chr]$chrom.start:
            DT.chr.uniq[chromosome == chr]$chrom.end), findInterval(
              c(DT.chr.uniq[chromosome == chr]$chrom.start:
                  DT.chr.uniq[chromosome == chr]$chrom.end), breakpts))
      #Create data.table of coordinates based on the intervals created
      dt.intervals <- do.call(rbind, by.default(dt, dt$V2, function(g){
        data.table::data.table(Start = min(g$V1), End = max(g$V1))}))
      rm(dt)
      #Count number of repeats in intervals
      ls.count <- lapply(
        X = seq(nrow(dt.intervals)), FUN = function(r){
          nrow(DT.pericentr.repeat[
            chromosome == chr & distance >= dt.intervals[r, ]$Start &
              distance <= dt.intervals[r, ]$End])
        })
      #Add count to dt intervals
      dt.intervals[, repeat.count := ls.count]
      #Keep intervals where count is strictly above 3 (density of more than 3 repeat by i)
      dt.intervals <- dt.intervals[repeat.count > 3]
      #Check if all remaining intervals are within the proximal half of chromosome
      bool.p <- all(
        dt.intervals$Start >= chr.proximal.half[chromosomes == chr]$p.prox.half)
      bool.q <- all(
        dt.intervals$End <= chr.proximal.half[chromosomes == chr]$q.prox.half)
      if(bool.p & bool.q){
        dt.intervals[, c("interval.size", "chromosome", "in.prox.half") := .(
          i, chr, TRUE)]
        dt.res <- dt.intervals
      } else {
        dt.res <- data.table(
          "Start" = 0, "End" = 0, repeat.count = 0, "interval.size" = i,
          "chromosome" = chr, "in.prox.half" = FALSE)
      }
      dt.res
    })
  rbindlist(ls.interval.res)
})

#Density of satellites is not the same from one chromosome to another
# -> differences in terms of elongation activity ?

#Keep the density cut-off returning the most intervals for each chromosome
DT.all.intervals <- rbindlist(ls.all.intervals)
DT.all.intervals <- DT.all.intervals[in.prox.half == TRUE]
DT.all.intervals[, N.intervals := .N, by = c("interval.size", "chromosome")]
# DT.all.intervals[
#   N.intervals == max(N.intervals), by = c("interval.size", "chromosome")]
DT.all.intervals <- DT.all.intervals[DT.all.intervals[, .I[
  N.intervals == max(N.intervals)], by = chromosome]$V1]
#If multiple interval sizes return maximum interval number, select the largest
# interval size.
DT.all.intervals <- DT.all.intervals[DT.all.intervals[, .I[
  interval.size == max(interval.size)],by = chromosome]$V1]

# Retrieve satellites within these intervals
ls.selected.satellites <- lapply(
  X = seq(nrow(DT.all.intervals)), FUN = function(i){
    DT.pericentr.repeat[chromosome == DT.all.intervals[i, ]$chromosome &
                          DT.all.intervals[i, ]$Start < distance &
                          distance < DT.all.intervals[i, ]$End]
  })
DT.selected.satellites <- rbindlist(ls.selected.satellites)
#Beta satellites seems actually more specific to pericentromeric regions than 
# alpha satellites, but at a larger distance from centromeres.
#HSATII even more pericentromeric than both beta and alpha.

#Save selected satellites
saveRDS(object = DT.selected.satellites,
        file = "~/DT.selected.satellites_without_beta.RDS")
