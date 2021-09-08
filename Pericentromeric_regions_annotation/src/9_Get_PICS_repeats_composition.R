#!/usr/bin/env Rscript
library("optparse")

opt_parser = OptionParser(description = "
##_Get Pericentromeric Interspersed Composite Sequence (PICS) repeats composition_##
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
Imports = c("GenomicRanges", "data.table")
invisible(lapply(Imports, library, character.only = T))

##PARAMETERS
#Load RepeatMasker annotation
hg19_repeats <- fread("~/rmsk_2020-03-22.txt.gz", nThread = 7)
#Load selected.PIR4
selected.PIR4 <- fread(file = "~/hg19_PIR4.bed")

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

#Convert repeats and PIR4 to GRanges
GR_repeats <- makeGRangesFromDataFrame(
  df = hg19_repeats, keep.extra.columns = TRUE, ignore.strand = FALSE,
  seqnames.field = "genoName", start.field = "genoStart", end.field = "genoEnd",
  strand.field = "strand", starts.in.df.are.0based = TRUE)
GR_repeats <- sort(GR_repeats)
GR_PIR4 <- makeGRangesFromDataFrame(
  df = selected.PIR4, keep.extra.columns = TRUE, ignore.strand = TRUE,
  seqnames.field = "chromosome", start.field = "start", end.field = "end")
GR_PIR4 <- sort(GR_PIR4)
seqlevels(GR_PIR4) <- gsub(
  pattern = "^chr", replacement = "", x = seqlevels(GR_PIR4))
#Get Repeats overlapping PIR4 sequences
DT.PIR4_repeats <- mergeByOverlaps(
  query = GR_PIR4, subject = GR_repeats, ignore.strand = TRUE)
DT.PIR4_repeats <- as.data.table(DT.PIR4_repeats)[, -c(
  "GR_PIR4.strand", "name", "alignment.score", "chromosome.cutoff.score",
  "frame", "repName", "repClass", "repFamily"), ]

#Retrieve the chain of repeats (contiguous repeats with same names are kept as 1)
DT.PIR4_repeats_pos <- DT.PIR4_repeats[GR_PIR4.frame == "+"][order(
  GR_PIR4.seqnames, GR_PIR4.start, GR_PIR4.end, GR_repeats.seqnames,
  GR_repeats.start, GR_repeats.end)][, repeat.chain := paste(GR_repeats.repName[
    c(TRUE, !GR_repeats.repName[
      -length(GR_repeats.repName)] == GR_repeats.repName[-1])],
    collapse = ", "), by = c("GR_PIR4.seqnames", "GR_PIR4.start", "GR_PIR4.end")
  ]
DT.PIR4_repeats_neg <- DT.PIR4_repeats[GR_PIR4.frame == "-"][order(
  GR_PIR4.seqnames, GR_PIR4.start, GR_PIR4.end, GR_repeats.seqnames,
  GR_repeats.start, GR_repeats.end)][, repeat.chain := paste(
    rev(GR_repeats.repName)[c(TRUE, !rev(GR_repeats.repName)[
      -length(rev(GR_repeats.repName))] == rev(GR_repeats.repName)[-1])],
    collapse = ", "), by = c("GR_PIR4.seqnames", "GR_PIR4.start", "GR_PIR4.end")
  ]
DT.PIR4_repeats <- rbind(DT.PIR4_repeats_pos, DT.PIR4_repeats_neg)[order(
  GR_PIR4.seqnames, GR_PIR4.start, GR_PIR4.end)]
#Extract chains
PIR4A_repeat_chain <- unique(DT.PIR4_repeats, by = c(
  "GR_PIR4.seqnames","GR_PIR4.start","GR_PIR4.end"))[GR_PIR4.name == "PIR4A"][
    order(-GR_PIR4.alignment.score, GR_repeats.start, GR_repeats.end)
  ]$repeat.chain
PIR4B_repeat_chain <- unique(DT.PIR4_repeats, by = c(
  "GR_PIR4.seqnames","GR_PIR4.start","GR_PIR4.end"))[GR_PIR4.name == "PIR4B"][
    order(-GR_PIR4.alignment.score, GR_repeats.start, GR_repeats.end)
  ]$repeat.chain
PIR4C_repeat_chain <- unique(DT.PIR4_repeats, by = c(
  "GR_PIR4.seqnames","GR_PIR4.start","GR_PIR4.end"))[GR_PIR4.name == "PIR4C"][
    order(-GR_PIR4.alignment.score, GR_repeats.start, GR_repeats.end)
  ]$repeat.chain
PIR4D_repeat_chain <- unique(DT.PIR4_repeats, by = c(
  "GR_PIR4.seqnames","GR_PIR4.start","GR_PIR4.end"))[GR_PIR4.name == "PIR4D"][
    order(-GR_PIR4.alignment.score, GR_repeats.start, GR_repeats.end)
  ]$repeat.chain

max_length <- max(c(length(PIR4A_repeat_chain), length(PIR4B_repeat_chain),
                    length(PIR4C_repeat_chain), length(PIR4D_repeat_chain))) 
#Save chains
dt.PIR4_repeat_chains <- data.table(
  PIR4A_repeat_chain = c(PIR4A_repeat_chain,
                         rep(NA, max_length - length(PIR4A_repeat_chain))),
  PIR4B_repeat_chain = c(PIR4B_repeat_chain,
                         rep(NA, max_length - length(PIR4B_repeat_chain))),
  PIR4C_repeat_chain = c(PIR4C_repeat_chain,
                         rep(NA, max_length - length(PIR4C_repeat_chain))),
  PIR4D_repeat_chain = c(PIR4D_repeat_chain,
                         rep(NA, max_length - length(PIR4D_repeat_chain))))
fwrite(x = dt.PIR4_repeat_chains, file = "~/PIR4_repeats_chains.tsv",
       sep = "\t")