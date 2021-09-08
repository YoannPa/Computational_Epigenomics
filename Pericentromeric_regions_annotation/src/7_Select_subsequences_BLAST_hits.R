#!/usr/bin/env Rscript
library("optparse")

opt_parser = OptionParser(description = "
##_Select subsequences BLAST hits_##############################################
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
Imports = c("data.table", "seqinr", "parallel")
invisible(lapply(Imports, library, character.only = T))

##PARAMETERS
#Load DT.subseq.blast
DT.subseq.blast <- fread(
  "~/PIR4(ABCD)_BLAST_results.csv")
#Load hg19 FASTA file 
hg19.fasta <- read.fasta(
  file = "~/hg19.fa.gz",
  as.string = TRUE, forceDNAtolower = FALSE, whole.header = TRUE)
#Load centromeres
centromeres <- fread("~/hg19_centromeres.bed")
#Keep only known chromosomes in hg19 fasta
hg19.fasta <- hg19.fasta[paste0("chr", c(1:22, "X", "Y", "M"))]

##ANALYSIS
#Create chromosomes length table
chr.length <- data.table("chromosomes" = names(hg19.fasta), "length" = vapply(
  X = hg19.fasta, FUN = nchar, FUN.VALUE = integer(length = 1L)))
fwrite(x = chr.length,
       file = "~/hg19_chromosomes.tsv",
       sep = "\t")
rm(hg19.fasta)

#Only keep hits on GRCh37.p13 Primary Assembly
DT.subseq.blast <- DT.subseq.blast[grepl(
  pattern = "[0-9XY], GRCh37.p13 Primary Assembly$", x = hits_def)]

#Get centers coordinate of all hg19 chromosome
centromeres[, width := end - start]
centromeres[, centers := start + (width/2)]
centers <- centromeres[, c(1,5), ]
#Save hg19 chromosome centers
fwrite(x = centers,
       file = "~/hg19_chromosomes_centers.bed",
       sep = "\t")
#Get coordinates of proximal half of both chromosome arms
chr.proximal.half <- merge(x = centers, y = chr.length, by = "chromosomes",
                           all.x = TRUE)
chr.proximal.half <- chr.proximal.half[order(match(
  chromosomes, chr.length$chromosomes[-25]))]
chr.proximal.half[, c("p.prox.half", "q.prox.half") := .(
  centers/2, centers + ((length - centers)/2))]
chr.proximal.half <- chr.proximal.half[, c(1, 4, 5), ]
#Save hg19 chromosome proximal halves
fwrite(x = chr.proximal.half, file = "~/hg19_chromosomes_proximal_half.bed",
       sep = "\t")

#Convert hits_def as factor
DT.subseq.blast[, hits_def := as.factor(hits_def)]
setattr(x = DT.subseq.blast$hits_def, name = "levels", value = paste0(
  "chr", gsub(
    pattern = "^Homo\\ssapiens\\schromosome\\s(\\d\\d*||X||Y),\\sGRCh37\\.p13\\sPrimary\\sAssembly$",
    replacement = "\\1", x = levels(DT.subseq.blast$hits_def))))
#Reorder levels
DT.subseq.blast[, hits_def := factor(hits_def, levels = levels(hits_def)[
  order(match(levels(hits_def), chr.length$chromosomes))])]
#For each subsequences keep top N hits with the smallest score possible as
# long as they remain within the proximal half of both chromosome arms
ls.selected.PIR4 <- lapply(X = levels(DT.subseq.blast$hits_def),FUN = function(chr){
  cat(paste0(chr, "\n"))
  ls.filename <- lapply(
    X = unique(DT.subseq.blast[hits_def == chr]$file_name), FUN = function(s){
      cat(paste0("\t", s, "\n"))
      ls.dt.subset <- mclapply(X = sort(unique(
        DT.subseq.blast[hits_def == chr & file_name == s]$score)), mc.cores = 7,
        FUN = function(i){
          #Get subset 
          dt.subset <- DT.subseq.blast[
            hits_def == chr & file_name == s & score >= i]
          dt.subset[, is.within.prox.half := .(
            subject_start >= chr.proximal.half[chromosomes == chr]$p.prox.half &
              subject_start <= chr.proximal.half[
                chromosomes == chr]$q.prox.half &
              subject_end <= chr.proximal.half[chromosomes == chr]$q.prox.half &
              subject_end >= chr.proximal.half[chromosomes == chr]$p.prox.half)]
          if(all(dt.subset$is.within.prox.half, na.rm = TRUE)){
            dt.subset[, cut.off.score := i]
            dt.subset <- dt.subset[order(subject_start, subject_end)]
            dt.subset
          }
        })
      ls.dt.subset <- ls.dt.subset[!sapply(ls.dt.subset, is.null)]
      if(length(ls.dt.subset) != 0){ ls.dt.subset[[1]] }
    })
  names(ls.filename) <- unique(DT.subseq.blast[hits_def == chr]$file_name)
  rbindlist(ls.filename)
})
names(ls.selected.PIR4) <- levels(DT.subseq.blast$hits_def)
#Create Bed annotation of PIR4 sub-sequences
selected.PIR4 <- rbindlist(ls.selected.PIR4)
selected.PIR4 <- selected.PIR4[, .(
  chromosome = hits_def, start = subject_start, end = subject_end,
  width = alignment_length, frame = subject_frame, name = file_name,
  alignment.score = score, chromosome.cutoff.score = cut.off.score)]
#Reverse start & end for alignment with negative frame
selected.PIR4[frame == -1, c("start", "end"):= .(end, start)]
#Change frame integer to signs
selected.PIR4[, frame := as.character(frame)]
selected.PIR4[frame == "-1", frame := "-"]
selected.PIR4[frame == "1", frame := "+"]
#Extract names from files
selected.PIR4[, name := gsub(pattern = "\\/1.xml", replacement = "", x = name)]
#Save annotation
fwrite(x = selected.PIR4, file = "~/hg19_PIR4.bed", sep = "\t")