#!/usr/bin/env Rscript
library("optparse")

opt_parser = OptionParser(description = "
##_Extract the 4 subsequences and BLAST them to hg19_###########################
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
Imports = c("data.table", "NCBI.BLAST2DT")
invisible(lapply(Imports, library, character.only = T))

##PARAMETERS
#load PIR4 query alignments
PIR4.subseq <- fread(
  file = "~/PIR4_query_alignments.csv")
#Get splitted queries from AC073318 (positions 71,401 to 120,576)
x <- data.frame("GB.access" = "AC073318", "Start" = 71401, "End" = 120576)
ls.seq <- split_queries(x = x, by = 7025, ncores = 7)

##ANALYSIS
#Extract the 4 sub-sequences
PIR4.subseq[, seq_name := gsub(pattern = "\\/1\\.xml$", replacement = "",
                               x = file_name)]
PIR4.subseq[, query_alignments := unlist(
  lapply(X = seq(nrow(PIR4.subseq)), function(r){
    substr(x = ls.seq[[PIR4.subseq[r,]$seq_name]],
           start = PIR4.subseq[r,]$min.query_start,
           stop = PIR4.subseq[r,]$max.query_end)}))]
#Merge continous sequences into the 4 sub-sequences
PIR4A <- paste(PIR4.subseq[file_name %in% c(
  "AC073318:71401-78425/1.xml", "AC073318:78426-85450/1.xml")]$query_alignments,
  collapse = "")
PIR4B <- paste(PIR4.subseq[file_name %in% c(
  "AC073318:85451-92475/1.xml", "AC073318:92476-99500/1.xml")]$query_alignments,
  collapse = rep(x = "n", 14)) #14bp gap filled with "n"
PIR4C <- PIR4.subseq[
  file_name == "AC073318:99501-106525/1.xml"]$query_alignments
PIR4D <- paste(PIR4.subseq[file_name %in% c(
  "AC073318:106526-113550/1.xml", "AC073318:113551-120576/1.xml")
]$query_alignments, collapse = "")
#Save sub-sequences in a FASTA file
ls.sub.seq <- list("PIR4A" = PIR4A, "PIR4B" = PIR4B, "PIR4C" = PIR4C,
                   "PIR4D" = PIR4D)
write.fasta(
  sequences = ls.sub.seq, names = names(ls.sub.seq), nbchar = 60,
  as.string = TRUE,
  file.out =
    "~/PIR4_subsequences_A_B_C_D.fasta")
#Submit PIR4 sub-sequence to NCBI BLAST servers
DT.subseq.blast <- get.NCBI.BLAST2DT(
  sequences = ls.sub.seq, db = "genomic/9606/GCF_000001405.25",
  res.dir = "~/result_BLAST", ncores = 7, auto.rm.dir = FALSE,
  email = "yoann.pageaud@gmail.com")
#Save DT.subseq.blast
fwrite(x = DT.subseq.blast,
       file = "~/PIR4(ABCD)_BLAST_results.csv")