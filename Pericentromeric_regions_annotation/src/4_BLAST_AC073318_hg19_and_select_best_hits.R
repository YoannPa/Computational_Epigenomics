#!/usr/bin/env Rscript
library("optparse")

opt_parser = OptionParser(description = "
##_BLAST AC073318 to hg19 and select best hits_#################################
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
Imports = c("data.table", "NCBI.BLAST2DT", "BiocompR")
invisible(lapply(Imports, library, character.only = T))

##ANALYSIS
#Get sequence from AC073318 (positions 71,401 to 120,576)
x <- data.frame("GB.access" = "AC073318", "Start" = 71401, "End" = 120576)
#The BLAST request for the sequence overload the NCBI servers: Split the query
# in subqueries before NCBI BLAST submission
ls.seq <- split_queries(x = x, by = 7025, ncores = 7)
#Submit sub-queries to NCBI BLAST
DT.PIR4 <- get.NCBI.BLAST2DT(
  sequences = ls.seq, db = "genomic/9606/GCF_000001405.25",
  res.dir = "~/result_BLAST", ncores = 7, auto.rm.dir = FALSE,
  email = "yoann.pageaud@gmail.com")
#Save DT.PIR4
fwrite(x = DT.PIR4, file = "~/PIR4_BLAST_results_7kb_fullseq.csv")

#Only keep hits on GRCh37.p13 Primary Assembly
DT.PIR4 <- DT.PIR4[grepl(pattern = "[0-9XY], GRCh37.p13 Primary Assembly$", x = hits_def)]

#Order by chromosome, start and end
DT.PIR4 <- DT.PIR4[order(hits_def, subject_start, subject_end)]
#TODO: add option to pass hlines & vlines to the function instead of guessing it
hist.e_val <- fancy.hist(
  x = DT.PIR4$`E-value`, nbreaks = 100, xmax = 0.5, ncores = 7,
  bin.col = "darkgreen") +
  xlab("Distribution of PIR4 BLAST hits E-values") +
  geom_vline(xintercept = 2, linetype = "dashed", color = "red", size = 1)
hist.score <- fancy.hist(x = DT.PIR4$score, nbreaks = 200, xmax = 1000,
                         ncores = 7, bin.col = "darkorange") +
  xlab("Distribution of PIR4 BLAST hits scores") + 
  geom_vline(xintercept = 120, linetype = "dashed", color = "red", size = 1)
hist.alignment <- fancy.hist(x = DT.PIR4$alignment_length, nbreaks = 200,
                             xmax = 1000, ncores = 7) +
  xlab("Distribution of PIR4 BLAST hits alignments length") + 
  geom_vline(xintercept = 140, linetype = "dashed", color = "red", size = 1)
ggsave(hist.e_val, filename = "~/Hist_PIR4_BLAST_E-values.pdf")
ggsave(hist.score, filename = "~/Hist_PIR4_BLAST_scores.pdf")
ggsave(hist.alignment, filename = "~/Hist_PIR4_BLAST_alignments.pdf")

#Keep hits with alignment length above 700 with E-value below 0.01,
# with score above 600
selected.PIR4 <- DT.PIR4[
  alignment_length > 700 & `E-value` < 0.01 & score > 600]
#Save selected.PIR4 hits
fwrite(x = selected.PIR4, file = "~/selected_PIR4_hits.csv")
