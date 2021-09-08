#!/usr/bin/env Rscript
library("optparse")

opt_parser = OptionParser(description = "
##_Extract PIR4 following BLAST hits_###########################################
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
setwd("/media/yoann/Disque Dur 4/PhD/Brors_Lab/")
Imports = c("data.table", "parallel")
invisible(lapply(Imports, library, character.only = T))

##PARAMETERS
#Look for hits happening in a specific order
right.order <- c("AC073318:71401-78425/1.xml", "AC073318:78426-85450/1.xml",
                 "AC073318:85451-92475/1.xml", "AC073318:92476-99500/1.xml",
                 "AC073318:99501-106525/1.xml", "AC073318:106526-113550/1.xml",
                 "AC073318:113551-120576/1.xml")
#Load selected.PIR4
selected.PIR4 <- fread(
  "~/selected_PIR4_hits.csv")

##ANALYSIS
#Extract sequences following each others, from one end or the other
ls.followseq <- mclapply(
  X = unique(selected.PIR4$hits_def),  mc.cores = 7, FUN = function(chr){
    dt.chr <- selected.PIR4[hits_def == chr]
    ls.chr <- lapply(X = seq(nrow(dt.chr)), FUN = function(i){
      #First element
      if(dt.chr[i,]$file_name == right.order[1]){
        #Positive frame
        if(i != nrow(dt.chr) & dt.chr[i,]$subject_frame == 1){
          if(i != nrow(dt.chr)-1 & dt.chr[i+1,]$file_name == right.order[2] &
             dt.chr[i+1,]$subject_frame == 1){
            if(i != nrow(dt.chr)-2 & dt.chr[i+2,]$file_name == right.order[3] &
               dt.chr[i+2,]$subject_frame == 1){
              if(i != nrow(dt.chr)-3 & dt.chr[i+3,]$file_name == right.order[4] &
                 dt.chr[i+3,]$subject_frame == 1){
                if(i != nrow(dt.chr)-4 &
                   dt.chr[i+4,]$file_name == right.order[5] &
                   dt.chr[i+4,]$subject_frame == 1){
                  if(i != nrow(dt.chr)-5 &
                     dt.chr[i+5,]$file_name == right.order[6] &
                     dt.chr[i+5,]$subject_frame == 1){
                    if(dt.chr[i+6,]$file_name == right.order[7] &
                       dt.chr[i+6,]$subject_frame == 1){
                      dt.chr[i:(i+6)]
                    } else { dt.chr[i:(i+5)] }
                  } else { dt.chr[i:(i+4)] }
                } else { dt.chr[i:(i+3)] }
              } else { dt.chr[i:(i+2)] }
            } else { dt.chr[i:(i+1)] }
          }
        } else if(i != 1 & dt.chr[i,]$subject_frame == -1){
          if(i != 2 & dt.chr[i-1,]$file_name == right.order[2] &
             dt.chr[i-1,]$subject_frame == -1){
            if(i != 3 & dt.chr[i-2,]$file_name == right.order[3] &
               dt.chr[i-2,]$subject_frame == -1){
              if(i != 4 & dt.chr[i-3,]$file_name == right.order[4] &
                 dt.chr[i-3,]$subject_frame == -1){
                if(i != 5 & dt.chr[i-4,]$file_name == right.order[5] &
                   dt.chr[i-4,]$subject_frame == -1){
                  if(i != 6 & dt.chr[i-5,]$file_name == right.order[6] &
                     dt.chr[i-5,]$subject_frame == -1){
                    if(dt.chr[i-6,]$file_name == right.order[7] &
                       dt.chr[i-6,]$subject_frame == -1){
                      dt.chr[(i-6):i]
                    } else { dt.chr[(i-5):i] }
                  } else { dt.chr[(i-4):i] }
                } else { dt.chr[(i-3):i] }
              } else { dt.chr[(i-2):i] }
            } else { dt.chr[(i-1):i] }
          }
        }
      } else if(dt.chr[i,]$file_name == right.order[7]){ #Last element
        if(i != 1 & dt.chr[i,]$subject_frame == 1){
          if(i != 2 & dt.chr[i-1,]$file_name == right.order[6] &
             dt.chr[i-1,]$subject_frame == 1){
            if(i != 3 & dt.chr[i-2,]$file_name == right.order[5] &
               dt.chr[i-2,]$subject_frame == 1){
              if(i != 4 & dt.chr[i-3,]$file_name == right.order[4] &
                 dt.chr[i-3,]$subject_frame == 1){
                if(i != 5 & dt.chr[i-4,]$file_name == right.order[3] &
                   dt.chr[i-4,]$subject_frame == 1){
                  if(i != 6 & dt.chr[i-5,]$file_name == right.order[2] &
                     dt.chr[i-5,]$subject_frame == 1){
                    if(dt.chr[i-6,]$file_name == right.order[1] &
                       dt.chr[i-6,]$subject_frame == 1){
                      dt.chr[(i-6):i]
                    } else { dt.chr[(i-5):i] }
                  } else { dt.chr[(i-4):i] }
                } else { dt.chr[(i-3):i] }
              } else { dt.chr[(i-2):i] }
            } else { dt.chr[(i-1):i] }
          }
        } else if(i != nrow(dt.chr) & dt.chr[i,]$subject_frame == -1){
          if(i != nrow(dt.chr)-1 & dt.chr[i+1,]$file_name == right.order[6] &
             dt.chr[i+1,]$subject_frame == -1){
            if(i != nrow(dt.chr)-2 & dt.chr[i+2,]$file_name == right.order[5] &
               dt.chr[i+2,]$subject_frame == -1){
              if(i != nrow(dt.chr)-3 & dt.chr[i+3,]$file_name == right.order[4] &
                 dt.chr[i+3,]$subject_frame == -1){
                if(i != nrow(dt.chr)-4 &
                   dt.chr[i+4,]$file_name == right.order[3] &
                   dt.chr[i+4,]$subject_frame == -1){
                  if(i != nrow(dt.chr)-5 &
                     dt.chr[i+5,]$file_name == right.order[2] &
                     dt.chr[i+5,]$subject_frame == -1){
                    if(dt.chr[i+6,]$file_name == right.order[1] &
                       dt.chr[i+6,]$subject_frame == -1){
                      dt.chr[i:(i+6)]
                    } else { dt.chr[i:(i+5)] }
                  } else { dt.chr[i:(i+4)] }
                } else { dt.chr[i:(i+3)] }
              } else { dt.chr[i:(i+2)] }
            } else { dt.chr[i:(i+1)] }
          }
        }
      } else if(dt.chr[i,]$file_name == right.order[4]){ #Middle element
        #Positive frame
        if(i != nrow(dt.chr) & dt.chr[i,]$subject_frame == 1){
          if(dt.chr[i+1,]$file_name == right.order[5] &
             dt.chr[i+1,]$subject_frame == 1){
            if(i != 1){
              if(dt.chr[i-1,]$file_name == right.order[3] &
                 dt.chr[i-1,]$subject_frame == 1){
                dt.chr[(i-1):(i+1)]
              } else { dt.chr[i:(i+1)] }
            } else { dt.chr[i:(i+1)] }
          }
        } else if(i != nrow(dt.chr) & dt.chr[i,]$subject_frame == -1){
          if(dt.chr[i+1,]$file_name == right.order[3] &
             dt.chr[i+1,]$subject_frame == -1){
            if(i != 1){
              if(dt.chr[i-1,]$file_name == right.order[5] &
                 dt.chr[i-1,]$subject_frame == -1){
                dt.chr[(i-1):(i+1)]
              } else { dt.chr[i:(i+1)] }
            } else { dt.chr[i:(i+1)] }
          }
        } 
        if(i != 1 & dt.chr[i,]$subject_frame == 1){
          if(dt.chr[i-1,]$file_name == right.order[3] &
             dt.chr[i-1,]$subject_frame == 1){
            if(i != nrow(dt.chr)){
              if(dt.chr[i+1,]$file_name == right.order[5] &
                 dt.chr[i+1,]$subject_frame == 1){
                dt.chr[(i-1):(i+1)]
              } else { dt.chr[(i-1):i] }
            } else { dt.chr[(i-1):i] }
          }
        } else if(i != 1 & dt.chr[i,]$subject_frame == -1){
          if(dt.chr[i-1,]$file_name == right.order[5] &
             dt.chr[i-1,]$subject_frame == -1){
            if(i != nrow(dt.chr)){
              if(dt.chr[i+1,]$file_name == right.order[3] &
                 dt.chr[i+1,]$subject_frame == -1){
                dt.chr[(i-1):(i+1)]
              } else { dt.chr[(i-1):i] }
            } else { dt.chr[(i-1):i] }
          }
        }
      }
    })
    ls.chr <- ls.chr[!sapply(ls.chr, is.null)]
  })
#Get largest coordinates of the aligned part of sequences
selected.PIR4 <- rbindlist(lapply(X = ls.followseq, FUN = rbindlist))
selected.PIR4[, c("min.query_start", "max.query_end") := .(
  min(query_start, na.rm = TRUE), max(query_end, na.rm = TRUE)), by = file_name]
PIR4.subseq <- unique(selected.PIR4[
  order(match(file_name, right.order)),
  c("file_name", "min.query_start", "max.query_end"), ], by = "file_name")
#There are 4 sub-sequences within PIR4!
#Save query coordinates of PIR4
fwrite(x = PIR4.subseq, file = "~/PIR4_query_alignments.csv")
