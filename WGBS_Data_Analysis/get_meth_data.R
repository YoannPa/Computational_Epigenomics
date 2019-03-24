
##IMPORTS
library(data.table)
library(GenomicRanges)
library(R.utils)
source("handle_directories.R")

##FUNCTIONS

# get.MC.data ##################################################################

#' Retrieve methylation calls from a given directory, automatically remove SNPs
#' from the data, format and save them in bed files and compute various
#' statistics on it.
#' (Warning: function can crash if the computer get disconnected from the
#' server. In such case:
#' 1- open the directory and copy the "temp.csv" file out,
#' 2- remove the latest sample folder generated,
#' 3- re-run the function.
#' 4- when the process will be completed you will have to merge manually
#' "temp.csv" and "MC_STATS.csv").
#' 
#' @param PID.dir       a \code{character} specifying the directory where all
#'                      the samples data are stored.
#' @param sample.names  a \code{list} specifying the names of samples to be
#'                      retrieved. If NULL, all samples in the directory will be
#'                      retrieved (Default: samples.names = NULL). 
#' @param out.dir       a \code{character} specifying the directory where the
#'                      data should be saved.
#' @param MC.name       a \code{character} to name the output folder
#'                      (Default: MC.name = "Meth_Calls").
#' @param GR.SNPs       a \code{character} specifying the path to the SNPs RDS
#'                      file as a genomic ranges (\code{GRanges}) object. If
#'                      there are no SNPs to remove, GR.SNPs must be NULL
#'                      (Default: GR.SNPs = NULL) and rm.CpGs.SNPs must be
#'                      FALSE.
#' @param rm.CpGs.SNPs  a \code{logical} to specify whether CpGs that are SNPs
#'                      should be removed or not (Default: rm.CpGs.SNPs = TRUE). 
#' @param chromosomes   a\code{character} vector specifying the chromosome to be
#'                      processed (Default: chromosomes = c(1:19,"MT","X","Y")).
#' @param write.files   a \code{logical} to specify whether the data retrieved
#'                      should be written out in the destination folder or not
#'                      (Default: writes.files = TRUE).
#' @param gz            a \code{logical} to specify whether the saved BED files
#'                      should be gzipped or not (Default: gz=TRUE).
#' @param ncores        a \code{integer} specifying the number of cores/threads
#'                      to be used for parallel processing
#'                      (Default: ncores = 1).
#' @param skip.amplicon a \code{logical} specifying whether to skip AmpliconSeq
#'                      samples or not (Default: skip.amplicon).
#' @param sample.subdir A \code{character} matching the name of the subdirectory
#'                      below the sample directory
#' @param paired        A \code{logical} specifying whether the WGBS data are
#'                      containing paired-end reads or not.
#' @param file.type     a \code{character} specifying the type of files to be
#'                      targeted. Currently only "merged" files are supported.
#' @return logs about the processing of methylation data, creates a folder
#'         containing the methylation calls, various statistics about the
#'         methylation data and the SNPs.
#' @author Yoann Pageaud.

get.MC.data<-function(PID.dir, sample.names=NULL, out.dir, MC.name="Meth_Calls",
                      GR.SNPs=NULL, rm.CpGs.SNPs=TRUE,
                      chromosomes=c(1:19,"MT","X","Y"), write.files=TRUE,
                      gz=TRUE, ncores=1, skip.amplicon = TRUE, sample.subdir,
                      paired=T, file.type="merged"){
  
  cat("Listing Samples in the directory...")
  all.samples<-list.files(PID.dir)
  cat("Done.\n")
  if(is.null(sample.names)){
    Sample.IDs <- all.samples
  } else {
    Sample.IDs <- all.samples[all.samples %in% sample.names]
  }
  #Create output dir 
  Output.Dir<- file.path(out.dir,MC.name)
  #List of supported methylation folders
  Supported.Mfolders<-c("methylationCalling","methylationCallingMetrics")
  #Set Working Directory
  setwd(PID.dir)
  
  if(rm.CpGs.SNPs){ #Load SNPs
    cat("Loading SNPs...")
    Gr_SNPs<-readRDS(GR.SNPs) 
    # Reduce Ranges and discard Metacolumns
    Gr_SNPs<-reduce(Gr_SNPs)
    cat("Done.\n")
  }
  #Create MC Dir if it does not already exists
  if(dir.exists(Output.Dir) == F){ dir.create(Output.Dir) }
  #Remove Amplicon-Seq 
  if(skip.amplicon){ Sample.IDs<-Sample.IDs[!(grepl(".+-1$",Sample.IDs))] }
  
  if(paired){
    reads.type<-"paired"
  } else {stop("reads type not supported yet.")}
  if(file.type == "merged"){data.type <-"merged-alignment"} else {
    stop("sample type not supported yet.")}
  
  #Create temp file
  file.create(file.path(Output.Dir,"temp.csv"))
  
  list.smpl.stat<-lapply(Sample.IDs, function(smpl_id){  
    #Create output path
    smpl.out<-file.path(Output.Dir,smpl_id)
    if(dir.exists(smpl.out) == F) {
      if(write.files){
        #Create output directories
        dir.create(smpl.out, showWarnings = FALSE)
      }
      #Define sample paths
      smplpath<-file.path(PID.dir,smpl_id,file.path(sample.subdir,reads.type),
                          data.type)
      #Make data dir
      data.dir<-c(file.path("methylation",file.type),
                  file.path("qualitycontrol",file.type))
      #Define methylation path
      methylation.path<-file.path(smplpath,data.dir[1])
      #Define QC path: QC.path<-file.path(smplpath,data.dir[2])
      
      ##Get Meta table
      if(dir.exists(smplpath)) {
        setwd(smplpath)
        #Check if directory contain any file or folder
        if(length(list.files()) != 0) {
          cat(paste0(smpl_id,"\n")) #Print Sample ID
          #Load Metadata Table
          smpl.meta.table<-read.csv(file = "metadataTable.tsv",header = T,
                                    sep = "\t",stringsAsFactors = F)
          
          if(length(unique(smpl.meta.table[,1]))==1){ #Check sample type
            Sample.type = unique(smpl.meta.table[,1])
          } else { stop("Multiple sample types found for this sample.") }
          
          ##Get Methylation Calling Data
          setwd(methylation.path) #Set Methylation Path
          #Get folders available
          if(dir.exists(Supported.Mfolders[1])) {
            cat("\tRetrieving Methylation Calling Data...")
            #Move in directory
            setwd(file.path(methylation.path,Supported.Mfolders[1]))
            #Set file prefix
            prefix<-paste(Sample.type,smpl_id,file.type, sep = "_")
            #Get list of files and remove index files
            list.MC.files<-list.files()[!(grepl(pattern = "bed.gz.tbi",
                                                x = list.files()))]
            #Split all file names following character "." as separator
            list.split<-strsplit(list.MC.files,split = "\\.")
            #Convert as matrix rows
            list.split<-lapply(list.split,matrix,nrow = 1)
            #Rbind the list
            df.MCfile<-as.data.frame(do.call(rbind, list.split),
                                     stringsAsFactors = F)[,-1]
            #Subset CG and CH Files and keep chromosomes of interest
            df.CG<-df.MCfile[which(df.MCfile$V5 == "CG"),]
            df.CG<-df.CG[df.CG$V4 %in% chromosomes,]
            # df.CH<-df.MCfile[which(df.MCfile$V5 == "CH"),]
            
            #Get MC and compute Stat for each chromosome
            list.chr.stat<-
              mclapply(seq(nrow(df.CG)), mc.cores = ncores,function(i){
                #Load file without chromosome and CG
                tbl.CG<-fread(paste(prefix,paste0(df.CG[i,],collapse = "."),
                                    sep = "."))[,c(2,3,5:7)]
                
                chr<-df.CG$V4[i]
                colnames(tbl.CG)<-c("pos","strand","SNPs","M","U")
                
                #Subset by forward and reverse strands
                df.forward<-tbl.CG[strand == "+",c(1,3:5)]
                colnames(df.forward)[1]<-"start"
                df.reverse<-tbl.CG[strand == "-",c(1,3:5)]
                colnames(df.reverse)[1]<-"end"
                
                #Create Granges of CpGs on the chromosome
                # (start and stop are 0-based)
                Gr_chr_CpG<-GRanges(seqnames = chr,
                                    ranges = IRanges(start = df.forward$start,
                                                     end = df.reverse$end+1))
                #Combine data
                DT.comb<-
                  data.table(pos = df.forward$start,
                             SNPs = rowMeans(cbind(df.forward$SNPs,
                                                   df.reverse$SNPs)),
                             M = rowSums(cbind(df.forward$M, df.reverse$M)),
                             U = rowSums(cbind(df.forward$U,df.reverse$U)))
                #Get Statistics
                DT.chr.stat<-
                  DT.comb[, .(Chr.AmountM = sum(M), Chr.AmountU = sum(U),
                              Chr.CpG.Total = DT.comb[,length(pos)],
                              Chr.CpG.Covered = DT.comb[M != 0 |
                                                          U !=0,length(pos)],
                              Chr.CpG.NotCovered =
                                DT.comb[M == 0 & U ==0,length(pos)],
                              Chr.CpG.SNP.Pred1 =
                                DT.comb[SNPs == 1, length(pos)],
                              Chr.CpG.SNP.Pred0.99to0.5 =
                                DT.comb[SNPs < 1 & SNPs >= 0.5, length(pos)],
                              Chr.CpG.SNP.Pred0.49to0.01 =
                                DT.comb[SNPs < 0.5 & SNPs > 0, length(pos)],
                              Chr.CpG.SNP.Pred0 = DT.comb[SNPs == 0,
                                                          length(pos)])]
                if(rm.CpGs.SNPs){
                  #Keep CpGs not overlapping SNPs
                  Gr_chr_CpG_NoSNPs<-subsetByOverlaps(Gr_chr_CpG,Gr_SNPs,
                                                      invert = T)
                  #Keep CpGs not being SNPs in 129P2_OlaHsd Mouse Strain
                  DT.comb<-DT.comb[pos %in% start(Gr_chr_CpG_NoSNPs)]
                }
                #Get number of CpGs predicted SNPs not removed
                Chr.CpG.SNP.Pred.NoRm<-
                  DT.comb[, .(Chr.CpG.SNP.Pred1.NoRm =
                                DT.comb[SNPs == 1, length(pos)],
                              Chr.CpG.SNP.Pred0.99to0.5.NoRm =
                                DT.comb[SNPs < 1 & SNPs >= 0.5, length(pos)],
                              Chr.CpG.SNP.Pred0.49to0.01.NoRm =
                                DT.comb[SNPs < 0.5 & SNPs > 0, length(pos)],
                              Chr.CpG.SNP.Pred0.NoRm =
                                DT.comb[SNPs == 0, length(pos)])]
                
                #Get percentage of CpG predicted SNPs removed
                Chr.CpG.SNP.Pred.Rm<-
                  round((DT.chr.stat[,c(6:9)] - Chr.CpG.SNP.Pred.NoRm)*100/
                          DT.chr.stat[,c(6:9)],2)
                colnames(Chr.CpG.SNP.Pred.Rm)<-
                  c("Pred1.%Rm","Pred0.99to0.5.%Rm","Pred0.49to0.01.%Rm",
                    "Pred0.%Rm")
                if(write.files){
                  #Print File
                  fwrite(x = DT.comb, file = file.path(smpl.out,
                                                       paste0(chr,".bed")),
                         sep = " ", col.names = FALSE)
                  if(gz){
                    #Gzip all files in parallel
                    gzip(file.path(smpl.out, paste0(chr,".bed")),
                         destname=file.path(smpl.out, paste0(chr,".bed.gz")))
                  }
                }
                #Save DT of Chromosome Stats in list
                cbind(DT.chr.stat,Chr.CpG.SNP.Pred.NoRm,Chr.CpG.SNP.Pred.Rm)
              })  #Close tasks for all chromosomes
            names(list.chr.stat)<-df.CG$V4
            #Rbind stat rows into a DT
            DT.chr.stat<-do.call(rbind, list.chr.stat)
            rownames(DT.chr.stat)<-names(list.chr.stat)
            #Apply sum() on col 1 to 13
            DT.chr.stat[, (1:13) := lapply(.SD, sum), .SDcols = 1:13]
            #Apply mean() on col 14 to 17
            DT.chr.stat[, (14:17) := lapply(.SD, mean, na.rm=T),.SDcols = 14:17]
            #Apply round() on col 14 to 17            
            DT.chr.stat[, (14:17) := lapply(.SD, round, digits = 2),
                        .SDcols = 14:17]
            #Keep first row only and re-order columns
            DT.chr.stat<-DT.chr.stat[1,c(3:5,1,2,6:17)]
            #Rename Columns
            colnames(DT.chr.stat)<-gsub(pattern = "Chr\\.",replacement = "",
                                        x = colnames(DT.chr.stat))
            #Save temporary in case of crash
            DF<-as.data.frame(DT.chr.stat)
            rownames(DF)<-smpl_id
            write.table(DF, file.path(Output.Dir,"temp.csv"), sep = ",",
                        row.names = T, col.names = F, append = T)
            #Save Sample Stat
            cat("Done.\n")
            DT.chr.stat
          } else {
            cat("No Methylation Calling Files.\n")
            if(dir.exists(Supported.Mfolders[2])) {
              cat("Retrieving Methylation Calling Metrics.\n")
              setwd(file.path(methylation.path,Supported.Mfolders[2]))
            } else { cat("No Methylation Calling Metrics.\n") }
          } #Close dir.exists() condition on the type of directory for a sample
        } else {
          cat("No File Found.\n")
        } #Close list.file() condition on the input directory
      } #Close dir.exists() condition on the input directory
    } else {
      warning("Sample data already found in the output directory")
      NULL
    }
  }) #Close tasks for all samples
  names(list.smpl.stat)<-Sample.IDs #Rename
  list.smpl.stat[sapply(list.smpl.stat, is.null)] <- NULL #Remove NULL elements
  #Rbind stat rows into a dataframe
  df.smpl.stat<-as.data.frame(do.call(rbind, list.smpl.stat))
  # rownames(df.smpl.stat)<-names(list.smpl.stat) #Rename
  #Set Working Dir back to data/ folder
  setwd(Output.Dir)
  #Load temporary csv file
  full.stat<-read.csv("temp.csv",header = F,row.names = "V1")
  colnames(full.stat)<-colnames(df.smpl.stat)
  #Print Methylation Data Table to CSV
  write.csv(x = full.stat, file = "MC_STATS.csv")
  cat(paste0('Statistics file "MC_STATS.csv" saved at:\n',getwd(),"\n"))
}
