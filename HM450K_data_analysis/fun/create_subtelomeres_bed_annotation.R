
#Create subtelomere bed annotation
ls.subtel <- lapply(X = chrs, function(i){
  #Get subtelomere from p arm
  p.dt <- subtelomere_blasts[
    chromosome == i & arm == "p" & subject_frame == 1,
    c("chromosome", "subject_start", "subject_end", "extremity")]
  #For distal part of the sequence
  if(nrow(p.dt[extremity == "distal"]) == 1){
    #Check if the distal position is before all proximal positions on the p arm
    if(p.dt[, .SD[which.min(subject_end)]]$extremity == "distal"){
      p.start <- p.dt[extremity == "distal"]$subject_start - 1
    } else { #Unexploitable alignment, get p arm telomere end position
      p.start <- telomeres[chromosomes == i][1]$end
    }
  } else if(nrow(p.dt[extremity == "distal"]) > 1){
    p.start <- min(p.dt[extremity == "distal"]$subject_start) - 1
  } else {
    if(nrow(telomeres[chromosomes == i]) == 0){
      #Chr17 in hg19 does not have any telomere defined
      p.start <- 0
    } else { p.start <- telomeres[chromosomes == i][1]$end }  
  }
  #For proximal part of the sequence
  if(nrow(p.dt[extremity == "proximal"]) == 1){
    p.end <- p.dt[extremity == "proximal"]$subject_end
  } else if(nrow(p.dt[extremity == "proximal"]) > 1){
    p.dt[extremity == "proximal", dist := abs(subject_end - p.start - 500000)]
    p.end <- p.dt[extremity == "proximal", .SD[which.min(dist)]]$subject_end
    # p.end <- max(p.dt[extremity == "proximal"]$subject_end)
  } else { p.end <- p.start + 500000 }
  #Get subtelomere from q arm
  q.dt <- subtelomere_blasts[
    chromosome == i & arm == "q" & subject_frame == -1,
    c("chromosome", "subject_start", "subject_end", "extremity")]
  #For distal part of the sequence
  if(nrow(q.dt[extremity == "distal"]) == 1){
    q.end <- q.dt[extremity == "distal"]$subject_start
  } else if(nrow(q.dt[extremity == "distal"]) > 1){
    q.end <- max(q.dt[extremity == "distal"]$subject_start)
  } else { 
    if(nrow(telomeres[chromosomes == i]) == 0){
      #Chr17 in hg19 does not have any telomere defined, set chr17 end position
      q.end <- 81195210
    } else { q.end <- telomeres[chromosomes == i][2]$start }
  }
  #For proximal part of the sequence
  if(nrow(q.dt[extremity == "proximal"]) == 1){
    q.start <- q.dt[extremity == "proximal"]$subject_end - 1
  } else if(nrow(q.dt[extremity == "proximal"]) > 1){
    q.dt[extremity == "proximal",
         dist := abs(q.end - (subject_end - 1) - 500000)]
    q.start <- q.dt[extremity == "proximal",
                    .SD[which.min(dist)]]$subject_end - 1
  } else { q.start <- q.end - 500000 }
  #Create chromosome subtelomeres annotation
  sub.dt <- data.table("chromosome" = i, "start" = c(p.start, q.start),
                       "end" = c(p.end, q.end))
  return(sub.dt)
})
#Get widths of subtelomeric regions
ls.subtel <- lapply(ls.subtel, function(i){i[, width := end - start]})
hg19.subtelomeric.bed <- rbindlist(ls.subtel)
fwrite(
  x = hg19.subtelomeric.bed,
  file = "PCAWG_Paper/data/hg19_annotations/hg19_subtelomeres.bed", sep = "\t",
  scipen = 2)

