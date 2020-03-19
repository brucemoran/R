exomeTumourMutationBurden <- function(GRplot){

  ##get exome for Illumina Nextera Rapid
  exomeBed <- fread("https://support.illumina.com/content/dam/illumina-support/documents/documentation/chemistry_documentation/samplepreps_nextera/nexterarapidcapture/nexterarapidcapture_exome_targetedregions_v1.2.bed", showProgress=FALSE, data.table=FALSE)

  ##triage chr tag, add header
  exomeBed[,1] <- gsub("chr","",exomeBed[,1])
  colnames(exomeBed) <- c("seqname", "start", "end")
  exomeGR <- makeGRangesFromDataFrame(exomeBed, ignore.strand=TRUE)
  exomeSize <- sum(width(exomeGR))/1000000

  ##overlap with input
  hits <- as.data.frame(findOverlaps(GRplot, exomeGR))
  ##output
  GRplotExome <- GRplot[hits$queryHits]
  TMB <- round(length(GRplotExome)/exomeSize,digits=1)
  return(TMB)
}
