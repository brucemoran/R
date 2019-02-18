#! R

##read VCFs into GrangesList, take overlap based on 'consensus' files (list)
##use for plotting
libs <- c("ensemblVEP", "org.Hs.eg.db", "customProDB", "tidyverse")
libsLoaded <- lapply(libs,function(l){suppressMessages(library(l, character.only = TRUE))})
strSplitFun <- function(input,sepn){
  lapply(input,function(f){strsplit(f,sepn)[[1]]})
}
strSplitVec <- function(inVec,sepn){
  sapply(seq_along(inVec),function(f){strsplit(inVec[f],sepn)[[1]]})
}

vcfVepAnnParseGR <- function(vcfIn){

  ##for single sample within a single VCF

  vcf <- readVcf(vcfIn)
  if(dim(vcf)[1] != 0){
    gr <- suppressWarnings(InputVcf(vcfIn))

    ##parse info
    infor <- info(header(vcf))

    ##VEP annotation naming
    annNames <- unlist(strSplitFun(infor[rownames(infor)=="ANN",]$Description,"\\|"))

    ##somatic
    ind <- gr[names(gr)]
    seqinfo(ind) <- seqinfo(vcf)[seqlevels(ind)]

    ##annotation by CANONICAL, and add to mcols
    indAnnDf <- t(as.data.frame(lapply(strSplitFun(ind$ANN,"\\|"),function(ff){
      if(ff[annNames=="CANONICAL"]=="YES"){
        if(is.null(ff)){ff<-rep("",length(annNames))}
        if(length(ff)!=length(annNames)){
          lengExtra <- length(annNames)-length(ff)
          ff<-c(ff,rep("",lengExtra))}
        return(ff)}
        else{
          return(rep("",length(annNames)))
        }
      })))
    colnames(indAnnDf) <- annNames

    if(sum(dim(indAnnDf)) != 0){
      values(ind) <- cbind(as.data.frame(mcols(ind)),indAnnDf)
      ind$ANN <- NULL
    }
    ind <-unique(ind)

    return(ind)
  }
  else{
    print("No variants found")
    return(GRanges())
  }
}

vcfVepAnnParseGRsoma <- function(vcfIn, germline=NULL){

  ##for somatic and germline samples within a single VCF
  ##this parses only the somatic, and takes 'germline' input
  ##defines the name of germline sample, must match header

  if(is.null(germline)){
    print("Require nominated germline sample, otherwise use vcfVepAnnParseGR, quitting...")
  }

  vcf <- readVcf(vcfIn)
  if(dim(vcf)[1] != 0){
    gr <- suppressWarnings(InputVcf(vcfIn))

    ##parse info
    infor <- info(header(vcf))

    ##VEP annotation naming
    annNames <- unlist(strSplitFun(infor[rownames(infor)=="ANN",]$Description,"\\|"))

    ##somatic
    somName <- names(gr)[names(gr)!=germline]
    print(paste0("Working on: ",somName))
    som <- gr[[somName]]
    seqinfo(som) <- seqinfo(vcf)[seqlevels(som)]

    ##annotation by CANONICAL, and add to mcols
    somAnnDf <- t(as.data.frame(lapply(strSplitFun(som$ANN,"\\|"),function(ff){
      if(ff[annNames=="CANONICAL"]=="YES"){
        if(is.null(ff)){ff<-rep("",length(annNames))}
        if(length(ff)!=length(annNames)){
          lengExtra <- length(annNames)-length(ff)
          ff<-c(ff,rep("",lengExtra))}
        return(ff)}
        else{
          return(rep("",length(annNames)))
        }
      })))
    colnames(somAnnDf) <- annNames

    if(sum(dim(somAnnDf)) != 0){
      values(som) <- cbind(as.data.frame(mcols(som)),somAnnDf)
      som$ANN <- NULL
    }
    som <-unique(som)

    return(som)
  }
  else{
    print("No variants found")
    return(GRanges())
  }
}

vcfVepAnnParseGRmultiMutect2 <- function(vcfIn){

  ##for multiple samples within a single VCF
  ##returns single mcol named per sample
  ##contains common REF, ALT, FILTER and sample GT:AD:AF:DP
  vcf <- readVcf(vcfIn)

  ##parse info
  infor <- info(header(vcf))

  ##read in all of VCF, list of GRanges format for multiple samples
  grVCFList <- suppressWarnings(InputVcf(vcfIn))

  ##make a GRList that we c together
  grList <- lapply(seq_along(grVCFList), function(VCF){

    ##individual VCF GRanges object
    ind <- grVCFList[[VCF]]
    seqinfo(ind) <- seqinfo(vcf)[seqlevels(ind)]

    ##parse out required values and rename by sample
    values(ind) <- cbind(mcols(ind)[c("REF", "ALT", "FILTER")],
                         set_names(data.frame(mcols(ind)[,"GT"],
                                              mcols(ind)[,"AD"],
                                              mcols(ind)[,"AD.1"],
                                              mcols(ind)[,"AF"]),
                                              paste0(names(grVCFList)[VCF],c("GT", "AD_REF", "AD_ALT", "AF"))))
    return(ind)
  })

  ##VEP annotation naming
  annNames <- unlist(strSplitFun(infor[rownames(infor)=="ANN",]$Description,"\\|"))

  ##annotation by CANONICAL, and add to mcols
  indAnnDf <- t(as.data.frame(lapply(strSplitFun(grVCFList[[1]]$ANN,"\\|"),function(ff){
    if(ff[annNames=="CANONICAL"]=="YES"){
      if(is.null(ff)){ff<-rep("",length(annNames))}
      if(length(ff)!=length(annNames)){
        lengExtra <- length(annNames)-length(ff)
        ff<-c(ff,rep("",lengExtra))}
      return(ff)}
      else{
        return(rep("",length(annNames)))
      }
    })))
  colnames(indAnnDf) <- annNames

  ##combine all
  grCombined <- grList[[1]]
  for(xx in 2:length(grList)){
    grCombined <- join_overlap_intersect(grCombined, grList[[xx]])
  }
  values(grCombined) <- c(mcols(grCombined), indAnnDf)
  return(grCombined)
}
