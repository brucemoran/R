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

vcfParseGR <- function(vcfIn, germline){

  vcf <- readVcf(vcfIn)
  gr <- suppressWarnings(InputVcf(vcfIn))

  ##parse info
  infor <- info(header(vcf))

  ##somatic
  if(!is.null(germline)){
    somName <- names(gr)[names(gr)!=germline]
  }
  if(is.null(germline)){
    somName <- names(gr)
  }
  print(paste0("Working on: ",somName))
  som <- gr[[somName]]
  ##ensure an AF is there, pisces has VF instead (thanks pisces dev=D)
  if(! "AF" %in% names(mcols(som))) {
    AD <- as.numeric(unlist(mcols(som)["AD"]))
    AD1 <- as.numeric(unlist(mcols(som)["AD.1"]))
    tot <- AD+AD1
    mcols(som)$AF <- AD1/tot
  }
  seqinfo(som) <- seqinfo(vcf)[seqlevels(som)]
  return(som)
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
    ind <- gr[names(gr)][[1]]
    seqinfo(ind) <- seqinfo(vcf)[seqlevels(ind)]

    ##annotation by CANONICAL, and add to mcols
    indAnnDf <- t(as.data.frame(lapply(ind$ANN,function(ff){
      ffu <- unique(unlist(ff))
      ffuret <- unlist(lapply(strSplitFun(ffu,"\\|"), function(fff){
        if(fff[annNames=="CANONICAL"]=="YES"){
          #print(ffu)
          if(length(fff)!=length(annNames)){
            lengExtra <- length(annNames)-length(fff)
            fff <-c(fff, rep("", lengExtra))
          }
          return(fff)
        }
      }))
      if(length(ffuret)>0){
        return(ffuret[1:43])
      }
      else{
        return(rep("", length(annNames)))
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
    annNames[1] <- gsub(" ","", rev(strsplit(annNames[1],":")[[1]])[1])

    ##somatic
    somName <- names(gr)[names(gr)!=germline]
    print(paste0("Working on: ",somName))
    som <- gr[[somName]]
    seqinfo(som) <- seqinfo(vcf)[seqlevels(som)]

    ##annotation by CANONICAL, and add to mcols
    somAnnDf <- t(as.data.frame(lapply(som$ANN,function(ff){
      ffu <- unique(unlist(ff))
      ffuret <- unlist(lapply(strSplitFun(ffu,"\\|"), function(fff){
        if(fff[annNames=="CANONICAL"]=="YES"){
          #print(ffu)
          if(length(fff)!=length(annNames)){
            lengExtra <- length(annNames)-length(fff)
            fff <-c(fff, rep("", lengExtra))
          }
          return(fff)
        }
      }))
      if(length(ffuret)>0){
        return(ffuret[1:43])
      }
      else{
        return(rep("", length(annNames)))
      }
    })))
    colnames(somAnnDf) <- annNames

    if(sum(dim(somAnnDf)) != 0){
      values(som) <- cbind(as.data.frame(mcols(som)),somAnnDf)
      som$ANN <- NULL
    }
    som <- unique(som)
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
  vcf <- suppressWarnings(readVcf(vcfIn))

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
  indAnnDf <- as.data.frame(indAnnDf)
  indAnnDf$HGVSp1 <- subHGVSp(indAnnDf$HGVSp)

  ##combine all
  grCombined <- grList[[1]]
  for(xx in 2:length(grList)){
    grCombined <- join_overlap_intersect(grCombined, grList[[xx]])
  }
  values(grCombined) <- c(mcols(grCombined), indAnnDf)
  return(grCombined)
}

##from: https://raw.githubusercontent.com/cnobles/gintools/master/R/unique_granges.R, pull requested but nitl granted use this
unique_granges <- function(sites, sum.counts = FALSE, counts.col = NULL, rmDupCols = NULL){
  # Checks and balance
  if(!class(sites) == "GRanges"){
    stop("Sites object is not a GRanges class.")}
  if(sum.counts & is.null(counts.col)){
    stop("Please specify the names of the column with count information.")}
  if(!is.null(counts.col)){
    if(!counts.col %in% names(GenomicRanges::mcols(sites))){
      stop("Could not find counts column name in sites object.")}}

  # Convert sites to a data.frame and remove duplicates
  if(!length(names(sites)) == length(unique(names(sites)))){
    message("Dropping rownames for data.frame conversion.")
    df <- GenomicRanges::as.data.frame(sites, row.names = NULL)
  }else{
    df <- GenomicRanges::as.data.frame(sites)
  }
  cols <- names(df)

  if(sum.counts){
    counts_pos <- match(counts.col, cols)}

  # Sum counts if needed
  if(!sum.counts){
    df <- dplyr::distinct(df)
  }else{
    df$counts <- df[,cols[counts_pos]]
    groups <- lapply(cols[-counts_pos], as.symbol)
    df <- dplyr::group_by_(df, .dots = groups) %>%
      dplyr::summarise(counts = sum(counts)) %>%
      dplyr::ungroup()
    names(df) <- c(cols[-counts_pos], cols[counts_pos])
  }

  if(!is.null(rmDupCols)){
    dupCols <- paste0(colnames(df),".1")
    dupCols <- dupCols[dupCols %in% colnames(df)]
    df <- dplyr::distinct(df) %>%
          dplyr::select(-!!dupCols)
  }
  # Rebuild GRanges object
  gr <- GenomicRanges::GRanges(
    seqnames = df$seqnames,
    ranges = IRanges::IRanges(start = df$start, end = df$end),
    strand = df$strand,
    seqinfo = GenomicRanges::seqinfo(sites)
  )

  GenomicRanges::mcols(gr) <- dplyr::select(df, 6:length(df)) %>% unique()
  gr
}

##create single-letter HGVS protein annotation (VEP outputs 3-letter)
##take vector, gsub out aa3 for aa1
subHGVSp <- function(inVec){
  lib <- c("bio3d")
  loadedLib <- lapply(lib,function(l){suppressMessages(library(l, character.only = TRUE))})

  aa1 <- bio3d::aa.table$aa1
  ##amino acid 3 letter to gsub HGVSp
  aa3 <- unlist(lapply(bio3d::aa.table$aa3,function(f){
    sp <- strsplit(f,"")[[1]];
    paste0(sp[1], tolower(sp[2]),tolower(sp[3]))
  }))

  ##include * for Ter
  aa1 <-c(aa1,"*")
  aa3 <- c(aa3, "Ter")

  unlist(lapply(inVec,function(f){
    #check matches (should be none or two)
    a3 <- aa3[!is.na(unlist(stringi::stri_match_all(f,regex=aa3)))]
    a1 <- aa1[!is.na(unlist(stringi::stri_match_all(f,regex=aa3)))]
    ##beauty:
    #https://stackoverflow.com/questions/19424709/r-gsub-pattern-vector-and-replacement-vector
    if(length(a3)>0){
      names(a1) <- a3
      str_replace_all(f,a1)
    }
    else{
      return("")
    }
  }))
}
