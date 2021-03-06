#! /usr/bin/R

##read VCFs into Granges

libs <- c("ensemblVEP", "customProDB", "tidyverse", "biomaRt", "magrittr", "reshape2")
libsLoaded <- lapply(libs,function(l){suppressMessages(library(l, character.only = TRUE))})

strSplitFun <- function(input,sepn){
  lapply(input,function(f){strsplit(f,sepn)[[1]]})
}
strSplitVec <- function(inVec,sepn){
  sapply(seq_along(inVec),function(f){strsplit(inVec[f],sepn)[[1]]})
}

##function to annotate GRanges from biomaRt
biomartAnno <- function(grIn=NULL, GRCh=NULL, version=NULL, chr=""){

  print("Reading BiomaRt data...")
  annoMart <- useEnsembl(biomart="ensembl",
                                   GRCh=GRCh,
                                   version=version,
                                   dataset="hsapiens_gene_ensembl")

  annoGenenameEnsEnt <- as_tibble(getBM(attributes=c('chromosome_name',
                                           'start_position',
                                           'end_position',
                                           'strand',
                                           'external_gene_name',
                                           'ensembl_gene_id',
                                           'entrezgene'),
                              mart = annoMart))
  colnames(annoGenenameEnsEnt) <- c("seqnames",
                                    "start",
                                    "end",
                                    "strand",
                                    "external_gene_name",
                                    "ensembl_gene_id",
                                    "entrezgene")

  annoGenenameEnsEnt %<>% dplyr::mutate(strand = unlist(lapply(strand, function(f){
                          if(f == 1){return("+")}
                          if(f == -1){return("-")}}
        )))
  annoGenenameEnsEnt$seqnames <- paste0(chr, annoGenenameEnsEnt$seqnames)

  grangesOut <- unique(GRanges(seqnames = annoGenenameEnsEnt$seqnames,
                       ranges = IRanges(annoGenenameEnsEnt$start,
                                        annoGenenameEnsEnt$end),
                       strand = annoGenenameEnsEnt$strand,
                       external_gene_name = annoGenenameEnsEnt$external_gene_name,
                       ensembl_gene_id = annoGenenameEnsEnt$ensembl_gene_id,
                       entrezgene = annoGenenameEnsEnt$entrezgene))

  if(is.null(grIn)){
      return(grangesOut)
  }
  if(!is.null(grIn)){
      ##UNTESTED
      grOut <- GenomicRanges::intersect(grangesOut, grIn)
      return(grOut)
  }
}

vcfParseAnnoGR <- function(vcfIn){

  ##read inputs
  print("Reading VCF input...")
  grVcf <- granges(readVcf(file=vcfIn))
  gr <- suppressWarnings(InputVcf(vcfIn))
  chrTest <- levels(seqnames(gr[[1]]))[1]
  if(length(grep("chr",chrTest))>0){
    chrTest <- "chr"
  }
  if(length(grep("chr",chrTest))==0){
    chrTest <- ""
  }
  ##annotate with biomart
  biomartAll <- biomartAnno(chr=chrTest)
  fol <- as_tibble(as.data.frame(findOverlaps(grVcf, biomartAll)))

  print("Annotating...")
  grVcfAnno <- do.call(rbind, lapply(unique(fol$queryHits), function(ff){
    shits <- fol[fol$queryHits==ff,]$subjectHits
     do.call(cbind,apply(as.data.frame(values(biomartAll[shits])),2,function(f){
       data.frame(paste(gsub(" ","",f),collapse=","))
     }))
   }
  ))
  names(grVcfAnno) <- names(values(biomartAll))

  ##take all grVcf in fol
  grVcf <- grVcf[unique(fol$queryHits)]
  values(grVcf) <- cbind(values(grVcf), DataFrame(grVcfAnno))

  ##apply samples GT into mcols of grVcf

  grT <- DataFrame(do.call(cbind, lapply(seq_along(names(gr)), function(f){
    grt <- gr[[names(gr)[f]]]
    grta <- grt[unique(fol$queryHits)]$GT
    as.data.frame(grta)
  })))
  names(grT) <- names(gr)
  values(grVcf) <- cbind(values(grVcf), grT)
  return(grVcf)
}

##return vector of mutations per annotated genes
variantsPerGene <- function(varVcf, rownameTag=NULL){

  if(is.null(rownameTag)){print("Require tag for rowname");break}

  genoVcf <- varVcf[,c(grep("external_gene",colnames(values(varVcf))),
                      grep(rownameTag,colnames(values(varVcf))))]
  geneVec <- unique(genoVcf$external_gene_name)
  varPerGene <- lapply(geneVec, function(gene){
                       length(genoVcf[genoVcf$external_gene_name==gene])
  })
  names(varPerGene) <- geneVec
  return(varPerGene)
}
