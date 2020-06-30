#! R

##function to bin GRanges object based on a specified region size (bp)
##this is primarily focussed on CNA data, to allow clustering of multiple samples
library(tidyverse)
library(GenomicRanges)
options(scipen=999)

##function
bin_granges <- function(gr, bin_size=NULL, name_in=NULL, mcol_in=NULL){

  ##check gr is A GRanges object
  if(!as.vector(class(gr)) %in% c("GRanges", "GRangesList", "list")){
    print("Input \'gr\' is not a GRanges, GRangesList nor list object, retry")
    break
  }

  ##if grList, set gr to a single entity therein for binning
  if(as.vector(class(gr)) %in% c("GRangesList", "list")){
    gri <- gr[[1]]
    if(!as.vector(class(gri)) == "GRanges"){
      print("Input \'gr\' list does not contain a GRanges object, retry")
    } else {
      gri <- gr
    }
  }

  ##set default
  if(is.null(bin_size)){
    print("Bin size is not specified, setting to 10MB")
    bin_size <- 10000000
  }

  ##make bin gr
  grb <- bin_maker(gri, bin_size)

  ##iterate over gr input; NB using seq_along so listify gr which is GRanges
  if(as.vector(class(gr)) == "GRanges"){
    gr <- list(gr)
  }

  ##set names as list_elem_x for gr if no names
  if(is.null(name_in) & is.null(names(gr))){
    name_in <- paste0("list_elem", seq_along(gr))
  } else {
    name_in <- names(gr)
  }

  ##which mcol value(s) to use from gr input (name or index)
  if(is.null(mcol_in)){
    mcol_in <- "Total_Copy_Number"
  }

  lapply(seq_along(gr), function(f){

    ##find hits and populate a named mcols(grb)
    grf <- gr[[f]]

    ##screen for diploid regions
    grfn <- grf[unlist(mcols(grf[,mcol_in]))!=2]

    ##findoverlaps
    hits <- as.data.frame(findOverlaps(grfn, grb, ignore.strand=TRUE))

    ##N.B. that hits can include multiple 'subjectHits', i.e. bins can be hit by more that one CNA
    ##here we test if: those are the same value (and include just that value)
    ##different value with one == 2, include other
    ##fail, reduce bin size below
    grbb$x <- mcols(grfn[hits$queryHits, mcol_in])
    gr$CGC_SYMBOLs <- "-"
    ##loop to collapse symbols per region
    for(x in 1:max(hits$queryHits)){
      hitsx <- as.vector(sort(unique(hits$SYMBOL[hits$queryHits==x])))
      hitsx <- hitsx[!is.na(hitsx)]
      if(length(hitsx)==0){gr$CGC_SYMBOL[x] <- NA; gr$CGC_SYMBOLs[x] <- 0}
      else{
        gr$CGC_SYMBOL[x] <- base::paste(hitsx[2:length(hitsx)], collapse=";")
        gr$CGC_SYMBOLs[x] <- length(hitsx)-1;
      }
    }
    })
}

##get input of a GRanges object, and use that seqinfo to bin by bin_size
bin_maker <- function(gri, bin_size){

  ##new GRanges across that bin based on that seqinfo
  sn <- seqnames(seqinfo(gri))
  sl <- seqlengths(seqinfo(gri))

  ##create IRanges for seqnames and seqlengths
  ##combine into bin GRange
  grbList <- lapply(seq_along(sn), function(f){
      st <- seq(1, sl[f], by=bin_size)
      nd <- c(seq(st[2], sl[f], by=bin_size), as.vector(sl[f])+1)-1
      GRanges(seqnames = sn[f], IRanges(start = st, end = nd), strand = "*")
  })
  grb <- unlist(as(grbList, "GRangesList"))
  mcols(grb) <- bin_size
  names(mcols(grb)) <- "bin_size"
  return(grb)
}

##multiple hits in the same bin are an issue, this reduces based on GRanges in a list
##iterate over the grList, find where multiple hits occur, set max_bin to the size of interval -1
##round to nearest kb

test_bin_reduce <- function(grList){

  if(!as.vector(class(grList))=="list"){
    print("Input 'grList' is not a list")
    break
  }

  ##make starter bin at 10MB
  gri <- grList[[1]]
  grb <- bin_maker(gri, bin_size=10000000)

  ##iterate over list, find multihits
  max_bins <- unlist(lapply(seq_along(grList), function(f){

    ##find hits and populate a named mcols(grb)
    grf <- grList[[f]]
    hits <- as.data.frame(findOverlaps(grf, grb, ignore.strand=TRUE))
    grfMulti <- as_tibble(hits) %>%
                dplyr::filter(subjectHits %in% as.numeric(names(tabshits[tabshits>1]))) %>%
                dplyr::select(queryHits) %>% unlist() %>% unique()

    max_bin <- 10000000
    if(length(grb) < length(hits$subjectHits)){
      ##multihits, find max_bin
      max_bin <- as_tibble(ranges(grf[grfMulti])) %>%
                 dplyr::select(width) %>%
                 unlist() %>% min() %>% signif(digits=2)
    }
    return(max_bin)
  }))
  return(max_bins)
}
