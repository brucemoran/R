library(tidyverse)
library(GenomicRanges)
options(scipen=999)

##function
reduce_intersect_granges <- function(grList, mcol_in="Total_Copy_Number"){

  ##check gr is A GRanges object
  if(!as.vector(class(grList)) %in% c("GRangesList", "list")){
    print("Input \'grList\' is not a GRangesList nor list object, retry")
    break
  }

  ##make reduced from grList
  grRed <- Reduce(GenomicRanges::intersect, grList)

  ##
  grf_mcols <- lapply(seq_along(grList), function(f){
    grf <- grList[[f]]
    hits <- as.data.frame(findOverlaps(grf, grRed, ignore.strand=TRUE))
    mcols(grf[,mcol_in][hits$queryHits])
  })

  mcols(grRed) <- grf_mcols
  names(mcols(grRed)) <- names(grList)
  return(grRed)
}
