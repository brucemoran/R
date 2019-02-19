#! R

##function to take DGEList, conditions data and a specific condition on which to colour MDS plot; option to printo ut =D
library("limma")

MDSplotter <- function(dge, conds, cond, printo=NULL){
    conds[] <- lapply(conds, factor)
    col.cond <- unlist(conds[cond])
    levels(col.cond) <- sample(rainbow(nlevels(col.cond)))
    col.cond <- as.character(col.cond)
    lcpm <- cpm(dge, log = TRUE, prior.count = 0.5, group = cond)
    naming <- unlist(conds[cond])
    condNames <- paste0(rownames(conds), "_", naming)
    plotMDS(lcpm, labels = condNames, col = col.cond)
    title(main = paste0("A. ", cond, " groups"))
    if(!is.null(printo)){
        pdf(printo)
        plotMDS(lcpm, labels=condNames, col=c(col.cond))
        title(main=paste0("A. ", cond, " groups"))
        dev.off()
    }
}
