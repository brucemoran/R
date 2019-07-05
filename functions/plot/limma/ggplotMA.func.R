##plotMA, but for ggplot

#plotWithHighlights(x, y, status = NULL, values = NULL, hl.pch = 16, hl.col = NULL, hl.cex = 1, legend = "topleft", bg.pch = 16, bg.col = "black", bg.cex = 0.3, pch = NULL, col = NULL, cex = NULL, ...)
library("tidyverse")
library("ggplot2")
library("ggrepel")

ggplotMA <- function(fit, coef = 1, anno=NULL, annocol=NULL, pval=0.1, title=NULL, topn=5, ...){

  ##standard execution of plotMA
  if (!is(fit, "MArrayLM")) stop("fit must be an MArrayLM")

  xlab <- "Average Mean Expression"
  ylab <- "Log Fold Change"

  ##table from which to plot; include P and adjusted
  plot.tbl <- tibble(ensembl_gene_id = rownames(fit$coef),
                     Amean = fit$Amean,
                     Coefs = fit$coef[,coef],
                     Lods = fit$lods[,coef],
                     Pval = fit$p.value[,coef]) %>%
              dplyr::mutate(mlog10Pval = -log10(Pval)) %>%
              dplyr::mutate(Adj.Pval = p.adjust(Pval, method="BH")) %>%
              na.omit()

  ##annotate with anno table if supplied
  if(! is.null(anno)){
    if("ensembl_gene_id" %in% colnames(anno)){
      plot.tbl <- left_join(plot.tbl, anno) %>% na.omit()
    }
  }

  ##which labels to use for annotation
  if(is.null(annocol)){
    annocol <- "ensembl_gene_id"
  }
  plot.tbl <- cbind(plot.tbl, annocol=plot.tbl[,grep(annocol, colnames(plot.tbl), value=F)])
  colnames(plot.tbl)[dim(plot.tbl)[2]] <- "annocol"
  plot.tbl <- plot.tbl %>% dplyr::mutate(Pval.Rank = rank(Pval))
  mapl <- ggplot(plot.tbl, aes(x=Amean, y=Coefs, label=external_gene_name)) +
          geom_point(size=0.5) +
          geom_point(data=subset(plot.tbl, Adj.Pval < pval), colour="red", size=0.5) +
          geom_text_repel(data=subset(plot.tbl, Adj.Pval < pval & Pval.Rank < 11), colour="red", segment.alpha=0.5, segment.size=0.2, segment.colour="red", fontface = "bold", cex=3) +
          labs(x=xlab, y=ylab, title=paste0("MA plot", subtitle=title))

  mapln <- ggplot(plot.tbl, aes(x=Amean, y=Coefs, label=external_gene_name)) +
          geom_point(size=0.5) +
          geom_point(data=subset(plot.tbl, Adj.Pval < pval), colour="red", size=0.5) +
          labs(x=xlab, y=ylab, title=paste0("MA plot", subtitle=title))

  return(list(plot.tbl, mapl, mapln))
}
