
##volcano plot in ggplot
library("tidyverse")
library("ggplot2")
library("ggrepel")

volcanoggplot <- function(fit, coef = 1, style = "p-value", highlight = 0, names = NULL, h1.col="blue", xlab = "Log2 Fold Change", ylab = NULL, pch=16, cex=0.35, anno=NULL, annocol=NULL, pval=0.1, title=NULL, topn=5, ...){

  ##standard execution of plotMA
  if (!is(fit, "MArrayLM")) stop("fit must be an MArrayLM")
  style <- match.arg(tolower(style), c("p-value", "b-statistic"))
  if (style == "p-value") {
    if (is.null(fit$p.value)) stop("No p-values found in linear model fit object")
    if (is.null(ylab)) ylab = "-log10(P-value)"
  }
  if (is.null(fit$lods)) stop("No B-statistics found in linear model fit object")
  if (is.null(ylab)) ylab = "Log Odds of Differential Expression"

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

  if(style == "p-value"){
    plot.tbl <- plot.tbl %>% dplyr::rename(y = mlog10Pval)
  }
  if(style != "p-value"){
    plot.tbl <- plot.tbl %>% dplyr::rename(y = Lods)
  }

  ##rank Pvalues for topn
  plot.tbl <- plot.tbl %>% dplyr::mutate(Pval.Rank = rank(Pval))

  ##plot object
  vcpl <- ggplot(plot.tbl, aes(x=Coefs, y=y, label=annocol)) +
          geom_point(size=0.5) +
          geom_point(data=subset(plot.tbl, Adj.Pval < pval), colour="red", size=0.5) +
          geom_text_repel(data=subset(plot.tbl, Adj.Pval < pval & Pval.Rank < 11), segment.alpha=0.5, segment.size=0.2, segment.colour="red", colour="red", cex=2, force=2) +
          labs(x=xlab, y=ylab, title=("Volcano Plot"), subtitle=title)
  vcpln <- ggplot(plot.tbl, aes(x=Coefs, y=y, label=annocol)) +
          geom_point(size=0.5) +
          geom_point(data=subset(plot.tbl, Adj.Pval < pval), colour="red", size=0.5) +
          labs(x=xlab, y=ylab, title=("Volcano Plot"), subtitle=title)
  return(list(plot.tbl, vcpl, vcpln))
}
