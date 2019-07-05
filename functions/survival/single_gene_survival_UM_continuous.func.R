##survival based on continuous expression of genes in all combinations
library(survival)

##genes must be a rowname of expruse
##rows of survuse must match cols of expruse
single_gene_survival_UM_continuous <- function(survuse, expruse, genes){

  ##object to be returned
  coxPHList <- fitList <- diffList <- medLinPredList <- as.list(1:length(genes))

  ##years, event
  event <- as.numeric(as.vector(survuse[,1]))
  years <- as.numeric(as.vector(survuse[,2]))

  for (c in 1:length(genes)){

    ##make input objects for survival
    currExprGenes <- t(expruse[rownames(expruse) %in% strsplit(genes[c],"\\.")[[1]],])

    ##median for one gene only used for univariate (tautology)
    if(c <= length(rownames(expruse))){
      medCurrExpr <- currExprGenes > median(currExprGenes)
      cox11 <- coxph(Surv(years, event) ~ medCurrExpr)
    }

    ##first cox model (multivariate to define ConcInd)
    cox1 <- coxph(Surv(years, event) ~ ., data=data.frame(currExprGenes))

    ##define linear predictor to allow single HR for multivariate analyses
    ##use as response term in coxPH, survfit to get dichotomised survival
    medLinPred <- as.numeric(cox1$linear.predictors > median(cox1$linear.predictors))
    medLinPredList[[c]] <- medLinPred

    cox2 <- coxph(Surv(years, event) ~ medLinPred)
    survfit2 <- survfit(Surv(years, event) ~ medLinPred)

    ##smuggle over what the time is based on (DMFS, RFS etc)
    ##this is just for naming plots(!)
    survfit2$type <- colnames(survuse)[2]

    ##make CI equal to cox1 (Steve Barrons request)
    ##this is accessed for 'cloud' plots vs HR
    cox1$concordance[4] <- as.numeric(summary(cox1)$concordance[1])
    cox2$concordance[4] <- as.numeric(summary(cox1)$concordance[1])
    names(cox1$concordance)[4] <- "as.numeric(summary(cox1))"
    names(cox2$concordance)[4] <- "as.numeric(summary(cox1))"

    if(c <= length(rownames(expruse))){
      ##this is like an easter egg!
      coxPHList[[c]] <- list(cox2,cox11)
      fitList[[c]] <- survfit2
    }
    if(c > length(rownames(expruse))){
      coxPHList[[c]] <- cox2
      fitList[[c]] <- survfit2
    }
  }
  names(coxPHList) <- names(fitList) <- genes

  survivalLists <- list(coxPHList, fitList, medLinPredList)
  return(survivalLists)
}
