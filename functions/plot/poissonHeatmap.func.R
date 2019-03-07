#! R

##takes input as table of features (row) by samples (col)
##NB if input tibble, specify col of features with featurecol argument
##optional to specify fontsize
##optional to input sleuth object, use soinput=TRUE
##this requires gene-mode scaled_reads_per_base
##allow annotation bars using 'ann' for metadata input table
##NB only samples in ann are used to make plot
##anncols tell which columns to use from ann

libs <- c("pheatmap", "PoiClaClu", "tidyverse", "NMF")
libsLoaded <- lapply(libs,function(l){
  suppressWarnings(suppressMessages(library(l, character.only = T)))
})

poissonHeatmap <- function(input, featurecol=NULL, fontsize=NULL, sounit=NULL, ann=NULL, anncols=NULL) {

  ##test if sleuth object input
  if(class(input)=="sleuth"){
    if(is.null(sounit)){
      sounit <- "scaled_reads_per_base"
    }
    inputclean <- dcast(input$obs_norm_filt,
                        target_id ~ sample,
                        value.var=sounit) %>%
                  column_to_rownames("target_id")
  }
  ##all other
  if(!class(input)=="sleuth"){
    if(is_tibble(input)){
      if(is.null(featurecol)){
        print("Please define \'featurecol\' variable to use tibble input")
        break
      }
      inputclean <- input %>%
                    column_to_rownames(featurecol)
    }
    if(!is_tibble(input)){
      inputclean <- input
    }
  }

  ##fontsize
  if(is.null(fontsize)){
    fontsized <- 6
  }
  if(!is.null(fontsize)){
    fontsized <- fontsize
  }

  if(!is.null(ann)){
    inputclean <- inputclean %>%
                  dplyr::select(unlist(ann$sample))

  }

  ##input should be OK now
  poisd <- PoissonDistance(t(inputclean))
  samplePoisDistMatrix <- as.matrix(poisd$dd)
  rownames(samplePoisDistMatrix) <- colnames(inputclean)
  colz <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
  aheatqc <- aheatmap(samplePoisDistMatrix,
                     annCol=ann[,anncols],
                     col=colz,
                     fontsize=fontsized)
  return(aheatqc)
}
