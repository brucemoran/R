#! R

##takes input as table of features (row) by samples (col)
##NB if input tibble, specify col of features with featurecol argument
##sleuth object input requires gene-mode scaled_reads_per_base
##allow annotation bars using 'ann' for metadata input table
##NB only samples in ann are used to make plot
##anncols tell which columns to use from ann
##vst transform data, and return is vstreturn is not null

libs <- c("tidyverse")
libsLoaded <- lapply(libs,function(l){
  suppressWarnings(suppressMessages(library(l, character.only = T)))
})

pcaPlot <- function(input, intgroup, featurecol=NULL, unitsinput=NULL, ann=NULL, annsamplecol=NULL, anncols=NULL, plottitle=NULL, vsttrans=TRUE, vstblind=TRUE, vstreturn=NULL) {

  ##test if sleuth object input
  if(class(input)=="sleuth"){
    if(is.null(unitsinput)){
      unitsinput <- "scaled_reads_per_base"
    }
    inputclean <- dcast(input$obs_norm_filt,
                        target_id ~ sample,
                        value.var=unitsinput) %>%
                  column_to_rownames("target_id")
  }
  ##all other
  if(!class(input)=="sleuth"){
    if(is.null(unitsinput)){
      unitsinput <- "unknown"
    }
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

  ##annotation data
  if(is.null(annsamplecol)){
    annsamplecol <- "sample"
  }
  if(!is.null(ann)){
    inputclean <- inputclean %>%
                  dplyr::select(ann[[annsamplecol]])
  }

  #plottitle
  if(is.null(plottitle)){
    plottitle <- ""
  }

  ##input should be VST transformed now
  inputvst <- varianceStabilizingTransformation(as.matrix(round(inputclean, 0)), blind=vstblind)

  #calculate PCA, proportion of variance attributable
  pca_res <- prcomp(t(inputvst), scale. = F, center = F)
  sum_pca_res <- summary(pca_res)
  percVar <- round(sum_pca_res$importance[2,1:2]*100,1)

  ##join with metadata
  pca_res_plot <- left_join(as_tibble((pca_res$x)[,1:2], rownames="sample"),
                          ann, by=annsamplecol) %>%
                mutate_if(is.character, as.factor) %>%
                dplyr::select(everything(), intgroup_ = one_of(intgroup))

  #plot object
  shape_vals <- seq(1,1+length(unique(pca_res_plot$intgroup_2)))
  pca_ploto <- ggplot(pca_res_plot,
                    aes(x = PC1, y = PC2, colour = intgroup_1, shape = intgroup_2, label = sample)) +
             geom_point(size=2) +
             geom_text(hjust=1, vjust=-0.4) +
             scale_shape_manual(values=shape_vals) +
             xlab(paste0("PC1: ", percVar[1], "% variance")) +
             ylab(paste0("PC2: ", percVar[2], "% variance")) +
             labs(title=plottitle,
                  subtitle=paste0("Units: ", unitsinput),
                  colour = intgroup[1],
                  shape = intgroup[2])
  if(!is.null(vstreturn)){
    return(list(pca_ploto, inputvst))
  }
  if(is.null(vstreturn)){
    return(pca_ploto)
  }
}
