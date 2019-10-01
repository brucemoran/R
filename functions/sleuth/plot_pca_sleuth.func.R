#https://groups.google.com/forum/#!topic/kallisto-sleuth-users/XrGXXWTUxsU
plot_pca_so <- function(so, group){
  mat <- sleuth:::spread_abundance_by(
         abund = so$obs_norm_filt,
         var = "scaled_reads_per_base",
         which_order = so$sample_to_covariates$sample)

  # Calculate PCA, log2 transfer and normalization between samples
  pc_x=1L
  pc_y=2L
  mat.pca <- prcomp(log2(t(mat)+1),center = TRUE,scale. = TRUE)

  #computation of variances
  eigenvalues <- (mat.pca$sdev) ^ 2
  var_explained <- eigenvalues * 100 / sum(eigenvalues)

  #set label names
  x_lab <- paste0('PC1 (', round(var_explained[1],digits=1))
  x_lab <- paste0(x_lab, '%)')
  y_lab <- paste0('PC2 (', round(var_explained[2],digits=1))
  y_lab <- paste0(y_lab, '%)')

  #Extract PC1 and PC2 to pcs, you can also change this PC1 and PC2 to others
  pcs <- sleuth:::as_df(mat.pca$x[, c(pc_x, pc_y)])
  pcs$sample <- rownames(pcs)
  rownames(pcs) <- NULL

  #add 'Group' information from experimental design
  pcs <- dplyr::left_join(pcs, so$sample_to_covariates,by = 'sample')

  #ggplot
  pc_x <- paste0('PC', pc_x)
  pc_y <- paste0('PC', pc_y)
  ggp <- ggplot(pcs, aes_string(pc_x, pc_y, colour = group))+
          geom_text(size = 4, alpha = 0.8, label=so$sample_to_covariates$sample) +
          ggtitle("PCA based on log2 normalized and filtered genes")+
          xlab(x_lab)+
          ylab(y_lab)+
          theme(text = element_text(size=10),axis.text = element_text(size=8),plot.title = element_text(size = 14))
  return(ggp)
}
