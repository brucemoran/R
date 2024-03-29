#! R

##script to run on LFQ data provided herein
##NB that data was parsed first using parse_BD2Han_LF!_141220.sh

##based on DEP: https://www.bioconductor.org/packages/release/bioc/vignettes/DEP/inst/doc/DEP.html

library("DEP")
library("tidyverse")
library("limma")
library("SummarizedExperiment")

#' run DEP on groups
#' @param project_name string to label outputs
#' @param sample_map_tb tibble with sample mapping info
#' @param prot_data_tb protein expression tibble, format is col_1 = protein IDs, col_2 ... col_n are samples annotated in sample_map_tb
#' @param sample_ID string of colname in sample_map_tb holding sample IDs in prot_data_tb tibble
#' @param group vector indicating which columns to use for DE analysis; multiple entries are combined in a new column which is used as a single new variable
#' @param row_data string indicating column of unique protein IDs
#' @param tech_reps column indicating technical replicate groupings
#' @param exclude_group string indicating group to exclude from analysis
#' @param convert_group vector length 2 indicating group to convert from, to
#' @param p_adj_val limma adjusted p-value threshold
#' @param plot_dep flag to allow not plotting of DEP QCs
#' @param out_dir directory to output to (created if not extant)

#' @return list of ggplot2 objects for printing (PCA, PCs, loadings)
#' @export

##format label, design replicate
run_expt_group <- function(project_name = "DEP_DE", sample_map_tb, sample_ID, group, prot_data_tb, row_data, tech_reps = NULL, exclude_group = NULL, convert_group = NULL, p_adj_val = 0.1, plot_dep = FALSE, out_dir = "DEP"){

  ##two or more groups to test on, create a single test column
  if(length(group)>1){
    sample_map_tb$new_group <- sample_map_tb[group] %>% tidyr::unite(col = new_group, sep = "_") %>% unlist()
    groupname <- paste(group, collapse = "_")
    group <- groupname
    colnames(sample_map_tb)[match("new_group", colnames(sample_map_tb))] <- groupname
  }

  groupname <- group
  if(!is.null(exclude_group)){
    groupname <- paste0(group, "_excl_", paste(exclude_group, collapse = "-"))
  }

  ##convert a group into another
  if(!is.null(convert_group)){
    if(!length(convert_group)==2){
      stop("Require a vector of length 2 in 'convert_group' to convert from, to...")
    } else {
      sample_map_tb[group][sample_map_tb[group] == convert_group[1]] <- convert_group[2]
      groupname <- paste0(group, "_conv", "_", convert_group[1])
    }
  }

  dir.create(paste0(out_dir, "/", groupname), showWarnings = FALSE, recursive = TRUE)

  if(!is.null(tech_reps)){
    reps <- as.numeric(factor(unlist(sample_map_tb[tech_reps])))
  } else {
    reps <- seq(1:dim(sample_map_tb)[1])
  }

  samples <- sample_map_tb[sample_ID]
  expdesign <- sample_map_tb %>%
               dplyr::mutate(replicate = reps) %>%
               dplyr::mutate(label = samples) %>%
               dplyr::select(label, condition = !!group, replicate) %>%
               dplyr::filter(!condition %in% grep(!!exclude_group, condition, value = TRUE)) %>%
               as.data.frame()
  expdesign$ID <- expdesign$label

  matched <- match(expdesign$label, colnames(prot_data_tb))
  prot_data_tb <- dplyr::select(.data = prot_data_tb, row_data, matched)
  colnames(prot_data_tb)[matched] <- expdesign$label
  prot_data_tbz <- prot_data_tb[, expdesign$label]
  row_dat <- data.frame(name = unlist(prot_data_tb[row_data]), ID = unlist(prot_data_tb[row_data]))
  rownames(row_dat) <- row_dat[,1]
  data_se <- SummarizedExperiment::SummarizedExperiment(assays = as.matrix(prot_data_tbz),
                                                        colData = expdesign,
                                                        rowData = row_dat)

  ##filter, normalise, impute, test_diffs
  # made_unique <- DEP::make_unique(proteins = prot_data_tb, names = name, ids = id)
  # matched <- match(expdesign$label, colnames(made_unique))
  # dep_se <- DEP::make_se(proteins_unique = made_unique, columns = matched, expdesign = expdesign)

  data_filt <- DEP::filter_missval(data_se, thr = 1)
  data_norm <- DEP::normalize_vsn(data_filt)
  data_imp_man <- DEP::impute(data_norm, fun = "man", shift = 1.8, scale = 0.3)
  data_diff_all_contrasts <- DEP::test_diff(se = data_imp_man,
                                       type = "all")
  dep <- DEP::add_rejections(data_diff_all_contrasts, alpha = p_adj_val)
  data_results <- DEP::get_results(dep)
  sig_results <- data_results %>% filter(significant)
  contrasts <- as_tibble(rowData(dep)) %>%
               dplyr::select(ends_with("_significant")) %>%
               colnames() %>%
               gsub("_significant","",.)

  if(plot_dep == TRUE){
    pdf(paste0("DEP/", groupname, "/", project_name, ".PCA_COR_HMs.", groupname, ".pdf"), onefile = TRUE)
      print(plot_pca(dep, x = 1, y = 2, n = dim(data_results)[1], point_size = 4))
      plot_cor(dep, significant = TRUE, lower = 0, upper = 1, pal = "Reds")
      plot_heatmap(dep, type = "centered", kmeans = TRUE,
                   col_limit = 4, show_row_names = FALSE,
                   indicate = c("condition", "replicate"))
      plot_heatmap(dep, type = "contrast", kmeans = TRUE,
                  col_limit = 10, show_row_names = FALSE)
      print(plot_numbers(data_filt))
      print(plot_coverage(data_filt))
      print(plot_normalization(data_filt, data_norm))
      print(plot_imputation(data_norm, data_imp_man))
    dev.off()

    pdf(paste0(out_dir, "/", groupname, "/", project_name, ".volcanos.", groupname, ".pdf"))
      lapply(contrasts, function(f){
        print(f)
        print(plot_volcano(dep, contrast = f, label_size = 2, add_names = TRUE))
      })
    dev.off()

  }

  ##run our own 'all' contrast limma to access t statistic
  se <- data_imp_man
  col_data <- colData(se)
  col_data$individual <- unlist(lapply(col_data$label, function(f){
                              paste(strsplit(f, "_")[[1]][1:2], collapse="_")
                         }))
  prot_data_tb <- assay(se)
  design <- model.matrix(formula(~0 + condition), data = col_data)
  colnames(design) <- gsub("condition", "", colnames(design))
  condits <- as.character(unique(col_data$condition))
  cntrst <- apply(utils::combn(condits, 2), 2, paste,
            collapse = " - ")

  ##duplicate correlation because of replicate individuals
  prot_data_tb.cor <- duplicateCorrelation(prot_data_tb,
                                  block = col_data$individual)
  print(paste0("Duplicate correlation: ", prot_data_tb.cor$consensus.correlation))
  fit <- lmFit(prot_data_tb,
               design = design,
               block = col_data$individual,
               correlation = prot_data_tb.cor$consensus.correlation)
  made_contrasts <- makeContrasts(contrasts = cntrst, levels = design)
  contrast_fit <- contrasts.fit(fit, made_contrasts)
  eB_fit <- eBayes(contrast_fit)
  retrieve_fun <- function(comp, fit = eB_fit) {
     res <- topTable(fit, sort.by = "t", coef = comp, number = Inf,
         confint = TRUE)
     res <- res[!is.na(res$t), ]
     fdr_res <- fdrtool::fdrtool(res$t, plot = FALSE, verbose = FALSE)
     res$qval <- fdr_res$qval
     res$lfdr <- fdr_res$lfdr
     res$comparison <- rep(comp, dim(res)[1])
     res <- rownames_to_column(res)
     return(res)
  }
  limma_res <- map_df(cntrst, retrieve_fun)

  ##make a per comparison set of sig results
  ##change naming to R friendly (no spaces, dashes...)
  limma_res$comparison <- unlist(lapply(limma_res$comparison, function(f){
      gsub(" - ", "_vs_", f)
    }))
  limma_res_list <- lapply(unique(limma_res$comparison), function(f){
      limma_res %>% dplyr::filter(adj.P.Val < 0.05,
                                  comparison %in% f) %>%
                    as_tibble()
    })
  names(limma_res_list) <- unique(limma_res$comparison)

  ##plot PCA
  plot_se <- se
  ggpcase <- BMplotPCAse(plot_se, intgroup = "condition", pc_limit = 10, pchz = c(16,17,3))
  pdf(paste0(out_dir, "/", groupname, "/", project_name, ".PCAse.", groupname, ".pdf"))
  print(ggpcase)
  dev.off()
  png(paste0(out_dir, "/", groupname, "/", project_name, ".PCAse.", groupname, ".png"), res = 1200)
  print(ggpcase)
  dev.off()

  volc_list <- lapply(unique(limma_res$comparison), function(f){
    lres <- limma_res[limma_res$comparison %in% f,]
    maxlogfc <- max(lres[lres$adj.P.Val < p_adj_val,]$logFC)+2
    minlogfc <- min(lres[lres$adj.P.Val < p_adj_val,]$logFC)
    maxp <- -log10(min(lres[lres$adj.P.Val < p_adj_val,]$adj.P.Val))
    ggp <- EnhancedVolcano::EnhancedVolcano(lres,
                                     lab = lres$rowname,
                                     x = "logFC",
                                     y = "adj.P.Val",
                                     title = f,
                                     xlim = c(minlogfc, maxlogfc),
                                     ylim = c(0, maxp),
                                     subtitle = paste0("Cut-off p.adj < ", p_adj_val),
                                     pCutoff = p_adj_val)
    return(ggp)
  })

  pdf(paste0(out_dir, "/", groupname, "/", project_name, ".enh_volcanos.", groupname, ".pdf"))
  print(volc_list)
  dev.off()
  png(paste0(out_dir, "/", groupname, "/", project_name, ".enh_volcanos.", groupname, ".png"), res = 1200)
  print(volc_list)
  dev.off()

  plots_list <- list(ggpcase, volc_list)
  save(dep, data_results, sig_results, data_se, plot_se, contrasts, limma_res, limma_res_list, plots_list, file = paste0(out_dir, "/", groupname, "/", project_name, ".", groupname, ".RData"))

  return(list(dep = dep, dep_results = data_results, dep_sig_results = sig_results, data_se = data_se, plot_se = plot_se, contrasts = contrasts, limma_res = limma_res, limma_res_list = limma_res_list, plots_list))
}

#' Plotting PCA function
#' @param se SummarizedExperiment object
#' @param intgroup which colname from colData(se) to be output on plot
#' @param pc_limit integer percent variance accounted by PC for inclusion in plots  (default: 10%)
#' @param colz colour vector for samples
#' @param pchz pch vector for samples
#' @return list of ggplot2 objects for printing (PCA, PCs, loadings)
#' @export

BMplotPCAse <- function(se, intgroup = NULL, pc_limit = 10, pchz = NULL) {

  ##https://stackoverflow.com/questions/8197559/emulate-ggplot2-default-color-palette
  gg_color_hue <- function(n) {
    hues = seq(15, 375, length = n + 1)
    hcl(h = hues, l = 65, c = 100)[1:n]
  }

  pca <- prcomp(t(SummarizedExperiment::assay(se)))
  sdPc <- apply(pca$x, 2, sd)
  percentVar <- sdPc^2/sum(sdPc^2)
  sdRot <- apply(pca$rotation, 2, sd)
  percentVarRot <- sdRot^2/sum(sdRot^2)
  if (!all(intgroup %in% names(SummarizedExperiment::colData(se)))) {
      stop("The argument 'intgroup' should specify columns of colData(xx)")
  }
  intgroup.df <- as.data.frame(SummarizedExperiment::colData(se)[, intgroup, drop = FALSE])
  group <- factor(apply(intgroup.df, 1, paste, collapse = " : "))

  if(is.null(pchz)){
    scale_shape_man_vals <- 1:dim(unique(intgroup.df))[1]
    scale_colour_man_vals <- gg_color_hue(dim(unique(intgroup.df))[1])
  } else {
    scale_shape_man_vals <- pchz
    scale_colour_man_vals <- gg_color_hue(length(pchz)+1)[c(1,2,4)]
    print(scale_colour_man_vals)
  }

  ##lapply to make all PC > pc_limit included
  pcv_u <- percentVar[percentVar > pc_limit/100]
  if(length(pcv_u)==1){
    pcv_u <- percentVar[1:2]
  }
  ggps <- lapply(2:length(pcv_u), function(f){

      d <- data.frame(PC1 = pca$x[, 1], PC2 = pca$x[, f], intgroup.df, names = colnames(se))
      colnames(d)[colnames(d) == "PC2"] <- paste0("PC", f)
      ggp <- pca_plot(d, group = intgroup, intgroup.df, pcv = percentVar[c(1,f)], scale_shape_man_vals, scale_colour_man_vals)
      return(ggp)

  })

  ##scree plot showing contributions of PCs
  pcv <- data.frame(PC = colnames(pca$x), percent_variance = percentVar)
  pcv <- pcv[order(pcv$percent_variance, decreasing = TRUE),]
  levels(pcv$PC) <- pcv$PC
  ggs <- ggplot2::ggplot(data = pcv, ggplot2::aes(x = PC, y = percent_variance )) +
         ggplot2::geom_col() +
         ggplot2::ggtitle("Proprotion of Variances of Principle Components") +
         ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90))

  ##loading top 5, bottom 5 from 10 top PCs
  top_10_pc <- pcv[1:10,]
  top_10_pc$percent_variance <- paste0(round(top_10_pc$percent_variance, digits = 3)*100, "%")
  pca_rot_10 <- pca$rotation[,top_10_pc$PC]
  pca_rot_10 <- pca_rot_10[,colnames(pca_rot_10)[!is.na(colnames(pca_rot_10))]]

  fives_list <- lapply(colnames(pca_rot_10), function(f){
    print(f)
    fc <- c(tail(sort(pca_rot_10[,f]),5), head(sort(pca_rot_10[,f]),5))
    # nfc <- dplyr::filter(.data = anno_tb, ensembl_gene_id %in% names(fc))

    return(data.frame(gene = names(fc), loading = fc, PC = f))
  })

  loadings_plot <- do.call(rbind,fives_list) %>%
                   dplyr::left_join(top_10_pc) %>%
                   dplyr::mutate(PC = paste0(PC, "\n(", percent_variance, ")"))

  ggl <- ggplot2::ggplot(data = loadings_plot,
                         ggplot2::aes_string(x = "PC", y = "loading", label = "gene")) +
         ggplot2::geom_point(ggplot2::aes_string(colour = "loading"),
                             size = 3) +
         ggplot2::scale_shape_discrete(solid = T) +
         ggrepel::geom_text_repel(colour = "black",
                                  size = 2,
                                  fontface = "bold") +
         ggplot2::ggtitle("Loadings plot, top 10 PCs by variance, top/bottom 5 genes") +
         ggplot2::xlab("PC (% variance accounted for)")
  return(list(ggps, ggs, ggl))
}

pca_plot <- function(d, group, intgroup.df, pcv, scale_shape_man_vals, scale_colour_man_vals){

  if(nlevels(group)>6){
    ggp <- ggplot2::ggplot(data = d,
                           ggplot2::aes_string(x = "PC1",
                                               y = colnames(d)[2],
                                               group = group)) +
           ggplot2::scale_shape_manual(values = scale_shape_man_vals) +
           ggplot2::scale_colour_manual(values = scale_colour_man_vals) +
           ggplot2::labs(paste0("PCA plot using ", colnames(intgroup.df)), x = paste0("PC1: ", round(pcv[1] *  100), "% variance"),y = paste0(colnames(d)[2], ": ", round(pcv[2] * 100), "% variance")) +
           ggrepel::geom_text_repel(label = rownames(d),
                                    colour = "black",
                                    size = 2,
                                    fontface = "bold") +
           ggplot2::geom_point(ggplot2::aes_string(shape = group,
                                                   colour = group,
                                                   fill = group),
                               size = 3) +
           ggplot2::ggtitle(paste0("PCA plot using ", colnames(intgroup.df)),
                            subtitle = paste0(colnames(d)[1], " vs. ", colnames(d)[2]))
    }

    if(nlevels(group)<=6){
    ggp <- ggplot2::ggplot(data = d, ggplot2::aes_string(x = "PC1",
                        y = colnames(d)[2],
                        group = group)) +
                        # ggplot2::aes(x = PC1, y = PC2, group = group, shape = group, colour = group)) +
           ggplot2::geom_point(ggplot2::aes_string(shape = group,
                                                   colour = group,
                                                   fill = group),
                               size = 3) +
           ggplot2::scale_shape_manual(values = scale_shape_man_vals) +
           ggplot2::scale_colour_manual(values = scale_colour_man_vals) +
           ggplot2::xlab(paste0("PC1: ", round(pcv[1] *  100), "% variance")) +
           ggplot2::ylab(paste0(colnames(d)[2], ": ", round(pcv[2] * 100), "% variance")) +
           ggrepel::geom_text_repel(label = rownames(d),
                                    colour = "black",
                                    size = 2,
                                    fontface = "bold") +
                                    # ggplot2::annotate("text",x=pca$x[,1], y = pca$x[,2]-0.4, label = colnames(x), cex = 1.6) +
           ggplot2::ggtitle(paste0("PCA plot using ", colnames(intgroup.df)),
                            subtitle = paste0(colnames(d)[1], " vs. ", colnames(d)[2]))
  }
  return(ggp)
}
