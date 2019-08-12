##function to take sleuth object, or simple matrix of values, and a list of genes
##produce expected plot, including tag if required
##so = sleuth object; else data.frame of tpm/srpb with genes as rownames
##conds = metadata with columns to be used
##condcol = name of column data to be used in plot colouring
##factors = factors in condcol to use in plot (two only!)
##genes = those to plot, must match in 'so' rownames
##mapgenes = tibble/df with external_gene_id and matching ensembl_gene_id of genes (reversible!)
##genemap = which 2 columns to map from and to (NULL <- ensembl_gene_id, external_gene_name)
plot_expression_so_srpb <- function(so, genes, mapgenes, conds, condcol, factors, tag, genemap=NULL){

    if(is.null(genemap)){
      genemap <- c("ensembl_gene_id", "external_gene_name")
    }

    if(class(so)=="sleuth"){
      log2_obs_norm_df <- dcast(so$obs_norm,
                                       target_id ~ sample,
                                       value.var="scaled_reads_per_base") %>%
                                dplyr::arrange(target_id) %>%
                                dplyr::mutate_if(is.numeric, log2) %>%
                                dplyr::filter(target_id %in% genes) %>%
                                as.data.frame() %>%
                                column_to_rownames(.,var="target_id")
    }
    if(class(so)!="sleuth"){
      log2_obs_norm_df <- so
    }
    log2srpb_genes <- log2_obs_norm_df

    #ensure genes required are in data, print out missing
    if(length(rownames(log2srpb_genes)) != length(genes)){
      print(paste0("Missing genes: ", genes[! genes %in% rownames(log2srpb_genes)]))
      genes <- genes[genes %in% rownames(log2tpm_genes)]
    }
    print(paste0("Working on: ", paste(rownames(log2srpb_genes), collapse=", ")))

    print(rownames(log2srpb_genes))
    if(dim(log2srpb_genes)[1]==0){
        print("No input found, check genes are in rownames of log2srpb input")
    }
    else{
        Factor <- rep(factors[1],dim(conds)[1])
        Factor[c(conds[condcol] == factors[2])]<-factors[2]

        mltgg <- data.frame(0,0,0)
        mltgg <- mltgg[-1,]
        genenames <- mapgenes %>% dplyr::filter(.data[[genemap[1]]] %in% genes) %>% dplyr::select(.data[[genemap[2]]])
        genenames <- unlist(c(genenames))

        for(xx in 1:length(genes)){
          mltgg <- rbind(mltgg, cbind(rep(genenames[xx], length(Factor)),
                                      melt(log2srpb_genes[rownames(log2srpb_genes) %in% genes[xx],]),
                                      Factor))
        }

        colnames(mltgg) <- c("Gene","sample","value","Factor")
        mltgg$Factor <- factor(mltgg$Factor,
                               levels = factors,
                               ordered = TRUE)
        ggp <- ggplot(data = mltgg, aes(x = Gene, y = value)) +
               geom_boxplot(aes(colour = Factor)) +
               geom_jitter(aes(colour = Factor), position=position_dodge(0.8)) +
               labs(y = "log2SRPB", title = tag, subtitle="Expression per Group per Gene") +
               scale_colour_manual(values = c("blue", "red")) +
               theme(axis.text.x = element_text(angle = 45, hjust = 1))
    }
    return(ggp)
}


plot_expression_so_tpm <- function(so, genes, mapgenes, conds, condcol, factors, tag, genemap=NULL){

    if(is.null(genemap)){
      genemap <- c("ensembl_gene_id", "external_gene_name")
    }

    if(class(so)=="sleuth"){
      log2_obs_norm_df <- dcast(so$obs_norm,
                                       target_id ~ sample,
                                       value.var="tpm") %>%
                                dplyr::arrange(target_id) %>%
                                dplyr::mutate_if(is.numeric, log2) %>%
                                dplyr::filter(target_id %in% genes) %>%
                                as.data.frame() %>%
                                column_to_rownames(.,var="target_id")
    }
    if(class(so)!="sleuth"){
      log2_obs_norm_df <- so
    }
    log2tpm_genes <- log2_obs_norm_df

    #ensure genes required are in data, print out missing
    if(length(rownames(log2tpm_genes)) != length(genes)){
      print(paste0("Missing genes: ", genes[! genes %in% rownames(log2tpm_genes)]))
      genes <- genes[genes %in% rownames(log2tpm_genes)]
    }
    print(paste0("Working on: ", paste(rownames(log2tpm_genes), collapse=", ")))

    if(dim(log2tpm_genes)[1]==0){
        print("No input found, check genes are in rownames of input")
    }
    else{
        Factor <- rep(factors[1],dim(conds)[1])
        Factor[c(conds[condcol] == factors[2])]<-factors[2]

        mltgg <- data.frame(0,0,0)
        mltgg <- mltgg[-1,]
        genenames <- mapgenes %>% dplyr::filter(.data[[genemap[1]]] %in% genes) %>% dplyr::select(.data[[genemap[2]]])
        genenames <- unlist(c(genenames))

        for(xx in 1:length(genes)){
          mltgg <- rbind(mltgg, cbind(rep(genenames[xx], length(Factor)),
                                      melt(log2tpm_genes[rownames(log2tpm_genes) %in% genes[xx],]),
                                      Factor))
        }

        colnames(mltgg) <- c("Gene","sample","value","Factor")
        mltgg$Factor <- factor(mltgg$Factor,
                               levels = factors,
                               ordered = TRUE)
        ggp <- ggplot(data = mltgg, aes(x = Gene, y = value)) +
               geom_boxplot(aes(colour = Factor)) +
               geom_jitter(aes(colour = Factor), position=position_dodge(0.8)) +
               labs(y = "log2TPM", title = tag, subtitle="Expression per Group per Gene") +
               scale_colour_manual(values = c("blue", "red")) +
               theme(axis.text.x = element_text(angle = 45, hjust = 1))
    }
    return(ggp)
}
