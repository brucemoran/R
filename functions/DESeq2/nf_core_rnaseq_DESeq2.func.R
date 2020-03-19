##used as a follow-on from aligning with STAR, featureCounts
##as per standard nf-core/rnaseq pipeline (v1.4 currently)

##run this from where 'nextflow run' and 'results' held
libs <- c("DESeq2", "tidyverse", "biomaRt", "ggplot2", "RColorBrewer", "gplots", "plyr")
libsLoaded <- lapply(libs,function(l){suppressMessages(library(l, character.only = TRUE))})
source("https://raw.githubusercontent.com/brucemoran/R/master/functions/biomaRt/biomaRt_anno_orth.func.R")

##takes input allowing DESeq2 analysis; further allows annotation if non-human, and writing of TPM from stringTie to annotated human-ortholog table
nf_core_rnaseq_DESeq2 <- function(TAG = NULL, SAMPLE_COND = NULL, COND_DESIGN = NULL, REF_COND = NULL, OUTDIR = NULL, ENS_VERSION = NULL, ORG_PREFIX = NULL, DELIM_SAMP = NULL){

  if(is.null(TAG)){
    print("Please provide a TAG to name your outputs")
  }

  if(is.null(SAMPLE_COND)){
    print("Please specify a SAMPLE_COND metadata/condition file in CSV format")
    break
  }

  if(is.null(COND_DESIGN)){
    print("Please specify a COND_DESIGN for DESeq2")
    print("NB this should be format: first + second + last")
    print("NBB that 'last' in design is condition of interest in DE analysis")
    print("NBBB that all elements of the design must be names of columns in SAMPLE_COND")
  }

  if(is.null(REF_COND)){
    print("Please specify a REF_COND reference condition for DESeq2")
    print("NB that this should be the control")
  }

  if(is.null(DELIM_SAMP) | is.na(DELIM_SAMP)){
    DELIM_SAMP = "\\."
  }

  ##read in merged featureCounts results
  ##parse names using "." for nicer colnames
  fc_merge <- read_tsv("results/featureCounts/merged_gene_counts.txt")
  coln_samples <- colnames(fc_merge[3:length(colnames(fc_merge))])
  n_coln_samples <- unlist(lapply(coln_samples, function(f){
    strsplit(f, DELIM_SAMP)[[1]][1]
  }))
  colnames(fc_merge) <- c("ensembl_gene_id",
                          "external_gene_name",
                          n_coln_samples)
  on_coln_samples <- n_coln_samples[order(n_coln_samples)]
  fc_merge_so <- fc_merge %>% dplyr::select(1, 2, on_coln_samples) %>%
           arrange(ensembl_gene_id)
  write_tsv(fc_merge_so, paste0(OUTDIR, "/", TAG, ".merged_gene_counts.fc.tsv"))

  fc_co <- fc_merge_so %>%
           dplyr::select(1, on_coln_samples) %>%
           as.data.frame() %>%
           column_to_rownames("ensembl_gene_id")

  ##remove zero-count genes
  nzero <- function(x){x[apply(x,1,sum)>0,]}
  nz_fc_co <- nzero(fc_co)

  ##read in condition data
  cond <- read_csv(SAMPLE_COND)
  cond_df <- as.data.frame(cond)
  rownames(cond_df) <- cond_df[,1]
  cond_df <- cond_df[,-1]

  ##break design into components
  design_vec <- unlist(lapply(COND_DESIGN, function(f){gsub(" ", "", strsplit(f, "\\+")[[1]])}))
  CONDITION <- rev(design_vec)[1]
  cond_df[,CONDITION] <- factor(cond_df[,CONDITION])
  cond_df[,CONDITION] <- relevel(cond_df[,CONDITION], ref = REF_COND)
  for(x in 2:dim(cond_df)[2]){
    cond_df[,rev(design_vec)[x]] <- factor(cond_df[,rev(design_vec)[x]])
  }

  ##DESeq2DataSet object
  dds <- DESeqDataSetFromMatrix(countData = nz_fc_co,
                                colData = cond_df,
                                design = formula(paste0("~ ", COND_DESIGN)))

  ##run DESeq2, make results
  ddseq <- DESeq(dds)

  ##make all contrasts of CONDITION, then set into named list
  combn_mat <- t(combn(levels(cond_df[,CONDITION]),2))
  combns <- apply(combn_mat, 1, function(f){paste(f, collapse="-")})
  res_list <- lapply(combns, function(f){
    ress <- na.omit(results(ddseq, contrast=c(CONDITION, strsplit(f, "-")[[1]])))
    ress_tf <- as_tibble(as.data.frame(ress), rownames="ensembl_gene_id")
    write_tsv(ress_tf, paste0(OUTDIR, "/", TAG, ".res.", f, ".DESeq2.tsv"))
    return(ress)
  })
  names(res_list) <- combns

  ##plots
  vsd <- vst(dds, blind=FALSE)
  sampleDists <- dist(t(assay(vsd)))
  sampleDistMatrix <- as.matrix(sampleDists)
  colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(256)
  hc <- hclust(sampleDists)

  pdf(paste0(OUTDIR, "/", TAG, ".heatmap.pdf"))
  heatmap.2(sampleDistMatrix, Rowv=as.dendrogram(hc),
            symm=TRUE, trace="none", col=colors,
            margins=c(2,10), labCol=FALSE, cexRow = 0.8 )
  dev.off()

  ##PCA
  bmpcaplot <- BMplotPCA(vsd, intgroup=c(CONDITION))
  ggsave(paste0(OUTDIR, "/", TAG, ".PCA.pdf"), bmpcaplot)

  ##annotate results if non-human
  if(!is.null(ORG_PREFIX)){
    biomart_ao <- biomaRt_anno_orth(GENOME1 = "hsapiens_gene_ensembl", GENOME2_PREFIX = ORG_PREFIX, VERSION = ENS_VERSION)
    str_write <- stringTie_res_join(ANNO_TABLE=biomart_ao)
  }
  if(is.null(ORG_PREFIX)){
    str_write <- stringTie_res_join()
  }
  write_tsv(str_write, path=paste0(OUTDIR, "/", TAG, ".stringTie.TPMs.tsv"))
}

BMplotPCA <- function(x, intgroup = NULL, ntop = 1500, returnData = FALSE) {
    rv <- rowVars(assay(x))
    select <- order(rv, decreasing = TRUE)[seq_len(min(ntop,
        length(rv)))]
    pca <- prcomp(t(assay(x)[select, ]))
    percentVar <- pca$sdev^2/sum(pca$sdev^2)
    if (!all(intgroup %in% names(colData(x)))) {
        stop("the argument 'intgroup' should specify columns of colData(dds)")
    }
    intgroup.df <- as.data.frame(colData(x)[, intgroup, drop = FALSE])
    group <- factor(apply(intgroup.df, 1, paste, collapse = " : "))
    d <- data.frame(PC1 = pca$x[, 1], PC2 = pca$x[, 2], group = group,
        intgroup.df, names = colnames(x))
    if (returnData) {
        attr(d, "percentVar") <- percentVar[1:2]
        return(d)
    }
    if(nlevels(group)>6){
      ggp <- ggplot(data=d, aes(x=PC1,y=PC2, group=group, colour=group, shape=group)) + scale_shape_manual(values=1:nlevels(group)) + labs(title=paste0("PCA plot using ",intgroup), x=paste0("PC1: ", round(percentVar[1] *  100), "% variance"),y=paste0("PC2: ", round(percentVar[2] * 100), "% variance")) + annotate("text",x=pca$x[,1],y=pca$x[,2],label=colnames(x),cex=1.6) + geom_point(size=3) + ggtitle(paste0("PCA plot using ",intgroup))
      }
      if(nlevels(group)<=6){
      ggp <- ggplot(data = d, aes(x = PC1, y = PC2, group=group, shape=group, colour=group)) + geom_point(size = 3) + scale_shape_discrete(solid=T) + xlab(paste0("PC1: ", round(percentVar[1] *  100), "% variance")) + ylab(paste0("PC2: ", round(percentVar[2] * 100), "% variance")) + annotate("text",x=pca$x[,1],y=pca$x[,2]-0.4,label=colnames(x),cex=1.6) + ggtitle(paste0("PCA plot using ",intgroup))
    }
    return(ggp)
}


stringTie_res_join <- function(ANNO_TABLE = NULL){

  ##FPKM/TPM from stringTie results
  st_res_tpm_list <- lapply(dir("results/stringtieFPKM", pattern="gene_abund.txt", full.names=TRUE), function(f){
    samp <- paste0(strsplit(strsplit(f,"\\/")[[1]][3], "\\.")[[1]][1], ".TPM")
    read_tsv(f) %>%
    dplyr::rename(ensembl_gene_id = 1) %>%
    dplyr::select(1, !!samp := TPM) %>%
    arrange(ensembl_gene_id)
  })
  st_res_tpm <- as_tibble(join_all(st_res_tpm_list, type="left"))

  if(is.null(ANNO_TABLE)){
    return(st_res_tpm)
  }
  if(!is.null(ANNO_TABLE)){
    coln <- colnames(ANNO_TABLE)[3]
    t1 <- st_res_tpm %>% dplyr::rename(!!coln := ensembl_gene_id)
    t2 <- left_join(ANNO_TABLE, t1)
    return(t2)
  }
}

genesGTF <- function(OUTDIR = OUTDR){

  ##get genes.gtf input
  genesgtf_file <- system("find work/stage -name genes.gtf | grep -v STAR", intern=T)
  system(paste0("cat ", genesgtf_file, " | perl -ane 'if(($F[1] eq \"ensembl\") && ($F[2] eq \"CDS\")){print $_;}' > ", OUTDIR, "/gene.ensembl.CDS.gtf"))
  cds_gtf_idnm <- read_table2(paste0(OUTDIR, "/gene.ensembl.CDS.gtf")) %>%
             dplyr::select(14,16)

  ##parse CDS from genes.gtf
  splitFun <- function(x){unlist(lapply(x, function(f){strsplit(f, '\\"')[[1]][2]}))}
  cds_gtf_nm <- tibble(ensembl_gene_id = splitFun(unlist(cds_gtf_idnm[,1])),
                       external_gene_name = splitFun(unlist(cds_gtf_idnm[,2]))) %>%
                distinct() %>%
                arrange(ensembl_gene_id)

  return(cds_gtf_nm)
}
