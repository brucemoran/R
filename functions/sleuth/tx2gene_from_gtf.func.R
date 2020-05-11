#! /usr/env/R

##take GTF file location as input, return a tx2gene object
##optional renaming of column names

library(rtracklayer)
tx2gene_from_gtf <- function(gtf, rn_gene_id = NULL, rn_gene_name = NULL, rn_transcript_id = NULL){

  my_obj <- import(gtf)

  as.data.frame(my_obj[my_obj$type == "transcript" & my_obj$source == "ensembl"]) %>%
  dplyr::select(!!rn_gene_id := gene_id, !!rn_gene_name := gene_name, !!rn_transcript_id := transcript_id)

}
