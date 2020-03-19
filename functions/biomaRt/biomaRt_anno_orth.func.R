##take two inputs of datasets in biomaRt::useMart
##return the tibble of the matched orthologs
library("biomaRt")
library("tidyverse")
biomaRt_anno_orth <- function(GENOME1 = "hsapiens_gene_ensembl", GENOME2_PREFIX = NULL, VERSION = NULL){

##need a way to access versions
##if you know you need a specific version for acces to a specific genome version
##then specify version; otehrwise defaults to latest by below
if(is.null(VERSION)){
  VERSION <- rev(strsplit(biomaRt::listMarts()[1,2], " ")[[1]])[1]
}

##specifies the HOST connect (URL)
HOST <- as_tibble(listEnsemblArchives()) %>%
        dplyr::filter(version %in% VERSION) %>%
        dplyr::select(url) %>% unlist()

##access the two genomes
mart1 <- biomaRt::useMart(biomart = "ensembl", dataset = GENOME1, host=HOST)
mart2 <- biomaRt::useMart(biomart = "ensembl", dataset = paste0(GENOME2_PREFIX, "_gene_ensembl"), host=HOST)

##first genome genes with homologs for second genome
ensid2gene1 <- as_tibble(biomaRt::getBM(attributes=c("ensembl_gene_id", "external_gene_name", paste0(GENOME2_PREFIX, "_homolog_ensembl_gene")), mart = mart1))

##second genome, with renaming to allow join with above
renames <- c(paste0(GENOME2_PREFIX,"_homolog_ensembl_gene"),
             paste0(GENOME2_PREFIX,"_homolog_external_name"))
ensid2gene2 <- as_tibble(biomaRt::getBM(attributes=c("ensembl_gene_id", "external_gene_name"), mart = mart2)) %>%
  dplyr::select(!!renames[1] := ensembl_gene_id, !!renames[2] := external_gene_name)

##join to get GENOME_1 (human?!) and homolog IDs
ensid2gene2orth <- left_join(ensid2gene1, ensid2gene2, by=renames[1]) %>%
                   na.omit()
return(ensid2gene2orth)
}
