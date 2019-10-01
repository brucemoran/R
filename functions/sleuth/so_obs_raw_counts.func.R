library("tidyverse")
library("sleuth")

##so = input sleuth object
##group_col = column on which to group and sum counts
so_obs_raw_counts <- function(so, group_col=NULL){

    ##function to sum all transcripts from a gene-mode sleuth object into counts of those genes
    ##this is then saved into the sleuth object in so$obs_norm$counts

    ##use first of the target_mapping if group_col NULL
    if(is.null(group_col)){
      group_col <- colnames(so$target_mapping)[!colnames(so$target_mapping) %in% "target_id"][1]
    }
    print(paste0("group_col -> ", group_col))

    ##cols to include as samples
    long_col <- so$sample_to_covariates$sample

    if(class(so) != "sleuth"){
        print("Require a sleuth class object as input")
        break
    }
    else{
      if(so$gene_mode != "TRUE"){
        print("Require a sleuth object in gene_mode")
        break
      }
      else{
          so_obs_raw_count_wide <- so$obs_raw %>%
                              pivot_wider(id_cols = target_id, names_from = sample, values_from = est_counts) %>%
                              inner_join(so$target_mapping, .) %>%
                              group_by_at(group_col) %>%
                              summarise_if(is.numeric, .funs = c(sum="sum")) %>%
                              ungroup() %>%
                              inner_join(so$target_mapping, .) %>%
                              rename_at(unlist(map(long_col, starts_with, vars = colnames(.))), list(~gsub("_sum", "", .)))  %>%
                              dplyr::select(-target_id) %>% distinct() %>%
                              arrange(ensembl_gene_id)

          #long format
          so_obs_raw_count_long <- so_obs_raw_count_wide %>%
                                   pivot_longer(cols = long_col,
                                   names_to = "sample",
                                   values_to = "counts") %>%
                                   dplyr::select(1,2,4,3)

          so$obs_raw_count$wide <- so_obs_raw_count_wide
          so$obs_raw_count$long <- so_obs_raw_count_long
          return(so)
      }
    }
}
