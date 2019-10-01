library("tidyverse")
library("sleuth")

##a situation arises where a single gene name maps to multiple gene identifiers
##e.g. ENSG00000271503, ENSG00000274233 both map to CCL5
##this function takes wide format input with at least 3 rows:
##1: gene name
##2: gene identifier
##3: expression value
##function returns a wide format object with highest expression name+id
##also returns full mapping info and mean expression value per name, id pairs

highest_exp_wide <- function(wide_object, name_col, id_col, value_col){

  ##create mean mapping of all, tally to count occurence
  mean_map <- wide_object %>%
              dplyr::select(name_col, id_col, value_col) %>%
              group_by_at(c(name_col, id_col)) %>%
              summarise_if(is.numeric, .funs = c(mean_value = mean)) %>%
              na.omit(mean_value) %>%
              add_tally()

  ##all multi-annotations
  mean_map_n <- mean_map %>% dplyr::filter(n > 1) %>%
                dplyr::filter(mean_value == max(mean_value)) %>%
                dplyr::select(-n) %>%
                add_tally()

  ##still multi-annotations
  mean_map_o <- mean_map_n %>% dplyr::filter(n > 1) %>%
                dplyr::filter(n > 1) %>%
                dplyr::select(-n)

  mean_map_o <- as_tibble(do.call(cbind, lapply(mean_map_o, function(f){ f[1] }))) %>% dplyr::mutate_at(3, as.numeric)

  ##format
  mean_map_1 <- mean_map %>% ungroup() %>%
                dplyr::filter(n == 1) %>%
                dplyr::select(-n)

  mean_map_2 <- mean_map_n %>% ungroup() %>%
                dplyr::filter(n == 1) %>%
                dplyr::select(-n)

  bind_rows(mean_map_1, mean_map_2, mean_map_o) %>%
  dplyr::rename(!!value_col := mean_value)
}
