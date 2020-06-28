library("tidyverse")
library("sleuth")

highest_exp_wide <- function(wide_object, name_col, id_col, value_col){

  ##a situation arises where a single gene name maps to multiple gene identifiers
  ##e.g. ENSG00000271503, ENSG00000274233 both map to CCL5
  ##this function takes wide format input with at least 3 rows:
  ##1: gene name
  ##2: gene identifier
  ##3: expression value(s) [NB can be index of columns]
  ##function returns a wide format object with highest expression name+id
  ##also returns full mapping info and mean expression value per name, id pairs

  ##create mean_value of all rows, then count name_col
  ##filter out all with name_col single mapping
  mean_map <- wide_object %>% dplyr::select(name_col, id_col, value_col) %>%
                              dplyr::mutate(mean_value = rowMeans(.[value_col])) %>%
                              group_by_at(vars(name_col)) %>%
                              na.omit(mean_value) %>%
                              add_count()%>%
                              dplyr::filter(n > 1) %>%
                              dplyr::filter(mean_value == max(mean_value)) %>%
                              ungroup() %>%
                              distinct_at(vars(-id_col), .keep_all = TRUE) %>%
                              dplyr::select(-n, -mean_value)
  ##remove those with multi mapping (in mean_map) from original input
  mean_map_o <- anti_join(wide_object, mean_map, by=name_col)

  ##return the two sets, bound
  bind_rows(mean_map_o, mean_map)
}
