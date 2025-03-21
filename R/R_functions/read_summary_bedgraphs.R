#' @title
#' Read Summary table bedGraph results.
#'
#' @description
#' Given a pipeline's pooled output summary.tsv, read all signal (bedGraph) 
#' files into a tibble with annotations for treatment, control, name, and 
#' spike fraction.
#'
#' @param summary_tsv Filename of pipelineTSV file to read from.
#'
#' @import dplyr
#' @import tidyr
#' @import magrittr
#' @import stringr
#' @import data.table
#'
#' @export

read_summary_bedgraphs <- function(summary_tsv){
  mutate  <- dplyr::mutate
  filter  <- dplyr::filter
  rename  <- dplyr::rename
  select  <- dplyr::select
  
  tb_fls  <-  fread(summary_tsv) %>% as_tibble
  
  # Get spikes.
  v_spikes <- tb_fls %>%
    filter(field == "tsv_spike") %>%
    mutate(name = gsub("_spike.tsv","",basename(value))) %>%
    rowwise %>%
    mutate(frac = {
      v_vals <- fread(value)$V2
      v_vals[2]/sum(v_vals)
    }) %>%
    vectify(frac,name)
  
  # Get bedGraphs.
  tb_fls %>%
    filter(field %in% c("bdgs_reps","bdgs_reps_ctrl","bdgs_pool","bdgs_pool_ctrl")) %>%
    mutate(name = gsub("_[treat|control].+","",basename(value)),
           set = ifelse(grepl("pool",field),"pool","reps"),
           cond = ifelse(grepl("ctrl$",field),"control","treat")) %>%
    select(-field) %>%
    mutate(spike = v_spikes[name],
           spike = ifelse(is.na(spike),mean(v_spikes),spike)) %>%
    rename(file = value) %>%
    return()
}
