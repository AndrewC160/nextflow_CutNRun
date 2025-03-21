#' @title
#' Read Summary table Peak Calling Results.
#'
#' @description
#' Given a pipeline's pooled output summary.tsv, read all peak files into a 
#' single GRanges object annotated with replicate/pool name and peak type 
#' (narrowPeaks, broadPeaks, and summits).
#'
#' @param summary_tsv Filename of pipelineTSV file to read from.
#'
#' @import dplyr
#' @import tidyr
#' @import magrittr
#' @import stringr
#' @import data.table
#' @import GenomicRanges
#'
#' @export

read_summary_peaks  <- function(summary_tsv){
  mutate  <- dplyr::mutate
  filter  <- dplyr::filter
  tb_lst  <- fread(summary_tsv) %>%
    filter(field %in% c("npks_reps","bpks_reps","npks_pool","bpks_pool","sums_pool")) %>%
    group_by(field) %>%
    mutate(rep = row_number()) %>%
    ungroup %>%
    mutate(name = 
      ifelse(grepl("^sums",field),
             gsub("_summits\\.bed$","",basename(value)),
             gsub("_peaks\\..+Peak$","",basename(value))),
      type = str_extract(value,"[:alpha:]+$")) %>%
    apply(1,function(row_in){
      if(row_in[['type']] == "bed"){
        tb_out <- fread(row_in[['value']],col.names=c("seqnames","start","end","name","nlog10q")) %>%
          mutate(
            score = 0,fc_summit=0,
            name = row_in[['name']],
            set = "summits")
      }else if(row_in[['type']] == "narrowPeak"){
        tb_out <- read_peak_file(row_in[['value']]) %>%
          mutate(name = row_in[['name']]) %>%
          makeGRangesFromDataFrame(keep.extra.columns = TRUE) %>%
          annotate_and_reduce_gr(
            retain_funcs = list(score = max,
                                fc_summit = max,
                                nlog10q = max)
          ) %>%
          as_tibble %>%
          select(-strand) %>%
          mutate(
            name=row_in[['name']],
            set = "narrowPeaks")
      }else{
        tb_out <- read_peak_file(row_in[['value']]) %>%
          select(-strand,-nlog10p) %>%
          mutate(
            name = row_in[['name']],
            set = "broadPeaks")
      }
      tb_out <- select(tb_out,seqnames,start,end,name,score,fc_summit,nlog10q,set)
      return(tb_out)
    }) %>% do.call(rbind,.) %>%
    makeGRangesFromDataFrame(keep.extra.columns = TRUE) %>%
    return()
}