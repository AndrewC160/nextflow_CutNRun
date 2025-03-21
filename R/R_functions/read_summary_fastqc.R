#' @title
#' Read Summary table FastQC results.
#'
#' @description
#' Given a pipeline's pooled output summary.tsv, read all FastQC files into a 
#' nested list of tibbles.
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

read_summary_fastqc <- function(summary_tsv){
  v_fls <- fread(summary_tsv) %>%
    dplyr::filter(grepl("^fastqc",field)) %>%
    dplyr::mutate(set = str_extract(field,"[:alnum:]+$")) %>% 
    dplyr::mutate(
      set = case_when(set == "trim1"~"fastqR1",
                      set == "trim2"~"fastqR2",
                      set == "filt"~"BAM")) %>%
    group_by(field) %>%
    dplyr::mutate(rep = row_number()) %>%
    ungroup %>%
    dplyr::mutate(name = paste0(set,"_",rep)) %>%
    vectify(value,name)
  
  lapply(names(v_fls), function(nm){
    tbs <- read_fastqc(v_fls[[nm]])
    lapply(tbs,function(tb){
      tb_out <- NULL
      if(inherits(tb,"data.frame")){
        tb_out <- mutate(tb,file = nm)
      }
      return(tb_out)
    })
  }) %>% setNames(names(v_fls)) %>%
    return()
}
