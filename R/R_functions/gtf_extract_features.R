#' @title GTF Extract Features
#'
#' @description
#' Given a GRange of genomic regions and a Tabix-indexed GTF file, read all
#' features in the specified region(s) and return as a GRanges object. Any
#' attributes from column 9 are processed into a table with NA values, and
#' these can be restricted using the <retain_attributes> argument.
#'
#' If a large number of features is being retrieved, <do_not_process_attributes>
#' can be used to skip regex steps. By default, if T/F are not provided to this
#' argument then processing is skipped with a warning if more than 10k features
#' are returned from the scanning step. If processing is skipped, column 9 is
#' returned as a semicolun-sparated list of attribute;"value" pairs.
#'
#' @param gr_window GRanges of regions to extract information from.
#' @param gtf_file Filename of Tabix-indexed GTF file to extract data from.
#' @param retain_attributes List of attributes to retain from column 9. For example, c("gene_name","transcript_id") will only keep gene names and transcript IDs while discarding the rest. Entries with one or both of these values missing will have NA values in their columns. No effect if <do_not_process_attributes> is TRUE.
#' @param do_not_process_attributes Should attribute column be returned unchanged (semi-colon separated list)? No default, in which case if more than 10000 features are returned processing is skipped with a warning.
#'
#' @import dplyr
#' @import tidyr
#' @import magrittr
#' @import GenomicRanges
#' @import Rsamtools
#' @import stringr
#' @import data.table
#'
#' @export

gtf_extract_features <- function(gr_window,gtf_file=NULL,retain_attributes=NULL,do_not_process_attributes){
  gtf_fl  <- gtf_file %||% system.file("extdata", "Homo_sapiens.GRCh38.104.chr.tabix.gtf.gz",package = "clugPac")
  gtf_fl  <- Rsamtools::TabixFile(gtf_fl)
  gr      <- gr_window
  seqlevelsStyle(gr) <- "Ensembl"
  gtf_txt <- Rsamtools::scanTabix(gtf_fl,param = gr)
  col_nms <- c("seqnames",
               "source",
               "feature",
               "start",
               "end",
               "score",
               "strand",
               "frame",
               "attribute")
  
  if(identical(gtf_txt[[1]],character(0))){
    tb_gtf <- setNames(rep(NA,length(col_nms)),col_nms) %>% rbind %>% as_tibble %>% filter(row_number() != 1)
  }else{
    tb_gtf <- gtf_txt %>%
      unlist %>%
      paste(collapse="\n") %>%
      fread(text=.,sep="\t",col.names=col_nms) %>%
      as_tibble
    if(missing(do_not_process_attributes)){
      if(nrow(tb_gtf) > 10000){
        do_not_process_attributes <- TRUE
        message("More than 10k features returned, skipping attribute processing. Force with <do_not_process_attributes>=TRUE.")
      }else{
        do_not_process_attributes <- FALSE
      }
    }
    if(!do_not_process_attributes){
      split_atts <- function(att_string=tb_gtf$attribute[1]){
        txt <- str_match_all(att_string,'([[:alnum:]_]+) "([:alnum:]+)"')
        setNames(txt[[1]][,3],txt[[1]][,2])
      }
      process_atts<- function(att_list=tb_gtf$attribute,keep_atts=NULL){
        atts_all  <- lapply(att_list,split_atts)
  
        att_nms   <- keep_atts %||% sapply(atts_all,names) %>% unlist %>% unique
  
        lapply(atts_all,function(att_v){
          setNames(att_v[att_nms],att_nms)
        }) %>% do.call(rbind,.) %>%
          as_tibble %>%
          return()
      }
  
      tb_atts <- process_atts(tb_gtf$attribute,keep_atts=retain_attributes)
      tb_gtf  <- tb_gtf %>%
        select(-attribute) %>%
        cbind(tb_atts) %>%
        makeGRangesFromDataFrame(keep.extra.columns = TRUE)
    }
  }
  return(tb_gtf)
}

# gtf_extract_features(gr_window,retain_attributes = "none")
# gtf_extract_features(gr_window,retain_attributes = c("gene_name","gene_source"))
# gtf_extract_features(gr_window = get_seqsizes(as_granges = TRUE)[1],retain_attributes = "none")


