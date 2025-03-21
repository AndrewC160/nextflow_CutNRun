#' @title
#' Get seqsizes
#'
#' @description
#' Retrieve chromosome sizes for hg38 genome. By default, uses included
#' hg38_seqsizes.tsv file from extdata. If the argument <coord_adjust> is given
#' (as a genomic range object), the adjusted start value of this range will be
#' returned as an integer. Otherwise, a vector of chromosome size values is
#' returned. If <as_granges> is true, this vector is converted into a GRanges
#' object with seq size adjustments in a meta column.
#'
#' @param coord_adjust GRanges coordinate to adjust; if provided, this coordinate will be returned with the appropriate adjustment made.
#' @param seqsize_file Defaults to NULL, in which case hg38_seqsizes.tsv is loaded from extdata. String of seqsizes TSV file to be used.
#' @param as_granges Boolean, should sizes be returned as a GRanges object containing one element per chromosome? Defaults to FALSE.
#'
#' @import dplyr
#' @import tidyr
#' @import magrittr
#' @import GenomicRanges
#'
#' @export

get_seqsizes    <- function(coord_adjust = NULL,seqsize_file=NULL,as_granges = FALSE){
  if(is.null(seqsize_file)){
    seqsize_file<- system.file("extdata","hg38_seqsizes.tsv",package="clugPac")
  }
  vec_out  <- read.table(seqsize_file,sep="\t",col.names=c("seqnames","bp")) %>%
    vectify(value_col = bp,name_col = seqnames)
  if(!is.null(coord_adjust)){
    return(setNames(start(coord_adjust) + vec_out[as.character(seqnames(coord_adjust))],NULL))
  }else if(as_granges){
    tibble(seqnames=names(vec_out),
           start = 1,end = vec_out) %>%
      mutate(start = as.double(start),
             end = as.double(end)) %>%
      mutate(adj_start = 1 + cumsum(dplyr::lag(end,default = 0))) %>%
      makeGRangesFromDataFrame(keep.extra.columns = TRUE) %>%
      return()
  }else{
    return(vec_out)
  }
}
