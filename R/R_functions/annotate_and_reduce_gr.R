#' @title Annotate and reduce GR
#'
#' @description
#' Given a GRanges object, flatten that object using GenomicRanges::reduce, but
#' specify columns and in-place functions via the <retain_funcs> list argument
#' to combine and retain column contents. For example, in the following list:
#'
#' list(score=max,padj=min,summits=list,annots=function(x) paste(unique(x),collapse=","))
#'
#' ...the "score" column would contain the maximum value from all entries,
#' "padj" would contain the minimum "padj" value, "summits" would contain a list
#' of all "summit" values, and "annots" would contain a string of comma-
#' separated unique values from "annots".
#'
#' Note that this function is considerably slower than a standard reduction, and
#' aggregation functions which involve string manipulation (such as the example
#' paste(unique(x)) above) are considerably slower than returning a list.
#'
#' @param gr_in GenomicRanges object to reduce.
#' @param retain_funcs List of functions to apply to retained column names, named as the columns to which they will apply. All column names must be present in the metadata columns of <gr_in> or an error will be thrown.
#' @param include_n Should a column 'n' tallying the total regions overlapping a reduced region be provided? Defaults to FALSE.
#' @param include_queryHits Should a column 'queryHits' containing a list of indices from the original GRange be included? Defaults to FALSE.
#'
#' @import GenomicRanges
#' @import tibble
#' @import dplyr
#' @import magrittr
#' @import tidyr
#' @import S4Vectors
#'
#' @export

annotate_and_reduce_gr <- function(gr_in,retain_funcs,include_n=FALSE,include_queryHits=FALSE,max_gap=1){
  mutate  <- dplyr::mutate
  arrange <- dplyr::arrange
  select  <- dplyr::select

  msng_nms <- setdiff(names(retain_funcs),colnames(mcols(gr_in)))
  if(length(msng_nms) > 0) stop("Column names '",paste(msng_nms,collapse=", "),"' not found in <gr_in>.")

  gr_red<- GenomicRanges::reduce(gr_in,min.gapwidth=max_gap)
  olaps <- findOverlaps(gr_in,gr_red)

  tb <- as_tibble(mcols(gr_in)[queryHits(olaps),names(retain_funcs),drop=FALSE])
  tb$subjectHits<- subjectHits(olaps)
  tb$queryHits  <- queryHits(olaps)
  tb <- group_by(tb,subjectHits)

  # Optional count and queryHits columns.
  if(include_n){
    tb <- mutate(tb,n=n())
  }
  if(include_queryHits){
    tb <- mutate(tb,queryHits=list(queryHits))
  }else{
    tb$queryHits <- NULL
  }

  for(col in names(retain_funcs)){
    tb <- mutate_at(tb,col,retain_funcs[[col]])
  }

  tb  <- tb %>%
    arrange(subjectHits) %>%
    distinct %>%
    ungroup %>%
    select(-subjectHits) %>%
    DataFrame()

  mcols(gr_red) <- tb
  return(gr_red)
}
