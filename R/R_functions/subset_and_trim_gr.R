#' @title Subset and truncate GRanges
#'
#' @description Given one GRanges object of features and another of length one
#' representing a region of interest, subset all features to only include those
#' that overlap the region of interest, and trim start/end locations to fit
#' within that window. If <note_truncated> is TRUE, append start/end_offset
#' columns with the number of basepairs removed from each object (i.e., if the
#' feature fell 100bp before the start of the window, add 100bp to the start
#' position and enter "100" in the start_offset column). Function can be run
#' recursively by submitting multiple gr_windows in a GRange, but if more than
#' 1,000 entries are provided an error will display. This can be disabled by
#' setting <check_length>=FALSE.
#'
#' An example usage would be:
#'
#' subset_and_trim(gr_input = gene_list_grange, gr_window = regions_of_interest_gr)
#'
#' ...Which would return a GRanges object for each window containing the
#' truncated genes within each window.
#'
#' @param gr_input GRanges object to subset.
#' @param gr_window GRanges object containing one or more windows to <gr_input> to.
#' @param note_truncated Defaults to FALSE; if TRUE a start_offs and end_offs column will be added containing the number of basepairs trimmed from the start/end, respectively.
#' @param check_length Should function check that fewer than 1,000 query windows have been submitted? Defaults to TRUE.
#'
#' @import GenomicRanges
#' @import data.table
#'
#' @export

subset_and_trim_gr  <- function(gr_input,gr_window,note_truncated=FALSE,check_length=TRUE){
  start <- GenomicRanges::start
  end   <- GenomicRanges::end
  ifelse<- data.table::fifelse

  if(length(gr_window) > 1){
    if(check_length & length(gr_window) > 1000) stop("More than 1,000 windows specified in gr_window, if this is intentional run function with 'check_length=FALSE'.")
    gr_out  <- vector(mode="list",length = length(gr_window))
    names(gr_out)   <- as.character(gr_window)
    for(i in 1:length(gr_out)){
      gr_out[[i]]   <- subset_and_trim_gr(gr_input,gr_window[i])
    }
  }else{
    gr_out    <- subsetByOverlaps(gr_input,gr_window)
    win_start <- start(gr_window)
    win_end   <- end(gr_window)
    start_offs<- win_start - start(gr_out)
    end_offs  <- end(gr_out) - win_end

    start(gr_out) <- ifelse(start_offs > 0,win_start,start(gr_out))
    end(gr_out)   <- ifelse(end_offs > 0, win_end,end(gr_out))

    if(note_truncated){
      mcols(gr_out)$start_offset  <- ifelse(start_offs > 0,start_offs,0)
      mcols(gr_out)$end_offset    <- ifelse(end_offs > 0, end_offs, 0)
    }
  }
  return(gr_out)
}
