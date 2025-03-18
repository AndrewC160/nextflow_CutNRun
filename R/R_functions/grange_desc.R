#' @title
#' GRanges description
#'
#' @description
#' Given a GRanges object, returns a simple annotation in the format:
#'
#' "Chr1: 1,000,000 - 2,000,000 (1MB)"
#'
#' @param grange_in GRange to describe, should be of length one (if not, only the first will be described.)
#' @param force_unit Enforce a desired unit (bp, kb, MB) for width readout.
#' @param digits Number of digits to display, defaults to 1.
#' @param sep Separation between width value and units, defaults to "".
#' @param append_ending Insert a string between the width printout and the closing parenthesis; intended for bin size if relevant (for instance, "(1MB, 2kb bin size)").
#'
#' @import stringr
#' @import GenomicRanges
#'
#' @export

grange_desc     <- function(grange_in,force_unit=NULL,digits=1,sep="",append_ending=""){
  gr <- grange_in[1]
  paste0(prettyTitle(as.character(seqnames(gr))),": ",
         prettyNum(GenomicRanges::start(gr),big.mark=",")," - ",
         prettyNum(GenomicRanges::end(gr),big.mark=",")," (",
         prettyBP(GenomicRanges::width(gr),force_unit = force_unit,digits = digits,sep = sep),append_ending,")") %>%
    return()
}
