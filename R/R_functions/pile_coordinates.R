#' @title
#' Pile coordinates
#'
#' @description
#' Given a vector of start/end locations and optional chromosome names, return
#' an integer number corresponding to a stack position for each read (if read 2
#' overlaps read one, it gets a 2; if read 3 overlaps 1 and 2, it gets a three;
#' if read 4 doesn't overlap 1 through 3 it gets a 1, etc.). Vector names
#' correspond to index of the original vectors (1:n). All coordinates should be
#' the same length.
#'
#' @param start_vals List or vector of start coordinates.
#' @param end_vals List or vector of end coordinates.
#' @param seqname_vals List or vector of chromosome names.
#' @param min_gap Minimum distance in basepairs between two items on the genome before they are incremented to another track. Defaults to NUll (1bp).
#'
#' @import GenomicRanges
#' @import dplyr
#' @import S4Vectors
#'
#' @export

pile_coords     <- function(start_vals,end_vals,seqname_vals = NULL,min_gap=NULL){
  if(length(start_vals) != length(end_vals)){
    stop("Input vector lengths differ (",length(start_vals)," starts, ",length(end_vals)," ends).")
  }
  if(is.null(seqname_vals)){
    gr_in <- GRanges(seqnames="chr1",ranges=IRanges(start_vals,end_vals))
  }else{
    gr_in <- GRanges(seqnames=seqname_vals,ranges=IRanges(start_vals,end_vals))
  }
  #"Pile": any contiguous region of reads. Split this way to minimize work within "for" loop.
  piles     <- GenomicRanges::reduce(gr_in)
  gr_in$pile<- subjectHits(findOverlaps(gr_in,piles))

  if(!is.null(min_gap)){
    gr_in <- resize(gr_in,width=width(gr_in) + min_gap,fix="center")
  }

  tb_lst <- gr_in %>%
    as_tibble %>%
    mutate(idx = row_number(),ymin=0) %>%
    arrange(seqnames,start,end) %>%
    group_split(pile)

  lapply(tb_lst, function(tb_in) {
    x <- as.data.frame(tb_in)
    x_maxs <- c(0)
    for(i in 1:nrow(x)){
      st_pos <- x$start[i]
      en_pos <- x$end[i]
      o_laps <- x_maxs > st_pos
      if(all(o_laps)){
        #If current range overlaps all previous layers, add a new layer.
        x_maxs  <- c(x_maxs,en_pos)
        y_min   <- length(x_maxs)
      }else{
        #Otherwise, set Y value to the lowest layer not overlapped, update that layers x-max.
        y_min   <- min(which(x_maxs <= st_pos))
        x_maxs[y_min] <- en_pos
      }
      x$ymin[i] <- y_min
    }
    return(x)
  }) %>%
    do.call(rbind,.) %>%
    arrange(idx) %>%
    vectify(ymin,idx) %>%
    return()
}
