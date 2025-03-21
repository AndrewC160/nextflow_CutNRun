#' @title
#' Randomize GRanges
#'
#' @description
#' Given an input GRange, randomly sample up to <max_sample> features within
#' this range and use regioneR's randomizeRegions() function to randomly
#' relocate them on the genome. Repeat this process <sample_num> times, and
#' return a single GRange with a column "set" labeling which iteration it was
#' produced in ("rand_3", for instance). Given that this function uses
#' regioneR, BSgenome.Hsapiens.UCSC.hg38 must be installed.
#'
#' @param gr_in GRange to randomize.
#' @param max_sample Maximum number of features in <gr_in> to use per iteration.
#' @param sample_num Number of random iterations to perform. Defaults to 3.
#' @param genome_in String containing the genome that specified regions correspond to. Defaults to "hg38".
#'
#' @import regioneR
#' @import BSgenome.Hsapiens.UCSC.hg38
#'
#' @export

randomize_grange<- function(gr_in,max_sample=10000,sample_num=3,genome_in="hg38"){
  if(length(gr_in) == 0){
    gr_out  <- gr_in
  }else{
    max_sample <- min(length(gr_in),max_sample)
    gr_out  <- sapply(1:sample_num, function(i) {
      gr <- randomizeRegions(sample(gr_in,size=max_sample,replace=FALSE),genome = genome_in,allow.overlaps = TRUE)
      gr$set = paste0("rand_",i)
      return(gr)
    }) %>%
      do.call(c,.)
  }
  return(gr_out)
}
