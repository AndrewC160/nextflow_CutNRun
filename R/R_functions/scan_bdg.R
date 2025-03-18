#' @title
#' Scan Tabix-indexed BedGraph file.
#'
#' @description
#' Scan one or more BGzipped and Tabix-indexed files within regions specified
#' in <gr_regions> GRange object. Function is run recursively if more than one
#' bgzipped file is specified, and results are concatenated into a single
#' table. Uses Rsamtools::TabixFile for scanning and returns a tibble. Note
#' that if the input gr_regions have names assigned, these will be substituted
#' into the region field. Otherwise this field will contain the as.character()
#' coercions of the gr_regions themselves.
#'
#' Function allows for "coarse downsampling" via the <coarse_downsampling_lines>
#' argument, as well. If the number of lines retrieved from a Tabix file is
#' greater than this threshold, the retrieved lines are randomly sampled to keep
#' approximately this number. This can ensure that large genomic regions can
#' still be represented while not retrieving tens of millions of lines of data
#' per sample.
#'
#' @param bgz_files List of filenames of tabix-indexed BDG files.
#' @param name_vec Vector of names for 'name' column. If missing, names can be retrieved from bgz_files, otherwise will use file base name. Must be the same length as bgz_files.
#' @param gr_regions GenomicRanges containing regions of interest.
#' @param col_names Column names of BedGraph file. Defaults to standard 5 columns (seq,start,end,score,name), but can be changed to include additional columns if necessary.
#' @param coarse_downsample_lines Threshold for corase downsampling. Defaults to Inf, in which case no coarse downsampling is done.
#'
#' @import Rsamtools
#' @import dplyr
#' @import tidyr
#' @import data.table
#' @import GenomicRanges
#'
#' @export

scan_bdg  <- function(bgz_files,name_vec,gr_regions,col_names,coarse_downsample_lines=Inf){
  if(missing(col_names)){
    col_names   <- c("seqnames","start","end","score","region")
  }
  if(length(bgz_files) > 1){
    if(!missing(name_vec)){
      if(length(name_vec) != length(bgz_files)){
        stop("Different number of files/names provided.")
      }
      names(bgz_files) <- name_vec
    }
    if(is.null(names(bgz_files))){
      names(bgz_files) <- paste0("bdg_",c(1:length(bgz_files)))
    }
    tb_out  <- lapply(names(bgz_files),function(nm) {
      scan_bdg(bgz_files = bgz_files[[nm]],
               gr_region = gr_regions,
               col_names=col_names,
               coarse_downsample_lines = coarse_downsample_lines) %>%
        as_tibble %>%
        mutate(name = nm)
    }) %>%
      do.call(rbind,.)
  }else{
    tb_file <- Rsamtools::TabixFile(file = bgz_files)
    txt_lst <- Rsamtools::scanTabix(tb_file,param=gr_regions)
    if(length(txt_lst) == 0){
      return(NULL)
    }
    txt_lst <- lapply(txt_lst, function(lst_in){
      if(length(lst_in) > coarse_downsample_lines){
        kp_lns  <- sample(1:length(lst_in),size=coarse_downsample_lines,replace=FALSE) %>% sort
        lst_out <- lst_in[kp_lns]
      }else{
        lst_out <- lst_in
      }
      return(lst_out)
    })

    if(length(txt_lst) > coarse_downsample_lines){
      kp_lns  <- sample(1:length(txt_lst),size=coarse_downsample_line,replace=FALSE) %>% sort
      txt_lst <- txt_lst[kp_lns]
    }
    if(missing(name_vec)){
      name_value  <- gsub(".bdg.gz","",basename(tb_file$path))
    }else{
      name_value  <- name_vec
    }
    tb_out  <- lapply(names(txt_lst), function(nm) {
      txt   <- txt_lst[[nm]]
      paste0(txt,"\t",nm)
    }) %>%
      do.call(c,.)
    if(length(tb_out) == 1){
      #Add newline in cases where only one value returned.
      tb_out <- paste0(tb_out,"\n")
    }
    tb_out  <- fread(text=tb_out,sep="\t",col.names = col_names) %>%
      as_tibble %>%
      mutate(name = name_value)

    #See if names are present in gr_regions.
    if(!is.null(names(gr_regions))){
      reg_nms<- gr_regions
      strand(reg_nms) <- "*"
      reg_nms<- as.character(reg_nms)

      gr_nms <- setNames(names(gr_regions),reg_nms)
      tb_out <- mutate(tb_out,region = gr_nms[region])
    }
  }
  return(tb_out)
}
