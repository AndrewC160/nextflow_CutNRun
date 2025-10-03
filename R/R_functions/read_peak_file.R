#'@title Read peak file
#'
#'@description Read in broad / narrowPeak files and append the correct column
#'names. Determines which file type using the file extension, so those must be
#'correct. Otherwise, is_narrowPeak can be set as TRUE or FALSE to override the
#'default detection behavior.
#'
#'A note re. MACS3 "score" values: these equate to int(-10*log10qvalue).
#'
#'@param file_name File to be read.
#'@param name_col If provided, a name column will be appended to the file. Defaults to "base_name," in which case the baseneame of the peaks file is used.
#'@param is_narrowPeak If provided as TRUE or FALSE, will override the detection component of the function to establish file as narrowPeak (TRUE) or broadPeak (FALSE).
#'
#'@import dplyr
#'@import data.table
#'@import tidyr
#'@import stringr
#'
#'@export

read_peak_file <- function(file_name,name_col=NULL,is_narrowPeak) {
  if(missing(is_narrowPeak)){
    ext <- tolower(str_extract(file_name,"\\.[:alnum:]+$"))
    if(ext %in% c(".narrowpeak",".np")){
      is_narrowPeak <- TRUE
    }else if(ext %in% c(".broadpeak",".bp")){
      is_narrowPeak <- FALSE
    }else{
      stop("File extension '",ext,"' not recognized.")
    }
  }
  # Suppress warnings in case fread() encounters empty peak files. Desired behavior is to just return
  # an empty tibble with appropriate column numbers/names (since this function frequently runs on
  # lists of files followed by rbind()).
  if(is_narrowPeak){
    col_nms<- c("seqnames","start","end","name","score","strand","fc_summit","nlog10p","nlog10q","summit_pos")
    tb_out <- fread(file_name,col.names = col_nms) %>% suppressWarnings
  }else{
    col_nms<- c("seqnames","start","end","name","score","strand","fc_summit","nlog10p","nlog10q")
    tb_out <- fread(file_name,col.names = col_nms) %>% suppressWarnings
  }
  if(nrow(tb_out) == 0){
    if(!is.null(name_col)) col_nms<- c(col_nms,"set")
    tb_out <- setNames(rep(NA,length(col_nms)),col_nms) %>% rbind %>% as_tibble %>% filter(row_number() != 1)
  }else{
    if(!is.null(name_col)){
      if(name_col == "base_name"){
        tb_out<- mutate(tb_out,set = basename(file_name))
      }else{
        tb_out<- mutate(tb_out,set = name_col)
      }
    }
  }
  return(tb_out)
}
