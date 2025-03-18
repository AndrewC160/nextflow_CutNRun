#' @title Parse FastQC zips.
#'
#' @description Given a FastQC output zip file, extract its content to a temp
#' directory and scan fastqc_data.txt. Remove temporary files, then parse data
#' file into constituent tables. If multiple Zip files are included, run
#' function recursively on each and append each into the same tables. Note that
#' tables with missing data are returned with the line "No data," and if run
#' recursively this may mean that some samples are not included in tables at
#' all ('No data' entries are omitted, so missing samples implicitly had no
#' data).
#'
#' Note that  some FastQC versions appear to insert additional lines before a
#' given table's header row; these are ignored (for instance, the sequence
#' duplication levels table may have a "#Total deduplication percentage" tag).
#'
#' @param zip_fastq FastQC Zip file(s) (with .zip ending) to be parsed.
#' @param stay_silent Should scan() report how many lines were read? Defaults to FALSE.
#'
#' @import tibble
#' @import dplyr
#' @import stringr
#' @import data.table
#' @import stringi
#'
#' @export

read_fastqc <- function(zip_fastqc,stay_silent=TRUE){
  if(length(zip_fastqc) > 1){
    dat_out <- lapply(zip_fastqc,read_fastqc,stay_silent=stay_silent)
    #ttl_sets<- names(dat_out[[1]])
    ttl_sets <- unique(unlist(sapply(dat_out,names,simplify = FALSE)))
    out_tbs <- lapply(ttl_sets,function(x){
      tb_out<- lapply(dat_out,function(y){
        y[[x]]
      })
      tb_out<- setdiff(tb_out,"No data") %>%
        do.call(rbind,.)
      return(tb_out)
    })
    names(out_tbs) <- ttl_sets
  }else{
    if(!grepl(".zip$",zip_fastqc)) stop("Input file should be a zipped FastQC file (with .zip ending).")
    base_nm   <- gsub(".zip","",basename(zip_fastqc))
    out_dir   <- file.path(tempdir(),"fastqc_out")
    unzip(zip_fastqc,exdir = out_dir,junkpaths = TRUE)
    txt_dat   <- scan(file.path(out_dir,"fastqc_data.txt"),what=character(),sep="\n",quiet = stay_silent)
    unlink(file.path(out_dir,"base_nm"))

    tb_summary<- tibble(start_row = grep(">>[^END]",txt_dat),
                        end_row = grep(">>END",txt_dat)) %>%
      mutate(section = gsub(">>","",txt_dat[start_row])) %>%
      separate(section,into=c("title","flag"),sep = "\t") %>%
      mutate(name = gsub("_fastqc$","",base_nm))

    out_tbs   <- apply(tb_summary,1,function(row_in) {
      st_row  <- as.integer(row_in['start_row']) + 1
      en_row  <- as.integer(row_in['end_row']) - 1
      if(st_row >= en_row){
        tb_out<- "No data"
      }else{
        txt <- txt_dat[st_row:en_row]
        #Some FastQC versions appear to arbitrarily stuff things before the table prefixed with a #; throw that out.
        header_row <- max(grep("^#",txt))
        txt <- txt[header_row:length(txt)]
        tb_out  <- fread(text=txt,sep = "\t") %>%
          mutate(name = gsub("_fastqc$","",base_nm)) %>%
          rename_all(function(x) gsub("#","",gsub(" ","_",gsub("'","",tolower(x)))))
      }
      return(tb_out)
    })
    out_tbs   <- c(list(tb_summary),out_tbs)
    names(out_tbs)  <- c("summary",gsub(" ","_",tolower(tb_summary$title)))
  }
  return(out_tbs)
}
