#' @title
#' GTF to genes
#'
#' @description
#' Given a GTF file, extract all genes and produce a GenomicRanges object. If a
#' cache file is provided, save a TSV copy in that location. If such a file
#' already exists and <overwrite_cache> == FALSE, load the existing TSV file.
#'
#' @param gtf_file Filename of GFF with gene annotations. If not provided, defaults to packaged "Homo_sapiens.GRCh38.104.chr.tabix.gtf.gz"
#' @param cache_file Filename of TSV file in which to store annotations. Defaults to NULL.
#' @param overwrite_cache Boolean, defaults to FALSE. Should TSV cache file be rebuilt from scratch?
#' @param as_grange Should data be returned as a GRange object? Defaults to TRUE.
#'
#' @import data.table
#' @import dplyr
#' @import tidyr
#' @import GenomicRanges
#' @import stringr
#'
#' @export

gtf_to_genes    <- function(gtf_file,cache_file=NULL,overwrite_cache=FALSE,as_grange=TRUE){
  select  <- dplyr::select
  mutate  <- dplyr::mutate
  filter  <- dplyr::filter
  rename  <- dplyr::rename
  group_rows <- dplyr::group_rows

  #If a cache file is provided, load that.
  tb_out  <- NULL
  if(!is.null(cache_file) & !overwrite_cache){
    if(file.exists(cache_file)){
      tb_out  <- fread(cache_file) %>% as_tibble
    }
  }

  #If no cached file, build the table fresh.
  if(is.null(tb_out)){
    # if(missing(gtf_file)){
    #   gtf_file <- system.file("extdata","Homo_sapiens.GRCh38.104.chr.tabix.gtf.gz",package = "clugPac")
    # }
    is_gzipped <- grepl(".gz$",gtf_file)

    tb_out <- paste(ifelse(is_gzipped,"zcat","cat"),
                gtf_file,"| awk '$3 == \"gene\"'") %>%
      fread(cmd = .,
            header=FALSE,
            col.names = c("seqnames",
                          "source",
                          "type",
                          "start",
                          "end",
                          "V6",
                          "strand",
                          "V8",
                          "vals")) %>%
      as_tibble %>%
      select(-V6,-V8) %>%
      mutate(seqnames= paste0("chr",seqnames)) %>%
      rowwise %>%
      mutate(vals = list(unlist(str_split(vals,pattern = ";")))) %>%
      ungroup %>%
      unnest(vals) %>%
      mutate(vals = trimws(vals),
             vals = gsub('"',"",vals),
             name = stringr::str_extract(vals,"^[[:alpha:]_]+"),
             val = stringr::str_extract(vals,"[^ ]+$")) %>%
      select(-vals) %>%
      filter(!is.na(name)) %>%
      pivot_wider(id_cols = c(seqnames,start,end,strand,type,source),
                  names_from = name,
                  values_from = val) %>%
      select(seqnames,start,end,strand,gene_name,gene_biotype,gene_id,gene_version,source,type) %>%
      mutate(gene_name = ifelse(is.na(gene_name),"",gene_name))
    #If a cached file was provided, save a copy there.
    if(!is.null(cache_file)){
      write.table(tb_out,cache_file,sep = "\t",quote = FALSE,col.names = TRUE,row.names = FALSE)
      message("Cached table to '",cache_file,"'.")
    }
  }
  tb_out <- mutate_at(tb_out,c("gene_biotype","source"),as.factor)
  if(as_grange){
    tb_out  <- makeGRangesFromDataFrame(tb_out,keep.extra.columns = TRUE)
  }
  return(tb_out)
}
