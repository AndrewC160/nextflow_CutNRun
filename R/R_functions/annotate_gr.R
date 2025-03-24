#' @title
#' Annotate GR
#'
#' @description
#' Given a GRanges object (query) and a second GRanges object to compare it to
#' (subject), identify all overlaps and annotate query by adding a new column
#' (cols_query) with values from a specified column on subject (cols_subject).
#' For instance, if cols_query is set as "peak_overlaps" and cols_subject
#' contains names of each peak in subject, function will return query with a
#' new column called "peak_overlaps" with the associated names from subject. If
#' no overlaps are found for a feature, it is annotated with the value of
#' <na_val>, which defaults to an empty string. If more than one overlap is
#' detected, all values are combined in a list, which can be collapsed into a
#' string separated by semicolons (collapse_as_string; can be slow) or reduced
#' to its first element (first_only). Function can also accept multiple query
#' and subject columns, provided these are the same length.
#'
#' If a null value is provided to gr_query, columns of NA values will be added
#' to the original query and retruned.
#'
#' @param gr_query GRanges object to be annotated.
#' @param gr_subject GRanges object to be overlapped
#' @param cols_query Column(s) to be created in output GRanges object.
#' @param cols_subject Column(s) containing values to be inserted into query columns.
#' @param max_gap Maximum distance between features; defaults to 0bp. Must be set mutually exclusive to min_overlap.
#' @param min_overlap Minimum overlap between features; defaults to 1bp. Must be set mutually exclusive to max_gap.
#' @param as_boolean Should annotations be TRUE/FALSE (overlapped/not overlapped)? Defaults to FALSE. If TRUE, column values are ignored in favor of boolean values..
#' @param na_vals List of value(s) to annotate with in the event no overlaps are detected. Defaults to an empty string. Can be length 1 to apply to all columns, or one entry per column. Note that it must be a list if different data types are used in each column.
#' @param collapse_as_string Should annotations from multiple values be collapsed into a single semicolon-separated string? Defaults to FALSE, in which case these are returned as lists.
#' @param first_only Should only the first annotation be returned? Defaults to false; does not occur if <collapse_as_string> is TRUE.
#'
#' @import stringr
#' @import GenomicRanges
#' @import S4Vectors
#' @import tibble
#' @import dplyr
#' @import magrittr
#' @import tidyr
#'
#' @export

annotate_gr <- function(gr_query,gr_subject,cols_query=NULL,cols_subject=NULL,max_gap=-1L,min_overlap=0L,as_boolean=FALSE,na_vals = "",collapse_as_string=FALSE,first_only=FALSE){
  if(is.null(cols_query) & is.null(cols_subject)) stop("No column names specified for gr_query and/or gr_subject.")
  cols_query<- cols_query %||% cols_subject # If no query columns provided, default to the subject column names.

  if(is.null(gr_subject)){
    mtx <- matrix(ncol=length(cols_query),nrow = length(gr_query))
    colnames(mtx) <- cols_query
    mcols(gr_query) <- cbind(mcols(gr_query),mtx)
  }else{
    olaps   <- findOverlaps(query = gr_query,subject=gr_subject,maxgap = max_gap,minoverlap = min_overlap)
    if(length(as_boolean) == 1){
      as_boolean <- rep(as_boolean,length(cols_query))
    }else if(length(as_boolean) != length(cols_query)){
      stop("If provided, length of 'as_boolean' (",length(as_boolean),
           ") should be either 1 to apply to all columns OR equal to the number of columns (",
           length(cols_query),").")
    }
    if(length(na_vals) == 1){
      na_vals <- rep(na_vals,length(cols_subject))
    }else if(length(na_vals) != length(cols_subject)){
      stop("If <na_vals> are specified, they must either be length 1 to apply to all columns OR the same length as <cols_subject>.")
    }

    for(i in 1:length(cols_query)){
      q_col <- cols_query[i]
      s_col <- cols_subject[i]
      na_val<- na_vals[i]
      a_bool<- as_boolean[i]

      if(a_bool){
        mcols(gr_query)[q_col] <- FALSE
        mcols(gr_query[queryHits(olaps)])[q_col] <- TRUE
      }else{
        if(is.null(cols_subject)){
          stop("If appended column is not boolean, at least one cols_subject value must be provided.")
        }
        if(length(cols_query) != length(cols_subject)){
          stop("Query and subject column names must be the same length.")
        }

        q_tb <- as_tibble(olaps) %>%
          mutate(subject_vals = mcols(gr_subject)[,s_col][subjectHits]) %>%
          group_by(queryHits) %>%
          summarize(subject_vals = list(subject_vals),.groups="drop")
        if(collapse_as_string){
          q_tb <- q_tb %>%
            rowwise %>%
            mutate(subject_vals = paste(subject_vals,collapse=";"))
        }else if(first_only){
          q_tb <- q_tb %>%
            mutate(subject_vals = sapply(subject_vals,function(x) x[[1]]))
        }
        q_tb  <- as.data.frame(q_tb)
        mcols(gr_query)[q_col] <- na_val
        mcols(gr_query)[[q_col]][q_tb$queryHits] <- q_tb$subject_vals
      }
    }
  }
  return(gr_query)
}
