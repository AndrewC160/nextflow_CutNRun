#' @title
#' Vectify
#'
#' @description
#' Starting with an input tibble, generate a named vector using one column for
#' names and another for values.
#'
#' @param table_in Input table.
#' @param value_col Name of column from which to derive values.
#' @param name_col Name of column from which to derive names.
#'
#' @import dplyr
#' @import stats
#'
#' @export

vectify         <- function(table_in,value_col,name_col){
  vls   <- select(table_in,!!as.name(deparse(substitute(value_col)))) %>%
    unlist(use.names=FALSE)
  nms   <- select(table_in,!!as.name(deparse(substitute(name_col)))) %>%
    unlist(use.names=FALSE) %>%
    as.character
  setNames(vls,nms)
}
