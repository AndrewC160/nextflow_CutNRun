#' @title Pretty numbers
#'
#' @description
#' Given a number value, determine if it should be reported in thousands
#' (k, 1e3), millions (M, 1e6), billions (B, 1e9), or trillions (T, 1e12). Round
#' off digits to the <digs> place. Force a unit with <force unit>, which can be
#' any of "T","B","M","k", or "" (for no units). Values can be in string format,
#' provided that they are coercible to double (values which do not convert will
#' be returned as empty strings), and commas will be ignored such that
#' "-223,442.05" can be read as "-223k" by default. NA values will be treated as
#' zero values.
#'
#' @param num_in Number to convert to a human-readable string.
#' @param drop_zero Should zero values be returned as ""? Defaults to FALSE, in which case "0" is returned.
#' @param force_unit Unit to convert numbers to, can be "T", "B", "M", "k" or "" for no units. Defaults to NULL, in which case output is decided automatically.
#' @param digs Digits to round numbers to. Defaults to 2.
#'
#' @import tibble
#' @import dplyr
#' @import magrittr
#' @import tidyr
#'
#' @export

prettyNumbers <- function(num_in,drop_zero=FALSE,force_unit=NULL,digs=2){
  if(length(num_in) > 1){
    out_p <- sapply(num_in,function(x) prettyNumbers(num_in=x,drop_zero=drop_zero,force_unit=force_unit,digs=digs))
  }else{
    n_val <- num_in %>%
      gsub(",","",.) %>%
      as.double %>%
      suppressWarnings
    if(is.na(n_val)) n_val <- 0
    if(n_val == 0){
      if(drop_zero){
        out_p <- ""
      }else{
        out_p <- "0"
      }
    }else{
      if(is.null(force_unit)){
        force_unit  <- case_when(
          abs(n_val) >= 1e12 ~ "T",
          abs(n_val) >= 1e9 ~ "B",
          abs(n_val) >= 1e6 ~ "M",
          abs(n_val) >= 1e3 ~ "k",
          TRUE ~ "")
      }
      out_p   <- case_when(
        force_unit == "T" ~ n_val/1e12,
        force_unit == "B" ~ n_val/1e9,
        force_unit == "M" ~ n_val/1e6,
        force_unit == "k" ~ n_val/1e3,
        TRUE ~ n_val) %>%
        round(digits=digs) %>%
        prettyNum(big.mark=",") %>%
        paste0(force_unit)
    }
  }
  return(out_p)
}
