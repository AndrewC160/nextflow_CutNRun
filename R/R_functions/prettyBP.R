#' @title
#' prettyBP
#'
#' @description
#' Express a number in basepairs with appropriate units.
#'
#' @param bp_val Numeric value of genomic distance in basepairs
#' @param force_unit If provided, unit will be forced. Can be "kb", "bp", "mb", or "MB" for kilobase, basepair, megabase, respectively.
#' @param digits Number of digits to be displayed; defaults to 1.
#' @param sep Separator between output number and unit, defaults to no space ("100MB").
#'
#' @import dplyr
#'
#' @export

prettyBP        <- function(bp_val,force_unit=NULL,digits=1,sep=""){
  if(length(bp_val) > 1){
    out_p <- sapply(bp_val,prettyBP,force_unit=force_unit,digits=digits,sep=sep)
  }else{
    if(!is.null(force_unit)){
      if(!force_unit %in% c("kb","bp","mb","MB")){
        stop("Unit '",force_unit,"' is unknown, can be kb, bp, or MB (case agnostic).")
      }else{
        unit_out<- force_unit
      }
    }else{
      unit_out<- case_when(abs(bp_val) < 1e3~"bp",
                           abs(bp_val) < 1e6~"kb",
                           TRUE ~ "MB")
    }
    bp_val <- as.double(bp_val)
    bp_val <-
      case_when(unit_out == "bp" ~ bp_val,
                unit_out == "kb" ~ bp_val / 1e3,
                unit_out == "mb" ~ bp_val / 1e6,
                unit_out == "MB" ~ bp_val / 1e6)
    is_int <- bp_val %% 1 == 0
    if(is.na(is_int)){ is_int <- FALSE } #Not sure what gives NA values when I run this within GGplot, but if that pops up change it to FALSE.
    if(is_int){
      bp_val  <- as.integer(bp_val)
    }else{
      bp_val  <- round(bp_val,digits=digits)
    }
    out_p <- bp_val %>%
      prettyNum(big.mark = ",") %>%
      paste0(sep,unit_out)
  }
  return(out_p)
}
