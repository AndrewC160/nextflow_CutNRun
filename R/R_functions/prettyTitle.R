#' @title
#' prettyTitle
#'
#' @description
#' Given a string, function capitalizes either the first or every word (separated
#' by 'sep'). Runs recursively if given a list or vector, and if the input is a
#' factor with levels and 'preserve_levels' is TRUE, existing levels are also
#' changed and then reapplied. Otherwise, vectors are returned as character type.
#'
#' @param str_in String or vector of strings to convert.
#' @param sep Separator between words. Can also be set to "camel" for CamelCase, in which case all capitol letters are treated as the first letter in a given word.
#' @param caps Which words should be capitalized? Defaults to "first" in which case only the first word, but can also be "all" for all words. To leave specific words in lower case, use override_words with gsub-compatible regex. For instance, when caps = "all" and override_words=c("a ","with"), "apples_with_a_pigeon" becomes "Apples with a Pigeon" (note the use of a space in "a ").
#' @param preserve_levels Boolean; should existing factor levels be updated and preserved (i.e., levels are also changed and order is maintained)? Defaults to TRUE. No effect if input is not a factor.
#' @param preserve_acronyms Should acronyms be detected and preserved? Treats subsequent capitol letters as acronyms.
#' @param override_words Vector of words to be over-ridden. All output words/levels are compared without case to words in override_words, and if they match they are over-ridden. For instance, "NcDNA" or "mirna" can be replaced with "ncDNA" and "miRNA."
#'
#' @import stringr
#' @import stringi
#'
#' @export

prettyTitle <- function(str_in,sep="_",caps="first",preserve_levels=TRUE,preserve_acronyms=TRUE,override_words=NULL){
  cap_func <- function(txt_in = "somatic"){
    paste0(toupper(substr(txt_in,1,1)),substr(txt_in,2,nchar(txt_in)))
  }
  camel_func <- function(txt_in = "CamelCase",save_acks=TRUE){
    #Parse CamelCase and snakeCase, and send to lowercase, but respect acronyms.
    c_txt   <- unlist(str_split(txt_in,pattern=""))

    #First letter is always boundary capitol (i.e, even in snakeCase).
    cap_pos <- c(1,grep("[A-Z]",unlist(str_split(txt_in,pattern = "")),ignore.case = FALSE)) %>% unique

    #Ignore subsequent capitols ("DNABaseLevel" becomes "DNA base level")
    is_first<- c(TRUE,!unlist(sapply(seq_along(cap_pos), function(i) cap_pos[i] == cap_pos[i-1] + 1)))
    is_last <- c(na.omit(unlist(sapply(seq_along(cap_pos), function(i) cap_pos[i+1] > cap_pos[i] + 1))),TRUE)
    is_term <- cap_pos == nchar(txt_in)
    #Boundary: Capitol OR first capitol in acronym.
    bounds  <- cap_pos[is_first | (is_last & !is_term)]

    #Do not lowercase acronyms (acronym: previous letter was capitol end of
    #string, UNLESS last capitol is not the last letter in the string).
    is_ack  <- c(na.omit(sapply(seq_along(cap_pos), function(i) !is_first[i+1])),FALSE) | is_term
    acks    <- cap_pos[is_ack]
    if(any(acks) & save_acks){
      c_txt[-acks] <- tolower(c_txt[-acks])
    }else{
      c_txt <- tolower(c_txt)
    }

    sapply(seq_along(bounds), function(i) {
      w_strt <- as.integer(bounds[i])
      w_end  <- as.integer(bounds[i + 1]) - 1
      if(is.na(w_end)){
        w_end<- nchar(txt_in)
      }
      return(paste(c_txt[w_strt:w_end],collapse=""))
    })
  }

  if(length(str_in) > 1){
    txt <- sapply(str_in,prettyTitle,sep=sep,caps=caps,preserve_levels=FALSE,override_words=override_words)#,preserve_levels=preserve_levels)
    if(is.factor(str_in) & preserve_levels){
      lv_nms <- prettyTitle(levels(str_in),sep=sep,caps=caps,override_words=override_words)
      txt <- factor(txt,levels=lv_nms)
    }
  }else{
    if(sep=="camel"){
      txt <- camel_func(str_in,save_acks = preserve_acronyms)
    }else{
      txt <- unlist(str_split(str_in,pattern=sep,simplify = FALSE))
    }
    if(caps == "first"){
      txt <- paste(cap_func(txt[1]),paste(txt[-1],collapse=" "),sep=" ")
    }else if(caps == "all"){
      txt <- paste(sapply(txt,cap_func),collapse=" ")
    }else{
      stop("Argument 'caps' can be either 'first' or 'all'.")
    }
    if(!is.null(override_words)){
      for(wrd in override_words){
        txt <- gsub(wrd,wrd,txt,ignore.case=TRUE)
      }
    }
    txt <- trimws(txt)
  }
  return(txt)
}
