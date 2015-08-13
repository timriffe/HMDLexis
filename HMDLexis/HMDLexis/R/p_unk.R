#' @title p_unk distribute popualtion counts of unknown ages.
#' 
#' @param Pop the standard internal Population data.frame, probably late down the production chain.
#' 
#' @return Pop rows containing UNK are removed, and these population counts have been distributed proportionally over other deaths of known age.
#' 
#' @importFrom compiler cmpfun
#' 
#' @export
#' 

p_unk <- function(Pop){
  Pop  <- Pop[!Pop$Age == "TOT", ]
  # this function can be run on any deaths data, whether it has UNK or not
  if (!any(Pop$Age == "UNK")){
    return(Pop)
  }
  p_unkYrSex <- cmpfun(function(PopYearSex){
      if (any(PopYearSex$Age == "UNK")){
        UNK                   <- PopYearSex$Population[PopYearSex$Age == "UNK"] 
        PopYearSex            <- PopYearSex[PopYearSex$Age != "UNK", ]
        PopYearSex$Population <- PopYearSex$Population + (PopYearSex$Population / sum(PopYearSex$Population)) * UNK
        PopYearSex            <- assignNoteCode(PopYearSex, "p_unk")
      }
      invisible(PopYearSex)
    })
  Pop <- do.call(rbind, lapply(X = split(Pop, list(Pop$Year, Pop$Sex)), FUN = p_unkYrSex))
  rownames(Pop) <- NULL
  invisible(resortPops(Pop))
}

