
#' @title d_unk distribute deaths of unknown ages
#' 
#' @param Deaths the standard internal Deaths data.frame, probably late down the production chain.
#' 
#' @return Deaths rows containing UNK are removed, and these deaths have been distributed proportionally over other deaths of known age.
#' 
#' @importFrom compiler cmpfun
#' @importFrom data.table data.table
#' 
#' @export
#' 

d_unk <- function(Deaths){
  # this function can be run on any deaths data, whether it has UNK or not
  if (!any(Deaths$Age == "UNK")){
    return(Deaths)
  }
  d_unkYrSex <- cmpfun(function(DeathsYearSex){
    if (any(DeathsYearSex$Age == "UNK")){
      UNKi          <- DeathsYearSex$Age == "UNK"
      UNK           <- sum(DeathsYearSex$Deaths[UNKi], na.rm = TRUE )
      DeathsYearSex <- DeathsYearSex[!UNKi, ]
      DeathsYearSex$Deaths <-  DeathsYearSex$Deaths + 
                     (DeathsYearSex$Deaths / sum(DeathsYearSex$Deaths, na.rm = TRUE)) * UNK
      DeathsYearSex <- assignNoteCode(DeathsYearSex, "d_unk()")
    }
    invisible(DeathsYearSex)
  })
  #Deaths <- do.call(rbind, lapply(X = split(Deaths, list(Deaths$Year, Deaths$Sex)), FUN = d_unkYrSex))
  
  Deaths <- data.table(Deaths)
  
  # it won't convert integer to double automatically... need to precoerce
  Deaths$Deaths <- as.numeric(Deaths$Deaths)
  
  # much faster than split apply combine
  D.Out <- as.data.frame(Deaths[,d_unkYrSex(.SD),by = list(Year,Sex)])
  
  rownames(D.Out) <- NULL
  
  D.Out <- resortDeaths(D.Out)
  
  invisible(resortDeaths(D.Out))
  
 }
