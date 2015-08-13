

#' @title d_unk_sex distribute deaths of unknown sex
#' 
#' @param Deaths the standard internal Deaths data.frame, early in the production chain.
#' 
#' @return Deaths rows containing UNK are removed, and these deaths have been distributed proportionally over other deaths of the same Year, Age, Lexis and AgeInterval and known Sex. If this is not possible for the given ID combination, we split the deaths equally between the sexes. 
#' 
#' @importFrom compiler cmpfun
#' @importFrom data.table data.table
#' 
#' @export
#' 
# Deaths <- data.frame(Year = c(rep(2010,20)))
d_unk_sex <- function(Deaths){
  # this function can be run on any deaths data, whether it has UNK or not
  if (!any(Deaths$Sex == "UNK")){
    return(Deaths)
  }
  # dim(Deaths)
  Deaths <- Deaths[Deaths$Age != "TOT", ]
 
  # a very conservative inner function. Better to name than to make anonymous
  d_unkYrAgeLex <- cmpfun(function(DeathsYearAgeLex){
      if (any(DeathsYearAgeLex$Sex == "UNK")){
        UNKi          <- DeathsYearAgeLex$Sex == "UNK"

        # in case there are several of the same age/year/lexis/interval combo with UNK sex
        # this is probably a rather useless contingency just for the sake
        # of robustness... Might actually be used in historical Japanese data,
        # however for Sex = "UNK" and Age = "UNK", fo delayed registrations...
        if (sum(UNKi) > 1){
          UnkRow           <- DeathsYearAgeLex[UNKi,][1,]
          UnkRow$Deaths    <- sum(DeathsYearAgeLex$Deaths[UNKi])
          DeathsYearAgeLex <- rbind(DeathsYearAgeLex[!UNKi, ], UnkRow)
          # recreate this
          UNKi             <- DeathsYearAgeLex$Sex == "UNK"
        }
        UNK               <- sum(DeathsYearAgeLex$Deaths[UNKi], na.rm = TRUE )
        KnownSex          <- DeathsYearAgeLex[!UNKi, ]
        
        # in case it's a unique shape, split 50-50
        if (nrow(KnownSex) == 0){
          KnownSex        <- rbind(DeathsYearAgeLex[UNKi,],DeathsYearAgeLex[UNKi,])
          KnownSex$Deaths <- UNK / 2 # if no observed sex for age/Lexis, split 50 / 50
          KnownSex$Sex    <- c("m", "f")
        } else {
          # special case 2: only one sex observed, still split 50-50
          if (nrow(KnownSex) == 1) {
            sx              <- KnownSex$Sex
            KnownSex        <- rbind(KnownSex, KnownSex)
            KnownSex$Sex    <- c("m", "f")
            KnownSex$Deaths[KnownSex$Sex != sx] <- 0
            KnownSex$Deaths <- KnownSex$Deaths + UNK / 2
          } else {
          # distribute proportionally
          KnownSex$Deaths <-  KnownSex$Deaths + 
            (KnownSex$Deaths / sum(KnownSex$Deaths, na.rm = TRUE)) * UNK
        }
        KnownSex <- assignNoteCode(KnownSex, "d_unk_sex()")
      }
        return(invisible(KnownSex))
      } 
      invisible(DeathsYearAgeLex)
    })
  
  Deaths <- data.table(Deaths)
  
  # it won't convert integer to double automatically... need to precoerce
  Deaths$Deaths <- as.numeric(Deaths$Deaths)
  
  # much faster than split apply combine
  D.Out <- as.data.frame(Deaths[,d_unkYrAgeLex(.SD),by = list(Year,Age,Lexis,AgeInterval)])

  rownames(D.Out) <- NULL
  
  D.Out <- resortDeaths(D.Out)

  invisible(D.Out)
  
}
