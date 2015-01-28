
#'
#' @title d_long adds 0s up to age 130 for years of deaths with no open age
#' 
#' @description Lexis triangles are created for the upper ages even if lower ages are RR or VV or something else. There may be overlap if the last age is not TL or TU, but it will be a 0 overlap, so it's innocuous. Intermediate gaps are also filled in. Year-Sex combos with open ages are left untouched.
#' 
#' @param Deaths the standard internal Deaths data.frame, probably early in its processing, but not necessarily.
#' 
#' @return Deaths the same deaths data.frame, but expanded to include ages up to 130 except for Year-Sex combos with pen ages, which are untouched. 
#' 
#' @importFrom compiler cmpfun
#' 
#' @export
#' 

d_long <- function(Deaths){
  
  # an internal function
  d_long_SexYr <- cmpfun(function(DeathsSexYr){
      if (any(DeathsSexYr$AgeInterval == "+", na.rm = TRUE)){
        return(DeathsSexYr)
      } else {
       
        Ages    <- rep(0:130, each = 2)
        Lexis   <- rep(c("TL","TU"),length(Ages) / 2)
        
        Dnew    <- as.data.frame(matrix(nrow = length(Ages), 
                              ncol = ncol(DeathsSexYr),
                              dimnames = list(NULL, colnames(DeathsSexYr))))
        
        Dnew$Year         <- unique(DeathsSexYr$Year)   
        Dnew$Age          <- Dnew$Agei <- Ages       
        Dnew$Lexis        <- rep(c("TL","TU"),length(Ages) / 2)
        Dnew$LDB          <- Dnew$AgeInterval  <-  Dnew$YearIntervali  <-  1
        Dnew$Access       <- "O"
        Dnew$Area         <- unique(DeathsSexYr$Area)
        Dnew$NoteCode1    <- "d_long()"
        Dnew$Deaths       <- 0
        Dnew$Sex          <- unique(DeathsSexYr$Sex)
        Dnew$PopName      <- unique(DeathsSexYr$PopName)
        
        rmID <- !(with(Dnew, paste(Year, Age, Lexis, sep = "-")) %in% 
                 with(DeathsSexYr, paste(Year, Age, Lexis, sep = "-")))
        Dnew <- Dnew[rmID, ]
        
        return(rbind(DeathsSexYr, Dnew))
      }
    }) # end internal function definition
  Deaths <- Deaths[Deaths$Age != "TOT", ]
  # remove UNK
  UNKTF <- any(Deaths$Age == "UNK")
  if (UNKTF){
    UNK           <- Deaths[Deaths$Age == "UNK", ]
    Deaths        <- Deaths[Deaths$Age != "UNK", ]
  }

# break apart, operate, stick back together
  Deaths <- do.call(rbind,
    lapply(
      split(Deaths, list(Deaths$Sex,Deaths$Year)), 
      d_long_SexYr
        )
      )
  
  if (UNKTF){
    Deaths <- rbind(Deaths, UNK)
  }
  Deaths   <- resortDeaths(Deaths)
  rownames(Deaths) <- NULL
  
  invisible(Deaths)
  
}

