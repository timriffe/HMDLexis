
#'
#' @title p_long adds 0s up to age 130 for years of populations with no open age
#' 
#' @description Intermediate gaps are also filled in, if they exist. This functions returns the object unmodified and with warning if there are still open intervals.
#' 
#' @param Pop the standard internal Population \texttt{data.frame}, probably late in its processing, but not necessarily. Intermediate gaps are also filled in. Year-Sex combos with open ages are left untouched.
#' 
#' @return Pop the same deaths data.frame, but expanded to include ages up to 130 except for Year-Sex combos with pen ages, which are untouched. 
#' 
#' @export
#' 

p_long <- function(Pop, OPENAGE = 130){
  
  # an internal function
  p_long_SexYr <- function(PopSexYr){
    if (any(PopSexYr$AgeInterval == "+", na.rm = TRUE)){
      warning("Pop still has open intervals!")
      return(PopSexYr)
    } else {
      
      # insert
      Ages    <- 0:OPENAGE
     
      Pnew    <- as.data.frame(matrix(nrow = length(Ages), 
                      ncol = ncol(PopSexYr),
                      dimnames = list(NULL, colnames(PopSexYr))))
      
      Pnew$Year         <- unique(PopSexYr$Year)   
      Pnew$Age          <- Pnew$Agei <- Ages       
      Pnew$LDB          <- Pnew$AgeIntervali <- Pnew$AgeInterval  <-  Pnew$YearInterval  <-  1
      Pnew$Access       <- "O"
      # TR: this area assignment seems OK to me. Safe to assume 
      # we'd never have more than one NoteCode for a given year & Sex?
      Pnew$Area         <- unique(PopSexYr$Area)
      Pnew              <- assignNoteCode(Pnew, "p_long")
      Pnew$Deaths       <- 0
      Pnew$Sex          <- unique(PopSexYr$Sex)
      Pnew$PopName      <- unique(PopSexYr$PopName)
      
      rmID <- !(with(Pnew, paste(Year, Age, sep = "-")) %in% 
          with(PopSexYr, paste(Year, Age, sep = "-")))
      Pnew <- Pnew[rmID, ]
      
      Pout <- resortPops(rbind(PopSexYr, Pnew))
      invisible(Pout)
    }
  } # end internal function definition
  Pop <- Pop[Pop$Age != "TOT", ]
  # remove UNK
  UNKTF <- any(Pop$Age == "UNK")
  if (UNKTF){
    UNK           <- Pop[Pop$Age == "UNK", ]
    Pop        <- Pop[Pop$Age != "UNK", ]
  }
  
# break apart, operate, stick back together. Would be easier and faster with data.table
  Pop <- do.call(rbind,
              lapply(
                split(Pop, list(Pop$Sex,Pop$Year)), 
                p_long_SexYr
              )
            )
  

  if (UNKTF){
    Pop <- rbind(Pop, UNK)
  }
  Pop   <- resortPops(Pop)
  rownames(Pop) <- NULL
  
  invisible(Pop)
}


