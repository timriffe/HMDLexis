
#'
#' @title p_ey2ny Move Dec 31 to Jan 1 of following year
#' 
#' @description This function can be harmlessley called on each country.
#' 
#' @param Pop the standard internal HMD population data.frame, probably early in processing.
#' 
#' @return Pop the same data.frame that was given as an argument, with this paricular date change undertaken where necessary
#' 
#' @export
#' 


p_ey2ny <- function(Pop){
  
  DM                <- paste(Pop$Day, Pop$Month, sep = "-")
  eyInd             <- DM == "31-12"
  
  Pop$Day[eyInd]    <- 1
  Pop$Month[eyInd]  <- 1
  Pop$Year[eyInd]   <- Pop$Year[eyInd] + 1
  
  # untested! this might break it!
  if (any(eyInd)){
    Pop[eyInd, ] <- assignNoteCode(Pop[eyInd, ], "p_ey2ny()")
  }
  
  invisible(Pop)
}
