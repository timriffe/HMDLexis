#' @title \code{PC2AC} convert a period-cohort (VV) matrix to an age-cohort matrix, wherein entries still represent the VV shape.
#'
#' @description
#' \code{PC2AC} is an auxiliary function available to call by any of the HMD function family.
#' 
#' @param PCmatrix a matrix (period by cohort) of quantities in VV shape. May work for other shapes too, but you'll need to experiment with Lexis = 1 or 2.
#' @param Lexis either 1 or 2. Per HMD internal convention, 1 indicates lower triangles and 2 indicates upper triangles. For population counts, 1 indicates January 1st and 2 indicates December 31st. Default is 2, so watch out- if you choose wrong, then column labels in the AP matrix output will be off by 1. If you're not sure, draw a diagram. VV gets a 2 by default!
#' 
#' @return a \code{ACmatrix} a matrix in age-cohort format, with correctly labeled dimensions. It will also have an attribute \code{Lexis}, used for automatic Lexis argument detection in the case that the matrix is placed back into \code{AC2AP()}. Note that if the original contents were the VV shape that age is the lower left corner, and that the shape spans two ages. Draw a picture if necessary.
#' 
#' @author Tim Riffe \email{triffe@@demog.berkeley.edu}
#' 
#' @export
#' 
#' @importFrom reshape2 melt
#' @importFrom reshape2 acast
#' @importFrom compiler cmpfun
#' 
PC2AC <- cmpfun(function(PCmatrix, Lexis){
  if (missing(Lexis)){
    Lexis         <- ifelse(is.null(attr(PCmatrix, "Lexis")), 1, attr(PCmatrix, "Lexis"))
  }
  longform        <- melt(PCmatrix, varnames = c("Year", "Cohort"), value.name = "value")
  head(longform)
  # assume Year is t+x+1. where 1 is the year.offset. Could also be 0. Test if not sure
  longform$Age    <- longform$Year - longform$Cohort + ifelse(Lexis == 2, -1, 0)
  longform        <- longform[!is.na(longform$value) & longform$Age > 0, ]
  # cast back to Age x Year matrix
  ACmatrix        <- acast(longform, Age ~ Cohort, value.var = "value")
  attr(ACmatrix, "Lexis") <- Lexis
  ACmatrix
})
