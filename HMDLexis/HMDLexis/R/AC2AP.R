#' @title \code{AC2AP} convert an age-cohort matrix to an age-period matrix
#'
#' @description
#' \code{AC2AP} is an auxiliary function called by \code{Exposures_per()} (and likely by LexisDB programs, once these are written). It is also useful for plain vanilla demography in a live R session...
#'
#' @param ACmatrix a matrix (age by cohort) of quantities (Deaths, Births, Exposures) in triangles, either upper or lower. The matrix must either represent lower triangles or upper triangles, not both. Whether triangles are upper or lower is indicated with the other argument, \code{Lexis}. This function can also be used to shift population counts from cohort to period matrices, as long as the population is January 1st or December 31st. Matrices must in any case have row names and column names, where rows are ages (integers, no '+'!) and columns are cohorts (also simple integers).
#' @param Lexis either 1 or 2. Per HMD internal convention, 1 indicates lower triangles and 2 indicates upper triangles. For population counts, 1 indicates January 1st and 2 indicates December 31st. Default is 2, so watch out- if you choose wrong, then column labels in the AP matrix output will be off by 1. If you're not sure, draw a diagram :-).
#' 
#' @return a \code{APmatrix} a matrix in age-period format, with correctly labeled dimensions. It will also have an attribute \code{Lexis}, used for automatic Lexis argument detection in the case that the matrix is placed back into \code{AP2AC()}.
#' 
#' @details This is a convenience function and it might not do what you think. Best to do trial and error to get things lined up right.
#' 
#' @author Tim Riffe \email{triffe@@demog.berkeley.edu}
#' 
#' @importFrom reshape2 melt
#' @importFrom reshape2 acast
#' @importFrom compiler cmpfun
#' 
#' @export
# Author: triffe
###############################################################################
AC2AP <- cmpfun(function(ACmatrix, Lexis){
  if (missing(Lexis)){
    Lexis         <- ifelse(is.null(attr(ACmatrix, "Lexis")), 2, attr(ACmatrix, "Lexis"))
  }
  # back to LDB format, but lacking Year
  longform        <- melt(ACmatrix, varnames = c("Age", "Cohort"), value.name = "value")
  # assume Year is t+x+1. where 1 is the year.offset. Could also be 0. Test if not sure
  longform$Year   <- longform$Cohort + longform$Age + ifelse(Lexis == 2, 1, 0)
  longform        <- longform[!is.na(longform$value), ]
  # cast back to Age x Year matrix
  APmatrix        <- acast(longform, Age ~ Year, value.var = "value")
  attr(APmatrix, "Lexis") <- Lexis
  APmatrix
})
