

#' Logical utility functions
#'
#' @aliases logic HMDutils
#' 
#' @description These logical functions are like the usual ones, but \code{NA} values are treated as \code{FALSE} by default. This is not an exhaustive list, but these are the ones that speed our coding, and reduce code clutter.
#' 
#' @param x,y any two vector that can be logically compared.
#' @name HMDlogic
#' 
#' @examples
#' \dontrun{
#' c(1,2,NA,4,5) == c(1,NA,3,4,NA)
#' # compare
#' c(1,2,NA,4,5) %==% c(1,NA,3,4,NA)
#' }
NULL
#' 

#' @rdname HMDlogic
'%==%' <- function(x,y){
  x == y & !is.na(x) & !is.na(y)
}

#' @rdname HMDlogic
'%!=%' <- function(x,y){
  x != y & !is.na(x) & !is.na(y)
}

#' @rdname HMDlogic
'%>%' <- function(x,y){
  x > y & !is.na(x) & !is.na(y)
}

#' @rdname HMDlogic
'%<%' <- function(x,y){
  x < y & !is.na(x) & !is.na(y)
}

#' @rdname HMDlogic
'%>=%' <- function(x,y){
  x >= y & !is.na(x) & !is.na(y)
}

#' @rdname HMDlogic
'%<=%' <- function(x,y){
  x <= y & !is.na(x) & !is.na(y)
}















