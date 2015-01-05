
#' @title KannistoSurvMin and KannistoSurv two functions for fitting and predicting the Kannisto survival function
#' 
#' @description These two functions are called by \code{d_KannistoSplitYrSex()}.
#' 
#' @param pars or \code{ab} a vector or list of the starting values for the a and b parameters of the Kannisto function. Default value of \code{c(a = .1, b = .1)}.
#' @param P the synthetic survival function, as passed in from \code{d_KannistoSplitYrSex()}.
#' @param oa the open age, as passed in from \code{d_KannistoSplitYrSex()}.
#' @param N the number of data points (same as \code{length(SyntheticSurvival)}).
#' @param x a vector of the ages of the data points used in \code{SyntheticSurvival}.
#'  
#' @importFrom compiler cmpfun
#' 
#' @return \code{KannistoSurvMin()} returns a residual to be minimized in \code{optim()} or similar. \code{KannistoSurv()} returns the predicted survival values (which may likely need to be rescaled later) for a given a,b and age vector.
#' 

KannistoSurvMin <- cmpfun(function(pars, P = SyntheticSurvival, oa = OA, N = N, x = Ages4Kannisto){
    Sx <- ((1 + pars["a"]) / (1 + pars["a"] * exp(pars["b"] * (x - oa + N)))) ^(1 / pars["b"])
    sum((log(Sx) - log(P))[P != 0]^2)
  })

KannistoSurv <- cmpfun(function(ab, x, oa = OA, N = 20){
    ((1 + ab["a"]) / (1 + ab["a"] * exp(ab["b"] * (x - oa + N)))) ^ (1 / ab["b"])
  })