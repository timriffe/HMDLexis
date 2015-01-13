#'
#' @title p_addCohortColumn adds Cohort column to Population data.frame
#' 
#' @description The cohort column will only be added to those entries where the age interval is equal to 1. Sometimes it will be easier to match populations using cohorts for certain LDB operations. This function is an auxiliary function available to all others. We assume Populations are January 1 for this function.
#' 
#' @param Population the standard internal Population object, at any point in processing
#' 
#' @return Population same data.frame with an extra column. 
#' 
#' @importFrom compiler cmpfun

p_addCohortColumn <- cmpfun(function(Pop){
    if (! "Cohort" %in% colnames(Pop)){
      # only valid for single-year age groups
      Interval1i        <- Pop$AgeInterval == "1" & !is.na( Pop$AgeIntervali)
      
      # this may overwrite
      Pop$Cohort <- NA
      # this takes care of TL and VH
      Pop$Cohort[Interval1i]   <- Pop$Year[Interval1i] - Pop$Agei[Interval1i] - 1
    }
    invisible(Pop)
  })

