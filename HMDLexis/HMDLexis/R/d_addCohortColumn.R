#'
#' @title d_addCohortColumn adds Cohort column to Deaths object
#' 
#' @description This will assign an integer cohort to TL, TU, VV, VH and RV entries. Optionally, it will only do so for triangles. RR and NA Lexis shapes get NA for cohort. 
#' 
#' @param Deaths the standard Deaths object, at any stage of production
#' @param TLTUonly logical, Default = \code{TRUE}. Should Cohorts only be given to triangles?
#' 
#' @importFrom compiler cmpfun
#' 

d_addCohortColumn <- function(Deaths, TLTUonly = TRUE){
    
    if (!"Cohort" %in% colnames(Deaths)){
      # this takes care of TL and VH
      Deaths$Cohort         <- Deaths$Year - Deaths$Agei
      
      NAi                   <- is.na(Deaths$Lexis)
      # this takes care of TU and VV 
      # testing use of %==% operator. See if stuff breaks
      TUVVi                 <- (Deaths$Lexis %==% "TU" | Deaths$Lexis %==% "VV") 
      Deaths$Cohort[TUVVi]  <- Deaths$Cohort[TUVVi] - 1
      
      # RV only for infants anyway, the rest are VV
      RVi                   <- Deaths$Lexis %==% "RV"
      Deaths$Cohort[RVi]    <- Deaths$Cohort[RVi] - 2
      
      # never assign Cohort to RR
      Deaths$Cohort[Deaths$Lexis %==% "RR"]      <- NA
      Deaths$Cohort[is.na(Deaths$Lexis)]         <- NA
      
      # in case we strictly want triangles
      if (TLTUonly){
        TLTUi <- !Deaths$Lexis %in% c("TL","TU") & !NAi
        Deaths$Cohort[TLTUi] <- NA
      }
    }
    invisible(Deaths)
}
