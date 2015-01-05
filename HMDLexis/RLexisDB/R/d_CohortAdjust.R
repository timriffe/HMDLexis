#'
#' @title d_CohortAdjust a function to adjust the redistributed open age group deaths according to cohort size differences detected in earlier ages
#' 
#' @description This is described in Appendix C, section 2. We take a very different programming strategy here than that employed in the legacy matlab scripts.
#' 
#' @param Dnew the interim data.frame, as passed in by \code{d_KannistoSplitYrSex()}. It has the same columns as the standard Deaths object, and draft death counts for the open age group as split into triangles upt to age 130.
#' @param OA the open age for this particular year/sex
#' @param VVmaster a matrix as passed in from \code{d_soainew()} via \code{d_KannistoSplitYrSex()}. This is a PC matrix (VV data in period-cohort format).
#' @param OAL the open age Lexis shape, \code{"RR"} or \code{"VV"}. If the original data were \code{"TL"} or \code{"TU"}, these will have been handled earlier and converted appropriately.
#' @param DOA deaths in the open age group, as passed in by \code{d_KannistoSplitYrSex()}. Used for scaling.
#' 
#' @return Dnew Same as received, but with open age deaths adjusted for cohort size.
#' 
#' @importFrom compiler cmpfun
#' 
#' @details \code{PC2AC()} is also called in a sort of novel way. PC data (VV shape) are put into the AC shape for easy named reference.


d_CohortAdjust <- cmpfun(function(Dnew, OA, ACVVmatrix, OAL, DOA, reproduce.matlab = TRUE){
  origCols          <- colnames(Dnew)
# adjustment happens by cohort, so let's assign a cohort to each triangle
  Dnew$Cohort       <- Dnew$Year - Dnew$Agei
  Dnew$Cohort[Dnew$Lexis == "TU"] <- Dnew$Cohort[Dnew$Lexis == "TU"] - 1

# get that as a free-standing vector too
  Cohorts           <- Dnew$Cohort
  CohortsU          <- unique(Cohorts)
  Ncoh              <- length(CohortsU)
# now add 2 cohorts to left and right for full selection
   
# need to space this just right. My own take on this might be confused...
# Dima: this shifts the age scale. Makes no difference for fitting
  Bump <- ifelse(OAL == "VV",
            ifelse(reproduce.matlab, 0, 0),
            ifelse(reproduce.matlab, 0, 1))  
# experimental function, possibly ill-conceived.
# to be clear, these are VV (PC) deaths organized by age-cohort!
# the as.character() bit picks out the relevant ages
  ACVVmatrix           <- ACVVmatrix[as.character((OA - 5 - Bump):(OA - 1 - Bump)),]
  
# get first ratios
  D1                   <- ACVVmatrix[, as.character(CohortsU - 2)]     # left 2
  D2                   <- ACVVmatrix[, as.character(CohortsU - 1)]     # left 1
  D3                   <- ACVVmatrix[, as.character(CohortsU + 1)]     # right 1
  D4                   <- ACVVmatrix[, as.character(CohortsU + 2)]     # right 2

  
  E1                   <- is.na(D1) # these will help us take an avg
  E2                   <- is.na(D2)
  E3                   <- is.na(D3)
  E4                   <- is.na(D4)

# this is the denominator for the avg
  DenomMostly4         <- (!E1) + (!E2) + (!E3) + (!E4)

# now impute 0s for NAs for numerator of avg
  D1[E1]               <- 0 
  D2[E2]               <- 0 
  D3[E3]               <- 0 
  D4[E4]               <- 0 
  
  Numerator            <- D1 + D2 + D3 + D4
  
  Rmat                 <- ACVVmatrix[, as.character(CohortsU)] / (Numerator / DenomMostly4)
  Rmat[is.na(Rmat)]    <- 1
  
  rvec                 <- colMeans(Rmat, na.rm = TRUE)
  
  # matlab does the same thing, but inclusive of the cohort in question
    if (reproduce.matlab){
      Dah       <- ACVVmatrix[, as.character(CohortsU)]
      Eah       <- is.na(Dah)
      Denomah   <- DenomMostly4 + (!Eah)
      Dah[Eah]  <- 0
      Numerah   <- Numerator + Dah
      Rmatah    <- ACVVmatrix[, as.character(CohortsU)] / (Numerah / Denomah)
      Rmatah[is.na(Rmatah)]    <- 1
      rvec      <- colMeans(Rmatah, na.rm = TRUE)
    }
# these are the known adjustment ratios...
# could also impute NAs with 1, would give different result, but same idea

  # sum(DeathsScaled1)
  DeathsScaled1        <- rvec[as.character(Cohorts)] * Dnew$Deaths
  Dnew$Deaths          <- DOA * (DeathsScaled1 / sum(DeathsScaled1))
  
  Dnew[, origCols]
})
# round(Dnew$Deaths ,4)