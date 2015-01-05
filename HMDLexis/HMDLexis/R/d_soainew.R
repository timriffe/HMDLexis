
#'
#' @title d_s0ainew a function to redistribute deaths in the open age group over higher age triangles.
#' 
#' @description This function takes deaths in the open age group and redistributes them over higher ages according to the Kannisto survival function as fit to a synthetic pseudo-survival function composed of the backward cumulative sum of old-age deaths preceding the open age. The matlab implementation of this contained many aspects that were out-of-protocol, and these have been reproduced, optionally. One can also toggle this off and get Tim's interpretation of the MP. This function depends on several other functions, most notably \code{d_KannistoSplitYrSex()}, which does the whole process for a given year and sex. The present function is more of a wrappers, though some prep work is done here. Two other functions are important, \code{d_remfluct()} and \code{d_CohortAdjust()}, which outsource ugly aspects of the redistribution method. \code{d_remfluct()} contains some out-of-MP aspects that have not been re-interpreted, whereas \code{d_CohortAdjust()} will toggle between the matlab and Tim's proposal. Other minor dependent functions include \code{d_KannistoSurvMin()}, \code{d_KannistoSurv()} (stored inside the former), and \code{PC2AC()}. This function should be called AFTER deaths have been standardized to triangles, and 0s imputed where necessary below the open age group. The open age group itself can either be RR or VV. If the open age group is declared as TL or TU, these are overwritten to RR and VV, respectively, for the time being.
#' 
#' @param Deaths the standard internal Deaths data.frame object. It MUST include only TL and TU entries, except the open age, which is ideally either RR or VV. 
#' @param reproduce.matlab logical. Default = \code{TRUE}. Should the quirks of the matlab implementation be reproduced?
#' 
#' @importFrom reshape2 acast
#' 
#' @return Deaths the standard internal Deaths data.frame, with open age entries gone and replaced with death triangles up to age 130.
#' 
#' @export

d_soainew <- function(Deaths, reproduce.matlab = TRUE){
 
# first, check that only TL and TU are in Lexis (except open age, which can be VV or RR.
# UNK redistributed, single year single age.
  stopifnot(!any(c("RR","VV","VH","RV") %in% Deaths$Lexis[Deaths$AgeInterval != "+"], na.rm = TRUE))
  stopifnot(!any(Deaths$AgeIntervali > 1, na.rm = TRUE))
  stopifnot(!any(Deaths$YearInterval > 1, na.rm = TRUE))

# -------------------------------------------------------------------  
# some temporary hacks:
# TOT is never necessary:
  Deaths <- Deaths[Deaths$Age != "TOT", ]
  
# while we're at it, slice off UNK, rbind back on later:
  UNKTF <- any(Deaths$Age == "UNK")
  if (UNKTF){
    UNK           <- Deaths[Deaths$Age == "UNK", ]
    Deaths        <- Deaths[Deaths$Age != "UNK", ]
  }
  
# these terminal Lexis shapes are curious
  Deaths$Lexis[with(Deaths, Lexis == "TL" & AgeInterval == "+" & !is.na(AgeInterval))] <- "RR"
  Deaths$Lexis[with(Deaths, Lexis == "TU" & AgeInterval == "+" & !is.na(AgeInterval))] <- "VV"

  Deaths         <- resortDeaths(Deaths)
# -------------------------------------------------------------------
# strategy to avoid nasty triangle looping over APC surface is to reshape
# to VV. 
  OrigCols       <- colnames(Deaths) 
  
  # who cares about open cohorts anyway?
  Deaths         <- d_addCohortColumn(Deaths, TLTUonly = TRUE)
  # head(Deaths[Deaths$Year == 1980 & Deaths$Agei == 99,])
# Cohort in columns, Year in rows
  VVmasterM      <- acast(Deaths[Deaths$AgeInterval != "+" & Deaths$Sex == "m", ], 
                                   Year ~ Cohort, 
                                   sum, 
                                   value.var = "Deaths")
  VVmasterF      <- acast(Deaths[Deaths$AgeInterval != "+" & Deaths$Sex == "f", ], 
                                   Year ~ Cohort, 
                                   sum, 
                                   value.var = "Deaths")                               
# these are VV shape, but organized in an AC matrix, meaning that each VV shape overlaps 2 ages
# the reference age is always the bottom one, but bear in mind that these parallelograms peek
# into the next age up (but don't overlap). We do this because selection happens by cohort over
# a given range of ages. Exotic at first, but very efficient.
  ACVVmatrixF             <- PC2AC(VVmasterF, 1) # 1 or 2? grrr as long as we know what's being refered to
  ACVVmatrixM             <- PC2AC(VVmasterM, 1)
# to handle edge cases, we pad out another 70 cohorts to the left (i.e. to account for an
# open age as low as 60 in the first year of data, where it is NOT advisable to run this function
# anyway. No need to add cohorts to the right- or is there?
  minCohM                 <- min(as.integer(colnames(ACVVmatrixM)))
  CohortsLeftM            <- (minCohM-71):(minCohM-1)
  LeftChunkM              <- matrix(ncol = length(CohortsLeftM), 
                                    nrow = nrow(ACVVmatrixM),
                                    dimnames= list(rownames(ACVVmatrixM), CohortsLeftM))
  ACVVmatrixM             <- cbind(LeftChunkM, ACVVmatrixM)
  
  minCohF                 <- min(as.integer(colnames(ACVVmatrixF)))
  CohortsLeftF            <- (minCohF-71):(minCohF-1)
  LeftChunkF              <- matrix(ncol = length(CohortsLeftF), 
                                    nrow = nrow(ACVVmatrixF),
                                    dimnames= list(rownames(ACVVmatrixF), CohortsLeftF))
  ACVVmatrixF             <- cbind(LeftChunkF, ACVVmatrixF)
# Tadj would happen on VVmaster in this function, namely of the Year rows
# that may need to happen live in d_CohortAdjust, though, because the coef 
# will depend on where you're standing (the reference year is always a 1)
# d_CohortAdjust will need to figure out internally how to use VVmaster
# -------------------------------------------------------------------
# on Deaths, for each year-sex combo 
# DeathsYearSex <- split(Deaths, list(Deaths$Year, Deaths$Sex))[[1]]
# DeathsYearSex <-Deaths[Deaths$Year == 1980 & Deaths$Sex == "f",]
  Deaths <- do.call(rbind,
                lapply(
                   split(Deaths, list(Deaths$Year, Deaths$Sex)), 
                   function(DeathsYearSex, ACVVmatrixF., ACVVmatrixM., reproduce.matlab.){
                     Sexit             <- unique(DeathsYearSex$Sex)
                     if (Sexit == "f"){
                       ACVVmatrix <- ACVVmatrixF.
                     } else {
                       ACVVmatrix <- ACVVmatrixM.
                     }
                     # avoiding a lengthy anonymous function here :-)
                     # DeathsYearSex <- Deaths[with(Deaths, Sex == "m" & Year == 1971),]
                     # this function calls d_remfluct() and d_CohortAdjust()
                     d_KannistoSplitYrSex(DeathsYearSex, 
                                          N = 20, 
                                          ACVVmatrix = ACVVmatrix, 
                                          reproduce.matlab = reproduce.matlab.)
                                   
      }, ACVVmatrixF. = ACVVmatrixF, ACVVmatrixM. = ACVVmatrixM, reproduce.matlab. = reproduce.matlab))

  # rbind deaths by Year and Sex together. That'll all be in a massive list of data.frames...
  Deaths         <- Deaths[, OrigCols] # this removes Cohort column
  
  if (UNKTF){
    Deaths       <- rbind(Deaths, UNK)
  }
  rownames(Deaths) <- NULL
  invisible(resortDeaths(Deaths))
}

