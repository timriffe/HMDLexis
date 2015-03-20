#'
#' @title p_ecm a simple extinct cohort method function, for modularity
#' 
#' @description This implementation seeks to follow the MP. Deaths must be finalized. Populations must be in single ages and all years. This function covers area B from Figure 6 'Methods used for population estimates'. The actual work is odone by \code{p_ecm_inner()}
#' 
#' @param Pop The standard internal Population data.frame, as
#' @param Deaths after all processing is done. Completed triangles.
#' @param a lowest age for which EC estimates should be given. Default 80.
#' @param reproduce.matlab logical. default is \code{FALSE}. This affects only the border cohort between EC and SR.
#' 
#' @return Population with old ages either filled in or imputed using the extinct cohort method
#' 
#' @details Conceivably, a population could have a long series with low open ages in the early part, but decent population data in the recent part. If we want a value of \code{A} lower than 80 for only part of the series, this can be achieved by subsetting \code{Dsex}, and making two calls to the function, then \code{rbind()}ing back together.
#' 
#' @export
#' 

p_ecm <- function(Pop, Deaths, a = 80, reproduce.matlab = FALSE){
 
  ColnamesKeep <- colnames(Pop)
  
  Pop <- Pop[Pop$Age != "TOT", ]
  
  # add Cohort Column to Deaths
  # MP uses omega, but we use the most recent extinct cohort, since
  # we mark Cohorts and can use them to select.
  Deaths        <- d_addCohortColumn(Deaths)
  Pop           <- p_addCohortColumn(Pop)
  
  # slice off UNK, rbind back on later:
  UNKi          <- Pop$Age == "UNK"
  UNKTF         <- any(UNKi)
  if (UNKTF){
    UNK           <- Pop[UNKi, ]
    Pop           <- Pop[!UNKi, ]
  }
  
  # Sex <- "f"; Sex <- "m"
  PopMF <- list()
  for (Sex in c("m","f")){
    Dsex              <- Deaths[Deaths$Sex == Sex, ]
    Psex              <- Pop[Pop$Sex == Sex, ]
    # 1) determine which cohorts are extinct. 
    #For now assume that this is separate for each sex
    omega             <- p_ecm_findOmega(Dsex, l = 5, threshold = 0.5)
    # p_ecm_inner() defined below, in same script. not quirky, but it is useful
    # to have these steps be modular, so that it can be called elsewhere
    ECpop <- p_ecm_inner(Dsex = Dsex, a = a, omega = omega, reproduce.matlab = reproduce.matlab)
    # Dima: should only be done for srecm, but not for ec.
    #    # an apparent quirk in the matlab code
    #    if (reproduce.matlab){
    #      ECpop <- rbind(ECpop, ECpop[nrow(ECpop), ])
    #      ECpop[nrow(ECpop), "Population"] <- 0
    #      ECpop[nrow(ECpop), "Age"]        <- ECpop[nrow(ECpop), "Agei"] <- omega["omega"]
    #      ECpop[nrow(ECpop), "Year"]       <- max(Psex$Year)
    #      ECpop$Population[ECpop$Cohort == omega["Cohmax"]] <- 
    #        ECpop$Population[ECpop$Cohort == omega["Cohmax"]] + omega["Ds"]
    #    }
    # append
#  LexisMap(acast(Psex, Agei~Year, value.var = "Population"),log=FALSE)
#  dev.new()
#  LexisMap(acast(ECpop, Agei~Year, value.var = "Population"),log=FALSE)
    # is.na(Cohort) picks out the open age groups, which we hope are all 80+

  Pold <- Psex[!with(Psex, (Year - Agei - 1) <= omega["Cohmax"] & Agei >= a), ColnamesKeep]
#  dev.new()
#  LexisMap(acast(Pold, Agei~Year, value.var = "Population"),log=FALSE)
    PopMF[[Sex]]      <- rbind(Pold, ECpop[, ColnamesKeep])
  }
  if (UNKTF){
    PopMF[["UNK"]]    <- UNK
  }
  # stick together males and females, resort and return
  Pop                 <- resortPops(do.call(rbind, PopMF))
  rownames(Pop)       <- NULL
  invisible(Pop)
}

#'
#' @title p_ecm_inner does the work of \code{p_ecm()} 
#' 
#' @description This implementation seeks to follow the MP, we have outsourced the work of \code{p_ecm()} to here for the sake of modularity. This function may also be called by the \code{p_ic()} function ecosystem in order to extend population counts in census 1 and census 2 prior to interpolating. Deaths must be finalized. This function covers area B from Figure 6 'Methods used for population estimates'. 
#' 
#' @param Dsex Deaths after all processing is done, subset of a single sex.
#' @param a lowest age for which EC estimates should be given. Default 80.
#' @param omega is the object returned by \code{p_ecm_findOmega()}.
#' @param reproduce.matlab logical. default is \code{FALSE}. This affects only the border cohort between EC and SR.
#' 
#' @return Population with old ages either filled in or imputed using the extinct cohort method
#' 
#' @details Conceivably, a population could have a long series with low open ages in the early part, but decent population data in the recent part. If we want a value of \code{A} lower than 80 for only part of the series, this can be achieved by subsetting \code{Dsex}, and making two calls to the function, then \code{rbind()}ing back together.
#' 
#' @importFrom reshape2 acast
#' 
#' @export
#' 

p_ecm_inner <- function(Dsex, a = 80, omega = NULL, reproduce.matlab = FALSE){
  # instead of feeding in Pop argument, just assume standard R LexisDB format of Pop DF.
  # this is what we do anyway by filling out columns below.
  Pop <- structure(list(PopName = character(0), Area = integer(0), Sex = character(0), 
      Age = character(0), AgeInterval = character(0), Type = character(0), 
      Day = integer(0), Month = integer(0), Year = integer(0), 
      RefCode = integer(0), Access = character(0), Population = numeric(0), 
      NoteCode1 = integer(0), NoteCode2 = integer(0), NoteCode3 = integer(0), 
      LDB = integer(0), Agei = integer(0), AgeIntervali = integer(0)), .Names = c("PopName", 
      "Area", "Sex", "Age", "AgeInterval", "Type", "Day", "Month", 
      "Year", "RefCode", "Access", "Population", "NoteCode1", "NoteCode2", 
      "NoteCode3", "LDB", "Agei", "AgeIntervali"), row.names = integer(0), class = "data.frame")
  
  if (!"Cohort" %in% colnames(Deaths)){
    Deaths            <- d_addCohortColumn(Deaths)
  }
  if (is.null(omega)){
    omega             <- p_ecm_findOmega(Dsex, l = 5, threshold = 0.5)
  }
  
  # 2) reshape deaths by period-cohort (PC = VV shape)
  VV                <- acast(Dsex[with(Dsex, Agei >= a & Cohort <= omega["Cohmax"]), ], 
                             Year ~ Cohort, sum, value.var = "Deaths")
  VV                <- rbind(VV, 0)   
  rownames(VV)[nrow(VV)] <- max(Dsex$Year) + 1
  if (reproduce.matlab){
    VV[nrow(VV),ncol(VV)] <- VV[nrow(VV),ncol(VV)] + omega["Ds"]
  }
  cohorts           <- as.integer(colnames(VV))
  years             <- as.integer(rownames(VV))
  
  # all possible ages given these years and cohorts (same dim as VV)
  AllAges           <- outer(years, cohorts, "-") - 1
  # years repeated, same dim as VV, for selection
  Allyears          <- replicate(length(cohorts), years)
  Allcohorts        <- t(replicate(length(years), cohorts))
  # the selection mask
  Mask              <- AllAges >= a & AllAges <= 130
  
  # 3) Equation 38
  ECpops            <- apply(VV[nrow(VV):1, ], 2, cumsum)[nrow(VV):1, ][Mask]
  
  # now select the corresponding ages and years
  ECAges            <- AllAges[Mask]
  ECYears           <- Allyears[Mask]
  ECCohorts         <- Allcohorts[Mask]
  
  # stick together and add remaining columns 
  ECpop             <- as.data.frame(
                          matrix(nrow=length(ECAges), 
                                 ncol = ncol(Pop), 
                                 dimnames = list(NULL, colnames(Pop))))
  ECpop$Age         <- ECpop$Agei            <- ECAges
  ECpop$Year        <- ECYears
  ECpop$Population  <- ECpops
  ECpop$Sex         <- unique(Dsex$Sex)
  ECpop$AgeInterval <- ECpop$AgeIntervali    <- 1
  ECpop$PopName     <- unique(Dsex$PopName)
  ECpop$Day         <- 1
  ECpop$Month       <- 1
  ECpop$LDB         <- 1
  ECpop$Cohort      <- ECCohorts
  ECpop$Access      <- "O" # presumably any data we invent is open access, even if the origin data are not?
  ECpop$NoteCode1   <- "p_ecm()" 
  ECpop[,colnames(Pop)]
}

