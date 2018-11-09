
# This R file contains the whole d_soainew() function family, which is admittedly overly complex and ad-hockish in nature.
# The current implementation only manages to approximate the matlab one, and apologies are due for
# an occassionally opaque and idiosyncratic R-implementation

# Contents:
# d_soainew()                top level wrapper
# d_KannistoSplitYrSex()     does the work for a single year and sex
# KannistoSurvMin()          Kannisto survival residual function
# KannistoSurv()             Kannisto survival given a, b and ages
# d_remfluct()               This one is ad-hickish, and a likely point 
#                            of difference with matlab
# d_CohortAdjust             Also not currently verified to match matlab output

# TR: in my opinion this family of functions should be revisited to simplify the 
#     method, and also to ensure that the implementation and MP description coincide.
#

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

d_soainew <- function(Deaths, sa=NULL, empy=NULL, reproduce.matlab = TRUE){
 
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
  ACVVmatrixF             <- PC2AC(VVmasterF, 1) # 1 or 2? grrr as long as we know what's being referred to
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
  f_split_deaths_call <-   function(DeathsYearSex, ACVVmatrixF., ACVVmatrixM., reproduce.matlab.){
    Sexit             <- unique(DeathsYearSex$Sex)
    if (Sexit == "f"){
      ACVVmatrix <- ACVVmatrixF.
    } else {
      ACVVmatrix <- ACVVmatrixM.
    }
    # avoiding a lengthy anonymous function here :-)
    
    # this function calls d_remfluct() and d_CohortAdjust()
    d_KannistoSplitYrSex(DeathsYearSex, 
                         N = 20, 
                         ACVVmatrix = ACVVmatrix, 
                         reproduce.matlab = reproduce.matlab.)
    
  }
  
  Deaths.split <- split(Deaths, list(Deaths$Year, Deaths$Sex))
  Deaths.lapply <- lapply(Deaths.split, f_split_deaths_call,
                          ACVVmatrixF. = ACVVmatrixF, ACVVmatrixM. = ACVVmatrixM, reproduce.matlab. = reproduce.matlab)
                        
  Deaths <- do.call(rbind, Deaths.lapply)
                

  # rbind deaths by Year and Sex together. That'll all be in a massive list of data.frames...
  Deaths         <- Deaths[, OrigCols] # this removes Cohort column
  
  if (UNKTF){
    Deaths       <- rbind(Deaths, UNK)
  }
  rownames(Deaths) <- NULL
  invisible(resortDeaths(Deaths))
}


#'
#' @title d_KannistoSplitYrSex() a function that redistributes the open age group for a particular Year-Sex.
#' 
#' @description This function actually does most of the work for \code{d_soainew()}, and it only exists as a separate function to avoid making a long anonymous function. The central construct inside \code{d_soainew()} will be the same, a long list of Year-Sex data.frames that later need to be \code{do.call("rbind",BigList)}ed together. The function outsources the ugly aspects of the open-age redistribution to just focus on the actual Kannisto fitting. The removal of unusual fluctuations happens in \code{d_remfluct()}, and the adjustment for cohort size happens in \code{d_CohortAdjust()}.
#' 
#' @param DeathsYearSex A Year-Sex chunk of the standard Deaths object, removed of UNK and TOT, TLTU only.
#' @param N the number of year data points to use when fitting the Kannisto survival function
#' @param VVmaster a matrix as passed in from \code{d_soainew()}. This is a PC matrix (VV data in period-cohort format). Passed onward to \code{d_CohortAdjust()} and not used directly in this function.
#' 
#' @return The same Year-Sex data.frame chunk, but expanded out to age 130 IFF there was an open age group.
#' 
#' @details The function calls \code{d_remfluct()} and \code{\d_CohortAdjust()}, which also depends on the \code{PC2AC()} utility function, as well as \code{d_KannistoSurvMin()} and \code{d_KannistoSurv()}, which are used directly.

d_KannistoSplitYrSex <- function(DeathsYearSex, N = 20, ACVVmatrix, reproduce.matlab){
  
  # return if no open age group (padded with 0s elsewhere)
  if (!any(is.na(DeathsYearSex$AgeIntervali))){
    return(DeathsYearSex)
  }
  
  ## CB: TODO:  robustify so that when Lexis is NA, the right thing happens or at least
  ## an error occurs with a proper messages.  Add check for multiple open age intervals
  ## in a given (Year,Sex) split.
  
  
  # if the open age group is a VV, then we need to 
  # 1) add the TL below it to it
  # 2) save the value so we can put it back later
  # open age parameters
  DOA               <- DeathsYearSex$Deaths[DeathsYearSex$AgeInterval == "+"]
  OA                <- DeathsYearSex$Agei[DeathsYearSex$AgeInterval == "+"]
  OAL               <- DeathsYearSex$Lexis[DeathsYearSex$AgeInterval == "+"]
  # this is that little triangle
  OATL              <- ifelse(OAL == "VV",DeathsYearSex$Deaths[with(DeathsYearSex, Age == OA & Lexis == "TL")], 0)
  # we don't redistribute this triangle, but we need it for this next step, which is in RR format.
  
  # get deaths back to RR temporarily, no other shapes or intervals
  deathsRR          <- DeathsYearSex[with(DeathsYearSex, Lexis == "TL" & !is.na(Lexis) & Agei < OA),]
  deathsRR$Deaths   <- deathsRR$Deaths + DeathsYearSex$Deaths[with(DeathsYearSex, Lexis == "TU" & !is.na(Lexis) & Agei < OA)]
  deathsRR$Lexis    <- "RR" # just for consistency, no big deal
  
  # note, we've removed the open age group itself for this interim object
  # this returns the last N years of deaths, with potential fluctuations removed, sort of.
  # if no fluctuations were detected according to that quirky function, we still get back the 
  # right object. Note, it's only ages and deaths and another indicator column to note which
  # ones were imputed. 
  
  dths.fit       <- d_remfluct(deathsRR, N = 25, spar1 = 0.9, spar2 = 0.3)
  
  # now do a reverse cumsum to get a synthetic survival on the last N-1 (19) years of it:
  # [this N is different from the one used as an argument in d_remfluct()...]
  
  # - N + 2 gives 19, whereas 1 gives 20, for a grand total of 21 data points.
  
  # Dima: probably a programming convenience in matlab: we take 21 points, including
  # the open age group, but of course ignore the opena ge group for fitting, which
  # leaves us with the required 20 points.
  ml21              <- ifelse(reproduce.matlab, 0, 1)
  
  dths.fit          <- dths.fit[(nrow(dths.fit) - N + 1 + ml21):nrow(dths.fit), ]
  
  # now the vectors are 20 long
  Deaths4Kannisto   <- c(dths.fit$DeathsForFitting, DOA + OATL)
  Ages4Kannisto     <- c(dths.fit$Age, OA)
  
  # get synthetic survival, where age OA-N = 1
  TOT               <- sum(Deaths4Kannisto)
  SyntheticSurvival <- rev(cumsum(rev(Deaths4Kannisto))) / TOT 
  # plot(SyntheticSurvival)
  # a = .1 and b = .1 are the starting values
  # KannistoSurvMin() is the OLS function for Kannisto survival ?optim
  ab                <- try(optim(c(a = .1, b = .1), 
      fn = KannistoSurvMin, 
      P = SyntheticSurvival, 
      oa = OA, 
      N = N, 
      x = Ages4Kannisto)$par)
  if (class(ab) == "try-error"){
    ab             <- optim(c(a = .1, b = .1), 
      fn = KannistoSurvMin, 
      method = "BFGS",
      P = SyntheticSurvival, 
      oa = OA, 
      N = N, 
      x = Ages4Kannisto)$par
  }
  
  # and this is just the formula itself for Kannisto survival, predicting out the tail
  KannSurvTail      <- KannistoSurv(ab = ab, x = seq(OA, 131, by = .5), oa = OA, N = N)
  #seq(OA, 131, by = .5) - OA + N
  # we don't want to impute the known TL, if it existed (taken care of?)
  if (OAL == "VV"){
    KannSurvTail <- KannSurvTail[-1]
  }
  #ddr <- (DOA / KannistoSurv(ab = ab, x = OA, oa = OA, N = N)) * -diff(KannSurvTail) 
  
  
  # get implied deaths distribution
  diffs             <- -diff(KannSurvTail) 
  
  # scale to sum to open age deaths
  Distr.draft       <- (diffs / sum(diffs)) * DOA 
  
  # per MP, remove counts less than .25, then rescale
  Distr.draft[Distr.draft < .25] <- 0
  Distr.draft       <- Distr.draft / sum(Distr.draft) * DOA
  
  # now expand back out
  Lexis             <- rep(c("TL", "TU"), ceiling(length(Distr.draft) / 2))
  Ages              <- rep(OA:130, each = 2)
  
  # again, account for possible known TL below open age
  if (OAL == "VV"){
    Lexis    <- Lexis[-1]
    Ages     <- Ages[-1]
  }
  
  # slap the newly distributed open age group into a standard deaths object
  Dnew              <- as.data.frame(
    matrix(nrow = length(Ages), ncol = ncol(DeathsYearSex), 
      dimnames = list(NULL, 
        colnames(DeathsYearSex)))
  )
  Dnew$PopName      <- unique(DeathsYearSex$PopName)
  Dnew$Area         <- unique(DeathsYearSex$Area)
  Dnew$YearInterval <- unique(DeathsYearSex$YearInterval)
  Dnew$Sex          <- unique(DeathsYearSex$Sex)
  Dnew$Year         <- unique(DeathsYearSex$Year)
  Dnew$Lexis        <- Lexis
  Dnew$Age          <- Dnew$Agei         <- Ages
  Dnew$Deaths       <- Distr.draft
  Dnew$AgeInterval  <- Dnew$AgeIntervali <- 1
  
  # The adjusts for relative cohort size using VVmaster, as passed in from d_soainew()
  Dnew              <- d_CohortAdjust(Dnew = Dnew, 
    OA = OA, 
    ACVVmatrix = ACVVmatrix, 
    OAL = OAL, 
    DOA = DOA, 
    reproduce.matlab = reproduce.matlab)
  
  # this can only happen if DOA is 0, but we'd get NaNs in Deaths due to 0 
  # in denom in various locations. Just overwrite
  Dnew$Deaths[is.nan(Dnew$Deaths)] <- 0
  Dnew              <- assignNoteCode(Dnew, "d_soainew()")
  DeathsYearSex     <- rbind(DeathsYearSex[DeathsYearSex$AgeInterval != "+", ], Dnew)
  
  invisible(resortDeaths(DeathsYearSex))
}


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

KannistoSurvMin <- function(pars, P = SyntheticSurvival, oa = OA, N = N, x = Ages4Kannisto){
    Sx <- ((1 + pars["a"]) / (1 + pars["a"] * exp(pars["b"] * (x - oa + N)))) ^(1 / pars["b"])
    sum((log(Sx) - log(P))[P != 0]^2)
  }

KannistoSurv <- function(ab, x, oa = OA, N = 20){
    ((1 + ab["a"]) / (1 + ab["a"] * exp(ab["b"] * (x - oa + N)))) ^ (1 / ab["b"])
  }


#'
#' @title d_remfluct a function called by \code{d_soainew()} to remove outliers prior to fitting the Kannisto survival function
#' 
#' @description In order to redistribute deaths in the open age group, one first fits the Kannisto survival function to an imaginary survival function built backward from death counts, similar to the extinct cohort method. Prior to backward summing the synthetic survival curve, we look at deaths to check for a particular variety of unusual fluctuation: that caused by death dearths, as happens when small cohorts pass through old ages. If such a fluctuation is found, deaths are imputed into the affected ages using a cubic spline. 
#' 
#' @details Some comments: The matlab implementation of this whole process is not in sync with MPv5, but attempt here to reproduce it. The following quirks have been identified at this time: 1) the function only removes death troughs, but not spikes, 2) fluctuations of longer than 7 years or of only one year are not treated, 3) 25 years are taken as reference rather than the 30 stated in the MP, 4) the first year in the window cannot be an outlier, 5) only smoothed values that differ from the original values by more than 10\% are kept. These rules may have been reasonable for the three test-countries used, but there is no reason to think that they are valid for the whole HMD. Further diagnostics are not produced on a regular basis, and should be part of the country specialist checklist. The particular spline used in matlab was \code{csaps()}, and the values of the smoothing parameter given in the MP, 0.0005, applies to that particular function. One can replicate \code{csaps()} in R by using \code{pspline::smooth.Pspline()}, with smoothing parameter \code{spar} set to $(1-p)/p$. One can come very close to replication by using the base routine \code{smooth.spline()} with \code{spar} set to 0.9 and 0.3.
#' 
#' @param deathsRR the universal deaths object AFTER all ages below the open age group have been split exclusively to triangles.
#' @param N the number of years to take when checking first differences for smoothness. The MP states 30, but the matlab used 25.
#' @param spar1 the smoothing parameter for the first cubic spline over first differences. In the MP and matlab as 0.0005.
#' @param spar2 the smoothing parameter for the second cubic spline over all ages of death counts, as seen in matlab as 0.9.
#'  
#' @importFrom compiler cmpfun
#' 
#' @return deaths.out a simple data.frame with 3 columns: Age, DeathsForFitting, and Original. These are just for the purpose of fitting the Kannisto survival function, so everything is in RR format. Ages are from OA up to OA-N to OA, Deaths are RR, with detected gaps imputed with smooth counts, and Original is an unused indicator of which counts were kept vs imputed. 
#' 

d_remfluct <- function(deathsRR, N = 25, spar1 = 0.9, spar2 = 0.3){
  # deathsRR are RR deaths including all columns, from a particular year and sex
   
  # already removed of open age, no need here
  NN                 <- nrow(deathsRR)
  DTHS               <- deathsRR$Deaths
  dths               <- DTHS[(NN - N):NN]
  Diffs              <- diff(dths)
  ages               <- deathsRR$Agei[(NN - N):NN ]
  agesd              <- ages[2:length(ages)]
  # this arg combo invokes the same spline procedure as the matlab function csaps() 
   ## fit.Diffs          <- c(predict(smooth.Pspline(
   ##                        x = agesd,
   ##                        y = Diffs, 
   ##                        norder = 2, 
   ##                        method = 1,
   ##                        spar = (1 - p2) / p2     # p given in MP and matlab as 0.0005
   ##                       ),agesd))

  ## CB: replace smooth.Pspline with smooth.spline(), which 
  ## is more robust and reduces a package dependency
  
  cat(paste("\n Year, Sex =", unique(deathsRR$Year), unique(deathsRR$Sex) ))
  
  fit.Diffs          <- predict(smooth.spline(
                          x = agesd,
                          y = Diffs, 
                          all.knots=TRUE,
                          spar = spar1     # p given in MP and matlab as 0.0005 for use with csaps()
                         ), agesd)$y     
  Errors             <- fit.Diffs - Diffs
  sigma              <- sd(Errors)
  upper              <- fit.Diffs + 1.8 * sigma # 1.8 is stated in MP and used in matlab
  lower              <- fit.Diffs - 1.8 * sigma
  
  outliers <- Diffs > upper | Diffs < lower
  
  # hack to reproduce matlab: can't have the very first diff be an outlier:
  outliers[1] <- FALSE
  
  # if there's only 1 point, or if there are more than 4 points, forget about it
  if (sum(outliers) <= 1 | sum(outliers) > 4){ # in matlab, not in MP!
    return(as.data.frame(cbind(Age = ages, DeathsForFitting = dths, Original = 1)))
  }
 
  # we want 2 outliers only, one of each sign. If there are 2 or more of the same sign, need to reduce to 1
  Signs             <- sign(Errors[outliers])
  
  # not likely, but if both outliers were of the same sign, we'd reduce to one.

  # two largest outliers of opposite sign
#  if (method == 1){
#    if (length(table(Signs)) !=2){
#      return(as.data.frame(cbind(Age = ages, DeathsForFitting = dths, Original = 1)))
#    }
#    outliers          <- as.logical(rowSums(sapply(c(-1,1), function(sgn, Errors, outliers, Signs){
#                                Errors == max(abs(Errors[outliers][Signs == sgn])) * sgn
#                         }, Errors = Errors, outliers = outliers, Signs = Signs)))
#  }
#  
  # two largest outliers overall, same direction or not
#  if (method == 2){
#    outlind  <- which(outliers)
#    outliers <- abs(Errors) %in% rev(sort(abs(Errors[outlind])))[1:2]
#  }
  # two largest outliers of opposite sign, IFF it's a trough

# TODO: only use imputed model values IFF the total number of outlier deviations is less than 5

  #if (method == 3){
    if (length(table(Signs)) !=2){
      return(as.data.frame(cbind(Age = ages, DeathsForFitting = dths, Original = 1)))
    }
    outliers          <- as.logical(rowSums(sapply(c(-1,1), function(sgn, Errors, outliers, Signs){
                                Errors == max(abs(Errors[outliers][Signs == sgn])) * sgn
                         }, Errors = Errors, outliers = outliers, Signs = Signs)))
                   
    if (sign(Errors[outliers])[1] != 1 | sum(outliers) == 1){
      return(as.data.frame(cbind(Age = ages, DeathsForFitting = dths, Original = 1)))
    }
  #}
  Inds              <- which(outliers)
  
  # not in MP, but found in matlab:
  # if we're thinking of removing data that span more than 7 years, forget it:
  if (abs(diff(Inds)) > 7){
    return(as.data.frame(cbind(Age = ages, DeathsForFitting = dths, Original = 1)))
  }
  
  # this is NOT identical to the matlab routine, but may give the same results. Needs checking
  RemoveInd1        <- Inds[which.max(Errors[outliers])] # pick out max (NOT ABSOLUTE DEV!)
  RemoveInd2        <- Inds[which.min(Errors[outliers])]
  RMind             <- sort(RemoveInd1:RemoveInd2) + 1
  RA                <- ages[RMind]
 
  # now fit spline to deaths, but with outlier ages (and intermediate ages) excluded
  # the matlab code uses deaths in ALL ages:
  agesall <- as.integer(deathsRR$Age)
  RMall   <- agesall %in% RA
  
  D.fill.in         <- predict(smooth.spline(
                               x = agesall[!RMall],
                               y = DTHS[!RMall], 
                               all.knots=TRUE,
                               spar = spar2    # CB/TR smoother param .9 as seen in matlab for use with csaps()
                             ), RA)$y          # this fit, for use with imputation, appears to be undocumented in the MP...
  dths.interp              <- dths
  dths.interp[RMind]       <- D.fill.in
  
  # not in MP!
  # per matlab code, only take it if it makes a greater than 10% difference:
  if (max(abs(dths.interp - dths) / (dths + 0.0001)) < .1){
    return(as.data.frame(cbind(Age = ages, DeathsForFitting = dths, Original = 1)))
  }
  
  # all we need for fitting in d_soainew() are ages and imputed (or not) deaths.
  deaths.out                 <- as.data.frame(cbind(Age = ages, DeathsForFitting = dths.interp, Original = 1))
  deaths.out$Original[RMind] <- 0
  
  return(deaths.out)
  # now spit back deaths for Kannisto fitting. Just give back the RR deaths, all ages.
}



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


d_CohortAdjust <- function(Dnew, OA, ACVVmatrix, OAL, DOA, reproduce.matlab = TRUE){
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
  }
# 