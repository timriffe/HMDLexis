
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
    dths.fit          <- d_remfluct(deathsRR, N = 25, p1 = 0.0005, p2 = 0.9)
    
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
    ab                <- try(optim(c(a=.1,b=.1), 
                               fn = KannistoSurvMin, 
                               P = SyntheticSurvival, 
                               oa = OA, 
                               N = N, 
                               x = Ages4Kannisto)$par)
    if (class(ab) == "try-error"){
      ab             <- optim(c(a=.1,b=.1), 
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
    #
    Dnew              <- d_CohortAdjust(Dnew = Dnew, 
                                        OA = OA, 
                                        ACVVmatrix = ACVVmatrix, 
                                        OAL = OAL, 
                                        DOA = DOA, 
                                        reproduce.matlab = reproduce.matlab)

    # this can only happen if DOA is 0, but we'd get NaNs in Deaths due to 0 in denom in various locations. Just overwrite
    Dnew$Deaths[is.nan(Dnew$Deaths)] <- 0
    
    DeathsYearSex     <- rbind(DeathsYearSex[DeathsYearSex$AgeInterval != "+", ], Dnew)
    
    invisible(resortDeaths(DeathsYearSex))
  }
# head(Dnew)
# head(d_CohortAdjust(Dnew = Dnew, OA = OA, ACVVmatrix = ACVVmatrix, OAL = OAL, DOA = DOA))
