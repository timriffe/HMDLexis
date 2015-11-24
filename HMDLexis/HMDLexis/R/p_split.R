

#'
#' @title split age intervals based on a future census
#' 
#' @description This function implements the MP verbal description of splitting population counts in age groups. Basically, we add deaths going back in time from census 2 (right-side), which will give an estimate for the single-age populations in the interval \code{[a,a+n]}. We then rescale these single age estimates to sum to the C1 population count in the interval.
#' 
#' @param C1 the standard \code{Pop} object, selected only for rows containing the first (left-side) census.
#' @param C2 the standard \code{Pop} object, selected only for rows containing the second (right-side) census.
#' @param Deaths the standard \code{Deaths} IDB object, after computations finalized.
#' @param Births the standard \code{Births} IDB object.
#' @param a the age where extinct cohort populations will take over. No need to split intervals above this age.
#' @param reproduce.matlab logical. Default \code{FALSE}. Should we reproduce all aspects of the matlab code? 
#' 
#' @details The only matlab oddity that is potentially affected by \code{reproduce.matlab} is the handling of dates. Matlab does not handle leap years. This function makes calls to two other LexisDB functions, such as \code{ypart()} and \code{p_addCohortColumn()}.
#' 
#' @return C2 in single ages
#' 
#' @export 
#' 
#' @importFrom reshape2 acast
#' @importFrom reshape2 melt

# General strategy:
# 1) adjust C2 as in ic, but not C1...
# 2) derive C1hat using C2 + VVcum
# 3) adjust C1hat using f1 (as per at beginning of ic for C1)
# 4) then rescale to sum to population in interval of C1.

# will need some reference objects, or clever matrix way to do this
#
#

p_split_inner <- function(C1, C2, Deaths, Births, a = 80, reproduce.matlab = FALSE){
  Deaths  <- d_addCohortColumn(Deaths)
  # cut off UNK if necessary
  # UNK should be redistributed prior to running p_ic()
  stopifnot(all(C1$Age != "UNK"))
  stopifnot(all(C2$Age != "UNK"))
  # there may be a need to extend open age group here, but it needs 
  # to happen before this function is called.
  # get rid of open age group:
  OAI      <- C1[C1$AgeInterval == "+", ]
  C1       <- C1[C1$AgeInterval != "+", ]
  
  Interval <- unique(C2$Year) - unique(C1$Year)
  aM       <- max(C1$Agei) + C1$AgeIntervali[which.max(C1$Agei)] - 1
  # Pop <- C2
  if (max(C2$Agei[C2$AgeIntervali %==% 1]) < (aM + Interval)){
    C2   <- p_soai(C2, Deaths, a = a, reproduce.matlab = reproduce.matlab)
    
    # need to pad in case of reproduce.matlab = TRUE...
    maxM <- max(C2$Agei[C2$Sex == "m"])
    maxF <- max(C2$Agei[C2$Sex == "f"])
    maxA <- min(c(maxM,maxF))
    if (maxA < 130){
      padA             <- (maxA + 1):130
      CHUNK            <- C2[1:length(padA), ]
      CHUNK$NoteCode3  <- "p_split()"
      CHUNK$Age        <- CHUNK$Agei <- padA
      CHUNK$AgeInterval <- CHUNK$AgeIntervali <- 1
      CHUNK$Population <- 0
      CHUNKF           <- CHUNK
      CHUNKF$Sex       <- "f"
      CHUNK$Sex        <- "m"
      CHUNKF           <- CHUNKF[CHUNKF$Agei > maxF, ]
      CHUNK            <- CHUNK [CHUNK$Agei > maxM, ]
      C2               <- resortPops(rbind(C2, CHUNK, CHUNKF))
    }
    
  } else {
    C2       <- C2[C2$AgeInterval != "+", ]
  }
  
  # Check for single ages in C2
  stopifnot(all(C2$AgeIntervali == 1))
  # upper age checks happen in the sex loop
  C1singles  <- C1[C1$AgeIntervali %==% 1, ]
  C1         <- C1[C1$AgeIntervali %!=% 1, ]
  # one Sex at a time.
  sexes  <- unique(C1$Sex)
  PopOut <- list()
  PopOut[["Singles"]] <- C1singles
  for (Sex in sexes){ # Sex <- "f"
    C1s  <- C1[C1$Sex == Sex, ]
    C2s  <- C2[C2$Sex == Sex, ]
    Dsex <- Deaths[Deaths$Sex == Sex, ]
    
    ####################################################################
    # Get general parameters sorted out                                #
    ####################################################################
    
    ##############
    # 1) f1, f2. #
    ##############
    # this lengthy code chunk might need to be externalized?
    # note, we can't just use the reproduce.matlab argument of ypart
    # because matlab oddly just takes the first element, which opens up the
    # possibility of there being multiple dates present for a particular year
    # this DOES occur in the InputDB... and I don't know whether it ever affects
    # results.
    if (!reproduce.matlab){
      f1 <- ypart(Year = unique(C1$Year), 
        Month = unique(C1$Month), 
        Day = unique(C1$Day), 
        reproduce.matlab = reproduce.matlab)
      f2 <- ypart(Year = unique(C2$Year), 
        Month = unique(C2$Month), 
        Day = unique(C2$Day), 
        reproduce.matlab = reproduce.matlab)
      
      if (length(f1) > 1 | length(f2) > 1){
        stop("multiple dates in use in the same year, makes intercensals tricky.\\Time to do some digging.")
      }
    } else {
      f1 <- ypart(Year = C1$Year[1], 
        Month = C1$Month[1], 
        Day = C1$Day[1], 
        reproduce.matlab = reproduce.matlab)
      f2 <- ypart(Year = C2$Year[1], 
        Month = C2$Month[1], 
        Day = C2$Day[1], 
        reproduce.matlab = reproduce.matlab)
    }
    
    ######################################################################
    # 2) determine the years we need for selecting data                  #
    ######################################################################
    
    years  <- unique(C1s$Year):unique(C2s$Year)  # so 'years' is useful in all cases.
    N      <- length(years) 
    
    ######################################################################
    # 3) adjust right-side census pops to single birth cohorts  #
    ######################################################################
    
    C2s$Population[-1] <- (1 - f2) * C2s$Population[-nrow(C2s)] + f2 * C2s$Population[-1]
    # no need to remove last row of C2, since it has been overwritten now.
    # the first row (infants) remain untouched.
    
    # now we add birth cohort column to C2:
    C2s        <- p_addCohortColumn(C2s) 
    
    # this accounts for the above shifting that we just did 
    # (not in the standard p_addCohortColumn() method)
    C2s$Cohort <- C2s$Cohort + 1
    
    #######################################################################
    # 4) make sure we have what we need in C2, then cut down C1 as needed #
    #######################################################################
    
    # This is a good place for a stop condition
    # this may involve *not* bothering to split C1 ages 80+.

    C1SingleAgesImplied     <- with(C1s,min(Agei):((max(Agei) + AgeIntervali[which.max(Agei)]) - 1))
    C1grouped               <- c(unlist(apply(C1s[, c("Agei","AgeIntervali")], 1, 
        function(x){
          rep(x[[1]],x[[2]])
        })))
    # the C1 cohorts we actually have represented 
    # TODO: this may be off by one, if we need to shift UP one to account for f1 blending
    # which we haven't yet done... Figure this out when we get so far.
    C1Cohorts               <- years[1] - C1SingleAgesImplied - 1
    
    CutThese                <- C1grouped[!C1Cohorts %in% C2s$Cohort]
    C1SliceOffAge           <- ifelse(length(CutThese)>= 1, min(CutThese),131)
    
    # cut down C1s to just needed ages
    C1s                     <- C1s[C1s$Agei < C1SliceOffAge, ]
    # cut down auxiliary-vectors too
    C1grouped               <- C1grouped[C1grouped < C1SliceOffAge]
    C1SingleAgesImplied     <- C1SingleAgesImplied[C1SingleAgesImplied < C1SliceOffAge]
    # these are the cohorts we'll work with (not including infant birth cohort)
    CohortsConnect          <- years[1] - C1SingleAgesImplied - 1
    
    # a simple check
    stopifnot(max(C1SingleAgesImplied) >= (a - 1))
    
    # NOTE TO SELF: figuring out first and last cohorts. Infant or not... 
    # Earlier cutoff needs to be thought through
    # figure out which cohorts to use for what (two slices)   
    Ccoh                    <- CohortsConnect
    IcohB                   <- years[1]
    ######################################################################
    # 5) get death matrices, triangles and PC parallelograms             #
    ######################################################################
    # OK, now for estimation components
    DL      <- acast(Dsex[with(Dsex, Year %in% years & Lexis == "TL"), ], 
      Year ~ Cohort, 
      sum, 
      value.var = "Deaths", 
      fill = NA_real_, drop = FALSE)
    DU      <- acast(Dsex[with(Dsex, Year %in% years & Lexis == "TU"), ], 
      Year ~ Cohort, 
      sum, 
      value.var = "Deaths", 
      fill = NA_real_, drop = FALSE) 
    VV      <- acast(Dsex[with(Dsex, Year %in% years), ], 
      Year ~ Cohort, 
      sum, 
      value.var = "Deaths", 
      fill = NA_real_, drop = FALSE)  
    
    # DL, DU only needed for Da,Db,Dc,Dd   : not used later
    # VV is used later to produce cumulative deaths
    
    ###########################################################################
    # 6) produce Da Db Dc Dd from MP equations, to adjust 1st and last years  #
    ###########################################################################
    
    # these elements are named as in MP equations: 
    # left-side adjustment for partial years
    Da      <- (1 - f1^2) * DL[1, ]
    Db      <- (1 - f1)^2 * DU[1, ]
    # right-side adjustment for partial years
    Dc      <- f2^2 * DL[nrow(DL), ]
    Dd      <- ((2 * f2) - f2^2) * DU[nrow(DU), ]
    
    ####################################################################
    # 7) Completed cohorts (cut objects down to conform dimensions)    #
    ####################################################################
    
    # select first and last year adjustments
    CDa              <- Da[as.character(Ccoh)]
    CDb              <- Db[as.character(Ccoh)]
    CDc              <- Dc[as.character(Ccoh)]
    CDd              <- Dd[as.character(Ccoh)]
    
    # CVV can now be used for cumulative deaths
    CVV              <- VV[, as.character(Ccoh)]
    CVV[1, ]         <- CDa + CDb
    CVV[nrow(CVV), ] <- CDc + CDd
    
    # get simple starting and ending vectors of population
    C2vec            <- C2s$Population[C2s$Cohort %in% Ccoh]
    # assuming no error / migration, we get:
    C1hat            <- C2vec + colSums(CVV)
    # this needs to be staggered to match f1 prior to forming pdf and rescaling!
    
#    # we want a pdf for each age group
#    CagesGroups      <- C1grouped[CohortsConnect %in% Ccoh]
#    Cpdf             <- unlist(tapply(C1hat,CagesGroups,function(x){x/sum(x)}))
#    # now multiply into C1 as needed
#    CPopRep          <- with(C1s, Population[Agei %in% CagesGroups])[as.character(CagesGroups)]
#    
#    Cpop             <- Cpdf * CPopRep
    
    ####################################################################
    # 8) Infant Cohort                                                 #
    ####################################################################
    
    # get IC starting pop for infant cohort
    # part is births:
    IB       <- Births$Births[Births$Year == IcohB & Births$Sex == Sex]
    C2B      <- C2s$Population[C2s$Cohort == IcohB]
    # we want a little segment, with a tricky formula. First make a
    # clean intercensal and get an estimate of the Jan 1 population of
    # year years[1] + 1.
    IBVV     <- VV[, as.character(IcohB)]
    IBVV[length(IBVV)] <- Dc[as.character(IcohB)] + Dd[as.character(IcohB)]
    IBerror  <- IB - (C2B + sum(IBVV))
    # now what proportion of this error gets us to age 0 Jan 1 years[1]+1?
    IBprop   <- .5  / (N - 3 / 2 + f2 )
    # our estimate for age zero next Jan1:
    H0est    <- IB - DL[1, as.character(IcohB)] + IBerror * IBprop
    C10_f1   <- (f1 * IB) - (f1^2) * (IB - H0est)
    
    # save C10_f1 for use in a moment...
    # now we have other pre-existing single-age cohorts to estimate for C1, similar to above:
    
    ####################################################################
    # 9) Stagger cohorts to match C1 single ages                       #
    ####################################################################
    
    C1hatAdj <- C1hat[1:(length(C1hat)-1)] * f1 + C1hat[2:length(C1hat)] * (1 - f1)
    C10hat   <- C1hat[1] * (1-f1) + C10_f1
    names(C10hat) <- years[1]
    
    # this is everything, but the cohort labels are meaningless, since this is mid year.
    # really we have a vector that starts at age zero and counts up in single ages... hat
    C1hat <- c(C10hat, C1hatAdj)
    
    ####################################################################
    # 10) rescale within age groups to pdf, then expand and multiply   #
    ####################################################################
    # we want a pdf for each age group
    
    C1hatpdf             <- unlist(tapply(C1hat[1:length(C1grouped)],C1grouped,function(x){
          x/sum(x)
        }
      ))
    # easy arithmetic to simply repeat each age group value N times, since C1grouped
    # was painstakingly made.
    PI               <- with(C1s, Population[Agei %in% C1grouped])
    names(PI)        <- with(C1s, Age[Agei %in% C1grouped])
    CPopRep          <- PI[as.character(C1grouped)]
    # now multiply into C1
    Cpop             <- C1hatpdf * CPopRep
    
    #################################################################### 
    # 11) now we need to deck the halls with all the other columns...  #
    ####################################################################
    Ps <- data.frame(PopName = unique(C1s$PopName),
      Area = unique(C1s$Area),
      Sex  = Sex,
      Age = C1SingleAgesImplied,
      AgeInterval = 1,
      Type = unique(C1s$Type),
      Day = unique(C1s$Day),
      Month = unique(C1s$Month),
      Year = unique(C1s$Year),
      RefCode = NA,
      Access = "O",
      Population = Cpop,
      NoteCode1 = "p_split()",
      NoteCode2 = NA,
      NoteCode3 = NA,
      LDB = 1,
      Agei = C1SingleAgesImplied,
      AgeIntervali = 1,
      stringsAsFactors = FALSE
    )
    PopOut[[Sex]] <- Ps
  }
  # exit Sex-loop, bind together, resort and return
  PopOut <- rbind(do.call(rbind,PopOut),OAI)
  PopOut <- resortPops(PopOut)
  rownames(PopOut) <- NULL
  invisible(PopOut)
}


#'
#' @title split age intervals based on single ages from the right
#' 
#' @description This function implements the MP verbal description of splitting population counts in age groups. Basically, we add deaths going back in time from census 2 (right-side), which will give an estimate for the single-age populations in the interval \code{[a,a+n]}. We then rescale these single age estimates to sum to the C1 population count in the interval.
#' 
#' @param Pop Standard \code{Population} IDB object. At minimum, the most recent year must be in single ages.
#' @param Deaths the standard \code{Deaths} IDB object, after computations finalized.
#' @param a the age where extinct cohort populations will take over. No need to split intervals above this age.
#' @param reproduce.matlab logical. Default \code{FALSE}. Should we reproduce all aspects of the matlab code? 
#' 
#' @details The only matlab oddity that is potentially affected by \code{reproduce.matlab} is the handling of dates. Matlab does not handle leap years. This function makes calls to two other LexisDB functions, such as \code{ypart()} and \code{p_addCohortColumn()}.
#' 
#' @return C2 in single ages
#' 
#' @export 
#' 
#' @importFrom reshape2 acast
#' @importFrom reshape2 melt

# in this data, we have 1960 and 1970 that are in 5-year age groups.
p_split <- function(Pop, Deaths, Births, a = 80, reproduce.matlab = FALSE){
  Pop     <- Pop[Pop$Age != "TOT", ]
  yearssp <- unique(Pop$Year[Pop$AgeIntervali %>% 1])
  yearssp <- sort(yearssp, decreasing = TRUE)
 
  years   <- sort(unique(Pop$Year))
  # a silly check
  if (max(yearssp) > max(Pop$Year)){
    stop("uh oh. We split grouped ages from the right,\n but your most recent year of data is grouped.\nTime to think... Function stopped")
  }

  for (yr in yearssp){ # yr <- 1970
    C1    <- Pop[Pop$Year == yr, ]
    C2    <- Pop[Pop$Year == years[years > yr][1], ]
    message("p_split(): year", yr)
    C1out <- p_split_inner(C1, C2, Deaths, Births, a = a, reproduce.matlab = reproduce.matlab )
    Pop   <- resortPops(rbind(Pop[Pop$Year != yr, ], C1out))
  }
 
  # return output...
  invisible(Pop)
}

