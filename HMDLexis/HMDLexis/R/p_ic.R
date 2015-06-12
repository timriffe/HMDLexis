# --------------------------------------------------------------------------
# p_ic() is a wrapper around p_ic_inner(), which is defined below in this same script.

#'
#' @title generalized intercensal estimation 
#' 
#' @description This function performs intercensal estimation for all intercensal periods containing at least one full calendar year. All age intervals must be single. Some checking is done to ensure upper ages are high enough. For C1 we must have a single-age estimate for population at age 79 at the lowest. C2 needs a single age estimate for at least age \code{79 + N}, where N is the number of years from C1 to C2, inclusive. If either census has open age groups that get in the way of this, we throw an error and stop the function. Open age groups must be dealt with prior to calling this function in any case. Likewise, age intervals greater than one must be split prior to calling this function.
#' 
#' @param Pop the standard \code{Pop} object, after age intervals and unknown ages have been dealt with.
#' @param Deaths the standard \code{Deaths} IDB object, after computations finalized.
#' @param Births the standard \code{Births} IDB object, after computations finalized.
#' @param reproduce.matlab logical. Default \code{FALSE}. Should we reproduce all aspects of the matlab code? 
#' 
#' @details The only Matlab oddity that is potentially affected by \code{reproduce.matlab} is the handling of dates. Matlab does not handle leap years. This function makes calls to two other LexisDB functions, such as \code{ypart()}, \code{d_addCohortColumn}, \code{p_addCohortColumn()} and \code{resortPops()}. Most work, and all MP formulas are found in the inner function, \code{p_ic_inner()}, which takes two clean censuses as its inputs.
#' 
#' @return Pop a standard population data.frame, in single ages, including all years in the intercensal period. \code{C1} and \code{C2} are only included in this object if they were January 1 estimates.
#' 
#' @export 
#' 
#' @importFrom reshape2 acast
#' @importFrom reshape2 melt

p_ic <- function(Pop, Deaths, Births, reproduce.matlab = FALSE){
  
  Deaths  <- d_addCohortColumn(Deaths)
  
  # remove TOT, as always
  Pop <- Pop[Pop$Age != "TOT", ]
  
  # UNK should be redistributed prior to running this function
  stopifnot(all(Pop$Age != "UNK"))
  
  Allyears  <- min(Pop$Year):max(Pop$Year)
  # we don't want to base intercensals on other intercensals, so we remove them 

  Cyears  <- sort(unique(Pop$Year))
  C1years <- Cyears[c(diff(Cyears) > 1, FALSE)]
  C2years <- rev(rev(Cyears)[c(diff(rev(Cyears)) < -1,FALSE)])
  
  # if this throws an error it means the detection strategy needs a rethink
  stopifnot(length(C1years) == length(C2years))
  
  # this is the total number of times we need to iterate...
  Y        <- length(C1years)
  
  # separate Population counts not involved in intercensals. 
  # Stick together later
  Pdontuse <- Pop[!Pop$Year %in% C2years, ]
  
  PIC <- list()
  
  # Y is the number of intercensal periods we need to deal with.
  # C1 <- C1i; C2 <- C2i
  for (i in 1:Y){ 
    C1i <- Pop[Pop$Year == C1years[i], ]
    C2i <- Pop[Pop$Year == C2years[i], ]
   
    PIC[[i]] <- p_ic_inner(C1 = C1i, 
                           C2 = C2i, 
                           Deaths = Deaths, 
                           Births = Births, 
                           reproduce.matlab = reproduce.matlab)
  }
  # maintain the unused populations
  PIC[[Y + 1]] <- Pdontuse
  
  Pout <- do.call(rbind, PIC)
  
  Pout <- resortPops(Pout)
  # cut off UNK if necessary
  invisible(Pout)
} 

# --------------------------------------------------------------------------
# p_ic_inner()

#'
#' @title generalized intercensal estimation between two censuses only
#' 
#' @description This function is called inside \code{p_ic()} on individual census pairs specified to be interpolated. We assume that all age intervals are 1, and we do light checking for upper ages. For C1 we must have a single-age estimate for population at age 79 at the lowest. C2 needs a single age estimate for at least age \code{79 + N}, where N is the number of years from C1 to C2, inclusive. If either census has open age groups that get in the way of this, we throw an error and stop the function. Open age groups must be dealt with prior to calling this function in any case. Likewise, age intervals greater than one must be split prior to calling this function.
#' 
#' @param C1 the standard \code{Pop} object, selected only for rows containing the first (left-side) census.
#' @param C2 the standard \code{Pop} object, selected only for rows containing the second (right-side) census.
#' @param Deaths the standard \code{Deaths} IDB object, after computations finalized.
#' @param Births the standard \code{Births} IDB object, after computations finalized.
#' @param reproduce.matlab logical. Default \code{FALSE}. Should we reproduce all aspects of the matlab code? 
#' 
#' @details The only matlab oddity that is potentially affected by \code{reproduce.matlab} is the handling of dates. Matlab does not handle leap years. This function makes calls to two other LexisDB functions, such as \code{ypart()}, \code{resortPops()}, and \code{p_addCohortColumn()}.
#' 
#' @return Pop a standard population data.frame, in single ages, including all years in the intercensal period. \code{C1} and \code{C2} are only included in this object if they were January 1 estimates.
#' 
#' @export 
#' 
#' @importFrom reshape2 acast
#' @importFrom reshape2 melt

p_ic_inner <- function(C1, C2, Deaths, Births, reproduce.matlab = FALSE){
  
  Deaths  <- d_addCohortColumn(Deaths)
  # cut off UNK if necessary 
  # UNK should be redistributed prior to running p_ic()
  stopifnot(all(C1$Age != "UNK"))
  stopifnot(all(C2$Age != "UNK"))
  # there may be a need to extend open age group here, but it needs 
  # to happen before this function is called.
  # get rid of open age group:
  C1  <- C1[C1$AgeInterval != "+", ]
  C2  <- C2[C2$AgeInterval != "+", ]
  
  # Check for single ages:
  stopifnot(all(c(all(C1$AgeIntervali == 1),all(C2$AgeIntervali == 1))))
  # upper age checks happen in the sex loop
  
  # one Sex at a time.
  sexes  <- c("m","f")
  PopOut <- list()
  for (Sex in sexes){ # Sex <- "m"
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
    # we ALWAYS produce an estimate for the last year
    # when f2 = 0, we reproduce C2.
    # when f1 = 0, we reproduce C1
    # if f1 > 0, we still need death data from C1 year, but
    #            but we won't re-estimate C1
    # if f2 > 0 we won't reproduce C2, but we produce an estimate for Jan1 of that year,
    # and still need deaths from that year...
    years  <- unique(C1s$Year):unique(C2s$Year)  # so 'years' is useful in all cases.
    N      <- length(years) 
    
    # let's now make sure ages go high enough:
    # anything higher than 80 will just be replaced by EC/SR
    stopifnot(max(C1s$Agei) >= 79)
    
    # C2 needs to connect cohorts that go back to age 80 at C1
    stopifnot(max(C2s$Agei) >= 79 + N)
    
    ######################################################################
    # 3) adjust starting and ending census pops to single birth cohorts  #
    ######################################################################
    
    # TODO: this method, currently in MP, could be improved given the monthly birth series,
    C1s$Population[-1] <- (1 - f1) * C1s$Population[-nrow(C1s)] + f1 * C1s$Population[-1]
    
    ###########################################################################
    # TR: sent email to Dima Jan 12, 2015 to clarify this matlab chunk. I assumed
    # that it's a bug, and I also proposed a fix in matlab....
    # ----------------------------------------------------------------
    # some old ages on the right-hand side of c2 are under-represented after adjustment.
    # this is the culprit line in matlab:
    #C1=(1-f1c)*c1(i-1)+f1c*c1(i);
    #if i-1+N+2>length(c2)-N   % condition catches too many ages....
    #    disp(i); % TEST
    #  C2=(1-f2c)*c2(i-1+N+1) % <- this is not right.
    #  else
    #    C2=(1-f2c)*c2(i-1+N+1)+f2c*c2(i-1+N+2); % this is not clear to me
    #  end
    ###########################################################################
    
    C2s$Population[-1] <- (1 - f2) * C2s$Population[-nrow(C2s)] + f2 * C2s$Population[-1]
    # no need to remove last row of C2, since it has been overwritten now.
    # the first row (infants) remain untouched.
    
    # now we add birth cohort column to C1 and C2:
    C1s        <- p_addCohortColumn(C1s) 
    C2s        <- p_addCohortColumn(C2s) 
    
    # this accounts for the above shifting that we just did 
    # (not in the standard p_addCohortColumn() method)
    C1s$Cohort <- C1s$Cohort + 1
    C2s$Cohort <- C2s$Cohort + 1
    
    ##############################################################
    # TODO: Tadj: either adjust C1 or C2, but not both, and be consistent with deaths
   ##############################################################
    
    ######################################################################
    # 4) get death matrices, triangles and PC parallelograms             #
    ######################################################################
    # OK, now for estimation components
    DL      <- acast(Dsex[with(Dsex, Year %in% years & Lexis == "TL"), ], 
      Year ~ Cohort, 
      sum, 
      value.var = "Deaths", 
      fill = NA_real_)
    DU      <- acast(Dsex[with(Dsex, Year %in% years & Lexis == "TU"), ], 
      Year ~ Cohort, 
      sum, 
      value.var = "Deaths", 
      fill = NA_real_) 
    VV      <- acast(Dsex[with(Dsex, Year %in% years), ], 
      Year ~ Cohort, 
      sum, 
      value.var = "Deaths", 
      fill = NA_real_)  
    
    ##############################################################
    # TODO: Tadj: either adjust each of the 3 matrices, but elementwise multiplication
    # with some sort of conformable Vx matrix, and consistent with the Pop adjustment
    # we just did
    ##############################################################
    
    
    
    # DL, DU only needed for Da,Db,Dc,Dd   : not used later
    # VV is used later to produce cumulative deaths
    
    ###########################################################################
    # 5) produce Da Db Dc Dd from MP equations, to adjust 1st and last years  #
    ###########################################################################
    
    # these elements are named as in MP equations: 
    # left-side adjustment for partial years
    Da      <- (1 - f1^2) * DL[1, ]
    Db      <- (1 - f1)^2 * DU[1, ]
    # right-side adjustment for partial years
    Dc      <- f2^2 * DL[nrow(DL), ]
    Dd      <- ((2 * f2) - f2^2) * DU[nrow(DU), ]
    
    ###########################################################################
    # 6) figure out which cohorts to use for what                             #
    ###########################################################################
    
    # which are complete?
    obs     <- colSums(!is.na(DL))
    Ni      <- as.integer(names(table(obs)))[which.max(table(obs))]
    
    # completed cohorts (C_):
    Ccoh    <- as.integer(colnames(DL)[obs == Ni])
    Ccoh    <- Ccoh[Ccoh < years[1]]
    Ccoh    <- Ccoh[Ccoh %in% C1s$Cohort & Ccoh %in% C2s$Cohort]
    # infant cohort (I_):
    Icoh    <- years[1]
    # newborn cohorts (N_):
    Ncoh    <- years[-c(1,length(years))]
    
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
    C1vec   <- C1s$Population[C1s$Cohort %in% Ccoh]
    C2vec   <- C2s$Population[C2s$Cohort %in% Ccoh]
    #C1s[C1s$Cohort %in% Ccoh,c("Age","Cohort","Population")]
    # assuming no error / migration, we get:
    C2hat   <- C1vec - rev(colSums(CVV))
    
    # remove recent year needs to happen AFTER C2hat, but before cumulative sums
    CVV     <- CVV[-nrow(CVV), ]
    # cumulative deaths (already discounted for Da Db, etc)
    CVVcum  <- t(apply(CVV, 2, cumsum))[rev(seq(ncol(CVV))),]
    
    # what is the total error (mostly migration we hope) by the end of the period?
    CDelta  <- C2vec - C2hat
    # full year time steps (potentially fragile)
    n       <- 0:(N - 2)
    # how to partition error (potentially fragile)
    prop    <- (1 - f1 + n) / (N - 1 - f1 + f2)
    # the cumulative error
    CDelta  <- outer(prop, CDelta, "*")
    #dim(CDelta);dim(CVV)
    
    # if dimensions do not conform, this breaks, and we investigate above
    Cpop    <- C1vec - CVVcum + t(CDelta)
    # length(C1vec);dim(CVVcum);dim(t(CDelta))
    # these are now a PC object, save to be augmented by the infant cohort
    colnames(Cpop) <- as.integer(colnames(Cpop)) + 1
    
    ####################################################################
    # 7) Infant cohort                                                 #
    ####################################################################
    
    # get IC starting pop for infant cohort
    # part is births:
    IB      <- Births$Births[Births$Year == Icoh & Births$Sex == Sex]
    # part is census:
    IC      <- C1s$Population[C1s$Cohort == Icoh]
    # combine using f1:
    I0      <- (1 - f1) * IB + f1 * IC # this is the starting pop
    
    # Da, Dc, Dd for infants (Db not needed)
    IDa     <- Da[as.character(Icoh)]
    IDc     <- Dc[as.character(Icoh)]
    IDd     <- Dd[as.character(Icoh)]
    
    # get deaths over cohort, and accumulate
    IVV     <- VV[, as.character(Icoh)]
    IVV[1]  <- IDa
    IVV[length(IVV)] <- IDc + IDd
    IVVcum  <- cumsum(IVV[-length(IVV)])
    # estimated size of cohort on arrival at C2
    Ihat    <- I0 - sum(IVV)
    # actual population at C2
    I2      <- C2s$Population[C2s$Cohort == Icoh]
    # the total error when based only on decrement
    IDelta  <- I2 - Ihat
    
    # determine time steps, proportions to distribute error over period
    n       <- 0:(N - 2)
    Iprop   <- ((1 - f1^2) / 2 + n)  / (N - 2 + (1 - f1^2) / 2 + f2 )
    
    # estimate population over period
    Ipop    <- I0 - IVVcum + Iprop * IDelta
    # increment years (deaths from prev map to next)
    names(Ipop) <- as.integer(names(Ipop)) + 1
    # put in matrix row form to easily combine later
    Ipop    <- matrix(Ipop, nrow=1, dimnames=list(Icoh,names(Ipop)),byrow=TRUE)
    
    ####################################################################
    # 8) Newborn cohorts                                               #
    ####################################################################
    
    # Births enter into LDB sorted, so no worries.
    NB      <- Births$Births[Births$Year %in% Ncoh & Births$Sex == Sex]
    # these go from C2 backward, so we want to rev them
    N2      <- rev(C2s$Population[C2s$Cohort %in% Ncoh])
    names(N2) <- Ncoh
    
    # right-side adjustment
    NDc     <- Dc[as.character(Ncoh)]
    NDd     <- Dd[as.character(Ncoh)]
    
    # fragile verify this 
    NVV     <- VV[-1, as.character(Ncoh)]
    NVV[nrow(NVV), ] <- NDc + NDd
    
    NVV[is.na(NVV)] <- 0
    
    # estimated pop size, assuming only death decrement
    N2hat  <- NB - colSums(NVV, na.rm = TRUE)
    
    NVV    <- NVV[-nrow(NVV), ]
    
    # total error on arrival at C2
    NDelta <- N2 - N2hat
    
    # getting prop to distribute error over, tricky
    Nprop  <- NVV * 0
    # K, k follow MP naming conventions
    for (coh in 1:length(Ncoh)){ # coh = 1
      K            <- N - 2 - coh
      k            <- 0:K
      prop         <- (2 * k + 1) / (2 * K  + 1 + 2 * f2)
      Nprop[coh, ] <- c(rep(0,coh-1), prop)
    }
    # Total error
    NDelta  <- Nprop * NDelta
    
    # cumualtive deaths
    NVVcum              <- t(apply(NVV, 2, cumsum))
    
    # estimate population over period
    Npop    <- NB - NVVcum + NDelta
    Npop[NVVcum == 0] <- NA
    
    colnames(Npop) <- as.integer(colnames(Npop)) + 1
    ####################################################################
    # 9) combine estimates                                             #
    ####################################################################   
    
    ####################################################################
    # TODO: us Vx matrix to adjust back to original reference population
    # Note there is no need to change the guts of what's happening above.
    # we transform once before, and then once again afterwards. Note also
    # that we're like 100% sure that this is how it needs to happen, and I
    # rate a high likelihood of some matlab idiosyncracy coming to the fore
    ####################################################################
    
    # move to long format 
    Ps     <- rbind(
      melt(Cpop, value.name = "Population", varnames = c("Cohort", "Year")),
      melt(Ipop, value.name = "Population", varnames = c("Cohort", "Year")),
      melt(Npop, value.name = "Population", varnames = c("Cohort", "Year"))
    )
    Ps             <- Ps[!is.na(Ps$Population), ]
    
    # need AP in the end...
    Ps$Agei        <- Ps$Age <- Ps$Year - Ps$Cohort - 1
    
    # now add on the remaining columns...
    Ps$Sex         <- Sex
    Ps$PopName     <- unique(C1s$PopName)
    Ps$Area        <- unique(C1s$Area)
    Ps$LDB         <- Ps$Month <- Ps$Day <- Ps$AgeInterval <- Ps$AgeIntervali <- 1
    Ps$Type        <- "IC"
    Ps$Access      <- "O"
    Ps$NoteCode1   <- "p_ic()"
    Ps$RefCode     <- Ps$NoteCode2 <- Ps$NoteCode3 <- NA
    
    #sort(colnames(Ps));sort(colnames(C1s))
    Ps             <- Ps[,colnames(C1s)]
    
    # swap in original C1 and C2 if necessary, no need to overwrite
#    if (f1 == 0){
#      Ps <- Ps[Ps$Year != years[1], ]
#      Ps <- rbind(p_addCohortColumn(C1[C1$Sex == Sex, ]), Ps)
#    }
    if (f2 == 0){
      Ps <- Ps[Ps$Year != years[N], ]
      Ps <- rbind(p_addCohortColumn(C2[C2$Sex == Sex, ]), Ps)
    }
    # add to list
    PopOut[[Sex]] <- Ps
  }
  PopOut        <- do.call(rbind, PopOut)
  PopOut$Cohort <- NULL
  
  PopOut        <- resortPops(PopOut)
  invisible(PopOut)
}











