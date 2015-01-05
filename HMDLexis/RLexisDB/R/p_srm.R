#'
#' @title p_srm a simple survivor ratio method function, for modularity.
#' 
#' @description This implementation seeks to follow the MP. LDB \code{Deaths} must be finalized. Populations must be in single ages and all years, but not necessarily up to age 130 in all years. An argument \code{reproduce.matlab} controls whether an unnecessary matlab kludge is included. This function calls various other LexisDB functions at the moment, such as \code{d_addCohortColumn()}, \code{p_addCohortColumn()}, \code{p_ecm_findOmega()}, as well as its own auxilliary function, \code{p_srm_updatePopGivenCmid()}, which serves to speed up a loop for serial dependency. 
#' 
#' @param Pop The standard internal Population data.frame, *after* running p_ecm().
#' @param Deaths after all processing is done. Completed triangles.
#' @param k the parameter 'k' from various equations in the section on SRM.
#' @param m the parameter 'm' from various equations in the section on SRM.
#' @param l the parameter 'l' from various equations in the section on SRM.
#' @param MaxAge controls whether we do SR 90+ or SR 85+ (or something else). Default 90. (85 is untested so far)
#' @param maxit maximum number of iterations to optimize the improvement coefficient, 'c'.
#' @param reproduce.matlab logical. Do we include the legacy matlab kludge? Default \code{FALSE}. If \code{TRUE} be sure to have run \code{p_ecm()} with this argument set to \code{TRUE} first.
#' 
#' @return Pop with old ages adjusted using the survivor ratio method
#' 
#' @note This very function will be generalized to account for Tadj eventually. Further testing is needed to check whether the parameter \code{MaxAge} adjusts things properly for use with the SR85+ method.
#' 
#' @importFrom compiler cmpfun
#' @importFrom reshape2 acast
#' @importFrom reshape2 melt
#' 
#' @export 
#' 

# k = 5;m = 5;l=5;a=80;MaxAge = 90;maxit = 100
p_srm <- cmpfun(function(Pop, 
    Deaths, 
    k = 5, 
    m = 5, 
    l = 5,
    a = 80,
    MaxAge = 90, 
    maxit = 100,
    reproduce.matlab = FALSE){
    
    # if totals are there, get rid of them
    Pop             <- Pop[Pop$Age != "TOT", ]
    # we make interim columns, namely
    ColnamesKeep    <- colnames(Pop)  
    
    # slice off UNK, rbind back on later:
    UNKi            <- Pop$Age == "UNK"
    UNKTF           <- any(UNKi)
    if (UNKTF){
      UNK           <- Pop[UNKi, ]
      Pop           <- Pop[!UNKi, ]
    }
    
    # add cohort column to make our lives easier
    Deaths          <- d_addCohortColumn(Deaths)
    Pop             <- p_addCohortColumn(Pop)
    
    # some globals that can be shared for males and females
    this.year       <- max(Pop$Year)
    this.yearc      <- as.character(this.year)
    
    # to be filled later
    PopMF           <- list()
    # sex loop Sex <- "m"
    # allows for single or both sexes
    sexes <- unique(Pop$Sex)
    sexes <- sexes[sexes %in% c("m", "f")]
    for (Sex in sexes){
      
      Psex            <- Pop[Pop$Sex == Sex, ]
      SUM             <- sum(with(Psex,Population[Year == this.year & Agei >= MaxAge]))
      YearsKeep       <- unique(Psex$Year)
      Psex            <- Psex[Psex$AgeInterval != "+", ]
      Dsex            <- Deaths[Deaths$Sex == Sex & Deaths$Year <= this.year  , ]
      
      # omega is typically different for males and females, must be in the loop:
      omega           <- p_ecm_findOmega(Dsex, l = l, threshold = 0.5)
      w               <- omega["omega"] # 
      
      # This presumes that reproduce.matlab was TRUE for p_ecm() as well
      # Dima: this is because we needed a value for the first non-extinct cohort 
      # value to cumulate backward from. It was added in the ec script IF reproduce.matlab was TRUE
      # if it wasn't TRUE in that script, then it's not necessary to remove here.
      # BUT still need to recheck indexing of omega.
      if (reproduce.matlab){
        SUM <- SUM - omega["Ds"]
      }
      # some helpful vectors for selecting and iterating:
      srages          <- MaxAge:(w - 1)                       # have tried w, w-1,w-2
      srcohorts       <- this.year - 1 - srages               # this is exact
      sryears         <- (this.year - (w - a - 1)):this.year # this is not particularly key
      
      # for name-based matrix indexing:
      sragesc         <- as.character(srages)
      srcohortsc      <- as.character(srcohorts)
      sryearsc        <- as.character(sryears)
      
      # more indexing (recycled in aux function)
      yr.denom        <- as.character((this.year - 1 - k):(this.year - m - k))
      yr.num          <- as.character((this.year - 1):(this.year - m))
      
      # remove population estimates that are located in the SR zone:
      Psex            <- Psex[!(Psex$Cohort %in% srcohorts & Psex$Year %in% sryears & Psex$Agei >= a), ] 
      # cleans out SR zone
      
      # AP population (contains 0s in SR zone)
     
      pop             <- acast(Psex, Agei ~ Year, fun.aggregate = sum,  value.var = "Population")
      
      # cut down to approximate SR zone (need lower ages for R* estimation)
      # need to do this in case only 1 year given in Pop
      # this is an intermediate processing step to handle the odd case
      # of when a single year of data is given, as can happen if called by p_ic()
      yrskeep         <- as.character((min(sryears) - 5):this.year)
      ageskeep        <- as.character(70:130)
      popBox          <- matrix(nrow = length(ageskeep), 
                                ncol = length(yrskeep),
                                dimnames = list(ageskeep,
                                                yrskeep))
      yrshave         <- intersect(yrskeep, colnames(pop))
      ageshave        <- intersect(ageskeep,rownames(pop))
      popBox[ageshave, yrshave] <- pop[ageshave, yrshave]
     
      pop             <- popBox
        #Psex[Psex$Year == 1990, ] # check oddity with HUN 1990 females. why skip ages 80-89?
      # TODO: fill in NAs with left-right linear interpolation. Border NAs
      #       receive constant values. This will patch case where ic() not yet called
      #       given that only needed for Rstar, this should be pretty robust. In case
      #       only a single year of population is given, it will assume constant cmid = 1
      
      
      # get EC pops for SR zone omega
      VV              <- acast(Dsex[with(Dsex, Year %in% sryears & Cohort %in% srcohorts), ], 
                               Year ~ Cohort, sum, value.var = "Deaths")
                            
      # 2) take cumulative sums to get extinct cohort populations, then convert it to an Age x Year matrix, slicing
      # the rev is to keep age oriented correctly
      VVcumsums       <- apply(VV, 2, function(x){
                                rev(cumsum(rev(x)))
                       })
      # ------------------------------------------------------
      # get it into temporary matlab-imitating shape.
      # ------------------------------------
      # make D_dot, a recycled vector.
      # aggregate 5x1 cohorts counting back from most recent year. 
      D_dot           <- acast(Dsex[with(Dsex, Year > (this.year - 6)), ], 
                               . ~ Cohort, sum, value.var = "Deaths")[, , drop = TRUE][srcohortsc]
      # starting values give as a c of 1
      c_l             <- 0
      c_r             <- 10
    # z <- 1
      res             <- rep(0, maxit)
      for (z in 1:maxit){   
        # each iteration we calibrate cmid and re-fill the entire area.
        cmid    <- (c_l + c_r) / 2
        # p_srm_updatePopGivenCmid() defined below, in same R script
        # it's quirky...
#        pop <- p_srm_updatePopGivenCmid(pop = pop, 
#                        VVcumsums = VVcumsums, 
#                        cmid = cmid, 
#                        D_dot = D_dot, 
#                        srages = srages, 
#                        sragesc = sragesc, 
#                        srcohortsc = srcohortsc, 
#                        this.yearc = this.yearc, 
#                        yr.denom = yr.denom, 
#                        yr.num = yr.num, 
#                        k = k, 
#                        w = w,
#                        a = a)
        for (i in length(srages):1){ # i <- 22
                        
                        # first make R_star, the Ratio, given 'c'
                        denom          <- sum(pop[as.character(srages[i] - k), yr.denom ], na.rm = TRUE)
                        num            <- sum(pop[sragesc[i],  yr.num], na.rm = TRUE)
                        R_star         <- num / denom
                        R_star[is.nan(R_star) | is.infinite(R_star)] <- 0
                        
                        # impute current age-c-iteration of 'this.year' pop
                        pop.this.yr.i  <- cmid * D_dot[i] * R_star / (1 - R_star)
                        srAPVVi        <- PC2AP(VVcumsums[, srcohortsc[i], drop = FALSE] + pop.this.yr.i, 
                          agemin = a, agemax = (w - 1), Lexis = 1) # have tried w-2,w-1,w, 
                        NAindi         <- is.na(srAPVVi) # need this
                        
                        # for imputation if using matlab population matrix, which we won't...
                        sr.ages        <- rownames(srAPVVi)[row(srAPVVi)][!NAindi]
                        sr.yrs         <- colnames(srAPVVi)[col(srAPVVi)][!NAindi]
                        
                        pop[cbind(sr.ages, sr.yrs)]  <- srAPVVi[!NAindi]
                        pop[sragesc[i], this.yearc]  <- pop.this.yr.i
        }
        # now adjust c
        Sumt      <- sum(pop[as.integer(rownames(pop)) >= MaxAge, this.yearc], na.rm = TRUE)
        res[z]    <- abs(Sumt - SUM)
        # stop optimizing cmid if we've balanced pop age 90+
        if (res[z] < 0.1){
          break
        }
        # ---------------------------
        # TR: this is clever of Dima, I admit:
        if (Sumt > SUM){
          if (c_r == c_l){
            c_l <- c_l / 2
          } else {
            c_r <- cmid
          }
        } else {
          if (c_r == c_l){
            c_r <- c_r * 2
          } else {
            c_l <- cmid
          }
        }
        # ---------------------------
        
      } # 
      cat("\nOptimizing", ifelse(Sex == "m", "males", "females"), "old-age rate of mortality improvement\n")
      cat("\nc =", cmid, "\n")
      
      Psr        <- melt(pop, varnames = c("Age", "Year"), value.name = "Population")
      # otherwise this could extend left of the data...
      Psr        <- Psr[Psr$Year %in% YearsKeep, ]
      # ADD COHORT, then select only sr cohorts and years. Don't redo EC data
      Psr$Cohort <- Psr$Year - Psr$Age - 1
      Psr        <- Psr[with(Psr, Cohort %in% srcohorts & Year %in% sryears), ]
      Psr        <- Psr[!is.na(Psr$Population), ]
      # now we have things cut down.
      
      # Make selection IDs
      PsrID      <- paste(Psr$Year, Psr$Age, sep = "-")
      PsexID     <- paste(Psex$Year, Psex$Age, sep = "-")
      rmID       <- !PsexID %in% PsrID
      
      SRpop      <- as.data.frame(matrix(NA, 
                      nrow = nrow(Psr), 
                      ncol = length(ColnamesKeep),
                      dimnames = list(NULL, ColnamesKeep)))
      
      SRpop$Age         <- SRpop$Agei            <- as.integer(Psr$Age)
      SRpop$Year        <- as.integer(Psr$Year)
      SRpop$Population  <- Psr$Population
      SRpop$Sex         <- Sex
      SRpop$AgeInterval <- SRpop$AgeIntervali    <- 1
      SRpop$PopName     <- unique(Psex$PopName)
      SRpop$Day         <- 1
      SRpop$Month       <- 1
      SRpop$LDB         <- 1
      SRpop$Access      <- "O" # presumably any data we invent is open access, even if the origin data are not?
      SRpop$NoteCode1   <- "p_srm()" # might as well leave this as a note
      SRpop$Cohort      <- Psr$Cohort
      
      # this is an unexplained feature of the matlab implementation...
#      if (reproduce.matlab){
#        SRpop$Population[SRpop$Cohort == Cohmax & SRpop$Age == w] <- 
#                             SRpop$Population[SRpop$Cohort == Cohmax & SRpop$Age == w] + omega["Ds"]
#      }
      
      PopMF[[Sex]] <- rbind(Psex[rmID, ColnamesKeep],SRpop[,ColnamesKeep])
    } # end Sex loop
    if (UNKTF){
      PopMF[["UNK"]]    <- UNK[,ColnamesKeep]
    }
    Pop                 <- resortPops(do.call(rbind, PopMF))[, ColnamesKeep]
    rownames(Pop)       <- NULL
    invisible(Pop)
    
  }
)

#'
#' @title p_srm_updatePopGivenCmid an auxiliary function called by \code{p_srm()}.
#' 
#' @description This functionalizes an inner loop-- that used for the optimization of 'c', an approximate rate of old-age improvement in mortality. All these weird arguments are created in the body of \code{p_srm} and should be examined there for more detail. This function is separated and compiled to increase speed. Calls aux function \code{PC2AP()}.
#' 
#' @param pop an AP matrix of running population counts in the SR zone, going as low as age 75 for purposes of estimation.
#' @param VVcumsums a PC matrix of EC death counts for the SR cohorts, added sequentially over the running population estimates for the cohorts of the most recent year.
#' @param cmid the given estimate of the improvement coefficient, 'c'.
#' @param D_dot 5x1 VV death counts for the most recent cohorts just left of the SR zone (explained in-line in the MP)
#' @param srages integer vector of ages for the SR zone, used for shifting.
#' @param sragesc character vector of ages for the SR zone, used for selecting.
#' @param srcohortsc character vector of SR cohorts, for selecting.
#' @param this.yearc character, the most recent year in the Population data.
#' @param yr.denom character vector of years to be selected for the denominator of R*.
#' @param yr.num character vector of years to be selected for the numerator of R*.
#' @param k the parameter 'k' from various equations in the section on SRM.
#' @param m the parameter 'm' from various equations in the section on SRM.
#' @param w omega.
#' @param a lower age bound in fitting (not necessarily overwritten). Should be lower than \code{AgeMax}
#' 
#' @importFrom compiler cmpfun
#' 
#' @export 
#' 


p_srm_updatePopGivenCmid <- cmpfun(function(pop, 
    VVcumsums, 
    cmid, 
    D_dot, 
    srages, 
    sragesc, 
    srcohortsc, 
    this.yearc, 
    yr.denom, 
    yr.num, 
    k, 
    w,
    a){
    # This procedure has serial-dependance, must be done in a loop... bummer
    for (i in length(srages):1){ # i <- 22
      
      # first make R_star, the Ratio, given 'c'
      denom          <- sum(pop[as.character(srages[i] - k), yr.denom ], na.rm = TRUE)
      num            <- sum(pop[sragesc[i],  yr.num], na.rm = TRUE)
      R_star         <- num / denom
      R_star[is.nan(R_star) | is.infinite(R_star)] <- 0
      
      # impute current age-c-iteration of 'this.year' pop
      pop.this.yr.i  <- cmid * D_dot[i] * R_star / (1 - R_star)
      srAPVVi        <- PC2AP(VVcumsums[, srcohortsc[i], drop = FALSE] + pop.this.yr.i, 
        agemin = a, agemax = (w - 1), Lexis = 1) # have tried w-2,w-1,w, 
      NAindi         <- is.na(srAPVVi) # need this
      
      # for imputation if using matlab population matrix, which we won't...
      sr.ages        <- rownames(srAPVVi)[row(srAPVVi)][!NAindi]
      sr.yrs         <- colnames(srAPVVi)[col(srAPVVi)][!NAindi]
      
      pop[cbind(sr.ages, sr.yrs)]  <- srAPVVi[!NAindi]
      pop[sragesc[i], this.yearc]  <- pop.this.yr.i
    }
    pop
  })

