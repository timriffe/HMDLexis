

#'
#' @title p_sra a simple and flexible survivor ratio method function, for modularity.
#' 
#' @description This implementation seeks to follow the MP. LDB \code{Deaths} must be finalized. Population data is only needed to get an estimate for population in an artifical open age group, age 90\+ by default, though this may be lowered or raised. An argument \code{reproduce.matlab} controls whether an unnecessary matlab kludge is included. This function calls various other LexisDB functions at the moment, such as \code{d_addCohortColumn()}, \code{p_addCohortColumn()}, and \code{p_ecm_findOmega()}. The parameter \code{bordercoh} is \code{FALSE} by defauls, unless function called by \code{p_soai()}.
#' 
#' @param Pop The standard internal Population data.frame, *after* running p_ecm().
#' @param Deaths after all processing is done. Completed triangles.
#' @param k the parameter 'k' from various equations in the section on SRM.
#' @param m the parameter 'm' from various equations in the section on SRM.
#' @param l the parameter 'l' from various equations in the section on SRM.
#' @param a the lower bound for population counts, almost always 80 (as in EC), but could be changed
#' @param A controls whether we do SR 90+ or SR 85+ (or something else). Default 90.
#' @param maxit maximum number of iterations to optimize the improvement coefficient, 'c'.
#' @param reproduce.matlab logical. Do we include the legacy matlab kludge? Default \code{FALSE}. This effects results very little, and is only used for exact matching of output.
#' @param bordercoh logical. Should be keep the border cohort between SR and EC zones?
#' 
#' @return Pop standard object, with recent old ages adjusted using the survivor ratio method
#' 
#' @note This very function will be generalized to account for Tadj eventually. Further testing is needed to check whether the parameter \code{A} adjusts things properly for use with the SR85+ method.
#' 
#' @details Years all data present in SR zone could be kept, as long as it is contained within the range of years present in \code{Pop}. If \code{Pop} is specified as 

#' @importFrom reshape2 acast
#' @importFrom reshape2 melt
#' @importFrom compiler cmpfun
#' 
#' @export 
#' 
#  tail(Pop)
# k = 5;m = 5;l=5;a=80;A = 90;maxit = 100;reproduce.matlab = TRUE
# A <- 85
# Pop <- P1970; 
# Deaths <- Deaths[Deaths$Year<=1969,]; library(reshape2)


p_sra <- function(Pop, 
    Deaths, 
    k = 5, 
    m = 5, 
    l = 5,
    a = 80,
    A = 90, 
    maxit = 100,
    reproduce.matlab = FALSE,
    bordercoh = FALSE){
    
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
    
    # added for when SRA used in years other than most recent year.
    Deaths <- Deaths[Deaths$Year < this.year, ]
    
    YearsFit        <- min(Deaths$Year):this.year
    YearsKeep       <- min(Pop$Year):this.year
    # to be filled later
    PopMF           <- list()
    # sex loop Sex <- "m"; Sex <- "f"
    # allows for single or both sexes
    sexes <- unique(Pop$Sex)
    sexes <- sexes[sexes %in% c("m", "f")]
    for (Sex in sexes){
      
      Psex            <- Pop[Pop$Sex == Sex, ]
      SUM             <- sum(with(Psex,Population[Year == this.year & Agei >= A]))
      Psex            <- Psex[Psex$AgeInterval != "+", ]
      Dsex            <- Deaths[Deaths$Sex == Sex, ]
      
      # omega is typically different for males and females, must be in the loop:
     
      t.l <- l ;
      omega           <- p_ecm_findOmega(Dsex, l = t.l, threshold = 0.5)
      #
      w               <- omega["omega"] 
     
      # some helpful vectors for selecting and iterating:
      srages          <- A:(w - 1)                       # have tried w, w-1,w-2
      srcohorts       <- this.year - 1 - srages               # this is exact
      sryears         <- (this.year - (w - a - 1)):(this.year-1)  # this is not particularly key
      # in case of p_soai(), the left side of this may need to tuck in.
      sryears         <- sryears[sryears >= min(YearsFit)]
      
      # for name-based matrix indexing:
      sragesc         <- as.character(srages)
      srcohortsc      <- as.character(srcohorts)
      sryearsc        <- as.character(sryears)
      
      # more indexing, recycled
      yr.denom        <- as.character((this.year - 1 - k):(this.year - m - k))
      yr.num          <- as.character((this.year - 1):(this.year - m))
      
      # remove population estimates that are located in the SR zone:
      Psex            <- Psex[!(Psex$Cohort %in% srcohorts & Psex$Year %in% sryears & Psex$Agei >= a), ] 
      # cleans out SR zone
      
      # ------------------------------------
      # make D_dot, a recycled vector.
      # aggregate 5x1 cohorts counting back from most recent year. 
      # create early to remove redundancy in repeated testing
      D_dot           <- acast(Dsex[with(Dsex, Year > (this.year - 6) & 
                                               #Year <= this.year & # this condition makes no difference for right-side
                                                                   # SRA application, but allows us flexibility for
                                                                   # p_soai() [mid-surface operations]
                                               Cohort > omega["Cohmax"]), ], 
                              . ~ Cohort, sum, value.var = "Deaths", fill = 0)[, , drop = TRUE][srcohortsc]
     #sort(unique(Dsex$Cohort))
      # VVcumsums are the same as EC, but for the SR region only. Recycled in later loop 
      # as a constant.
      VV              <- acast(Dsex[with(Dsex, Year %in% sryears & Cohort %in% srcohorts), ], 
                               Year ~ Cohort, 
                              sum, 
                              value.var = "Deaths")
      
      # VVcumsums is like ECpop, except for non-extinct cohorts. We just to need
      # to add in year t populations for ages maxAge+ to add backwards from...
      VVcumsums <- apply(VV, 2, function(x){
          rev(cumsum(rev(x)))
        })
      
      # pop is a generic holder, in AP format (the ratios follow A and P...)
      #yrskeep         <- as.character((min(sryears) - k):this.year)
      ageskeep        <- as.character((a - 10):130)
      pop             <- matrix(nrow = length(ageskeep), 
                                ncol = length(YearsFit),
                                dimnames = list(ageskeep,
                                  YearsFit)
                                )
      # get EC pops for SR zone omega
      #
      ECcohs <- omega["Cohmax"]:(omega["Cohmax"] - 12)
      
      # sum deaths in PC cells (VV)
      ECVV   <- acast(Dsex[Dsex$Cohort %in% ECcohs, ], 
                      Year ~ Cohort, sum, value.var = "Deaths")
      # add vector of 0s for most recent year 
      ECVV   <- rbind(ECVV,0)  
      rownames(ECVV)[length(rownames(ECVV))] <- this.yearc
      
      # matlab implementations uses the fitted value, Ds, which is an average from a much larger
      # triangle, to start off the most recent extinct cohort... I think this is odd, and actually
      # this artifact makes the cohort neither EC nor SR, but something else (more similar to EC).
      if (reproduce.matlab){
        ECVV[nrow(ECVV), ncol(ECVV)] <- omega["Ds"]
      }
      
      # the backward cumulation of these are the extinct cohort populations, 
      # used to start off estimation of known survivor ratios...
      ECpop   <- apply(ECVV, 2, function(x){
                           rev(cumsum(rev(x))) 
                         })
      # move to AP format, to place into population results container, pop
      ECpopAP <- PC2AP(ECpop)
      
      rn      <- intersect(rownames(ECpopAP),rownames(pop))
      cn      <- intersect(colnames(ECpopAP),colnames(pop))
      pop[rn, cn]       <-   ECpopAP[rn, cn]       
   
      # -----------------------------------------------------------
      # starting values give as a c of 5, oddly high, although this is
      # a multiplier, Inf would be the counterpart to zero...
      c_l             <- 0
      c_r             <- 10
 
      cmid <- res             <- rep(0, maxit)
      for (z in 1:maxit){  
        cmid[z]    <- (c_l + c_r) / 2
        pop.this.yr <- R_star <- num <- denom <- rep(0,length(srages))
        for (i in length(srages):1){ # i <- 19 ; i <- 11                     
          # first make R_star, the Ratio, given 'c'
   
          denom[i]          <- sum(pop[as.character(srages[i] - k), yr.denom ], na.rm = TRUE)
          num[i]            <- sum(pop[sragesc[i],  yr.num], na.rm = TRUE)
          R_star[i]         <- num[i] / denom[i]
          R_star[is.nan(R_star) | is.infinite(R_star)] <- 0
          
          # impute current age-c-iteration of 'this.year' pop
          pop.this.yr[i]    <- cmid[z] * D_dot[i] * R_star[i] / (1 - R_star[i])
          srAPVVi           <- PC2AP(VVcumsums[, srcohortsc[i], drop = FALSE] + pop.this.yr[i] , 
                                     agemin = a, 
                                     agemax = w, 
                                     Lexis = 1) # have tried w-2,w-1,w, 
          NAindi            <- is.na(srAPVVi) # need this
        
          # for imputation if using AP matrix
          sr.ages           <- rownames(srAPVVi)[row(srAPVVi)][!NAindi]
          sr.yrs            <- colnames(srAPVVi)[col(srAPVVi)][!NAindi]
          
          pop[cbind(sr.ages, sr.yrs)]  <- srAPVVi[!NAindi]
          pop[sragesc[i], this.yearc]  <- pop.this.yr[i] 
        }
                      
        # now adjust c
        Sumt      <- sum(pop[as.integer(rownames(pop)) >= A, this.yearc], na.rm = TRUE)
        res[z]    <- abs(Sumt - SUM)
        # stop optimizing cmid if we've balanced pop age 90+
        if (res[z] < 0.1){ # this res copied from matlab
          break
        }
        # ---------------------------
        # adjust bounds
        if (Sumt > SUM){
          c_r <- cmid[z]
        } 
        if (Sumt < SUM){
          c_l <- cmid[z]
        }    
      } # 
      
      # now reshape and filter the new population estimates
      Psr        <- melt(pop, varnames = c("Age", "Year"), value.name = "Population")
      # otherwise this could extend left of the data...
      Psr        <- Psr[Psr$Year %in% YearsKeep, ]
      # ADD COHORT, then select only sr cohorts and years. Don't redo EC data
      Psr$Cohort <- Psr$Year - Psr$Age - 1
      # this EXCLUDES the border cohort, started off with Ds value
      if(bordercoh){
        Psr        <- Psr[Psr$Cohort %in% c(srcohorts, omega["Cohmax"]) & Psr$Age >= a, ] 
      } else {
        Psr        <- Psr[Psr$Cohort %in% srcohorts & Psr$Age >= a, ] 
      }
    
     
      # this likely won't happen
      Psr$Population[is.na(Psr$Population)] <- 0
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
      SRpop$Access      <- "O"       # presumably any data we invent is open access, even if the origin data are not?
      SRpop             <- assignNoteCode(SRpop, "p_sra()") 
      SRpop$Cohort      <- Psr$Cohort
      
      # TR: added to carry over Area assignment in newly created data
      # this might be able to stay here even after Tadj has been incorporated
      SRpop             <- assignArea(SRpop, Psex)

      PopMF[[Sex]]      <- rbind(Psex[rmID, ColnamesKeep],SRpop[,ColnamesKeep])
    } # end Sex loop
    if (UNKTF){
      PopMF[["UNK"]]    <- UNK[,ColnamesKeep]
    }
    Pop                 <- resortPops(do.call(rbind, PopMF))[, ColnamesKeep]
    rownames(Pop)       <- NULL
    invisible(Pop)
    
  }













