
#' @title a function for postcensal population estimation
#' 
#' @description This function is in principle simple, because we assume no error. This function should be called prior to \code{p_sr()} or \code{p_srecm()} in all cases. All age groups in \code{Pop} must be single (no open ages or groups). This function calls a few other LexisDB functions, such as \code{d_addCohortColumn()}, \code{p_addCohortColumn()}, \code{ypart()}, and \code{resortPops()}.
#'
#' @param Pop the standard Pop object, ideally with open age split, e.g., after running \code{p_srecm()} on the data.
#' @param Deaths the standard Deaths object, after all processing is completed.
#' @param Births the standard Births object, only needs to be specified in \code{MPVERSION >= 7}.
#' @param MPVERSION numeric, default value is 5. This only toggles error distribution for newborn cohorts, as in description above.
#' @param reproduce.matlab logical, defaults to \code{FALSE}. This only affects date handling in the calculation of MP variable \code{f2}.
#'
#' @export 
#' 
#' @importFrom reshape2 acast
#' @importFrom reshape2 melt
#'

p_postcensal <- function(Pop, Deaths, Births, MPVERSION = 5, reproduce.matlab = FALSE){
  yr2    <- max(Deaths$Year) + 1 # because we shoot for Jan 1
  yr1    <- max(Pop$Year) 
  C1     <- Pop[Pop$Year == yr1, ]
  Deaths <- Deaths[Deaths$Year >= yr1, ]
  Deaths <- d_addCohortColumn(Deaths)
  C1     <- C1[C1$AgeInterval != "+", ]
  
  # ergo, we can't have UNK, do this earlier in processing.
  stopifnot(all(C1$AgeIntervali == 1))
  # note Pop date might not be Jan 1. Could be mid year.
  # even so, yr2 would be Jan 1 of NEXT year, so we're good here
  years  <- yr1:(yr2 - 1)
  N      <- length(years)
  
  PopOut   <- list()
  # i.e., if the currently right-most population counts are mid-year, we still keep them
  # in this function unique(Pop$Year)
  PopOut[["LeftSide"]] <- Pop[Pop$Year < yr2, ]
  # 
  for (Sex in c("f","m")){ # Sex <- "m"
    C1s  <- C1[C1$Sex == Sex, ]
    Dsex <- Deaths[Deaths$Sex == Sex, ]
    
    # get f1 only (no f2 in postcensals)
    if (!reproduce.matlab){
      f1 <- ypart(Year = unique(C1s$Year), 
        Month = unique(C1s$Month), 
        Day = unique(C1s$Day), 
        reproduce.matlab = reproduce.matlab,
        detect.mid.year = TRUE,
        detect.start.end = TRUE)
      
      if (length(f1) > 1){
        stop("multiple dates in use in the same year, makes intercensals tricky.\\Time to do some digging.")
      }
    } else {
      f1 <- ypart(Year = C1s$Year[1], 
        Month = C1s$Month[1], 
        Day = C1s$Day[1], 
        reproduce.matlab = reproduce.matlab,
        detect.mid.year = TRUE,
        detect.start.end = TRUE)
    }
    
    # adjust right-side population if necessary. 
    # Actually to have run ecm() we must have Jan 1 in pops below age a
    C1s <- p_addCohortColumn(C1s)
    
    # this may be problematic (lose top observation?
    if (f1 > 0){
      
      C1s$Population[-1] <- (1 - f1) * C1s$Population[-nrow(C1s)] + f1 * C1s$Population[-1]
      C1s$Cohort <- C1s$Cohort + 1
    }
    
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
    
    # right-side adjustment for partial years
    Da      <- (1 - f1^2) * DL[1, ]
    Db      <- (1 - f1)^2 * DU[1, ]
    
    # the cohorts present in the first pop year, up to age 130 by default (should be param)
    # infant cohort (I_):
    NewCohorts   <- yr1:(yr2-1)
    # infants treated separately
    Icoh         <- NewCohorts[1]
    # newborn cohorts (N_):
    if (length(NewCohorts[-1]) > 0){
      Ncoh <- NewCohorts[-1]
    } else {
      Ncoh <- NULL
    }
    ####################################################################
    # Complete cohorts                                                 #
    ####################################################################
    Ccoh    <- C1s$Cohort[!C1s$Cohort %in% c(Icoh, Ncoh)]
    
    CDa              <- Da[as.character(Ccoh)]
    CDb              <- Db[as.character(Ccoh)]
    CDb[is.na(CDb)]  <- 0
    
    # CVV can now be used for cumulative deaths
    CVV              <- VV[, as.character(Ccoh), drop = FALSE]
    CVV[is.na(CVV)]  <- 0
    
    # cumulative deaths
    CVVcum           <- apply(CVV, 2, cumsum)
    
    # if we're only projecting a year, then we might lose dimensions in apply...
    dim(CVVcum)      <- dim(CVV)
    dimnames(CVVcum) <- dimnames(CVV)
    # get simple starting and ending vectors of population
    C1vec            <- C1s$Population[C1s$Cohort %in% Ccoh]
    names(C1vec)     <- C1s$Cohort[C1s$Cohort %in% Ccoh]
    # assuming no error / migration, we get:
    # this is complete, but lacking cohorts for which we have births
    PPC              <- C1vec - t(CVVcum)[names(C1vec), ]
    dim(PPC)         <- dim(t(CVVcum)[names(C1vec), , drop = FALSE])
    dimnames(PPC)    <- dimnames(t(CVVcum)[names(C1vec), , drop = FALSE])
    
    ####################################################################
    # Infant cohort                                                    #
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
    
    # get deaths over cohort, and accumulate
    IVV     <- VV[, as.character(Icoh)]
    IVV[1]  <- IDa
    
    IVVcum  <- cumsum(IVV)
    # estimated size of cohort on Jan 1 on way to Jan 1, yr2
    Ipop    <- I0 - IVVcum
    
    # put in matrix row form to easily combine later
    Ipop    <- matrix(Ipop, nrow=1, dimnames=list(Icoh,years),byrow=TRUE)
    
    ####################################################################
    # Newborn cohorts                                                  #
    ####################################################################
    # TR: in the case of a mid-year census in year t, where we only need precensals
    # back to Jan 1 of year t, then we skip this part:
    if (!is.null(Ncoh)){
      # Births enter into LDB sorted, so no worries.
      NB      <- Births$Births[Births$Year %in% Ncoh & Births$Sex == Sex]
      # these go from C2 backward, so we want to rev them
      
      # fragile verify this 
      NVV     <- VV[, as.character(Ncoh), drop = FALSE]     
      NVV[is.na(NVV)] <- 0
      
      # cumualtive deaths
      NVVcum              <- t(apply(NVV, 2, cumsum))
      dim(NVVcum)         <- dim(NVV)
      dimnames(NVVcum)    <- dimnames(NVV)
      # estimate population over period
      Npop    <- NB - NVVcum 
      Npop[NVVcum == 0] <- NA
      
      colnames(Npop) <- as.integer(colnames(Npop)) + 1
    } else {
      # TR: just to make the code keep chugging
      Npop <- PPC * NA
    }
    
    Ps     <- rbind(
      melt(PPC, value.name = "Population", varnames = c("Cohort", "Year")),
      melt(Npop, value.name = "Population", varnames = c("Cohort", "Year")),
      melt(Ipop, value.name = "Population", varnames = c("Cohort", "Year"))
    )
    Ps             <- Ps[!is.na(Ps$Population), ]
    
    # add on other columns
    # need AP in the end...
    Ps$Year        <- Ps$Year + 1 #
    Ps$Agei        <- Ps$Age <- Ps$Year - Ps$Cohort - 1
    Ps             <- Ps[Ps$Age <= 130, ]
    
    # now add on the remaining columns...
    Ps$Sex         <- Sex
    Ps$PopName     <- unique(C1s$PopName)
    Ps$Area        <- unique(C1s$Area)
    Ps$LDB         <- Ps$Month <- Ps$Day <- Ps$AgeInterval <- Ps$AgeIntervali <- 1
    Ps$Type        <- "postcensal"
    Ps$Access      <- "O"
    Ps$NoteCode1   <- "p_postcensal()"
    Ps$RefCode     <- Ps$NoteCode2 <- Ps$NoteCode3 <- NA
    Ps             <- Ps[, colnames(C1s)]
    Ps$Cohort      <- NULL
    
    # TR: added 24.05.2016
    # negative numbers are possible in high ages. I'm going to make an executive decision
    # here. If it's above age 111, we can replace negatives with 0,
    # but if a negative occurs below age 111 then we'll abort with predjudice and
    # an informative warning message.
    
    if (any(Ps$Population < 0)){
      AgesNeg <- Ps$Agei[Ps$Population < 0]
      if (!all(AgesNeg >= 111)){
        cat("\npostcensal estimates have some negative population counts in ages under 115
We therefore abort this function. Either there was a negative pop count you need to deal with, 
or there are deaths being deducted from a zero or small pop count and dropping
below zero. If that is the case, and this is a minor data quality issue, then deal with
it prior to calling this function.")
cat(AgesNeg)
       return(Pop)
      }
     
      Ps$Population[Ps$Population < 0] <- 0
    }
    
    
    PopOut[[Sex]]  <- Ps
    
  }
  PopOut        <- do.call(rbind, PopOut)
  

  
  
  
  PopOut        <- resortPops(PopOut)
  invisible(PopOut)
}




