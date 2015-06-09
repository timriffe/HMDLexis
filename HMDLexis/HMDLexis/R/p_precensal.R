
#' @title a function for precensal population estimation
#' 
#' @description This function is in principle simple, because we assume no error, except possibly in the case of newborn cohorts. Whether or not we redistribute error over newborn cohorts in the precensal period is toggled with the \code{MPVERSION} argument. At the time of this writing, \code{MPVERSION >= 7} will distribute error over these cohorts, otherwise, we assume no error. This function should be called near the end of population estimation in the country script. All age groups in \code{Pop} must be single (no open ages or groups). So, typically, one would already have run \code{p_srecm()} or similar prior to calling this function. A little triangle of data in the top left lexis region may be left empty, and this can be filled by re-calling \code{p_srecm()}, for example. This function calls a few other LexisDB functions, such as \code{d_addCohortColumn()}, \code{p_addCohortColumn()}, \code{ypart()}, and \code{resortPops()}.
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


### TODO: We have a dimension problem if we're doing a 1-year precensal.
p_precensal <- function(Pop, Deaths, Births, MPVERSION = 5, reproduce.matlab = FALSE){
  yr1    <- min(Deaths$Year) 
  yr2    <- min(Pop$Year)
  C2     <- Pop[Pop$Year == yr2, ]
  Deaths <- Deaths[Deaths$Year <= yr2, ]
  Deaths <- d_addCohortColumn(Deaths)
  
  stopifnot(all(C2$AgeIntervali == 1))
  # note Pop date might not be Jan 1. Could be mid year.
  
  years  <- yr1:yr2 # possibly the source of the problem.
  N      <- length(years)
  
  PopOut   <- list()
  PopOut[["RightSide"]] <- Pop[Pop$Year > yr2, ]
  # 
  for (Sex in c("f","m")){ # Sex <- "f"
    C2s  <- C2[C2$Sex == Sex, ]
    Dsex <- Deaths[Deaths$Sex == Sex, ]
    # get f2 only (no f1 in precensals)
    if (!reproduce.matlab){
      f2 <- ypart(Year = unique(C2s$Year), 
        Month = unique(C2s$Month), 
        Day = unique(C2s$Day), 
        reproduce.matlab = reproduce.matlab)
      
      if (length(f2) > 1){
        stop("multiple dates in use in the same year, makes intercensals tricky.\\Time to do some digging.")
      }
    } else {
      f2 <- ypart(Year = C2s$Year[1], 
        Month = C2s$Month[1], 
        Day = C2s$Day[1], 
        reproduce.matlab = reproduce.matlab)
    }
    
    # adjust right-side population if necessary. 
    # Actually to have run ecm() we must have Jan 1 in pops below age a
    C2s <- p_addCohortColumn(C2s)
    
    # this may be problematic (lose top observation?
    if (f2 > 0){
      C2s$Population[-1] <- (1 - f2) * C2s$Population[-nrow(C2s)] + f2 * C2s$Population[-1]
      C2s$Cohort <- C2s$Cohort + 1
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
    Dc      <- f2^2 * DL[nrow(DL), ]
    Dd      <- ((2 * f2) - f2^2) * DU[nrow(DU), ]
    
    # the cohorts present in the first pop year, up to age 130 by default (should be param)
    Ncoh    <- yr1:(yr2 - 1)
    Ccoh    <- C2s$Cohort[!C2s$Cohort %in% Ncoh]
    

    CDc              <- Dc[as.character(Ccoh)]
    CDd              <- Dd[as.character(Ccoh)]
    CDc[is.na(CDc)]  <- 0

    # CVV can now be used for cumulative deaths
    CVV              <- VV[, as.character(Ccoh), drop = FALSE]
    CVV[is.na(CVV)]  <- 0
    CVV[nrow(CVV), ] <- CDc + CDd
    
    CVVcumsums <- apply(CVV,2,function(x){
        rev(cumsum(rev(x)))
      })
    # get simple starting and ending vectors of population
    C2vec   <- C2s$Population[C2s$Cohort %in% Ccoh]
    names(C2vec) <- C2s$Cohort[C2s$Cohort %in% Ccoh]
    #

    # this is complete, but lacking cohorts for whom we have births
    PPC   <- C2vec + t(CVVcumsums)[, names(C2vec)]
    
    ####################################
    # now newborn cohorts:             #
    ####################################
    # these go from C2 backward, so we want to rev them
    N2      <- rev(C2s$Population[C2s$Cohort %in% Ncoh])
    names(N2) <- Ncoh
    
    # right-side adjustment
    NDc     <- Dc[as.character(Ncoh)]
    NDd     <- Dd[as.character(Ncoh)]
    NDd[is.na(NDd)] <- 0
    # fragile verify this 
    NVV     <- VV[, as.character(Ncoh), drop = FALSE]
    NVV[nrow(NVV), ] <- NDc + NDd
    
    NVV[is.na(NVV)] <- 0
    
    # this is the TR-proposed amendment, TBD for MPv7
    if (MPVERSION >= 7){
      
      # TR: code copied and pasted, with minor modification, from p_ic_inner().
      # this redistributes error for infants, as we do in intercensal estimation.
      
      # Births enter into LDB sorted, so no worries.
      NB        <- Births$Births[Births$Year %in% Ncoh & Births$Sex == Sex]
      names(NB) <- Ncoh
      # these go from C2 backward, so we want to rev them
      
      # estimated pop size, assuming only death decrement
      N2hat     <- NB - colSums(NVV, na.rm = TRUE)
      # total error on arrival at C2
      NDelta    <- N2 - N2hat
      
      # getting prop to distribute error over, tricky
      # remove last row of NVV (note it's PC data), because the row labels refer to cohort of deaths,
      # but what we really want is Jan 1 of year t, not the last year, though, which could stick out past our data.
      NVV       <- NVV[-nrow(NVV), , drop = FALSE]
      Nprop     <- NVV * 0
      # K, k follow MP naming conventions
      for (coh in 1:length(Ncoh)){ # like in MP (note our N might be different)
        K            <- N - coh - 1
        k            <- 0:K
        prop         <- (2 * k + 1) / (2 * K  + 1 + 2 * f2)
        Nprop[coh, ] <- c(rep(0,coh-1), prop)
      }
      # Total error
      NDelta  <- Nprop * NDelta
      
      # cumulative deaths, for decrement *from* births
      NVVcum              <- t(apply(NVV, 2, cumsum))
      
      # estimate population over period
      Npop    <- NB - NVVcum + NDelta
      Npop[Nprop == 0] <- NA
      
      colnames(Npop) <- as.integer(colnames(Npop)) + 1
    }
    # This is the matlab error-free assumption
    if (MPVERSION < 7){
      # generate the same output assuming no error.
      # this was default behavior until TR noticed this.
      
      # we're always dealing with something close to a triangle in this zone, so this works
      NVV <- NVV * !upper.tri(NVV, TRUE)
      NVV <- NVV[-1, , drop = FALSE]
      
      # note the other version is a simple cumsum, but here we need 
      # rcumusum() because it's an increment from the right
      NCPcumsum <- t(apply(NVV,2,function(x){
          rev(cumsum(rev(x)))
        }))
      # it's this easy because we have no error...
      Npop <- N2 + NCPcumsum
      Npop[lower.tri(Npop)] <- NA
 
    }

      
      Ps     <- rbind(
        melt(PPC, value.name = "Population", varnames = c("Cohort", "Year")),
        melt(Npop, value.name = "Population", varnames = c("Cohort", "Year"))
      )
      Ps             <- Ps[!is.na(Ps$Population), ]
    
      # add on other columns
      # need AP in the end...
      Ps$Agei        <- Ps$Age <- Ps$Year - Ps$Cohort - 1
      Ps             <- Ps[Ps$Age <= 130, ]
      
      # now add on the remaining columns...
      Ps$Sex         <- Sex
      Ps$PopName     <- unique(C2s$PopName)
      Ps$Area        <- NA
      Ps$LDB         <- Ps$Month <- Ps$Day <- Ps$AgeInterval <- Ps$AgeIntervali <- 1
      Ps$Type        <- "precensal"
      Ps$Access      <- "O"
      Ps$NoteCode1   <- "p_precensal()"
      Ps$RefCode     <- Ps$NoteCode2 <- Ps$NoteCode3 <- NA
      Ps             <- Ps[, colnames(C2s)]
      Ps$Cohort      <- NULL
      PopOut[[Sex]] <- Ps
      
    }
    PopOut        <- do.call(rbind, PopOut)

    PopOut        <- resortPops(PopOut)
    invisible(PopOut)
  }














