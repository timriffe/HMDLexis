#'
#' @title p_ic_prexisting connects cohorts between two censuses, distributes error.
#' 
#' @description Census 1 and census 2 must already be Jan 1st. This function does not account for newly born cohorts, nor will it extend out to age beyond the open age in census 2. Territorial adjustments are not accounted for in this implementation. Censuses must be in single ages, UNK already redistributed.
#' 
#' @param Pop the standard HMD internal population data.frame. After unknowns are redistributed.
#' @param Deaths the standard HMD internal deaths data.frame, after all processing completed.
#' 
#' @return Pop with new lines for intermediate years.
#' 
#' @importFrom reshape2 acast
#' @importFrom reshape2 melt
#' @importFrom compiler cmpfun
#' 
#' @export 
#' 

p_ic_prexisting <- cmpfun(function(Pop, Deaths, reproduce.matlab = TRUE){
  
  ColKeep     <- colnames(Pop)
  # never keep TOT, just a useless filter. could be dealt with in readInputDB
  Pop         <- Pop[Pop$Age != "TOT", ]
  Pop         <- p_ey2ny(Pop)
  # UNK needs to be dealt with first.
  # this is a safety so we don't mess that up
  UNKTF       <- any(Pop$Age == "UNK")
  if (UNKTF){
    cat("p_unk() was necessary. You should do this first outside this function!\n")
    Pop       <- p_unk(Pop)
  }
  # slice off open ages, rbind back on later:
  OPTF        <- any(Pop$AgeInterval == "+")
  if (OPTF){
    OP        <- Pop[Pop$AgeInterval == "+", ]
    Pop       <- Pop[Pop$AgeInterval != "+", ]
  }
  # add Cohort Column to Pop and Deaths, used cohorts to select and organize
  Deaths      <- d_addCohortColumn(Deaths)
  Pop         <- p_addCohortColumn(Pop)
 
  # -------------------------------------------------
  # an object to catch newly created data
  PopOut      <- list()
  # outer loop over sex .Sex <- "f"; i <- 1
  for (.Sex in c("f", "m")){
    Psex      <- Pop[Pop$Sex == .Sex, ]
    Dsex      <- Deaths[Deaths$Sex == .Sex, ]
    
    # OK, detect censuses to connect
      
    Allyears  <- min(Psex$Year):max(Psex$Year)
    # we don't want to base intercensals on other intercensals, so we remove them 
    # during estimation
    FunkyDetectionStrategy <- diff(Allyears %in% Psex$Year[Psex$Type != "HMDintercensal"])
    # our 'censuses' we do not care about the Type (for now). 
    # Intercensals can be based on other estimates too..
    yrs       <- union(Allyears[which(FunkyDetectionStrategy == -1)],
                         Allyears[which(FunkyDetectionStrategy == 1) + 1])
    # cut down Pops to just these
    Cpop      <- Psex[Psex$Year %in% yrs & Psex$Type != "HMDintercensal", ]
   
    # inner loop over census gaps
    Censuses  <- list()
    for (i in 1:(length(yrs) - 1)){
      # pick out census 1 and census 2 i<-4
      C1      <- Cpop[Cpop$Year == yrs[i], ]
      C2      <- Cpop[Cpop$Year == yrs[i + 1], ]
      
      # fractions of years are being worked in little by little f2 <- 0
      # get census date fractions. If there is more than one date in a given year,
      # e.g. composite sources, this *ought* to break. Matlab bases everything on the 
      # Date in row 1...
      if (reproduce.matlab){
        f1 <- ypart(Year = unique(C1$Year), 
          Month = unique(C1$Month), 
          Day = unique(C1$Day), 
          reproduce.matlab = reproduce.matlab)
        f2 <- ypart(Year = unique(C2$Year), 
          Month = unique(C2$Month), 
          Day = unique(C2$Day), 
          reproduce.matlab = reproduce.matlab)
      } else {
        f1 <- ypart(Year = C1$Year[1], 
          Month = C1$Month[1], 
          Day = C1$Day[1], 
          reproduce.matlab = reproduce.matlab)
        f2 <- ypart(Year = C2$Year[1], 
          Month = C2$Month[1], 
          Day = C2$Day[1], 
          reproduce.matlab = reproduce.matlab)
        if (length(f1) > 1 | length(f2) > 1){
          stop("multiple dates in use in the same year, makes intercensals tricky.\\Time to do some digging.")
        }
      }
     
      # adjust census1 population counts for Jan 1st.
      Nages <- nrow(C1)
      C1$Population[2:Nages] <- (1 - f1) * C1$Population[1:(Nages - 1)] + f1 * C1$Population[2:Nages]
      # TODO: what happens to age 0 if C1 is midyear? do we throw it out? 
      #       do we instead count up from births, bypassing C1 all the way to C2?
      #       that's what I'll do here
      if (f1 > 0){
        C1 <- C1[C1$Aegi > 0, ]
      }
          
      # adjust census population counts for Jan 1st.
      Nages   <- nrow(C2)
      C2$Population[2:Nages] <- (1 - f2) * C2$Population[1:(Nages - 1)] + f2 * C2$Population[2:Nages]
      # young people in C2 are bypassed anyway here, no need to adjust
          
      # which cohorts can we connect?
      CohKeep <- intersect(C1$Cohort, C2$Cohort)
        
      # yr x cohort matrices
      TL      <- acast(Dsex[with(Dsex, Year >= yrs[i] & Year <= yrs[i + 1] & Cohort %in% CohKeep & Lexis == "TL"), ], 
                                Year ~ Cohort, sum, value.var = "Deaths", fill = NA_real_)
      TU      <- acast(Dsex[with(Dsex, Year >= yrs[i] & Year <= yrs[i + 1] & Cohort %in% CohKeep & Lexis == "TU"), ], 
                                Year ~ Cohort, sum, value.var = "Deaths", fill = NA_real_)   
      # if f1 = 0 = Jan 1st, this just multiplies by 1. Otherwise we discount
      TL[1, ]         <- (1 - f1^2) * TL[1, ] # Da eq 23
      TU[1, ]         <- (1 - f1)^2 * TU[1, ] # Db eq 24
          
      # likewise for the time of Census 2
      TL[nrow(TL), ]  <- f2^2 * TL[nrow(TL), ] # Dc eq 25
      TU[nrow(TU), ]  <- (2 * f2 - f2^2) * TU[nrow(TU), ] # Dd eq 26
          
      VV              <- TL + TU
   
      # cut down
      C1      <- C1[C1$Cohort %in% CohKeep, ]
      C2      <- C2[C2$Cohort %in% CohKeep, ]
      # inserts NAs for missing combos rather than 0s. More robust.
      #VV      <- acast(Dsex[with(Dsex, Year >= yrs[i] & Year < yrs[i + 1] & Cohort %in% CohKeep), ], 
      #                 Year ~ Cohort, sum, value.var = "Deaths", fill = NA_real_)
      # returns a matrix, which we reorient to match C1 and C2
      VVc     <- t(apply(VV, 2, cumsum))
      VVc     <- VVc[nrow(VVc):1, ]
      # check if you want:
      #all(C1$Cohort == rownames(VVc))
      #all(C2$Cohort == rownames(VVc))
      # ---------------------------------------
      # taking advantage of R's idiosyncractic behavior
      # matrix of total error, which we then distribute equally by dividing out the interval (ncol())
      
      # ------------------------------
      # I guess this is what you'd want to mess with Sebastian:
      # ------------------------------
      Error                 <- VVc * 0 + C1$Population - VVc[,ncol(VVc)] - C2$Population
      # avg migs / deaths unaccounted for per year

      # partition error using f1 and f2
      parts                 <- VVc * 0 + 1
      parts[,1]             <- 1 - f1
      parts[,ncol(parts)]   <- f2
      prop                  <- parts / rowSums(parts)
     
      Error                 <- Error * prop
      # now take the cumumulative sum of errors along each cohort
      CumError              <- t(apply(Error, 1, cumsum))
      # first census, minus cumulative deaths minus cumulative error, voila
      NewPops               <- C1$Population - VVc - CumError

      NewPops               <- NewPops[, 1:(ncol(NewPops) - 1)]
      # relabel Years in colnames
      colnames(NewPops)     <- as.integer(colnames(NewPops)) + 1
      
      # now NewPops will in no case include the year of C1, but it will include the year of
      # C2, sometimes redundantly with the census that was used. That same census will still
      # be used for posterior intercensals and newly born cohorts, but we will overwrite it
      # eventually, even if it was a Jan 1st census.
      
      # 0 remaining residual at census 2:
      NewPops               <- melt(NewPops, varnames = c("Cohort", "Year"), value.name = "Population")
      NewPops$Age           <- NewPops$Year - NewPops$Cohort - 1
      
      # we don't keep the Cohort column at this time
      # now prepare in standard format
      PopOutSx              <- as.data.frame(
                                matrix(nrow=nrow(NewPops), 
                                      ncol = ncol(Cpop), 
                                      dimnames = list(NULL, colnames(Cpop))))
      PopOutSx$Age         <- PopOutSx$Agei          <- NewPops$Age
      PopOutSx$Year        <- NewPops$Year
      PopOutSx$Population  <- NewPops$Population
      PopOutSx$Sex         <- .Sex
      PopOutSx$AgeInterval <- PopOutSx$AgeIntervali    <- 1
      PopOutSx$PopName     <- unique(Cpop$PopName)
      PopOutSx$Day         <- 1
      PopOutSx$Month       <- 1
      PopOutSx$LDB         <- 1
      PopOutSx$NoteCode1   <- "p_ic_prexisting()"
      PopOutSx$Type        <- "HMDintercensal"
      # what about Area?
      Censuses[[i]]         <- PopOutSx[, ColKeep] 
    }
    # combines intercensal estimates into a single long df for the given sex, puts into output list
    PopOut[[.Sex]]  <-  do.call(rbind, Censuses)
  }
  # of there was Open Age data, put it back on.
  if (OPTF){
    PopOut[["OP"]]  <- OP
  }
    
  # this just gathers the new estimates
  PopNew  <- do.call(rbind, PopOut)
 
  # Now, we don't want to overwrite data already used
  #PNids   <- with(PopNew, paste(Year, Sex, Age))
  #PIids   <- with(Pop, paste(Year, Sex, Age))
  # discard part if necessary 
  # we won't overwrite data- if so, then it needs to be marked with LDB = 0
  #PopNew  <- PopNew[!PNids %in% PIids, ]
  # stick on new data, sort, and return without printing to console

  # !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  # ALERT, this Pop output may contain redundancies. This is intentional.
  # you can opt NOT to base further intercensal estimates on these very data
  # by picking out pop data of Type "HMDintercensal" see above for year-selection strategy
  # and further filters
  Pop     <- resortPops(rbind(Pop[, ColKeep], PopNew[, ColKeep]))
  rownames(Pop) <- NULL
 
  invisible(Pop)
})
# PopOutSx[PopOutSx$Year == 1995 & PopOutSx$Agei == 30,"Population"]

