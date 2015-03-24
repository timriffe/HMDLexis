
#'
#' @title redistributes open age group population counts using SRA method (a single year and sex of data)
#' 
#' @description This function mirrors behavior of the matlab function, but not the present live function. The live function was found to have a few bugs, which were fixed for purposes of testing. Since, this function calls \code{p_sra()}, we pass the relevant arguments down the chain. There is one shortcoming with this implementation, namely that \code{p_sra() only works with Jan 1 in the current implementation. The resultant error is likely very small, however. This function is called by \code{p_soai()}. Further critiques are given in the code comments and description of the parent function, \code{p_soai()}.
#' 
#' @param PopYrSex A single year and sex of the standard population object. Single ages only. 
#' @param DeathsSex the standard Deaths \code{data.frame}, for a single sex, all years.
#' @param k the parameter 'k' from various equations in the section on SRM.
#' @param m the parameter 'm' from various equations in the section on SRM.
#' @param l the parameter 'l' from various equations in the section on SRM.
#' @param maxit parameter passed to \code{p_sra()}, which is iterative. 100 is more than enough.
#' @param reproduce.matlab logical. Do we include legacy matlab kludges? Default \code{FALSE}. This effects results very little, and is only used for exact matching of output.
#' 
#' @return PopYrSex, with open age group redistributed over higher ages, up to and including omega, as detected by \code{p_ecm_findOmega()}.
#' 
#' @export
#' 
# takes Pop data one year at a time.
p_soai_inner_matlab <- function(PopYrSex, DeathsSex, 
  k = 5, 
  l = 5,
  m = 5, 
  maxit = 100,
  reproduce.matlab = FALSE){
  
  # 
  do.this <- any(PopYrSex$AgeInterval == "+")
  if (!do.this){
    return(PopYrSex)
  }
  
  # now get on with it.
  A     <- with(PopYrSex, Agei[AgeInterval == "+"])
  
  Year  <- unique(PopYrSex$Year)
  Month <- unique(PopYrSex$Month)
  Day   <- unique(PopYrSex$Day)
  if (A < 70){
    cat("\nSR method is not advised for open ages below 70.\n", Year, "returned unmodified\n")
    return(PopYrSex)
  }
  # from matlab. Not sure what this is about
  # if length(unique(d(:,3)))<11 | oa < 70 
  # p(ind,8) = NaN;
  # pnew     = p;
  # return
  # end
  
  # what we can abort if the first year of deaths data isn't far enough before the first year of pop data...
  # but where should that threshold be? 
  # Looking in SRA, it appears 5 years of deaths data are required for the MP defaults
  DeathsSex   <- DeathsSex[DeathsSex$Year < Year, ]

  # TODO: in an ideal world, SR method could account for partial years, by use of f1.
  # this should be implemented and enter into MPv7 with a footnote perhaps (no need for
  # lengthy documentation, although really the whole SR section of the MP required revision.
  
  # this roughly mirrors the matlab procedure
  
  if (A >= 90){
    PopSR <- p_sra(PopYrSex, 
                   DeathsSex, 
                   A = 90, a = 80, k = k, m = m, l = l, maxit = maxit,
                   reproduce.matlab = reproduce.matlab, bordercoh = TRUE)
  } else {
    PopSR <- p_sra(Pop = PopYrSex, 
                   Deaths = DeathsSex, 
                   A = A, a = min(80, A), k = k, m = m, l = l, maxit = maxit,
                   reproduce.matlab = reproduce.matlab, bordercoh = TRUE)
  }

  # do not alter ages below open age
  Pout        <- PopYrSex[PopYrSex$AgeInterval != "+", ]
  POA         <- PopYrSex$Population[PopYrSex$AgeInterval == "+"]
  # slice off needed ages
  PopSR       <- PopSR[PopSR$Year == Year & PopSR$Agei >= A, ]
  
  f1          <- ypart( Year = Year, 
                        Month = Month,
                        Day = Day,
                        reproduce.matlab = reproduce.matlab)
  # stagger as in IC 
  # TR: I didn't think hard on this, just copied the matlab
  PopSR$Population <- f1 * c(PopSR$Population[-1], 0) + 
                      (1 - f1) * PopSR$Population
  # the first notecode is necessarily due to p_sra(), but now we've again altered the data.
  PopSR$NoteCode2  <- "p_soai()"
  # rescale to match
  PopSR$Population <- (PopSR$Population / sum(PopSR$Population)) * POA
   
  # match dates.
  PopSR$Month      <- Month
  PopSR$Day        <- Day
  PopSR$Type       <- unique(PopYrSex$Type)
  PopSR$Area       <- unique(PopYrSex$Area)
  # combine data

  Pout             <- resortPops(rbind(Pout, PopSR))
   
  invisible(Pout)
}

#' @title redistributes population counts in the open age group
#' 
#' @description There are two methods available in this function, toggled by the argument \code{reproduce.matlab}. Users should bear in mind a few things. In matlab, the user did not directly use this function. Instead it was called under strict assumptions within the \code{p_ic()} function ecosystem as needed. However, this use was not subject to scrutiny or review, and the function design went un-noticed for years, until the R translation project. Both the old method and a new suggestion are implemented. if \code{reproduce.matlab = TRUE}, then we use SR to extend all population counts that would be detected as the right side of an intercensal estimate. If \code{reproduce.matlab = FALSE}, we calculate the entire sr-ecm area, and redistribute all open-age population counts according to the pattern observated at the open age and above in each year of population data with an open age group. This differs in a few key ways 1) we rely on a mix of EC and SR estimates, depending on the year. 2) We select all years with open age groups, rather than just the right-side of intercensal periods. This is innocuous for long stretches of annual data and adds little computational overhead, but is also gains us the ability to extend open age groups for the first years of data, rather than relying on \code{p_ecm_area()}, which will likely deprecate. 3) We deal with mid-year estimations in a more thorough manner.
#' 
#' @param Pop The standard internal Population data.frame, *after* splitting age groups in ages lower than open and after redistributing unknown ages.
#' @param Deaths after all processing is done. Completed triangles.
#' @param k passed to \code{p_sra()}-- the parameter 'k' from various equations in the section on SRM.
#' @param l passed to \code{p_sra()}-- the parameter 'l' from various equations in the section on SRM.
#' @param m passed to \code{p_sra()}-- the parameter 'm' from various equations in the section on SRM.
#' @param a lower age bound. Default 80, according to MP. If population contains some years of data with open ages below 80, lower the value of this parameter.
#' @param A passed to \code{p_sra()}-- controls whether we do SR 90+ or SR 85+ (or something else). Default 90.
#' @param maxit parameter passed to \code{p_sra()}, which is iterative. 100 is more than enough.
#' @param reproduce.matlab logical. Default \code{FALSE}. See description.
#' 
#' @export
#' 
#' @importFrom reshape2 acast
#' 

# k <- m <- l <- 5; reproduce.matlab = TRUE
# A <- 85; a <- 80; maxit = 100

p_soai <- function(
  Pop, 
  Deaths, 
  k = 5, 
  l = 5,
  m = 5, 
  A = 90,
  a = 80,
  maxit = 100, 
  reproduce.matlab = FALSE){
  
  # never use TOT
  Pop     <- Pop[Pop$Age != "TOT", ]
  
  Dyears  <- range(Deaths$Year)
  Pyears  <- sort(unique(Pop$Year))
  
  # explained below
  if (reproduce.matlab){
    # we can only do this for pops with at least k + m years of deaths prior to year of pop.
    # sooo,
    MinYear <- min(Dyears) + m + k + 1
    
    # first, do any of Years have open ages?
    Yrs     <- with(Pop, unique(Year[AgeInterval == "+"]))
    Yrs     <- Yrs[Yrs > MinYear]
    # move right to left
    Yrs     <- sort(Yrs, decreasing = TRUE)
    
    # these are years on the right side of gaps > 1 year
    C2years <- rev(rev(Pyears)[c(diff(rev(Pyears)) < -1,FALSE)])
    Yrs     <- Yrs[Yrs %in% C2years]
    # at a minimum, the left 11 years of data...
    NoTouch <- Pop[!Pop$Year %in% Yrs, ] 
    Pop     <- Pop[Pop$Year %in% Yrs, ]
    #
    
    # not efficient, but easy to manage pieces this way...
    Pout <- list()
    Pout[["NoTouch"]] <- NoTouch
    # sex <- "f"
    # yr <- 1992
    
    for (sex in c("m","f")){
      DeathsSex <- Deaths[Deaths$Sex == sex, ]
      Psex      <- Pop[Pop$Sex == sex, ]
      for (yr in Yrs){
        Pout[[paste0(yr, sex)]] <- 
                    p_soai_inner_matlab(PopYrSex = Psex[Psex$Year == yr, ], 
                                        DeathsSex = DeathsSex, 
                                        k = k,
                                        l = l,
                                        m = m,
                                        maxit = maxit,
                                        reproduce.matlab = reproduce.matlab)
      }
    }
  }
    # TR: SR was used in matlab to split open age groups for population counts,
    # which in HMD input data exclusively arises in an IC situation. So, for that reason
    # it was previously burried, and also not subject to group review. SR is not a good choice.
    # because
    # 1) SR is intended for almost-extinct cohorts, as a second-best approximation
    #    to results that would be obtained from EC. To use SR in EC-regions of the Lexis
    #    surface is therefore unnecessary, and not our best approximation to the truth.
    # 2.1) SR relies on death data for the 11 years *prior* to the population year in question.
    #    Even if using SR, we would prefer to use death data *centered* on the population year.
    # 2.2) we therefore cannot handle open age groups in the first 11 years of data...   
    # 3) SR is computationally intensive (slow), and also suffers its own drawbacks, such as 
    #    recursive reliance on starting values, which in the matlab implementation are ill-chosen.
  
    # I'm going to implement something similar in spirit, but more methodologically sound:
    # step 1) calculate the entire srecm region that will eventually replace pop estimates anyway.
    # step 2) appropriately use f1 to select mid-year values, as needed (along cohorts)
    # step 3) rescale to match in open age group (or 90+ if OA > 90).
    # step 4) only return needed values- discard rest of srecm region.
    # *after this step, we are ready for p_ic()
  if (!reproduce.matlab) {
    # k = 5; l = 5; m = 5; A = 5; a = 80; maxit = 100; A = 85
  # reproduce.matlab = FALSE
    PSRECM <- p_srecm(Pop,
      Deaths, 
      k = k, 
      l = l, 
      m = m, 
      a = a, 
      A = A, 
      maxit = maxit, 
      reproduce.matlab = reproduce.matlab)
    # PSRECM is a reference object, used for scaling.
    PSRECM  <- PSRECM[PSRECM$Agei >= a, ]           
    
    # now start the business. 
    Pout    <- list()
    Yrs     <- with(Pop, unique(Year[AgeInterval == "+"]))
    # maybe not do this for the extreme right-side year:
    Yrs     <- Yrs[Yrs < max(Pop$Year)]
    NoTouch <- Pop[!Pop$Year %in% Yrs, ] 
    Pout[["NoTouch"]] <- NoTouch
    Pop     <- Pop[Pop$Year %in% Yrs, ]
    
    for (sex in c("f","m")){ # sex <- "m"
      RefSex   <- acast(PSRECM[PSRECM$Sex == sex, ], Agei ~ Year, sum, value.var = "Population")
      Psex     <- Pop[Pop$Sex == sex, ]
      for (yr in Yrs){ # yr <- 1981
        PSY      <- Psex[Psex$Year == yr, ]
        OAi      <- PSY$AgeInterval == "+"
        OA       <- PSY$Agei[OAi]
        if (OA > 90){
          OA     <- 90
          OAP    <- sum(PSY$Population[ PSY$Agei >= 90])
        } else {
          OAP    <- PSY$Population[OAi]
        }
        
        
        # the ages to extract
        Ages             <- OA:130
        f1               <- ypart(Year = yr, 
                              Month = unique(Psex$Month[Psex$Year == yr]), 
                              Day = unique(Psex$Day[Psex$Year == yr]), 
                              detect.mid.year = TRUE,
                              reproduce.matlab = reproduce.matlab)
        
        Pleft            <- RefSex[as.character(Ages),as.character(yr)]
        Pright           <- RefSex[as.character(Ages),as.character(yr + 1)]
        Pref             <- Pleft * (1 - f1) + f1 * c(Pright[-1], 0)
        
        Pest             <- OAP * Pref / sum(Pref)
        
        Pouti            <- PSY[1:length(Ages), ]
        Pouti$Agei       <- Pouti$Age <- Ages
        Pouti$Population <- Pest
        Pouti$NoteCode1  <- "p_soai()"
        
        Pouti            <- rbind(PSY[PSY$Agei < OA, ], Pouti)
        Pout[[paste0(yr, sex)]] <- Pouti
      }
      
    }
  }
  # end alternative loop
  Pout <- do.call(rbind, Pout)
  rownames(Pout) <- NULL
  Pout <- resortPops(Pout)
  
  invisible(Pout)
}

# end
Pout[Pout$Year == 1981, ]

