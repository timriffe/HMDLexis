# TODO: TR : 1) for future MP revision: there is no reason this function shouldn't follow the exact
# same methodology as IC()
# 2) also no reason to ignore monthly birth variation, which we now have. S0, two
# improvements are possible.

#' @title Internal work function of \code{p_movedata()}
#'
#' @description See \code{?p_movedata} for some details on the quirks of this function. 
#' 
#' @param PopL The left-side population object. A single year of data, including both males and females. Single ages only. referred to as C1 in MP.
#' @param PopR The right-side population object. A single year of data, including both males and females. Single ages only. referred to as C2 in MP. Must be adjacent years.
#' @param detect.mid.year logical. if \code{TRUE}, June 30 or July 1 will always return .5.
#' @param detect.start.end logical. default \code{TRUE}. Should Jan 1 always be 0 and Dec 31 always be 1?
#' @param reproduce.matlab logical. Default TRUE. Assume 365 days in a year.
#' @param OPENAGE pad with 0s out to this age, if necessary.
#' 
#' @export


p_movedata_inner <- function(
  PopL, 
  PopR, 
  detect.mid.year = TRUE, 
  detect.start.end = TRUE, 
  reproduce.matlab = FALSE){
  
  # 1) left and right need to be from consecutive years...
  # otherwise return PopR.
  if (unique(PopR$Year) - unique(PopL$Year) > 1){
    return(PopR)
  }
  # can just toss open age groups
  # if we're here, then we might be able to do interpolation.
  # must have all single years.
  stopifnot(all(PopL$AgeIntervali==1))
  stopifnot(all(PopR$AgeIntervali==1))
  
  
  PopL        <- PopL[order(PopL$Sex, PopL$Agei), ]
  PopR        <- PopR[order(PopR$Sex, PopR$Agei), ]
  # add date prop column to PopL, PopR
  if (!reproduce.matlab){
    f1 <- ypart(Year = unique(PopL$Year), 
                Month = unique(PopL$Month), 
                Day = unique(PopL$Day), 
                detect.mid.year = detect.mid.year,
                detect.start.end = detect.start.end,
                reproduce.matlab = reproduce.matlab)
    f2 <- ypart(Year = unique(PopR$Year), 
                Month = unique(PopR$Month), 
                Day = unique(PopR$Day), 
                detect.mid.year = detect.mid.year,
                detect.start.end = detect.start.end,
                reproduce.matlab = reproduce.matlab)
    
    if (length(f1) > 1 | length(f2) > 1){
      stop("multiple dates in use in the same year, makes intercensals tricky.\\Time to do some digging.")
    }
  } else {
    f1 <- ypart(Year = PopL$Year[1], 
                Month = PopL$Month[1], 
                Day = PopL$Day[1], 
                reproduce.matlab = reproduce.matlab)
    f2 <- ypart(Year = PopR$Year[1], 
                Month = PopR$Month[1], 
                Day = PopR$Day[1], 
                reproduce.matlab = reproduce.matlab)
  }
  
  totspan     	  <- 1 - f1 + f2
  propl           <- (1 - f1) / totspan
  propr           <- f2 / totspan
  
  # split sexes. no need for loop, pretty simple here
  ML              <- PopL[PopL$Sex == "m", ]
  FL              <- PopL[PopL$Sex == "f", ]
  MR    	        <- PopR[PopR$Sex == "m", ]
  FR              <- PopR[PopR$Sex == "f", ]
  
  
  # can only interpolate shared ages
  agesm           <- union(ML$Agei, MR$Agei)
  agesf           <- union(FL$Agei, FR$Agei)
  
  # another hopefully-redundant check:
  stopifnot(all(0:max(agesm)) %in% ML$Agei)
  stopifnot(all(0:max(agesm)) %in% MR$Agei)
  stopifnot(all(0:max(agesf)) %in% FL$Agei)
  stopifnot(all(0:max(agesf)) %in% FR$Agei)
  
  Mout            <- MR
  Fout            <- FR
  
  Mout$Month      <- 1
  Mout$Day        <- 1
  Fout$Month      <- 1
  Fout$Day        <- 1
  
  # TR the actual calc is super simple here. simple is good, but we can
  # do better than this 
  Mout$Population <- MR$Population * propr + ML$Population * propl 
  Fout$Population <- FR$Population * propr + FL$Population * propl 
  
  Pout <- rbind(Fout, Mout)
  
  Pout <- assignNoteCode(Pout, "p_movedata()")
   
  Pout
}

#' @title Move a series of midyear single-age populations estimates to January 1.
#'
#' @description This function uses linear interpolation, per the MP, to estimate January 1 population. Interpolation takes place between like-ages in adjacent years. Time is calculated either assuming 365 days per year (\code{reproduce.matlab = TRUE}) or with leap-year detection (\code{reproduce.matlab = FALSE}). This function covers the behavior of the matlab \code{p_movedata()} and \code{p_my()}, depending on the settings of \code{detect.mid.year} and \code{detect.start.end}. For \code{p_my()}, set these two arguments to \code{TRUE}, otherwise \code{FALSE}. It seems, however, that when a statistical office declare June 30 or July 1 that these are both intended to be midyear estimates, and so TR thinks it wise to just set the defaults to \code{TRUE} and not change them. \code{reproduce.matlab} only affects the leapyear aspect of estimating year proportions. See \code{?ypart} for more info on that. In general, TR thinks it would be best to replace this method with a) \code{p_ic()} for single-year spans or b) a version of \code{p_ic()} for single-year spans that accounts for within-cohort variation in the population distribution, e.g., from the monthly birth distribution of a given cohort. One quirk of this function is that Jan 1 of year t is considered adjacent to Dec 31 of year t+1. This means that it makes a difference if you run \code{p_ey2ny()} before this or not!
#' 
#' @param Pop The standard HMD population \code{data.frame}
#' @param detect.mid.year logical. if \code{TRUE}, June 30 or July 1 will always return .5.
#' @param detect.start.end logical. default \code{TRUE}. Should Jan 1 always be 0 and Dec 31 always be 1?
#' @param reproduce.matlab logical. Default TRUE. Assume 365 days in a year.
#' 
#' @export

p_movedata <- function(Pop, detect.mid.year = TRUE, detect.start.end = TRUE, reproduce.matlab = FALSE, OPENAGE = 130){
  
  # TR: added this line. This pads any years with no open age groups out to 130.
  Pop <- suppressWarnings(p_long(Pop, OPENAGE = OPENAGE))
  # assumes we have single-age population counts already.
  # should do 5-yr splitting beforehand...
  # move from right to left. We'll do this using inefficient programming,
  # for the sake of consistency...
  years      <- sort(unique(Pop$Year))
  
  # if you were going to split open ages, you should have done it first.
  # stopifnot(any(Pop$AgeIntervali %!=% 1))
   
  Pop <- Pop[Pop$AgeIntervali %==% 1, ]
  # calculate jumps in completed years.
  # this allows for use of this method to fill in a Jan 1, yr t+1 estimate based on dates as far apart
  # as Jan 1, year t, and Dec 31, year t+1, i.e. an actual 2-year gap.
  # we won't fill in for Jan 1, t+1 if our two dates are Jan 1, t and Jan1, t+2. 
  # That's the cutoff, that's our rule. So it makes a difference if we run
  # ey2ny() before doing this function...
  Pout <- list()
  
  # TR: unsure what would happen if years not continuous, but this works for now.
  for (yr in rev(years)){ #
    PopR <- Pop[Pop$Year == yr, ]
    if ((yr - 1) %in% years){
      PopL <- Pop[Pop$Year == (yr - 1), ]
      
      PopR <- p_movedata_inner( PopL, 
                                PopR, 
                                detect.mid.year = detect.mid.year, 
                                detect.start.end = detect.start.end, 
                                reproduce.matlab = reproduce.matlab)
    }
    Pout[[as.character(yr)]] <- PopR
  }
  
  Pout <- do.call(rbind, Pout)
  Pout <- resortPops(Pout)
  
  Pout
  
}


















