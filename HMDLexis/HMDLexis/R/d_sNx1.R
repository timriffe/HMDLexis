#'
#' @title split RR deaths in age groups to single ages
#' 
#' @description This function iterates over each year and sex of data, isolating RR data (grouping to RR in the case of infants) and then splitting age intervals greater than one to single ages. All data in other Lexis formats are unchanged, as are pre-existing single-age data in RR format. From an error standpoint, this function may be run either before or after \code{d_unk()}. The actual work is done by the inner function, \code{d_sNx1_inner()}, which itself calls \code{d_sNx1_spline()} once data are prepared.
#' 
#' @param Deaths the standard LexisDB Deaths object.
#' 
#' @return Deaths, with all RR entries in single ages.
#'
#' @export
#' 
#' @importFrom magrittr %>%
d_sNx1 <- function(Deaths){
  
  #
  if (!any(Deaths$AgeIntervali > 1 & Deaths$Lexis == "RR")){
    cat("\nNo need to call d_sNx1()\nReturned Deaths unchanged\n")
    return(Deaths)
  }
  # many scripts include this lines, since we simply never need to keep TOT on hand...
  Deaths <- Deaths[Deaths$Age != "TOT", ]
  
  UNKi <- Deaths$Age == "UNK"
  UNKTF <- any(UNKi)
  if (UNKTF){
    UNK    <- Deaths[UNKi, ]
    Deaths <- Deaths[!UNKi, ]
  }
  
  # tried this using a data.table .SD chunk operator and had scoping issues
  # found it not worth investigating, when this works fine. If speedup is desired,
  # look to data.table or dplyr in future.
#  Deaths <- do.call(
#              rbind,
#                lapply( # this lapply iterates over year and sex
#                  split(Deaths, f = list(Deaths$Year, Deaths$Sex)), 
#                  d_sNx1_inner)
#              )
  
  # magrittr pipes are slightly more legible than the above...
  Deaths <- split(Deaths,f = list(Deaths$Year, Deaths$Sex)) %>%
    lapply(., d_sNx1_inner) %>%
    do.call(rbind, .) 
  
  if (UNKTF){
    Deaths <- resortDeaths(rbind(UNK, Deaths))
  }          
  
  invisible(Deaths)
  
}


#'
#' @title a function that splits age groups for a single sex and year
#' 
#' @description This is the inner function that gets iterated inside \code{d_sNx1()}. Some matlab oddities are preserved and others are not, pending tests. This function can be fed any single-sex-year chunk of deaths data and is innocuous if there is nothing to be done. If there is RR data to be split, any pre-existing TL, TL, single-age RR or other Lexis shapes are preserved. For any RR data with age intervals greater than 1, we return the single age RR estimates. See MP for methods description, or check out the in-house spline function \code{d_sNx1_spline()}.
#' 
#' @param DSexYr a chunk of the standard Deaths object, a single year and sex of data.
#' 
#' @return The same data, with 
#' 
#' @export
#' 
#' @importFrom pracma pchip
#' 
#' @details The \code{pchip()} function of the \code{pracma} package is only imported to follow the matlab robustness procedures. In R, we have the option of a monotonic interpolating spline too, but this has not been implemented. I'd (TR) propose using base R splines as fallbacks from our in-house method rather than importing a particular brand of cubic hermite splines.

d_sNx1_inner <- function(DSexYr){
  if (nrow(DSexYr) == 0){
    return(DSexYr)
  }
  # A convenience. Now this function can be run over any HMD deaths data
  # and it will not affect data that it shouldn't affect.
  if (!any(DSexYr$AgeIntervali > 1 & DSexYr$Lexis == "RR") | is.na(any(DSexYr$AgeIntervali > 1))){
    return(DSexYr)
  }
  
  # remove anything we don't want:
  DRR     <- DSexYr[DSexYr$Lexis == "RR", ]
  DOther  <- DSexYr[DSexYr$Lexis != "RR", ]
  
  DRR4fit <- DRR
  
  # if infant age group > 1, we need a bottom condition for the spline...
  if (DRR4fit$AgeIntervali[which.min(DRR4fit$Agei)] > 1){
    LeftCondition             <- DRR4fit[which.min(DRR4fit$Agei), ]
    LeftCondition$Age         <- LeftCondition$Agei <- LeftCondition$Agei - 1
    LeftCondition$AgeInterval <- LeftCondition$AgeIntervali <- 1
    LeftCondition$Deaths      <- 0
    DRR4fit                   <- rbind(LeftCondition, DRR4fit)
  }
  
  # 2) if infants are in triangles, group to RR
  # AgeInterval must be 1 if we're talking about a triangle,
  # so the above block is innocuous this block
  if (any(c("TL","TU") %in% DRR4fit$Lexis[DRR4fit$Agei == 0])){
    ind0          <- DRR4fit$Agei == 0
    Age0          <- DRR4fit[ind0, ][1, ]
    Age0$Lexis    <- "RR"
    Age0$Deaths   <- sum(DRR4fit$Deaths[ind0 & DRR4fit$Lexis %in% c("TU","TL","RR")])
    DRR4fit       <- rbind(Age0, DRR4fit[!ind0, ])
  }
  # an overly conservative resort, just in case
  DRR4fit                <- resortDeaths(DRR4fit)
  
  # 4) final age in vector is either the highest non-open age, plus its interval,
  # or else 106
  DRR4fit$Age4fit        <- DRR4fit$Agei + DRR4fit$AgeIntervali
  # NA will only pick up "+" age group, since we killed UNK and TOT
  if (any(is.na(DRR4fit$Age4fit))){
    oai                  <- DRR4fit$AgeInterval == "+"
    DRR4fit$Age4fit[oai] <- max(106, DRR4fit$Agei[oai] + 5)
  }
  if (any(is.na(DRR4fit$Age4fit))){
    cat(unique(DRR4fit$Year))
    stop("investigate NA")
  }
  singleAges             <- min(DRR4fit$Age4fit):max(DRR4fit$Age4fit)
  DRRCumsum              <- cumsum(DRR4fit$Deaths)
  
  # there are no single or 5-year age groups, use equivalent of matlab pchip
  if (! any(c(1,5) %in% DRR4fit$Age4fit)){
    # unknown how often this is used:
    cat("Careful: used pchip()\n")
    DRRCumsum1x1         <- pracma::pchip(DRR4fit$Age4fit, DRRCumsum, singleAges)
    names(DRRCumsum1x1)  <- singleAges
  } else {
    # this is meant to be the standard case:
    DRRCumsum1x1         <- d_sNx1_spline(DRR4fit$Age4fit, DRRCumsum, singleAges)
  }
  
  # If a discrepancy is found, check if by some odd chance this function
  # was being called in matlab. I doubt it would ever be called, though.
  # then in matlab there is another condition that isn't so clear
  # I'm not replicating this for the time being.
  #  if isnan(d1(1))
  #        disp('Standad spline');
  #        D1=spline(A, Y, a1);
  #  end
  
  # now clean output:
  # TR: I think the splines produce fairly exact output, but this takes
  # care of machine precision issues, which don't really matter since 
  # single age deaths will be decimal anyway...
  DRRCumsum1x1[as.character(DRR4fit$Age4fit)] <- DRRCumsum
  # de-cumulate deaths:
  DRR1x1              <- diff(c(0, DRRCumsum1x1))
  # can't have negatives:
  DRR1x1[DRR1x1 < 0]  <- 0
  # give age names for nice indexing
  names(DRR1x1)       <- singleAges - 1 # remove 1, back to completed age
  # now figure out what to keep and what to throw:
  rmind <- DRR$AgeIntervali > 1 & !is.na(DRR$AgeIntervali)
  DRRreplace          <- DRR[rmind, ]
  # i.e. was used for fitting, but not to be replaced
  DRRappend           <- DRR[!rmind, ] # could be age 0 and open age...
  # if by some odd change, age 0 was triangles, these are preserved.
  # matlab did not preserve them, instead grouping to RR.
  
  # and what single ages does this chunk correspond to? (we'll only grab those from DRR1x1...)
  from                <- DRRreplace$Agei
  to                  <- from + DRRreplace$AgeIntervali - 1
  singlekeep          <- sort(c(unlist(mapply(seq, from, to))))
  
  # we also need to ensure that N-year age groups with 0 deaths have
  # zero deaths in their corresponding single ages
  agesrep             <- sort(c(unlist(mapply(rep, from, DRRreplace$AgeIntervali))))
  deathsrep           <- sort(c(unlist(mapply(rep, DRRreplace$Deaths, DRRreplace$AgeIntervali))))
  names(deathsrep)    <- agesrep # use to indicate zeros
  zeros               <- deathsrep == 0
  # use cheap character indexing...
  Dsingle             <- DRR1x1[as.character(singlekeep)]
  Dsingle[zeros]      <- 0
  Dnew                <- DRRreplace[1:length(Dsingle), ]
  
  # now the verbose assigning of all the co
  Dnew$Deaths         <- Dsingle
  Dnew$Age            <- Dnew$Agei  <- singlekeep
  Dnew$RefCode        <- DSexYr$RefCode[1] #cheap ...
  Dnew$PopName        <- DSexYr$PopName[1]
  Dnew$NoteCode1      <- "d_sNx1()"
  Dnew$LDB            <- 1
  Dnew$AgeInterval    <- Dnew$AgeIntervali <- 1
  Dnew$Lexis          <- "RR"
  Dnew$YearInterval   <- 1
  Dnew$Sex            <- DSexYr$Sex[1]
  Dnew$Year           <- DSexYr$Year[1]
  Dnew$Access         <- "O"
  
  # stick stuff together:
  D_out               <- rbind(Dnew, DRRappend, DOther)
  
  D_out               <- resortDeaths(D_out)
  D_out
}


#' @title d_s5x1_spline the spline function used for splitting Nx1 age-period death counts into single ages.
#' 
#' @description This function implements the spline method described in Appendix B of the Methods Protocol. It is called by \code{d_sNx1()}, which does the outer loop, data prep, and converts the output of this function into the standard LexisDB data format. This function preserves death counts for split intervals by virtue of passing through exact points on the cumulative deaths curve.
#' 
#' @param k completed age plus age interval
#' @param Y cumulative deaths
#' @param age the single ages into which counts are to be split.
#' 
#' @return a vector of the estimated single-age cumulative death counts.
#' 
#' @export
#' 
#' @importFrom compiler cmpfun
#' 

d_sNx1_spline <- cmpfun(function(k, Y, age){
    y                 <- age * 0
    K                 <- length(k)
    omega             <- k[K]    
    
    A                 <- outer(k[-K], k[-K], "-")^3
    A[upper.tri(A)]   <- 0
    A                 <- cbind(1, rbind(0, cbind(outer(k[-K], 1:3, "^"), A)))
    # This 'A' is from eq B3
    A     <- rbind(A, 
              c(1, omega, omega^2, omega^3, (omega - k[1:(K - 1)])^3), 
              c(0:3, rep(0, K - 1)),
              c(0, 1, 2 * omega, 3 * omega^2, 3 * (omega - k[1:(K - 1)])^2)
            )
    
    # Y stays the same, loop eliminated
    d     <- c(0, Y,
      ifelse(any(k == 1),
        (Y[k == 5] - Y[k == 1]) / 2, 
        Y[k == 5] / 3),
      0)
    
    CC               <- solve(A) %*% d # eq B3
    
    x                <- min(age):max(k) 
    # several loops eliminated
    Adiff            <- t(outer(x, k[-length(k)], "-") )
    Adiff[Adiff < 0] <- 0 
    Adiff            <- Adiff ^ 3 * CC[5:length(CC)]
    
    y <- CC[1] + 
         CC[2] * x + 
         CC[3] * x^2 + 
         CC[4] * x ^ 3 +
         colSums(Adiff)
    
    names(y) <- x
    y[as.character(age)] # select specified ages
    # diff(y) gives back the deaths, not guaranteed positive...
  })

#### test compare
#k = c(1,  5, 10, 15, 20, 25)
#Y = c(14,  29,  73, 149, 267, 556)
#age = 2:20
#res <- 1e-8
#all(abs(
#d_s5x1_spline(k,Y,age) - 
# # output from matlab on same inputs:
#1e2 * c(0.198644778481017,   0.233886075949374,   0.259684335443050,
#0.290000000000018,   0.336649873417745,   0.402876202531682,
#0.489777594936758,   0.598452658227910,   0.730000000000093,
#0.883284050633025,   1.048232531645716,   1.212538987341949,
#1.363896962025540,   1.490000000000250,   1.589413924050907,
#1.704193670886361,   1.887266455696526,   2.191559493671229,
#2.670000000000420)) < res)



