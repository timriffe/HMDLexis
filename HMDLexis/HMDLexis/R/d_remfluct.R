
#'
#' @title d_remfluct a function called by \code{d_soainew()} to remove outliers prior to fitting the Kannisto survival function
#' 
#' @description In order to redistribute deaths in the open age group, one first fits the Kannisto survival function to an imaginary survival function built backward from death counts, similar to the extinct cohort method. Prior to backward summing the synthetic survival curve, we look at deaths to check for a particular variety of unusual fluctuation: that caused by death dearths, as happens when small cohorts pass through old ages. If such a fluctuation is found, deaths are imputed into the affected ages using a cubic spline. 
#' 
#' @details Some comments: The matlab implementation of this whole process is not in sync with MPv5, but attempt here to reproduce it. The following quirks have been identified at this time: 1) the function only removes death troughs, but not spikes, 2) fluctuations of longer than 7 years or of only one year are not treated, 3) 25 years are taken as reference rather than the 30 stated in the MP, 4) the first year in the window cannot be an outlier, 5) only smoothed values that differ from the original values by more than 10\% are kept. These rules may have been reasonable for the three test-countries used, but there is no reason to think that they are valid for the whole HMD. Further diagnostics are not produced on a regular basis, and should be part of the country specialist checklist. The particular spline used in matlab was \code{csaps()}, and the values of the smoothing parameter given in the MP, 0.0005, applies to that particular function. One can replicate \code{csaps()} in R by using \code{pspline::smooth.Pspline()}, with smoothing parameter \code{spar} set to $(1-p)/p$.
#' 
#' @param deathsRR the universal deaths object AFTER all ages below the open age group have been split exclusively to triangles.
#' @param N the number of years to take when checking first differences for smoothness. The MP states 30, but the matlab used 25.
#' @param p1 the smoothing parameter for the first cubic spline over first differences. In the MP and matlab as 0.0005.
#' @param p2 the smoothing parameter for the second cubic spline over all ages of death counts, as seen in matlab as 0.9.
#'  
#' @importFrom pspline smooth.Pspline
#' @importFrom compiler cmpfun
#' 
#' @return deaths.out a simple data.frame with 3 columns: Age, DeathsForFitting, and Original. These are just for the purpose of fitting the Kannisto survival function, so everything is in RR format. Ages are from OA up to OA-N to OA, Deaths are RR, with detected gaps imputed with smooth counts, and Original is an unused indicator of which counts were kept vs imputed. 
#' 

d_remfluct <- cmpfun(function(deathsRR, N = 25, p1 = 0.0005, p2 = 0.9){
  # deathsRR are RR deaths including all columns, from a particular year and sex
  
  # already removed of open age, no need here
  NN                 <- nrow(deathsRR)
  DTHS               <- deathsRR$Deaths
  dths               <- DTHS[(NN - N):NN]
  Diffs              <- diff(dths)
  ages               <- deathsRR$Agei[(NN - N):NN ]
  agesd              <- ages[2:length(ages)]
  # this arg combo invokes the same spline procedure as the matlab function csaps() 
  fit.Diffs          <- c(predict(smooth.Pspline(
                          x = agesd,
                          y = Diffs, 
                          norder = 2, 
                          method = 1,
                          spar = (1 - p1) / p1     # p given in MP and matlab as 0.0005
                         ),agesd))
          
  Errors             <- fit.Diffs - Diffs
  sigma              <- sd(Errors)
  upper              <- fit.Diffs + 1.8 * sigma # 1.8 is stated in MP and used in matlab
  lower              <- fit.Diffs - 1.8 * sigma
  
  outliers <- Diffs > upper | Diffs < lower
  
  # hack to reproduce matlab: can't have the very first diff be an outlier:
  outliers[1] <- FALSE
  
  # if there's only 1 point, or if there are more than 4 points, forget about it
  if (sum(outliers) <= 1 | sum(outliers) > 4){ # in matlab, not in MP!
    return(as.data.frame(cbind(Age = ages, DeathsForFitting = dths, Original = 1)))
  }
 
  # we want 2 outliers only, one of each sign. If there are 2 or more of the same sign, need to reduce to 1
  Signs             <- sign(Errors[outliers])
  
  # not likely, but if both outliers were of the same sign, we'd reduce to one.

  # two largest outliers of opposite sign
#  if (method == 1){
#    if (length(table(Signs)) !=2){
#      return(as.data.frame(cbind(Age = ages, DeathsForFitting = dths, Original = 1)))
#    }
#    outliers          <- as.logical(rowSums(sapply(c(-1,1), function(sgn, Errors, outliers, Signs){
#                                Errors == max(abs(Errors[outliers][Signs == sgn])) * sgn
#                         }, Errors = Errors, outliers = outliers, Signs = Signs)))
#  }
#  
  # two largest outliers overall, same direction or not
#  if (method == 2){
#    outlind  <- which(outliers)
#    outliers <- abs(Errors) %in% rev(sort(abs(Errors[outlind])))[1:2]
#  }
  # two largest outliers of opposite sign, IFF it's a trough

# TODO: only use imputed model values IFF the total number of outlier deviations is less than 5

  #if (method == 3){
    if (length(table(Signs)) !=2){
      return(as.data.frame(cbind(Age = ages, DeathsForFitting = dths, Original = 1)))
    }
    outliers          <- as.logical(rowSums(sapply(c(-1,1), function(sgn, Errors, outliers, Signs){
                                Errors == max(abs(Errors[outliers][Signs == sgn])) * sgn
                         }, Errors = Errors, outliers = outliers, Signs = Signs)))
                   
    if (sign(Errors[outliers])[1] != 1 | sum(outliers) == 1){
      return(as.data.frame(cbind(Age = ages, DeathsForFitting = dths, Original = 1)))
    }
  #}
  Inds              <- which(outliers)
  
  # not in MP, but found in matlab:
  # if we're thinking of removing data that span more than 7 years, forget it:
  if (abs(diff(Inds)) > 7){
    return(as.data.frame(cbind(Age = ages, DeathsForFitting = dths, Original = 1)))
  }
  
  # this is NOT identical to the matlab routine, but may give the same results. Needs checking
  RemoveInd1        <- Inds[which.max(Errors[outliers])] # pick out max (NOT ABSOLUTE DEV!)
  RemoveInd2        <- Inds[which.min(Errors[outliers])]
  RMind             <- sort(RemoveInd1:RemoveInd2) + 1
  RA                <- ages[RMind]
 
  # now fit spline to deaths, but with outlier ages (and intermediate ages) excluded
  # the matlab code uses deaths in ALL ages:
  agesall <- as.integer(deathsRR$Age)
  RMall   <- agesall %in% RA
  
  D.fill.in         <- c(predict(smooth.Pspline(
                               x = agesall[!RMall],
                               y = DTHS[!RMall], 
                               norder = 2, 
                               method = 1,
                               spar = (1-p2) / p2    # smoother param .9 as seen in matlab, but not in MP.
                             ), RA))
  dths.interp              <- dths
  dths.interp[RMind]       <- D.fill.in
  
  # not in MP!
  # per matlab code, only take it if it makes a greater than 10% difference:
  if (max(abs(dths.interp - dths) / (dths + 0.0001)) < .1){
    return(as.data.frame(cbind(Age = ages, DeathsForFitting = dths, Original = 1)))
  }
  
  # all we need for fitting in d_soainew() are ages and imputed (or not) deaths.
  deaths.out                 <- as.data.frame(cbind(Age = ages, DeathsForFitting = dths.interp, Original = 1))
  deaths.out$Original[RMind] <- 0
  
  # finally, we're only going to do this if it makes a difference
  
  return(deaths.out)
  # now spit back deaths for Kannisto fitting. Just give back the RR deaths, all ages.
})

# for testing:
#source("/data/commons/triffe/LexisDB_R/testing/R/saveRdeaths2M.R")
#source("/data/commons/triffe/LexisDB_R/testing/R/loadMdeaths2R.R")


