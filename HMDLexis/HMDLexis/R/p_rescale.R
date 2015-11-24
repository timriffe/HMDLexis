
# p_agg() is an internal for p_rescale(), but it work work as a generic top-level function too, just going against pretty much everything the LexisDB is designed to do

#' @title a function to aggregate population counts to N-year age groups
#' 
#' @description This function is the opposite of what HMD usually does! It exists for the sake of rescaling HMD intercensal estimates to known 5-year age groups for the sake of US states. This is a peculiar data situation: we have our own intercensal estimates for Jan 1 in single year age groups, but we also have official intercensals in 5-year age groups for mid-year. We'd first need to move mid year to jan 1, but in order to do a good job of this, we first split to single ages. We don't want to just take these jan1 single ages as-is because we think the HMD \code{p_ic()} single-age distribution is better, as it uses info on both left and right, whereas \code{p_split()} only uses info from right... It's really a nerdy detail that only a demographer would care about. Anyway, in order to use the 'moved' official intercensals, we need to group back to 5-year age groups, and that's what this function is for. This function will NOT be exported, and it will be treated as an internal for \code{p_rescale()}.
#' 
#' @param Pop standard HMD Population \code{data.frame}, in single ages.
#' @param N the N-year age group you want.
#' 
#' @details Not that this function does not send things to the usual 'abridged' case. age 0 is not treated differently. The only N that make sense are 5 or 10, really. Not that official intercensals were for ages 0-4, etc, for US states, so this is consistent for this limited use.
#' 
#' @importFrom data.table data.table
#' @importFrom data.table setnames
# library(data.table)
p_agg <- function(Pop, N = 5){
  # Pop <- PopN
  # determine lower bound for age intervals, assigned to each
  # single age within interval, for grouping
  Pop$ageN          <- Pop$Agei - Pop$Agei %% N
  
  # operation much easier in data.table 
  Pop               <- data.table(Pop)
  
  # aggregate populations by new age groups
  Pop[,PN := sum(Population),by=list(Year, Sex, ageN)]
  
  # remove old population and age columns
  Pop$Population    <- NULL
  Pop$Agei          <- NULL
  
  # reassign new interval
  Pop$AgeIntervali  <- N
  
  # take care of open, just in case
  Pop$AgeIntervali[Pop$AgeInterval == "+"] <- "+"

  # reassign column names
  setnames(Pop,"PN","Population")
  setnames(Pop,"ageN","Agei")
  
  # select to remove redundancy, as data.frame again
  Pop               <- as.data.frame(Pop)  
  
  Pop               <- Pop[as.integer(Pop$Age) == Pop$Agei, ]
  
  #a ssign note code to mark that it has passed through this function
  Pop               <- assignNoteCode(Pop,"p_agg()")
  
  invisible(Pop)
}
# this function is needed for an interation of the US states LDB protocol.
# For a certain period of time, we have single age censal estimates (approx 1970 and 1980),
# and also have official intercensal estimates in 5-year age groups. In this case, we 
# presume that the USCB has taken better account of migration. We could either split the 5-year
# age data using proportions from the two censuses, or else do HMD intercensal methods, and then
# rescale to match the 5-year age group totals. We opt for the second choice, although the
# differences may be minor.

#' @title rescale one population using a second population that has grouped ages, but whose totals are more believable

#' @description Note, any dates that match will get rescaled. Any dates that don't match will not get rescaled! For US states, we therefore move move midyear official estimates to Jan 1. In order to do so they were first split into single ages. Here we need them back in 5-year ages. Therefore \code{p_agg()} gets run to create \code{P5} just before doing this.

#' @param P1 a standard HMD Population \code{data.frame} in single ages. This one gets rescaled.
#' @param PN a standard HMD Population \code{data.frame}, that is in N-year age groups. This one is used as the standard.
#' 
#' @export
#' 
#' @importFrom data.table data.table
#' 
p_rescale <- function(Pop, PopN, N = 5){
  # Pop is the single age data, and Standard Pop is the grouped data. If it's not pre-grouped,
  # we group it.  Pop <- P1960_1;PopN <- P1960_5
  
  # dates must match to rescale. Take care to do this before running this function
  PopN$Date         <- with(PopN, paste(Day, Month, Year, sep = "-"))
  Pop$Date          <- with(Pop, paste(Day, Month, Year, sep = "-"))
  DatesRescale      <- intersect(PopN$Date,Pop$Date )
  
  # just a simple check
  if (length(DatesRescale) == 0){
    cat("\nThere are no comparable dates in the two populations\nNO RESCALING HAPPENED!\nReturning Pop as given, but you need to check your script because this function had no effect.\n")
    return(Pop)
  }
  
  # continue if rescaling will happen
  # these are the 2 objects we actually work with.
  PopR              <- Pop[Pop$Date %in% DatesRescale, ]    # gets rescaled
  PopN              <- PopN[PopN$Date %in% DatesRescale, ]  # standard used to scale to
  
  # we tack onto this later.
  PopKeep           <- Pop[!Pop$Date %in%  DatesRescale, ]
  
  # remove unneeded date column
  PopKeep$Date      <- NULL
  # make sure the standard pop is actually in 5-year age groups...
  if (any(PopN$AgeIntervali %==% 1)){
    PopN            <- p_agg(PopN, N = N)
  }
  # TR: this only makes sense if the standard pop is not in single ages, because otherwise we'd just
  # take the standard as true...
  
  # add N-Year age group to PopR
  PopR$ageN         <- PopR$Agei - PopR$Agei %% N
  # these ages should coincide with PopN
  
  # now get Population sums by ageN
  # operation much easier in data.table 
  # library(data.table)
  PopR              <- data.table(PopR)
  
  # aggregate populations by new age groups
  PopR[ , PN := sum(Population), by = list(Year, Sex, ageN)]
  
  # now we do homebrew merging.
  PopR              <- as.data.frame(PopR)
  
  # make a key vector, and make the standard indexable by same kind of name
  keyR              <- with(PopR, paste(Year, Sex, ageN))
  rownames(PopN)    <- with(PopN, paste(Year, Sex, Agei))
  
  # copy over standard pop this automatically takes care of ordering
  PopR$PSN          <- PopN[keyR, "Population"]
  
  # any(is.na(PopR$PSN ))
  if (any(is.na(PopR$PSN))){
    if (min(PopR$Agei[is.na(PopR$PSN)]) < 80 ){
      cat("Looks like the population you're rescaling has ages higher\nthan those available in the standard population...\nWe keep counts for those ages as-is.\n")
    }
  }
  # no scaling for these ages
  PopR$PSN[is.na(PopR$PSN)] <- PopR$PN[is.na(PopR$PSN)]
  
  # and here we rescale:
  PopR$Population     <- PopR$Population * (PopR$PSN / PopR$PN)
 
  if (any(is.infinite( PopR$Population))){
    cat("Warning: there are some Inf Population counts now,\nsince we rescale using a coeficient that could theoretically have a zero in the denominator...\nAhhhh data quality\n")
  }
  
  # remove auxilliary columns
  PopR$Date           <- NULL
  PopR$ageN 	        <- NULL
  PopR$PN             <- NULL
  PopR$PSN            <- NULL
  
  # note we used this function
  PopR                <- assignNoteCode(PopR, "p_rescale()")
  
  # stick back together, sort, return invisibly
  Pout                <- rbind(PopKeep, PopR )
  Pout                <- resortPops(Pout)
  invisible(Pout)
}







