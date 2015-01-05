#'
#' @title ypart function to determine the proportion of a year passed as of a particular date
#' 
#' @description The fraction returned by this is used e.g. for intercensal estimates. Function uses 'lubridate' package to handle dates elegantly.
#' 
#' @param Year 4-digit year (string or integer)
#' @param Month month digits (string or integer, 1 or 2 characters)
#' @param Day Day of month digits (string or integer, 1 or 2 characters)
#' @param reproduce.matlab logical. Default TRUE. Assume 365 days in a year.
#'
#' @importFrom lubridate ymd
#' @importFrom lubridate floor_date
#' @importFrom lubridate yday
#' 
#' @export
#' 

ypart <- function(Year, Month, Day, reproduce.matlab = TRUE){
  
  if (reproduce.matlab){
    Day   <- as.integer(Day)
    Month <- as.integer(Month)
    if (is.na(Day) & is.na(Month)){
      return(.5)
    }
    if (Day == 1 & Month == 1){
      return(0)
    }
    
    M <- c(0,31,59,90,120,151,181,212,243,273,304,334)
    return((M[Month] + Day) / 365)
  }
  
  # get into date class
  Date  <- ymd(paste(Year, Month, Day, sep = "/"))
  # what was jan 1st?
  Jan1  <- floor_date(Date, "year") 
  # how many days have passed?
  Days  <- yday(Date) - 1
  # if we want to account for possible leap years, we get next year's Jan 1st (our Dec 31st)
  Dec31 <- floor_date(ymd(paste(Year + 1, Month, Day, sep = "/")), "year") 
  # and the difference is the year length
  Denom <- as.integer(Dec31 - Jan1)
  # only gives same as matlab on Jan 1st.
  Days / Denom
}
