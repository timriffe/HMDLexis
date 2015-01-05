#'
#' @title yint get interval as fraction of full years
#' 
#' @description Either assume 365 days in the year, or get the precise duration.
#' 
#' @param Day1 Day of first date
#' @param Month1 Month of first date
#' @param Year1 Year of first date
#' @param Day1 Day of second date
#' @param Month1 Month of second date
#' @param Year1 Year of second date
#' @param reproduce.matlab logical. default \code{TRUE}. Assume 365 days in all years?
#' 
#' @return decimal value of year fraction (can be greater than 1)
#' 
#' @importFrom lubridate decimal_date
#' @importFrom lubridate ymd
#' 
#' @export
#' 
yint <- function(Day1, Month1, Year1, Day2, Month2, Year2, reproduce.matlab = TRUE){
  if (reproduce.matlab){
    return(abs(Year1 - Year2 + (Day1 - Day2) / 365 +  (Month1 - Month2) / 12))
  }
  
  # we can be more exacting, if desired:
  abs(decimal_date(ymd(paste0(Year1, sprintf("%02d", Month1), sprintf("%02d", Day1)))) -
      decimal_date(ymd(paste0(Year2, sprintf("%02d", Month2), sprintf("%02d", Day2)))))
  
}
