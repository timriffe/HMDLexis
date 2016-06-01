

#' Logical utility functions
#'
#' @aliases logic HMDutils
#' 
#' @description These logical functions are like the usual ones, but \code{NA} values are treated as \code{FALSE} by default. This is not an exhaustive list, but these are the ones that speed our coding, and reduce code clutter.
#' 
#' @param x,y any two vector that can be logically compared.
#' @name HMDlogic
#' 
#' @examples
#' \dontrun{
#' c(1,2,NA,4,5) == c(1,NA,3,4,NA)
#' # compare
#' c(1,2,NA,4,5) %==% c(1,NA,3,4,NA)
#' }
NULL
#' 

#' @rdname HMDlogic
'%==%' <- function(x,y){
  x == y & !is.na(x) & !is.na(y)
}

#' @rdname HMDlogic
'%!=%' <- function(x,y){
  x != y & !is.na(x) & !is.na(y)
}

# note this is incompatible with magrittr!
#' @rdname HMDlogic
'%>%' <- function(x,y){
  x > y & !is.na(x) & !is.na(y)
}

#' @rdname HMDlogic
'%<%' <- function(x,y){
  x < y & !is.na(x) & !is.na(y)
}

#' @rdname HMDlogic
'%>=%' <- function(x,y){
  x >= y & !is.na(x) & !is.na(y)
}

#' @rdname HMDlogic
'%<=%' <- function(x,y){
  x <= y & !is.na(x) & !is.na(y)
}



#'
#' @title a function to determine whether a year is a leap year. 
#' 
#' @description In order to remove lubridate dependency, we self-detect leap years and adjust February accordingly.
#' 
#' @param Year integer of year to query
#' 
#' @return logical is the Year a leap year or not
#' 
#' @export

isLeapYear <- function (Year){      # CB: mostly good algorithm from wikipedia
  ifelse(
    ( (Year %% 4) == 0  &  (Year %% 100) != 0   ) | ( (Year %% 400) == 0 ),
    TRUE, FALSE )
}
#years <- 1700:2200
#all(lubridate::leap_year(years) == isLeapYear(years))

#'
#' @title ypart function to determine the proportion of a year passed as of a particular date
#' 
#' @description The fraction returned by this is used e.g. for intercensal estimates. Function uses 'lubridate' package to handle dates elegantly.
#' 
#' @param Year 4-digit year (string or integer)
#' @param Month month digits (string or integer, 1 or 2 characters)
#' @param Day Day of month digits (string or integer, 1 or 2 characters)
#' @param reproduce.matlab logical. Default TRUE. Assume 365 days in a year.
#' @param detect.mid.year logical. if \code{TRUE}, June 30 or July 1 will always return .5.
#' @param detect.start.end logical. default \code{TRUE}. Should Jan 1 always be 0 and Dec 31 always be 1?
#' 
#' @export
#' 

ypart <- function(Year, Month, Day, reproduce.matlab = TRUE, detect.mid.year = FALSE, detect.start.end = TRUE){
  M <- c(0, 31, 59, 90, 120, 151, 181, 212, 243, 273, 304, 334)
  if (reproduce.matlab){
    Day   <- as.integer(Day)
    Month <- as.integer(Month)
    if (is.na(Day) & is.na(Month)){
      return(.5)
    }
    if (Day == 1 & Month == 1){
      return(0)
    }
    return((M[Month] + Day) / 365)
  }
  # this chunk written just to avoiding writing p_my()
  if (detect.mid.year){
    .d <- as.integer(Day)
    .m <- as.integer(Month)
    if ((.d == 30 & .m == 6) | (.d == 1 & .m == 7)){
      return(.5)
    }
  }
  
  if (detect.start.end){
    .d <- as.integer(Day)
    .m <- as.integer(Month)
    if (.d == 1 & .m == 1){
      return(0)
    }
    if(.d == 31 & .m == 12){
      return(1)
    }
  }
  
  monthdur    <- diff(c(M,365))
  monthdur[2] <- monthdur[2] + isLeapYear(Year)
  M           <- cumsum(monthdur) - 31
  return((M[Month] + Day) / sum(monthdur))
  
  # get into date class
#  Date  <- ymd(paste(Year, Month, Day, sep = "/"))
#  # what was jan 1st?
#  Jan1  <- floor_date(Date, "year") 
#  # how many days have passed?
#  Days  <- yday(Date) - 1
#  # if we want to account for possible leap years, we get next year's Jan 1st (our Dec 31st)
#  Dec31 <- floor_date(ymd(paste(Year + 1, Month, Day, sep = "/")), "year") 
#  # and the difference is the year length
#  Denom <- as.integer(Dec31 - Jan1)
#  # only gives same as matlab on Jan 1st.
#  Days / Denom
}


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
#' @param detect.mid.year logical. default \code{FALSE}. Should June 30 and July 1 be considered .5?
#' @param detect.start.end logical. default \code{TRUE}. Should Jan 1 always be 0 and Dec 31 always be 1?
#' 
#' @return decimal value of year fraction (can be greater than 1)
#' 
#' @export
#' 
yint <- function(Day1, Month1, Year1, Day2, Month2, Year2, reproduce.matlab = TRUE, 
  detect.mid.year = FALSE, detect.start.end = TRUE){
  if (reproduce.matlab){
    return(abs(Year1 - Year2 + (Day1 - Day2) / 365 +  (Month1 - Month2) / 12))
  }
  
  # we can be more exacting, if desired:
#  abs(decimal_date(ymd(paste0(Year1, sprintf("%02d", Month1), sprintf("%02d", Day1)))) -
#      decimal_date(ymd(paste0(Year2, sprintf("%02d", Month2), sprintf("%02d", Day2)))))
#  
  Ypart1 <- ypart(Year = Year1, 
    Month = Month1, 
    Day = Day1, 
    reproduce.matlab = reproduce.matlab, 
    detect.mid.year = detect.mid.year, 
    detect.start.end = detect.start.end)
  Ypart2 <- ypart(Year = Year2, 
    Month = Month2, 
    Day = Day2, 
    reproduce.matlab = reproduce.matlab, 
    detect.mid.year = detect.mid.year, 
    detect.start.end = detect.start.end)
#  
  (1 - Ypart1) + abs(Year2 - Year1) + Ypart2
  
}


#' @title a cheap way to choose which column to assign a NoteCode to
#' 
#' @description One property of the LexisDB scripts that might be useful for downstream checks is the ability to trace which functions have modified a given data object. These can go into NoteCode slots. This function writes \code{code} to the first unoccupied \code{NoteCode} column. If all three \code{NoteCode} columns are occupied, it concatenates the end of the third column. This way we preserve a full history. Unfortunately it gets split between columns. Oh well. Good for eyeballing. This function written for the sake of modularity.
#' 
#' @param X the HMD data object that presumably has three \code{NoteCode} columns
#' @param code character string to assign to the column, typically the name of the function operating on \code{X}.
#' 
#' @export
#' 

assignNoteCode <- function(X, code){
 
  Free <- colSums(is.na(as.matrix(X[,c("NoteCode1","NoteCode2","NoteCode3")]))) > 0
  if (any(Free)){
    NoteCol <- paste0("NoteCode",min(which(Free)))
    X[[NoteCol]] <- code
  } else {
    X$NoteCode3 <- paste(X$NoteCode3, code, sep = " ")
  }
  X
}


#' @title assign the correct Area code to each year
#' 
#' @description Some newly created Population data, e.g., made by \texttt{p_ecm()} needs to have Area assigned to the new data. This function just modulates the way of accomplishing this.
#' 
#' @param Pnew an interim Population object created by some other function
#' @param PopRef the original Population object from which Area codes are to be carried over
#' 
#' @export
#' 

assignArea <- function(Pnew, PopRef){
  
  Years <- unique(PopRef$Year)
  Areas <- sapply(Years, function(yr,P){
      code <- unique(P$Area[P$Year == yr], na.rm = TRUE)
      stopifnot(length(code) == 1)
      code
    }, P = PopRef)
  names(Areas) <- Years
  Pnew$Area   <- Areas[as.character(Pnew$Year)]
  Pnew
}
# Pnew <- ECpop
# PopRef <- Psex


#'
#' @title make a Date class column for Population
#' @description The need for this is evident in cases where we have more than one date per year. This can happen, for example, if the most recent Input data year is a mid-year estimate. In that case, if \code{p_movedata()} is run before \code{p_postcensal}, there will be two dates in the most recent year, and stuff can misbehave.
#' 
#' @param Pop the standard Pop object, ideally with open age split, e.g., after running \code{p_srecm()} on the data.
#' @return Pop the same object, but with a date class column added.
#' 
#' @export
p_Date <- function(Pop){
  Pop$Date <- as.Date(with(Pop,paste(Year,Month,Day,sep="-")))
  Pop
}




