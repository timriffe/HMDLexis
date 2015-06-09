
#' @title a function to convert Celeste's RDC death tabulations to the standard HMD R internal format
#' 
#' @description of course, if Celeste's format changes, then this will need to change. For instance, a capitalization anomaly in the column name PopName is dealt with. This function may or may not be phased out in the future, and at this time it is used specifically in the RDC environment.
#' 
#' @param DeathsTab deaths, as read into R in Celeste's format
#' 
#' @return Deaths the standard HMD R internal deaths \code{data.frame}.
#' 
#' @export
#' 
Tab2Standard <- function(DeathsTab){
  
  # only used here, internal function, defined internally just so it doesn't get lost!
  redistYrSex <- function(DSY){
    TOT         <- DSY$Deaths[DSY$Age == "ALL"]
    DSY         <- DSY[DSY$Age != "ALL", ]
    SUM         <- sum(DSY$Deaths, na.rm=TRUE)
    Redist      <- TOT - SUM
    NAi         <- is.na(DSY$Deaths)
    NNA         <- sum(NAi)
    RedistVal   <- Redist / NNA
    if (is.nan(RedistVal) | is.infinite(RedistVal)){
      RedistVal <- 0
    }
    DSY$Deaths[NAi]    <- RedistVal
    DSY$NoteCode2[NAi] <- "redistYrSex()"
    DSY
  }
  
  # correct PopName caps
  colnames(DeathsTab)[1]      <- "PopName"
  # correct Sex to lower
  DeathsTab$Sex               <- tolower(DeathsTab$Sex)
  # coerce to numeric/integer as case may be
  DeathsTab$Agei              <- suppressWarnings(as.integer(DeathsTab$Age))
  DeathsTab$Year              <- as.integer(DeathsTab$Year)
  DeathsTab$Deaths            <- as.numeric(DeathsTab$Deaths)
  # in case some deaths were marked NA by using negative (danger danger)
  # this is innocuous if tabs are done all in-RDC and there are no small nums
  DeathsTab$Deaths[DeathsTab$Deaths < 0] <- NA
  
  # redistribute as needed (innocuous if it isn't needed)
  if (any(is.na(DeathsTab$Deaths))){
    DSYL      <- split(DeathsTab, list(DeathsTab$Year, DeathsTab$Sex))
    DeathsTab <- do.call(rbind, lapply(DSYL, redistYrSex))
    rownames(DeathsTab)       <- NULL
  }
  # create standard columns which may or may not be used
  DeathsTab$YearReg           <- DeathsTab$Year
  DeathsTab$AgeIntervali      <- DeathsTab$AgeInterval      <- 1
  NAi                         <- is.na(DeathsTab$Age)
  DeathsTab$AgeIntervali[NAi] <- DeathsTab$AgeInterval[NAi] <- NA
  
  DeathsTab$NoteCode1         <- "Tab2Standard()"
  DeathsTab$NoteCode2         <- DeathsTab$NoteCode3        <- 
                                     DeathsTab$Area         <- NA
  DeathsTab$RefCode           <- "RDC Tabulation"
  DeathsTab$LDB               <- DeathsTab$YearInterval     <- 
                                   DeathsTab$Area           <- 1

  DeathsTab$Access            <- "O"
  DeathsTab$Access[DeathsTab$Deaths %<% 5 ]                 <- "C"
#  dput(colnames(Deaths))
   cnames <- c("PopName", "Area", "Year", "YearReg", "YearInterval", "Sex", 
      "Age", "AgeInterval", "Lexis", "RefCode", "Access", "Deaths", 
      "NoteCode1", "NoteCode2", "NoteCode3", "LDB", "Agei", "AgeIntervali"
   )
   colnames(DeathsTab)
  DeathsTab[, cnames]
}

#head(Tab2Standard(DeathsTab))
#
#
#










