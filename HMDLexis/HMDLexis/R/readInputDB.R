
# Author: triffe
###############################################################################

#' @title read in the InputDB, in its entirety: XXXpop, XXXdeath, XXXtadj, XXXbirth, XXXmonthly
#'
#' @description \code{readInputDB()} is the function that starts off the LexisDB programs by getting everything into R in memory, doing several checks, cleaning, etc. This should standardize all inputs for use downstream.
#'
#' @details this function is essentially version independent. It will read \code{XXXmonthly} births (v6) if present, and return them in the list, but does not need these. Dame with \code{XXXtadj}, as not all countries have undergone a territorial adjustment,
#'
#' @param WORKING the file path the country working folder. The last part of this path is typically the country abbreviation
#' @param XXX the standard HMD country abbreviation. If absent, \code{XXX} is inferred from the file path as its last segment. This is key for concatenating file names.
#' @param log.file log message are sent to the console by default (\code{""}). If you want these message written to a log file instead, specify a full path, including filename and extension.
#' @param InputDB.name the name of the InputDB folder, default \code{"InputDB"}. Must minimally contain files \code{XXXpop.txt}, \code{XXXdeath.txt}, \code{XXXbirth.txt}.
#' @param save.bin logical, should a copy of the function output be saved to the \code{XXX/Rbin/} directory?
#' @param verbose. logical. Default = \code{TRUE}. Should informative messages print to the console. If \code{log.file} is specified, most messages don't go to the console anyway.
#' 
#' @return a \code{list} of all InputDB items as R \code{data.frames}
#' 
#' @export

readInputDB <- function(WORKING = "/data/commons/hmd/HMDWORK/DNK", 
          XXX = NULL, 
          log.file = "", 
          InputDB.name = "InputDB",
          save.bin = FALSE,
          verbose = TRUE,
          strict = TRUE){
  # cheap verbosity toggle
  log.empty <- log.file == ""
  if (log.empty & !verbose){
    log.file <- tempfile()
  }
  # -------------------------------------------------------------------------------------------------------------
  # ensure paths valid, extract XXX if not given
  {
  if (!file.exists(WORKING)){
    stop("WORKING must be a valid file path")
  }
  InputDB.path      <- file.path(WORKING, InputDB.name)
  if (!file.exists(InputDB.path)){
    stop("InputDB.name must exist in WORKING")
  }
  if (is.null(XXX)){
    parts           <- rev(unlist(strsplit(WORKING, split = .Platform$file.sep)))
    XXX             <- parts[!parts == ""][1] 
  }
  }
  # make the LDB functions verbose! 
 cat("\nReading in InputDB items for", XXX, "...\n", "run at: ",  as.character(Sys.time()),"\n\n", file = log.file)
  # -------------------------------------------------------------------------------------------------------------
  # sort out which files to read in
  
  {
  # what files could there be? 
  potential.file.names <- paste0(XXX, c("pop.txt", "death.txt", "birth.txt", "monthly.txt", "tadj.txt"))
  names(potential.file.names) <- c("pop", "death", "birth", "monthly", "tadj")
  # what files are there?
  files.in.folder <- list.files(InputDB.path)
  
  # relevant ones:
  files.we.want   <- files.in.folder[files.in.folder %in% potential.file.names]
  
  # make sure the minimum three files are present
  if (! all(potential.file.names[1:3] %in% files.we.want)){
    stop("\nDirectory", InputDB.path, "\nmust contain each of: ", paste(potential.file.names[1:3], collapse = ", "), "\n")
  }
  }
  # -------------------------------------------------------------------------------------------------------------
  # quick TRUE / FALSE for tadj and monthly
  tadjTF      <- potential.file.names["tadj"] %in% files.we.want
  monthlyTF   <- potential.file.names["monthly"] %in% files.we.want
  # -------------------------------------------------------------------------------------------------------------
  # read in files.
  {
  Deaths     <- read.table(file.path(InputDB.path, grep(pattern = "death", x = files.we.want, value = TRUE)), 
                  sep = ",", na.strings = ".", header = TRUE, stringsAsFactors = FALSE, strip.white = TRUE)
  Pop        <- read.table(file.path(InputDB.path, grep(pattern = "pop", x = files.we.want, value = TRUE)), 
                  sep = ",", na.strings = ".", header = TRUE, stringsAsFactors = FALSE, strip.white = TRUE)
  Births     <- read.table(file.path(InputDB.path, grep(pattern = "birth", x = files.we.want, value = TRUE)), 
                  sep = ",", na.strings = ".", header = TRUE, stringsAsFactors = FALSE, strip.white = TRUE)
  if (tadjTF){
    Tadj        <- read.table(file.path(InputDB.path, grep(pattern = "tadj", x = files.we.want, value = TRUE)), 
                     sep = ",", na.strings = ".", header = TRUE, stringsAsFactors = FALSE, strip.white = TRUE)
  }
  if (monthlyTF){
    Monthly     <- read.table(file.path(InputDB.path, grep(pattern = "monthly", x = files.we.want, value = TRUE)), 
                     sep = ",", na.strings = ".", header = TRUE, stringsAsFactors = FALSE, strip.white = TRUE)
  }
  }
  # -------------------------------------------------------------------------------------------------------------
  # expected file headers:
  {
    header.expect.death     <- c("PopName", "Area", "Year", "YearReg", "YearInterval", "Sex", 
      "Age", "AgeInterval", "Lexis", "RefCode", "Access", "Deaths", 
      "NoteCode1", "NoteCode2", "NoteCode3", "LDB")
    header.expect.pop       <- c("PopName", "Area", "Sex", "Age", "AgeInterval", "Type", "Day", 
      "Month", "Year", "RefCode", "Access", "Population", "NoteCode1", 
      "NoteCode2", "NoteCode3", "LDB")
    header.expect.birth     <- c("PopName", "Area", "Sex", "Year", "YearReg", "RefCode", "Access", 
      "Births", "NoteCode1", "NoteCode2", "NoteCode3", "LDB")
    header.expect.tadj      <- c("PopName", "Year", "Age", "Area1", "Area2", "Sex", "RefCode", 
      "Access", "Type", "Value", "NoteCode1", "NoteCode2", "NoteCode3", "LDB")
    header.expect.monthly   <- c("PopName", "Area", "Year", "YearReg","Month","Vital","Births",
      "Access", "NoteCode1", "NoteCode2", "NoteCode3", "RefCode", "LDB")
  }

  # -------------------------------------------------------------------------------------------------------------
  # do headers match properly? If not, check after toupper() - if toupper() does it, then warn and relabel, otherwise stop
  {
    # Deaths
    if (!all(colnames(Deaths) == header.expect.death)){
      if (all(toupper(colnames(Deaths)) == toupper(header.expect.death))){
        colnames(Deaths) <- header.expect.death
        cat("\nThere was a case error in the column names of ", 
          potential.file.names["death"],"\nnames reassigned correctly, but you should change this!\n\n", file = log.file, append = TRUE)
      } else {
        stop("\nEither the spelling or the order of the column names is off in ", 
          potential.file.names["death"],"\nCorrect this before continuing. The correct columns should be:\n\n", 
          paste(header.expect.death, collapse = ", "),"\n")
      }
    }
    # Births
    if (!all(colnames(Births) == header.expect.birth)){
      if (all(toupper(colnames(Births)) == toupper(header.expect.birth))){
        colnames(Births) <- header.expect.birth
        cat("\nThere was a case error in the column names of ", potential.file.names["birth"], 
          "\nnames reassigned correctly, but you should change this!\n\n", file = log.file, append = TRUE)
      } else {
        stop("\nEither the spelling or the order of the column names is off in ", 
          potential.file.names["birth"], ".\nCorrect this before continuing. The correct columns should be:\n\n", 
          paste(header.expect.birth, collapse = ", "),"\n")
      }
    }
    # Population
    if (!all(colnames(Pop) == header.expect.pop)){
      if (all(toupper(colnames(Pop)) == toupper(header.expect.pop))){
        colnames(Pop) <- header.expect.pop
        cat("\nThere was a case error in the column names of ", 
          potential.file.names["pop"], "\nnames reassigned correctly, but you should change this!\n\n", file = log.file, append = TRUE)
      } else {
        stop("\nEither the spelling or the order of the column names is off in ", 
          potential.file.names["pop"], ".\nCorrect this before continuing. The correct columns should be:\n\n", 
          paste(header.expect.pop, collapse = ", "),"\n")
      }
    }
    # Territorial adjustment
    if (tadjTF){
      if (ncol(Tadj) == length(header.expect.tadj)){
        if (!all(colnames(Tadj) == header.expect.tadj)){
          if (all(toupper(colnames(Tadj)) == toupper(header.expect.tadj))){
            colnames(Tadj) <- header.expect.tadj
            cat("\nThere was a case error in the column names of ", 
              potential.file.names["tadj"], "\nnames reassigned correctly, but you should change this!\n\n", 
              file = log.file, append = TRUE)
          } else {
            stop("\nEither the spelling or the order of the column names is off in ", 
              potential.file.names["tadj"], ".\nCorrect this before continuing. The correct columns should be:\n\n", 
              paste(header.expect.tadj, collapse = ", "),"\n")
          }
        } 
      } else {
        # special case for tadj, since technically it is not documented to have an LDB column:
        if (all(toupper(colnames(Tadj)) == toupper(header.expect.tadj[-length(header.expect.tadj)]))){
          Tadj$LDB <- 1
          colnames(Tadj) <- header.expect.tadj
          cat("\nThere was no LDB column in ", potential.file.names["tadj"], 
            ".\nPreviously, for some reason, the tadj file did not have an LDB column, but it is supposed to now.\nThe column has been added with an assumed value of 1 in all rows.\nPlease add this to the original csv version of the file, as it will not be written out\n\n",
            file = log.file, append = TRUE)
        } else {
          stop("\nThe columns of ", 
            potential.file.names["tadj"], ".\nneed to be fixed\nCorrect this before continuing. The correct columns should be:\n\n", 
            paste(header.expect.tadj, collapse = ", "),"\n")
        }
      } 
    }
    # Monthly Births
    if (monthlyTF){
      if (!all(colnames(Monthly) == header.expect.monthly)){
        if (all(toupper(colnames(Monthly)) == toupper(header.expect.monthly))){
          colnames(Monthly) <- header.expect.monthly
          cat("\nThere was a case error in the column names of ", 
            potential.file.names["monthly"], "\nnames reassigned correctly, but you should change this!\n\n", file = log.file, append = TRUE)
          
        } else {
          header.expect.monthly2 <- gsub(header.expect.monthly, pattern = "NoteCode", replacement = "Note")
          if (all(toupper(colnames(Monthly)) == toupper(header.expect.monthly2))){
            colnames(Monthly) <- header.expect.monthly
            cat("\nThere was an error in the column names of", 
              potential.file.names["monthly"], "\n'Note' changed to 'NoteCode', plus there may have been other case errors\nnames reassigned correctly, but you should change this!\n\n", file = log.file, append = TRUE)
          } else {
            stop("\nEither the spelling or the order of the column names is off in ", 
              potential.file.names["monthly"], ".\nCorrect this before continuing. The correct columns should be:\n\n", 
              paste(header.expect.monthly, collapse = ", "),"\n")
          }
        }
      }
    }
  }
  # -------------------------------------------------------------------------------------------------------------
  # only work with LDB = 1 (if there are any errors in rows where LDB = 0, this isn't the place to catch them
  {
  Pop          <- Pop[Pop$LDB == 1, ]
  Deaths       <- Deaths[Deaths$LDB == 1, ]
  Births       <- Births[Births$LDB == 1, ]
  if (tadjTF){
      Tadj    <- Tadj[Tadj$LDB == 1, ]
  }
  if (monthlyTF){
    Monthly   <- Monthly[Monthly$LDB == 1, ]
  }
  }
  # -------------------------------------------------------------------------------------------------------------
  # check that PopName = XXX, try toupper(), stop if fails
  {
  # Deaths
  if (!all(Deaths$PopName == XXX)){
    if (all(toupper(Deaths$PopName) == XXX)){
      Deaths$PopName <- XXX
      cat("\nCase error in ", potential.file.names["death"], " in the PopName column.\nThis was corrected to ", 
        XXX, " but should be changed in the original file\n", file = log.file, append = TRUE)
    } else {
      stop("Problem in ", potential.file.names["death"], " in the PopName column\nWas expecting all values to be ", XXX, "\n")
    }
  }
  # Births unique(Births$PopName)
  if (!all(Births$PopName == XXX)){
    if (all(toupper(Births$PopName) == XXX)){
      Births$PopName <- XXX
      cat("\nCase error in ", potential.file.names["birth"], " in the PopName column.\nThis was corrected to ", 
        XXX, " but should be changed in the original file\n\n", file = log.file, append = TRUE)
    } else {
      stop("Problem in ", potential.file.names["birth"], " in the PopName column\nWas expecting all values to be ", XXX, "\n")
    }
  }
  # Population 
  if (!all(Pop$PopName == XXX)){
    if (all(toupper(Pop$PopName) == XXX)){
      Pop$PopName <- XXX
      cat("\nCase error in ", potential.file.names["pop"], " in the PopName column.\nThis was corrected to ", 
        XXX, " but should be changed in the original file\n\n", file = log.file, append = TRUE)
    } else {
      stop("Problem in ", potential.file.names["pop"], " in the PopName column\nWas expecting all values to be ", XXX, "\n")
    }
  }
  # Territorial adjustment
  if (tadjTF){
    if (!all(Tadj$PopName == XXX)){
      if (all(toupper(Tadj$PopName) == XXX)){
        Tadj$PopName <- XXX
        cat("\nCase error in ", potential.file.names["tadj"], " in the PopName column.\nThis was corrected to ", 
          XXX, " but should be changed in the original file\n\n", file = log.file, append = TRUE)
      } else {
        stop("Problem in ", potential.file.names["tadj"], " in the PopName column\nWas expecting all values to be ", XXX, "\n")
      }
    }
  }
  # Monthly Births
  if (monthlyTF){
    if (!all(Monthly$PopName == XXX)){
      if (all(toupper(Monthly$PopName) == XXX)){
        Monthly$PopName <- XXX
        cat("\nCase error in ", potential.file.names["monthly"], " in the PopName column.\nThis was corrected to ", 
          XXX, " but should be changed in the original file\n\n", file = log.file, append = TRUE)
      } else {
        stop("Problem in ", potential.file.names["monthly"], " in the PopName column\nWas expecting all values to be ", XXX, "\n")
      }
    }
  }
  }
  # -------------------------------------------------------------------------------------------------------------
  # check value columns (Deaths, Births, etc) - if not numeric it's because they can't be without introducing NAs
  {
  # Deaths
  # special case for AUS "*" = 1.5 (values of 1 or 2 were replaced with "*"
  if (XXX == "AUS"){
    Deaths$Deaths[Deaths$Deaths == "*"] <- 1.5
    Deaths$Deaths <- as.numeric(Deaths$Deaths)
  }
  
  if (!is.numeric(Deaths$Deaths)){
    cat("\nThe 'Deaths' column in", potential.file.names["death"], 
      "did not read in as numeric\nwhich means that coercing the column to numeric will cause NAs... bad.\n", file = log.file, append = TRUE)
    cat("\nCheck out the following rows:\n", which(is.na(as.numeric(Deaths$Deaths))), "\n", file = log.file, append = TRUE)
    cat("\nWhich had the following values:\n", Deaths$Deaths[is.na(as.numeric(Deaths$Deaths))], "\n", file = log.file, append = TRUE)
    stop("Take care of this issue first.")
  }
  # Population
  if (!is.numeric(Pop$Population)){
    cat("\nThe 'Population' column in", potential.file.names["pop"], 
      "did not read in as numeric\nwhich means that coercing the column to numeric will cause NAs... bad.\n", file = log.file, append = TRUE)
    cat("\nCheck out the following rows:\n", which(is.na(as.numeric(Pop$Population))), "\n", file = log.file, append = TRUE)
    cat("\nWhich had the following values:\n", Pop$Population[is.na(as.numeric(Pop$Population))], "\n", file = log.file, append = TRUE)
    stop("Take care of this issue first.")
  }
  # Births
  if (!is.numeric(Births$Births)){
    cat("\nThe 'Births' column in", potential.file.names["birth"], 
      "did not read in as numeric\nwhich means that coercing the column to numeric will cause NAs... bad.\n", file = log.file, append = TRUE)
    cat("\nCheck out the following rows:\n", which(is.na(as.numeric(Births$Births))), "\n", file = log.file, append = TRUE)
    cat("\nWhich had the following values:\n", Births$Births[is.na(as.numeric(Births$Births))], "\n", file = log.file, append = TRUE)
    stop("Take care of this issue first.")
  }
  # Tadj
  if (tadjTF){
    if (!is.numeric(Tadj$Value)){
      cat("\nThe 'Value' column in", potential.file.names["tadj"], 
        "did not read in as numeric\nwhich means that coercing the column to numeric will cause NAs... bad.\n", file = log.file, append = TRUE)
      cat("\nCheck out the following rows:\n", which(is.na(as.numeric(Tadj$Value))), "\n", file = log.file, append = TRUE)
      cat("\nWhich had the following values:\n", Tadj$Value[is.na(as.numeric(Tadj$Value))], "\n", file = log.file, append = TRUE)
      stop("Take care of this issue first.")
    }
  }
  # Monthly births
  if (monthlyTF){
    if (!is.numeric(Monthly$Births)){
      cat("\nThe 'Birthss' column in", potential.file.names["births"], 
        "did not read in as numeric\nwhich means that coercing the column to numeric will cause NAs... bad.\n", file = log.file, append = TRUE)
      cat("\nCheck out the following rows:\n", which(is.na(as.numeric(Monthly$Births))), "\n", file = log.file, append = TRUE)
      cat("\nWhich had the following values:\n", Monthly$Births[is.na(as.numeric(Monthly$Births))], "\n", file = log.file, append = TRUE)
      stop("Take care of this issue first.")
    }
  }
  }
  # -------------------------------------------------------------------------------------------------------------
  # check Year, Month, Day, Lexis, Type (Tadj) columns for appropriate ranges
  {
  possible.years  <- 1700:(as.integer(as.character(format(Sys.time(), "%Y"))))
  possible.months <- c(1:12, "TOT", NA) # check if can contain UNK
  possible.days   <- as.character(c(1:31), NA)
  possible.lexis  <- c("RR", "TL", "TU", "VV", "VH", "RV", NA) # NA for Age = "UNK"
  possible.tadj.types <- c("Rb", "Vx")
  # Deathss
  if (!all(Deaths$Year %in% possible.years)){
    cat("\nThe 'Year' column in", potential.file.names["death"], 
      "either has values earlier than 1700 or later than last year.\n", file = log.file, append = TRUE)
    stop("Take care of this issue first.")
  }
  Deaths[is.na(Deaths$Lexis), ]
  if (!all(Deaths$Lexis %in% possible.lexis)){
    cat("\nThe 'Lexis' column in", potential.file.names["death"], 
      "either has values earlier than 1700 or later than last year.\n", file = log.file, append = TRUE)
    stop("Take care of this issue first.")
  }
  # Population
  if (!all(Pop$Year %in% possible.years)){
    cat("\nThe 'Year' column in", potential.file.names["pop"], 
      "either has values earlier than 1700 or later than last year.\n", file = log.file, append = TRUE)
    stop("Take care of this issue first.")
  }
  if (!all(Pop$Month %in% possible.months)){
    cat("\nThe 'Month' column in", potential.file.names["pop"], 
      "has values outside the set:(1-12, 'TOT', NA)\n", file = log.file, append = TRUE)
    stop("Take care of this issue first.")
  }
  if (!all(Pop$Day %in% possible.days)){
    cat("\nThe 'Day' column in", potential.file.names["pop"], 
      "has values outside the set: (1-31, NA)\n", file = log.file, append = TRUE)
    stop("Take care of this issue first.")
  }
  # Birthss
  if (!all(Births$Year %in% possible.years)){
    cat("\nThe 'Year' column in", potential.file.names["birth"], 
      "either has values earlier than 1700 or later than last year.\n", file = log.file, append = TRUE)
    stop("Take care of this issue first.")
  }
  # Tadj
  if (tadjTF){
    if (!all(Tadj$Year %in% possible.years)){
      cat("\nThe 'Year' column in", potential.file.names["tadj"], 
        "either has values earlier than 1700 or later than last year.\n", file = log.file, append = TRUE)
      stop("Take care of this issue first.")
    }
    if (!all(Tadj$Type %in% possible.tadj.types)){
      cat("\nThe 'Type' column in", potential.file.names["tadj"], 
        "contains at least one value not in the set ('Rb','Vx').\n", file = log.file, append = TRUE)
      stop("Take care of this issue first.")
    }
  }
  # Monthly births
  if (monthlyTF){
    if (!all(Births$Year %in% possible.years)){
      cat("\nThe 'Year' column in", potential.file.names["births"], 
        "either has values earlier than 1700 or later than last year.\n", file = log.file, append = TRUE)
      stop("Take care of this issue first.")
    }
    if (!all(Births$Month %in% possible.months)){
      cat("\nThe 'Month' column in", potential.file.names["monthly"], 
        "has values outside the set:(1-12, 'TOT', NA)\n", file = log.file, append = TRUE)
      stop("Take care of this issue first.")
    }
  }
  }
  # -------------------------------------------------------------------------------------------------------------
  # check Age and AgeInterval values where appropriate - let's say no higher than 130?
  {
    if (strict){
      possible.ages       <- c(0:130, "UNK", "TOT")
    } else {
      possible.ages       <- c(0:150, "UNK", "TOT")
    }
  #possible.ages       <- c(0:130, "UNK", "TOT")
  possible.intervals  <- c(1:10, "+", NA) # NA possible if Age is "UNK"
  # Deaths
  if (!all(Deaths$Age %in% possible.ages)){
    cat("\nThe 'Age' column in", potential.file.names["death"], 
      "has values outside the set (0-130, 'UNK', 'TOT').\n", file = log.file, append = TRUE)
    stop("Take care of this issue first.")
  }
  if (!all(Deaths$AgeInterval %in% possible.intervals)){
    cat("\nThe 'AgeInterval' column in", potential.file.names["death"], 
      "has values outside the set (0-10, '+', NA).\n", file = log.file, append = TRUE)
    stop("Take care of this issue first.")
  }
  # Population
  if (!all(Pop$Age %in% possible.ages)){
    cat("\nThe 'Age' column in", potential.file.names["pop"], 
      "has values outside the set (0-130, 'UNK', 'TOT').\n", file = log.file, append = TRUE)
    stop("Take care of this issue first.")
  }
  if (!all(Pop$AgeInterval %in% possible.intervals)){
    cat("\nThe 'AgeInterval' column in", potential.file.names["pop"], 
      "has values outside the set (0-10, '+', NA).\n", file = log.file, append = TRUE)
    stop("Take care of this issue first.")
  }
  # Tadj
  if (strict){
    possible.tadj.ages <- c(NA, 0:130) 
  } else {
    possible.tadj.ages <- c(NA, 0:150) 
  }
  #possible.tadj.ages <- c(NA, 0:130) 
  if (tadjTF){
    if (!all(Tadj$Age %in% possible.tadj.ages)){
      cat("\nThe 'Age' column in", potential.file.names["tadj"], 
        "has values outside the set (0-130, NA).\n", file = log.file, append = TRUE)
      stop("Take care of this issue first.")
    }
  }
  }
  # -------------------------------------------------------------------------------------------------------------
  # check Sex- can only be "m", "f" or "b". if not, coerce to lower, warn. if still not valid throw error
  {
  possible.sexes <- c("m", "f", "b")
  # Deaths
  if (!all(Deaths$Sex %in% possible.sexes)){
    if (all(tolower(Deaths$Sex) %in% possible.sexes)){
      Deaths$Sex <- tolower(Deaths$Sex)
      cat("\nCase error in ", potential.file.names["death"], 
        " in the Sex column.\nCorrected to lower case but should be changed in the original file\n", file = log.file, append = TRUE)
    } else {
      stop("Problem in ", potential.file.names["death"], " in the Sex column\nonly values in ('f','m') accepted\n")
    }
  }
  # Births
  if (!all(Births$Sex %in% possible.sexes)){
    if (all(tolower(Births$Sex) %in% possible.sexes)){
      Births$Sex <- tolower(Births$Sex)
      cat("\nCase error in ", potential.file.names["birth"], 
        " in the Sex column.\nCorrected to lower case but should be changed in the original file\n", file = log.file, append = TRUE)
    } else {
      stop("Problem in ", potential.file.names["birth"], " in the Sex column\nonly values in ('f','m') accepted\n")
    }
  }
# Population
  if (!all(Pop$Sex %in% possible.sexes)){
    if (all(tolower(Pop$Sex) %in% possible.sexes)){
      Pop$Sex <- tolower(Pop$Sex)
      cat("\nCase error in ", potential.file.names["pop"], 
        " in the Sex column.\nCorrected to lower case but should be changed in the original file\n", file = log.file, append = TRUE)
    } else {
      stop("Problem in ", potential.file.names["pop"], " in the Sex column\nonly values in ('f','m') accepted\n")
    }
  }
# Territorial adjustment
  if (tadjTF){
    if (!all(Tadj$Sex %in% possible.sexes)){
      if (all(tolower(Tadj$Sex) %in% possible.sexes)){
        Tadj$Sex <- tolower(Tadj$Sex)
        cat("\nCase error in ", potential.file.names["tadj"], 
          " in the Sex column.\nCorrected to lower case but should be changed in the original file\n", file = log.file, append = TRUE)
      } else {
        stop("Problem in ", potential.file.names["tadj"], " in the Sex column\nonly values in ('f','m') accepted\n")
      }
    }
  }
  }
  # -------------------------------------------------------------------------------------------------------------
  # check Tadj$Area codes in Pop$Area codes - this could be more sophisticated. Monthly has Area too, but not sure what for
  if (tadjTF){
    if (!all(Tadj$Area1 %in% Pop$Area) | !all(Tadj$Area2 %in% Pop$Area) | !all(Pop$Area %in% c(Tadj$Area1, Tadj$Area2))){ # there are Pop$Area 's not in Tadj (the most recent)
      stop(potential.file.names["tadj"], " Area1 and Area2 columns must contain values present in ", potential.file.names["pop"], " Area column")
    }
  }
  # -------------------------------------------------------------------------------------------------------------
  # check Population types unique(Pop$Type) sum(is.na(Pop$Type)) Pop[is.na(Pop$Type),]
  possible.pop.types <- c("R", "O", "C", "E", "B")
  if (!all(Pop$Type %in% possible.pop.types)){
    foreign.values <- unique(Pop$Type)[!unique(Pop$Type) %in% possible.pop.types]
    stop(potential.file.names["pop"], " column 'Type' contains at least 1 value not in ('R', 'O', 'C', 'E', 'B'):\n", foreign.values)
  }
  # -------------------------------------------------------------------------------------------------------------
  # sort all files.
  {
  # Deaths:
  Deathss.in                       <- Deaths
  Deaths$Age[Deaths$Age == "UNK"]   <- "888"
  Deaths$Age[Deaths$Age == "TOT"]   <- "999"
  Deaths$Age                       <- as.integer(Deaths$Age)
  if (!all(diff(with(Deaths, order(Year, Sex, Age, Lexis))) == 1)){
    cat(potential.file.names["death"],
      " rows were not ordered correctly. Do this!\nLexis within Age (UNK and TOT penultimate and final, respectively) within Sex within Year\nFunction continued..\n\n", file = log.file, append = TRUE)
  }
  Deaths                           <- Deaths[with(Deaths, order(Year, Sex, Age, Lexis)), ]
  Deaths$Age                       <- as.character(Deaths$Age)
  Deaths$Age[Deaths$Age == "999"]   <- "TOT"
  Deaths$Age[Deaths$Age == "888"]   <- "UNK"
  #---------------------------------------------------------
  # Population:
  Pop.in                          <- Pop
  Pop$Age[Pop$Age == "TOT"]       <- "999"
  Pop$Age[Pop$Age == "UNK"]       <- "888"
  Pop$Age                         <- as.integer(Pop$Age)
  if (!all(diff(with(Pop, order(Year, Month, Sex, Age))) == 1)){
    cat(potential.file.names["pop"],
      " rows were not ordered correctly. Do this!\nAge (UNK and TOT penultimate and final, respectively) within Sex within Month within Year\nFunction continued..\n\n", file = log.file, append = TRUE)
  }
  Pop                             <- Pop[with(Pop, order(Year, Month, Sex, Age)), ]
  Pop$Age                         <- as.character(Pop$Age)
  Pop$Age[Pop$Age == "888"]       <- "UNK"
  Pop$Age[Pop$Age == "999"]       <- "TOT"
  #---------------------------------------------------------
  # Births:
  if (!all(diff(with(Births, order(Year, Sex))) == 1)){
    cat(potential.file.names["birth"],
      " rows were not ordered correctly. Do this!\nSex within Year\nFunction continued..\n\n", file = log.file, append = TRUE)
  }
  Births                           <- Births[with(Births, order(Year, Sex)), ]
  #---------------------------------------------------------
  # Tadj:
  if (tadjTF){
    Tadj$Age[Tadj$Type == "Rb"]   <- -1
    if (!all(diff(with(Tadj, order(Year, Sex, Age))) == 1)){
      cat(potential.file.names["tadj"],
        " rows were not ordered correctly. Do this!\nAge within Sex within Year (imagining that Type = 'Rb' age is -1..)\nFunction continued..\n\n", file = log.file, append = TRUE)
    }
    Tadj                          <- Tadj[with(Tadj, order(Year, Sex, Age)), ]
    Tadj$Age[Tadj$Type == "Rb"]   <- NA
  }
  # Monthly:
  if (monthlyTF){
    Monthly$Month[Monthly$Month == "UNK"] <- "888"
    Monthly$Month[Monthly$Month == "TOT"] <- "999"
    Monthly$Month                         <- as.integer(Monthly$Month)
    if (!all(diff(with(Monthly, order(Year, Month))) == 1)){
      cat(potential.file.names["monthly"],
        " rows were not ordered correctly. Do this!\nMonth (UNK and TOT penultimate and final, respectively) within Year \nFunction continued..\n\n", file = log.file, append = TRUE)
    }
    Monthly                               <- Monthly[with(Monthly, order(Year, Month)), ]
    Monthly$Month                         <- as.character(Monthly$Month)
    Monthly$Month[Monthly$Month == "888"] <- "UNK"
    Monthly$Month[Monthly$Month == "999"] <- "TOT"
  }
  }
  # -------------------------------------------------------------------------------------------------------------
  # check for duplicates
  {
  # Deaths:
  d.ind     <- duplicated(apply(Deaths, 1, function(x){
                    paste(x, collapse ="-")
                  }))
  if (any(d.ind)){
    cat(potential.file.names["death"], 
      "contains duplicate rows! Duplicates were removed,\nbut you should find and remove these in the original file!\n\n", file = log.file, append = TRUE)
    Deaths   <- Deaths[!d.ind, ]
  }
  # Pop:
  d.ind     <- duplicated(apply(Pop, 1, function(x){
                  paste(x, collapse ="-")
                }))
  if (any(d.ind)){
    cat(potential.file.names["pop"], "contains duplicate rows! Duplicates were removed,\nbut you should find and remove these in the original file!\n\n", file = log.file, append = TRUE)
    Pop   <- Pop[!d.ind, ]
  }
  # Births:
  d.ind     <- duplicated(apply(Births, 1, function(x){
                  paste(x, collapse ="-")
                }))
  if (any(d.ind)){
    cat(potential.file.names["birth"], "contains duplicate rows! Duplicates were removed,\nbut you should find and remove these in the original file!\n\n", file = log.file, append = TRUE)
    Births   <- Births[!d.ind, ]
  }
  # Tadj:
  if (tadjTF){
    d.ind     <- duplicated(apply(Tadj, 1, function(x){
                    paste(x, collapse ="-")
                  }))
    if (any(d.ind)){
      cat(potential.file.names["tadj"], "contains duplicate rows! Duplicates were removed,\nbut you should find and remove these in the original file!\n\n", file = log.file, append = TRUE)
      Tadj   <- Tadj[!d.ind, ]
    }
  }
  # Monthly:
  if (monthlyTF){
    d.ind     <- duplicated(apply(Monthly, 1, function(x){
                    paste(x, collapse ="-")
                  }))
    if (any(d.ind)){
      cat(potential.file.names["monthly"], "contains duplicate rows! Duplicates were removed,\nbut you should find and remove these in the original file!\n\n", file = log.file, append = TRUE)
      Monthly   <- Monthly[!d.ind, ]
    }
  }
  }
  # -------------------------------------------------------------------------------------------------------------
  # Check if Monthly births sum to Total Births (from Births, not TOT)
  if (monthlyTF){
    TOTind <- Monthly$Month == "TOT"
    Btot <- tapply(Births$Births, Births$Year, sum)
    Mtot <- Monthly$Births[TOTind]
    names(Mtot) <- Monthly$Year[TOTind]
    # cut down to same years
    Btot <- Btot[names(Btot) %in% names(Mtot)]
    Mtot <- Mtot[names(Mtot) %in% names(Btot)]
    # max relative difference:
    mdiff <- max(abs((Btot - Mtot)/Btot))
    if (mdiff > 1e-6){
      yearsWithDiff <- names(Btot)[(Btot - Mtot) != 0]
      cat(potential.file.names["monthly"],"totals (TOT) do not match total births in",
        potential.file.names["birth"],"\nThe maximum relative difference was",mdiff,
        "\nThe following years showed differences:\n",paste(yearsWithDiff, collapse = " ,"),
        "\nFYI. This may not be fatal, but deserves to be investiagated\n\n", file = log.file, append = TRUE)
    }
  }
  # -------------------------------------------------------------------------------------------------------------
  # Add columns to Deaths where necesary
  {
  # rather major inclusion: add Integer columsn for Age, AgeInterval
  Deaths$Agei                                     <- Deaths$Age
  Deaths$AgeIntervali                             <- Deaths$AgeInterval
  # these should be the only non-integerable cases
  Deaths$Agei[Deaths$Agei == "UNK"]               <- NA 
  Deaths$Agei[Deaths$Agei == "TOT"]               <- NA
  Deaths$AgeIntervali[Deaths$AgeIntervali == "+"] <- NA
  # now coerce to integer
  Deaths$AgeIntervali                             <- as.integer(Deaths$AgeIntervali)
  Deaths$Agei                                     <- as.integer(Deaths$Agei)
  }
  # -------------------------------------------------------------------------------------------------------------
  # Add columns to Population where necesary
  {
  Pop$Agei                                         <- Pop$Age
  Pop$AgeIntervali                                 <- Pop$AgeInterval
    # these should be the only non-integerable cases
  Pop$Agei[Pop$Agei == "UNK"]               <- NA 
  Pop$Agei[Pop$Agei == "TOT"]               <- NA
  Pop$AgeIntervali[Pop$AgeIntervali == "+"]        <- NA
    # now coerce to integer
  Pop$AgeIntervali                                 <- as.integer(Pop$AgeIntervali)
  Pop$Agei                                         <- as.integer(Pop$Agei)
  }
  # -------------------------------------------------------------------------------------------------------------
  # Create generic Tadj file, wherein non-tadj years have 1s for vx, rb
  # this is tricky. probaby a better way to program it...
  {
  Popyrs      <- range(Pop$Year)
  Birthyrs    <- range(Births$Year)
  Dyrs        <- range(Deaths$Year)
  Myrs        <- ifelse(monthlyTF,range(Monthly$Year),Popyrs)
  Tyrs        <- ifelse(tadjTF,range(Tadj$Year),Popyrs)
  minYr <- min(c(Popyrs,Birthyrs,Myrs,Tyrs,Dyrs))
  maxYr <- max(c(Popyrs,Birthyrs,Myrs,Tyrs,Dyrs))
  all.yrs <- minYr:maxYr

  NN <- max(possible.tadj.ages,na.rm=TRUE)+2
  GenericTadj <- data.frame(PopName = XXX, 
    Year = rep(all.yrs,each = NN*2), 
    Age = rep(possible.tadj.ages,length(all.yrs)*2),
    Area1 = NA,
    Area2 = NA,
    Sex = rep(c(rep("f",NN),rep("m",NN)),length(all.yrs)),
    Type = rep(c("Rb",rep("Vx",NN-1)),length(all.yrs)*2),
    Value = 1,
    LDB = 1)
  if (tadjTF){
    tadjyrs <- unique(Tadj$Year)
    fullyrs <- c(all.yrs[1],tadjyrs,rev(all.yrs)[1])
    # assign Area1, Area2
    for (i in 1:length(tadjyrs)){
      GenericTadj$Area1[GenericTadj$Year %in% (fullyrs[i]:fullyrs[i + 1])]           <- unique(Tadj$Area1[Tadj$Year == tadjyrs[i]]) 
      GenericTadj$Area2[GenericTadj$Year %in% (fullyrs[i]:(fullyrs[i + 1] -1 ))]     <- unique(Tadj$Area1[Tadj$Year == tadjyrs[i]]) 
      
      GenericTadj$Area1[GenericTadj$Year %in% ((fullyrs[i + 1] + 1):fullyrs[i + 2])]   <- unique(Tadj$Area2[Tadj$Year == tadjyrs[i]]) 
      GenericTadj$Area2[GenericTadj$Year %in% ((fullyrs[i + 1]):fullyrs[i + 2])]   <- unique(Tadj$Area2[Tadj$Year == tadjyrs[i]]) 
    }
    # impute values where needed:
    AllIDs  <- with(GenericTadj,paste(Year,Sex,Age,Type,sep="-"))
    TadjIDs <- with(Tadj,paste(Year,Sex,Age,Type,sep="-"))
    GenericTadj$Value[AllIDs %in% TadjIDs] <- Tadj$Value[TadjIDs %in% AllIDs]
  } else {
    Area <- unique(Pop$Area)
    GenericTadj$Area1 <- Area
    GenericTadj$Area2 <- Area
  }
  }
  # TODO: discuss other checks that might be necessary or reasonable. Can ignore note codes
  # -------------------------------------------------------------------------------------------------------------
  # prepare output
  {
  output <- list(Deaths = Deaths, Population = Pop, Births = Births, Tadj = GenericTadj)
  
  if (monthlyTF){
    output[["Monthly"]] <- Monthly
  } else {
    output[["Monthly"]] <- NULL
  }
  }
  # save as
  if (save.bin){
    dir.path   <- file.path(WORKING, "Rbin")
    if(!file.exists(dir.path)) {
      dir.create(dir.path)
      Sys.chmod(dir.path, mode = "2775", use_umask = FALSE)
      chgrp <- paste0("chgrp hmdcalc ",dir.path)
      system(chgrp)
    }
    # prepare file  
    out.path   <- file.path(dir.path, "IDB_data.frames.Rdata")
    # saving happens here
    save(output, file = out.path)
    Sys.chmod(out.path, mode = "2775", use_umask = FALSE)
    chgrp <- paste0("chgrp hmdcalc ",out.path)
    system(chgrp)
    cat("\n\nR binary copy of function output saved to:\n", out.path, "\n\n", file = log.file, append = TRUE)
  }
  
  # unlink junk log if necessary  
  if (log.empty & !verbose){
    unlink(log.file)
  }
  
  # return invisibly
  cat("\nNo fatal errors, but do mind warnings returned, if any. InputDB read in successfully.\n\n", file = log.file, append = TRUE)
  if (!log.empty){
    browser.path <- file.path("http://www.mortality.org/Working", XXX, "LDBlog")
    cat("log file written to:", log.file, "\n\nCan also be viewed in browser here:\n", browser.path, "\n")
  }
  invisible(output)
} # end of function
