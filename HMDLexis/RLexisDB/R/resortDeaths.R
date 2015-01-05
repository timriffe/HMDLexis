#'
#' @title resortDeaths a utility function for sorting the standard Deaths data.frame
#' 
#' @description The sorting hierarchy is: Lexis within Age within Sex within Year.
#' 
#' @param Deaths the standard internal Deaths data.frame
#' 
#' @return Deaths the standard internal Deaths data.frame, same dimensions but sorted
#' 
#' @importFrom compiler cmpfun
#' 
#' @export

resortDeaths <- cmpfun(function(Deaths){
  Deaths$Age[Deaths$Age == "UNK"]   <- "888"
  Deaths$Age[Deaths$Age == "TOT"]   <- "999"
  Deaths$Age                        <- as.integer(Deaths$Age)
  Deaths$Lexis[Deaths$Age == "888"]   <- "ZZZ"
  Deaths                            <- Deaths[with(Deaths, order(Year, Sex, Age, Lexis)), ]
  Deaths$Age                        <- as.character(Deaths$Age)
  Deaths$Age[Deaths$Age == "999"]   <- "TOT"
  Deaths$Age[Deaths$Age == "888"]   <- "UNK"
  Deaths$Lexis[Deaths$Lexis == "ZZZ"] <- NA
  Deaths
})
