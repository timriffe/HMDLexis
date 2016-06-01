#'
#' @title resortPops a utility function for sorting the standard Population data.frame
#' 
#' @description The sorting hierarchy is: Age within Sex within Month within Year.
#' 
#' @param Pop the standard internal Population data.frame
#' 
#' @return Pop the standard internal Population data.frame, same dimensions but sorted
#' 
#' @importFrom compiler cmpfun
#' 
#' @export 
#' 

resortPops <- cmpfun(function(Pop){
  Pop                             <- p_Date(Pop)
  Pop$Age[Pop$Age == "TOT"]       <- "999"
  Pop$Age[Pop$Age == "UNK"]       <- "888"
  Pop$Age                         <- as.integer(Pop$Age)
  Pop                             <- Pop[with(Pop, order(Date, Sex, Age)), ]
  Pop$Age                         <- as.character(Pop$Age)
  Pop$Age[Pop$Age == "888"]       <- "UNK"
  Pop$Age[Pop$Age == "999"]       <- "TOT"
  Pop$Date                        <- NULL
  Pop
})
