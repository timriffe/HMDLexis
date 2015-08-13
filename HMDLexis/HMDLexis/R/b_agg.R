
#' @title a function to aggregate late-registered births into the appropriate years.

#' @description InputDB files typically contain a \code{"YearReg"} column, but it does not appear to be used either implicitly or explicitly in the LexisDB files. Further, oftentimes when a particular combination of Year-Sex is extracted from the data in the matlab code it is assumed that this will be a single number. In the case of additional late registered events, this operation will in fact return a vector, and could lead to unexpected results. This function is written to explicitly sum-in late registered events prior to further matlab operations, and can be considered a pre-processing step until such time as the matlab is phased out. In order to combine matlab and R in this way, we've decided to name the births file in the original raw state as \code{XYZbirth_pre.txt}, and the output of this function as the standard \code{XYZbirth.txt} so that the matlab programs will be able to find it. This is a temporary kludge for Japan, hopefully, although there are other countries where this procedure may be needed in the future, depending on how data come. This function may have other unexpected side-effects where redundant deaths have not been properly dealt with in the matlab code. For now, this function will be called on the CS-side, prior to submitting the InputDB, as this function needs to be used prior to completion of the LexisDB R programs.
#' 
#' @param Births.In This should be the standard Births data.frame, as typically returned by \code{readInputDB()}. 
#' 
#' @importFrom data.table data.table
#' 
#' @export 
#' 

b_agg <- function(Births.In){
  # remove TOT
  Births.In <- Births.In[!Births.In$Age == "TOT", ]
  
  # aggregate by all potentially relevant columns
  B.Agg            <- as.data.frame(
                        data.table(Births.In)[, sum(Births), 
                          by = c("PopName","Area","Year","YearInterval","Sex","LDB")])
  # V1 = Deaths
  colnames(B.Agg)[colnames(B.Agg) == "V1"] <- "Births"

  # fill in unused columns with NA
  B.Agg$NoteCode3  <- B.Agg$NoteCode2 <- B.Agg$RefCode <- 
    B.Agg$Access   <- B.Agg$YearReg <- NA
  
  # make note that this is no longer an input file
  B.Agg            <- assignNoteCode(B.Agg, "b_agg()")

  # resort columns
  B.Agg            <- B.Agg[,colnames(Births.In)]
  
  # return invisibly
  invisible(B.Agg)
}
