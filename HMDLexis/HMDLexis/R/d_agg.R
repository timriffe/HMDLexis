#' @title a function to aggregate redundant death counts
#' 
#' @description InputDB files typically contain a \code{"YearReg"} column, but it does not appear to be used either implicitly or explicitly in the LexisDB files. Further, oftentimes when a particular combination of Year-Sex-Age-Lexis is extracted from the data in the matlab code it is assumed that this will be a single number. In the case of additional late registered events, this operation will in fact return a vector, and could lead to unexpected results. This function is written to explicitly sum-in late registered events prior to further matlab operations, and can be considered a pre-processing step until such time as the matlab is phased out. In order to combine matlab and R in this way, we've decided to name the deaths file in the original raw state as \code{XYZdeath_pre.txt}, and the output of this function as the standard \code{XYZdeath.txt} so that the matlab programs will be able to find it. This is a temporary kludge for JPN and FIN, hopefully, although there are other countries where this procedure may be needed in the future, such as Finland or Australia, depending on how data come. This function may have other unexpected side-effects where redundant deaths have not been properly dealt with in the matlab code. Another side effect of this function is to add Lexis triangle death counts appropriately after splitting RR, VV or VH death counts in the case that other triangles already existed. This arises in cases where TL and TU are explicitly used to split death counts proportionally. That behavior is not document in the MP, but it is used in \code{d_rr2tltu()} and \code{d_svv()} at a minimum (possibly elsewhere too).
#' 
#' @param Deaths This should be the standard Deaths data.frame, as typically returned by \code{readInputDB()}. 
#' 
#' @importFrom data.table data.table
#' 
#' @export 
#' 


d_agg <- function(DeathsIn){

  
  # aggregate by all potentially relevant columns
  DAgg            <-  as.data.frame(
                         data.table(DeathsIn)[, Deaths := sum(Deaths), 
                          by = c("Year","YearInterval","Sex","Age","AgeInterval","Lexis")], 
                        stringsAsFactors == FALSE)
      
  throw    <- duplicated(DAgg[,c("Year","YearInterval","Sex","Age","AgeInterval","Lexis")])
  D.out    <- DAgg[!throw,]
  
 
  # resort columns
  D.Out            <- resortDeaths(D.out[,colnames(DeathsIn)])
  D.Out            <- assignNoteCode( D.Out, "d_agg()")
  # return invisibly
  invisible(D.Out)
}
