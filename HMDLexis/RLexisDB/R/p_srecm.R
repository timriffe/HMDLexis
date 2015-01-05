#' @title p_srecm extinct cohort followed by survivor ratio methods, wrapper.
#' 
#' @description This is a wrapper, strictly. \code{p_ecm()} is first called, then \code{p_srm()}, with UNK handled first.
#' 
#' @param Pop The standard internal Population data.frame, *after* running p_ecm().
#' @param Deaths after all processing is done. Completed triangles.
#' @param k passed to \code{p_sra()}-- the parameter 'k' from various equations in the section on SRM.
#' @param m passed to \code{p_sra()}-- the parameter 'm' from various equations in the section on SRM.
#' @param a passed to \code{p_sra()}-- the parameter 'm' from various equations in the section on SRM.
#' @param A passed to \code{p_sra()}-- controls whether we do SR 90+ or SR 85+ (or something else). Default 90. (85 is untested so far)
#' @param maxit passed to \code{p_srm()}-- maximum number of iterations to optimize the improvement coefficient, 'c'.
#' @param reproduce.matlab passed to \code{p_srm()} and \code{p_ecm()}-- logical. Do we include the legacy matlab kludge? Default \code{FALSE}.
#' 
#' @export
#' 

p_srecm <- function(Pop, Deaths, k = 5, m = 5, a = 80, A = 90, maxit = 100, reproduce.matlab = FALSE){
  
  # slice off UNK, rbind back on later:
  UNKi          <- Pop$Age == "UNK"
  UNKTF         <- any(UNKi)
  if (UNKTF){
    UNK           <- Pop[UNKi, ]
    Pop           <- Pop[!UNKi, ]
  }
  
  # just put one inside the other
  Pop               <- p_sra(p_ecm(Pop = Pop, 
                                   Deaths = Deaths, 
                                   a = a,
                                   reproduce.matlab = reproduce.matlab), 
                             Deaths, 
                             k = k, 
                             m = m, 
                             a = a,
                             A = A, 
                             maxit = maxit, 
                             reproduce.matlab = reproduce.matlab)
  
  if (UNKTF){
    Pop <- resortPops(rbind(Pop, UNK))
  }
 
  invisible(Pop)
  

}
