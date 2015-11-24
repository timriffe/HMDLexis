#' @title p_srecm extinct cohort followed by survivor ratio methods, wrapper.
#' 
#' @description This is a wrapper, strictly. \code{p_ecm()} is first called, then \code{p_srm()}, with UNK handled first.
#' 
#' @param Pop The standard internal Population data.frame, *after* running p_ecm().
#' @param Deaths after all processing is done. Completed triangles.
#' @param k passed to \code{p_sra()}-- the parameter 'k' from various equations in the section on SRM. Default 5, per MP.
#' @param l passed to \code{p_sra()}-- the parameter 'l' from various equations in the section on SRM. Default 5, per MP.
#' @param m passed to \code{p_sra()}-- the parameter 'm' from various equations in the section on SRM. Default 5, per MP.
#' @param a lower age bound. Default 80, according to MP..
#' @param A passed to \code{p_sra()}-- controls whether we do SR 90+ or SR 85+ (or something else). Default 90. (85 is untested so far)
#' @param maxit passed to \code{p_srm()}-- maximum number of iterations to optimize the improvement coefficient, 'c'.
#' @param reproduce.matlab passed to \code{p_srm()} and \code{p_ecm()}-- logical. Do we include the legacy matlab kludge? Default \code{FALSE}.
#' 
#' @export
#' 
## CAB: particularly nasty and hard to find bug....function def below used
## default arg of 'l=l', which resulted in the variable 'l' being held for
## valuation, unevaluated because of lazy evaluation, and then passed through to
## 'p_sra'.  This results in the error message "Error: promise already under
## evaluation: recursive default argument reference or earlier problems?"
## stemming from 'l' being undefined yet still passed as an argument.

p_srecm <- function(Pop, Deaths, k = 5, l = 5, m = 5, a = 80, A = 90, maxit = 100, reproduce.matlab = FALSE){
  
 
  
  # slice off UNK, rbind back on later:
  UNKi          <- Pop$Age == "UNK"
  UNKTF         <- any(UNKi)
  if (UNKTF){
    UNK           <- Pop[UNKi, ]
    Pop           <- Pop[!UNKi, ]
  }
  
  # just put one inside the other
  # NO, enough of this recursive nested crap.  Burns hours of my time in 
  # debugging hell
  
#   Pop               <- p_sra(p_ecm(Pop = Pop, 
#                                    Deaths = Deaths, 
#                                    a = a,
#                                    reproduce.matlab = reproduce.matlab), 
#                              Deaths, 
#                              k = k, 
#                              l = l,
#                              m = m, 
#                              a = a,
#                              A = A, 
#                              maxit = maxit, 
#                              reproduce.matlab = reproduce.matlab)
#   
  Pop <- p_ecm( Pop = Pop,
                Deaths = Deaths,
                a = a,
                reproduce.matlab = reproduce.matlab
  )
  
  Pop <- p_sra(Pop, 
               Deaths,
               k = k,
               l = l,
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
