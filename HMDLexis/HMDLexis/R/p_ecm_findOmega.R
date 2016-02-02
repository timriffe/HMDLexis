

#'
#' @title p_ecm_findOmega finds the omega neede in order to do the extinct cohort method
#' 
#' @description This function only takes Deaths (completed) as its input. The extinct cohort method requires the determination of an age by which cohorts are considered extinct. This function finds that age and returns it so that \code{p_ecm()} can do its job cleaner. Until checked more thoroughly, this could be off by one.
#' 
#' @param Deaths The standard Deaths data.frame, in its final unrounded state after all processing.
#' @param l how many recent cohorts will be checked to compute the averages that determine omega? Default = 5.
#' @param threshold how many deaths do we consider 'close to 0'? Per MP, we use 0.5, which is reasonable, since we average over 15 values.
#' 
#' @return omega The age at which cohorts will be considered extinct
#' 
#' @importFrom reshape2 acast
#' 

p_ecm_findOmega <- function(Deaths, l = 5, threshold = 0.5){
  
  
  # add Cohort column to Deaths (possibly to be made standard later on?
  Deaths      <- d_addCohortColumn(Deaths)
  
  # get VV deaths but in AC format..
  VV          <- acast(Deaths[Deaths$Agei >= 80,], Year ~ Cohort, sum, value.var = "Deaths")
  N           <- nrow(VV)
  cohorts     <- as.integer(colnames(VV))
  this.year   <- max(as.integer(rownames(VV)))

  # this is a selection matrix
  Blocki      <- lower.tri(matrix(0, nrow = 5, ncol = 5), diag = TRUE)
  Ncoh        <- length(cohorts)
  
  # last l rows
  VV          <- VV[(N - l + 1):N, ]
  # VV contains counts in the VV shape and in PC format, hence age is in the diagonal
  
  Dstilde <- Dtilde      <- rep(0, Ncoh - l)
  for (i in 1:(Ncoh - l)){
    Left       <- matrix(FALSE, nrow = l, ncol = (i - 1))
    Right      <- matrix(FALSE, nrow = l, ncol = (Ncoh - l - i + 1))
    Dtilde[i]  <- sum(VV[cbind(Left, Blocki, Right)])
  }
  
  Ds          <- Dtilde[(Dtilde / l) > threshold][1] / l # Ds is for reproducing matlab
  Cohmax      <- max(cohorts[1:(Ncoh - l)][Dtilde / l <= threshold]) + l + 1 # test: added one
  omega       <- this.year - Cohmax
  c(omega = omega, Cohmax = Cohmax, Ds = Ds)
}


