
##' Hermitean Interpolation Polynomials
##' 
##' Piecewise Cubic Hermitean Interpolation Polynomials.
##' 
##' \code{pchip} is a `shape-preserving' piecewise cubic Hermite polynomial
##' approach that apptempts to determine slopes such that function values do
##' not overshoot data values.  \code{pchipfun} is a wrapper around
##' \code{pchip} and returns a function.  Both \code{pchip} and the function
##' returned by \code{pchipfun} are vectorized.
##' 
##' \code{xi} and \code{yi} must be vectors of the same length greater or equal
##' 3 (for cubic interpolation to be possible), and \code{xi} must be sorted.
##' \code{pchip} can be applied to points outside \code{[min(xi), max(xi)]},
##' but the result does not make much sense outside this interval.
##' 

##' @param xi,yi x- and y-coordinates of supporting nodes.
##' @param x x-coordinates of interpolation points.
##' @return Values of interpolated data at points \code{x}.
##' 
##' All functions in this file have been copied out of the pracma
##' package by Hans Werner Borchers. CRAN version 1.8.3
##' HMD needed to limit package dependencies and so we've forked this functionality.
##' This function was needed in order to replicate the pchip() function of matlab,
##' which unfortunately is how deaths in age groups are split into sigle ages in 
##' certain circumstances. This method is in serious need of revision, as it has 
##' flaws. CB and TR have been thinking and working on this recently. The pchip()
##' solution was a patch-over.
##' CB went and found these functions in the pracma package and now they live here
##' as a snapshot

pchip <- function(xi, yi, x) {
  stopifnot(is.numeric(xi), is.numeric(yi), is.numeric(x))
  # xi <- c(xi); yi <- c(yi); x <- c(x)
  if (!is.sorted(xi))
    stop("Argument 'xi' must be a sorted vector of real numbers.")
  n <- length(xi);
  if (length(yi) != n)
    stop("Arguments 'xi', 'yi' must be vectors of equal length.")
  if (n <= 2)
    stop("At least three points needed for cubic interpolation.")
  
  # First derivatives
  h <- diff(xi)
  delta <- diff(yi) / h
  d <- .pchipslopes(h, delta)
  
  # Piecewise polynomial coefficients
  a <- (3*delta - 2*d[1:(n-1)] - d[2:n]) / h
  b <- (d[1:(n-1)] - 2*delta + d[2:n]) / h^2;
  
  # Find subinterval indices k so that xi[k] <= x < xi[k+1]
  k <- rep(1, length(x))
  for (j in 2:(n-1)) {
    k[xi[j] <= x] <- j
  }
  
  # Evaluate interpolant
  s <- x - xi[k]
  v <- yi[k] + s*(d[k] + s*(a[k] + s*b[k]))
  
  return(v)
}


.pchipslopes <- function(h, delta) {
  
  # Slopes at interior points
  n <- length(h) + 1
  d <- numeric(length(h))
  k <- which(sign(delta[1:(n-2)]) * sign(delta[2:(n-1)]) > 0) + 1
  w1 <- 2*h[k] + h[k-1]
  w2 <- h[k]+2*h[k-1]
  d[k] <- (w1+w2) / (w1/delta[k-1] + w2/delta[k])
  
  # Slopes at endpoints
  d[1] <- .pchipend(h[1], h[2], delta[1], delta[2])
  d[n] <- .pchipend(h[n-1], h[n-2], delta[n-1], delta[n-2])
  
  return(d)
}


.pchipend <- function(h1, h2, del1, del2) {
  # Noncentered, shape-preserving, three-point formula.
  d <- ((2*h1 + h2)*del1 - h1*del2) / (h1 + h2)
  if (sign(d) != sign(del1)) {
    d <- 0
  } else if ((sign(del1) != sign(del2)) && (abs(d) > abs(3*del1))) {
    d <- 3*del1
  }
  return(d)
}


pchipfun <- function(xi, yi) {
  stopifnot(is.numeric(xi), is.numeric(yi))
  # xi <- c(xi); yi <- c(yi)
  if (!is.sorted(xi))
    stop("Argument 'xi' must be a sorted vector of real numbers.")
  n <- length(xi);
  if (length(yi) != n)
    stop("Arguments 'xi', 'yi' must be vectors of equal length.")
  if (n <= 2)
    stop("At least three points needed for cubic interpolation.")
  
  function(x) pchip(xi, yi, x)
}

is.sorted <- function(a) !is.unsorted(a)