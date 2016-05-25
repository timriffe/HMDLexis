\name{AP2AC}
\alias{AP2AC}
\title{\code{AP2AC} convert an age-cohort matrix to an age-period matrix}
\usage{
  AP2AC(APmatrix, Lexis)
}
\arguments{
  \item{APmatrix}{a matrix (age by period) of quantities
  (Deaths, Births, Exposures) in triangles, either upper or
  lower. The matrix must either represent lower triangles
  or upper triangles, not both. Whether triangles are upper
  or lower is indicated with the other argument,
  \code{Lexis}. This function can also be used to shift
  population counts from cohort to period matrices, as long
  as the population is January 1st or December 31st.
  Matrices must in any case have row names and column
  names, where rows are ages (integers, no '+'!) and
  columns are cohorts (also simple integers).}

  \item{Lexis}{either 1 or 2. Per HMD internal convention,
  1 indicates lower triangles and 2 indicates upper
  triangles. For population counts, 1 indicates January 1st
  and 2 indicates December 31st. Default is 2, so watch
  out- if you choose wrong, then column labels in the AP
  matrix output will be off by 1. If you're not sure, draw
  a diagram :-).}
}
\value{
  a \code{ACmatrix} a matrix in age-cohort format, with
  correctly labeled dimensions. It will also have an
  attribute \code{Lexis}, used for automatic Lexis argument
  detection in the case that the matrix is placed back into
  \code{AC2AP()}.
}
\description{
  \code{AP2AC} is an auxiliary function called by
  \code{Exposures_per()} (and likely by LexisDB programs,
  once these are written). It is also useful for plain
  vanilla demography in a live R session...
}
\author{
  Tim Riffe \email{triffe@demog.berkeley.edu}
}

