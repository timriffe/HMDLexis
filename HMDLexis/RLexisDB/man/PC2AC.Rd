\name{PC2AC}
\alias{PC2AC}
\title{\code{PC2AC} convert a period-cohort (VV) matrix to an age-cohort matrix, wherein entries still represent the VV shape.}
\usage{
  PC2AC(PCmatrix, Lexis)
}
\arguments{
  \item{PCmatrix}{a matrix (period by cohort) of quantities
  in VV shape. May work for other shapes too, but you'll
  need to experiment with Lexis = 1 or 2.}

  \item{Lexis}{either 1 or 2. Per HMD internal convention,
  1 indicates lower triangles and 2 indicates upper
  triangles. For population counts, 1 indicates January 1st
  and 2 indicates December 31st. Default is 2, so watch
  out- if you choose wrong, then column labels in the AP
  matrix output will be off by 1. If you're not sure, draw
  a diagram. VV gets a 2 by default!}
}
\value{
  a \code{ACmatrix} a matrix in age-cohort format, with
  correctly labeled dimensions. It will also have an
  attribute \code{Lexis}, used for automatic Lexis argument
  detection in the case that the matrix is placed back into
  \code{AC2AP()}. Note that if the original contents were
  the VV shape that age is the lower left corner, and that
  the shape spans two ages. Draw a picture if necessary.
}
\description{
  \code{PC2AC} is an auxiliary function available to call
  by any of the HMD function family.
}
\author{
  Tim Riffe \email{triffe@demog.berkeley.edu}
}

