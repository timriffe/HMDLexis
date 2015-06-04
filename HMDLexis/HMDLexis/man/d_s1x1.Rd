\name{d_s1x1}
\alias{d_s1x1}
\title{d_s1x1, a function to split RR death counts into Lexis triangles based on regression coefficients.}
\usage{
  d_s1x1(Deaths, Births)
}
\arguments{
  \item{Deaths}{LexisDB internal Deaths \code{data.frame},
  long format, all columns. Format as given by
  \code{readInputDB()}.}

  \item{Births}{LexisDB internal Births \code{data.frame},
  as given by \code{readInputDB()}.}
}
\description{
  The HMD Methods Protocol gives a formula to split Lexis
  1x1 death counts (squares) into triangles based on the
  results of a regression. Function can be run innocuously
  on any Deaths data, even if no 1x1 RR is present.
}
\details{
  This function does not yet deal with territorial
  adjustments.
}

