\name{p_ecm_inner}
\alias{p_ecm_inner}
\title{p_ecm_inner does the work of \code{p_ecm()}}
\usage{
  p_ecm_inner(Dsex, a = 80, omega = NULL,
    reproduce.matlab = FALSE)
}
\arguments{
  \item{Dsex}{Deaths after all processing is done, subset
  of a single sex.}

  \item{a}{lowest age for which EC estimates should be
  given. Default 80.}

  \item{omega}{is the object returned by
  \code{p_ecm_findOmega()}.}

  \item{reproduce.matlab}{logical. default is \code{FALSE}.
  This affects only the border cohort between EC and SR.}
}
\value{
  Population with old ages either filled in or imputed
  using the extinct cohort method
}
\description{
  This implementation seeks to follow the MP, we have
  outsourced the work of \code{p_ecm()} to here for the
  sake of modularity. This function may also be called by
  the \code{p_ic()} function ecosystem in order to extend
  population counts in census 1 and census 2 prior to
  interpolating. Deaths must be finalized. This function
  covers area B from Figure 6 'Methods used for population
  estimates'.
}
\details{
  Conceivably, a population could have a long series with
  low open ages in the early part, but decent population
  data in the recent part. If we want a value of \code{A}
  lower than 80 for only part of the series, this can be
  achieved by subsetting \code{Dsex}, and making two calls
  to the function, then \code{rbind()}ing back together.
}

