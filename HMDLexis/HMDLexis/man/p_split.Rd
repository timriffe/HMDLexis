\name{p_split}
\alias{p_split}
\title{split age intervals based on single ages from the right}
\usage{
  p_split(Pop, Deaths, Births, a = 80,
    reproduce.matlab = FALSE)
}
\arguments{
  \item{Pop}{Standard \code{Population} IDB object. At
  minimum, the most recent year must be in single ages.}

  \item{Deaths}{the standard \code{Deaths} IDB object,
  after computations finalized.}

  \item{a}{the age where extinct cohort populations will
  take over. No need to split intervals above this age.}

  \item{reproduce.matlab}{logical. Default \code{FALSE}.
  Should we reproduce all aspects of the matlab code?}
}
\value{
  C2 in single ages
}
\description{
  This function implements the MP verbal description of
  splitting population counts in age groups. Basically, we
  add deaths going back in time from census 2 (right-side),
  which will give an estimate for the single-age
  populations in the interval \code{[a,a+n]}. We then
  rescale these single age estimates to sum to the C1
  population count in the interval.
}
\details{
  The only matlab oddity that is potentially affected by
  \code{reproduce.matlab} is the handling of dates. Matlab
  does not handle leap years. This function makes calls to
  two other LexisDB functions, such as \code{ypart()} and
  \code{p_addCohortColumn()}.
}

