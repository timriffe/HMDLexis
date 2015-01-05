\name{d_KannistoSplitYrSex}
\alias{d_KannistoSplitYrSex}
\title{d_KannistoSplitYrSex() a function that redistributes the open age group for a particular Year-Sex.}
\usage{
  d_KannistoSplitYrSex(DeathsYearSex, N = 20, ACVVmatrix,
    reproduce.matlab)
}
\arguments{
  \item{DeathsYearSex}{A Year-Sex chunk of the standard
  Deaths object, removed of UNK and TOT, TLTU only.}

  \item{N}{the number of year data points to use when
  fitting the Kannisto survival function}

  \item{VVmaster}{a matrix as passed in from
  \code{d_soainew()}. This is a PC matrix (VV data in
  period-cohort format). Passed onward to
  \code{d_CohortAdjust()} and not used directly in this
  function.}
}
\value{
  The same Year-Sex data.frame chunk, but expanded out to
  age 130 IFF there was an open age group.
}
\description{
  This function actually does most of the work for
  \code{d_soainew()}, and it only exists as a separate
  function to avoid making a long anonymous function. The
  central construct inside \code{d_soainew()} will be the
  same, a long list of Year-Sex data.frames that later need
  to be \code{do.call("rbind",BigList)}ed together. The
  function outsources the ugly aspects of the open-age
  redistribution to just focus on the actual Kannisto
  fitting. The removal of unusual fluctuations happens in
  \code{d_remfluct()}, and the adjustment for cohort size
  happens in \code{d_CohortAdjust()}.
}
\details{
  The function calls \code{d_remfluct()} and
  \code{\d_CohortAdjust()}, which also depends on the
  \code{PC2AC()} utility function, as well as
  \code{d_KannistoSurvMin()} and \code{d_KannistoSurv()},
  which are used directly.
}

