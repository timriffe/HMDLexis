\name{p_srm_updatePopGivenCmid}
\alias{p_srm_updatePopGivenCmid}
\title{p_srm_updatePopGivenCmid an auxiliary function called by \code{p_srm()}.}
\usage{
  p_srm_updatePopGivenCmid(pop, VVcumsums, cmid, D_dot,
    srages, sragesc, srcohortsc, this.yearc, yr.denom,
    yr.num, k, w)
}
\arguments{
  \item{pop}{an AP matrix of running population counts in
  the SR zone, going as low as age 75 for purposes of
  estimation.}

  \item{VVcumsums}{a PC matrix of EC death counts for the
  SR cohorts, added sequentially over the running
  population estimates for the cohorts of the most recent
  year.}

  \item{cmid}{the given estimate of the improvement
  coefficient, 'c'.}

  \item{D_dot}{5x1 VV death counts for the most recent
  cohorts just left of the SR zone (explained in-line in
  the MP)}

  \item{srages}{integer vector of ages for the SR zone,
  used for shifting.}

  \item{sragesc}{character vector of ages for the SR zone,
  used for selecting.}

  \item{srcohortsc}{character vector of SR cohorts, for
  selecting.}

  \item{this.yearc}{character, the most recent year in the
  Population data.}

  \item{yr.denom}{character vector of years to be selected
  for the denominator of R*.}

  \item{yr.num}{character vector of years to be selected
  for the numerator of R*.}

  \item{k}{the parameter 'k' from various equations in the
  section on SRM.}

  \item{m}{the parameter 'm' from various equations in the
  section on SRM.}

  \item{w}{omega.}
}
\description{
  This functionalizes an inner loop-- that used for the
  optimization of 'c', an approximate rate of old-age
  improvement in mortality. All these weird arguments are
  created in the body of \code{p_srm} and should be
  examined there for more detail. This function is
  separated and compiled to increase speed. Calls aux
  function \code{PC2AP()}.
}

