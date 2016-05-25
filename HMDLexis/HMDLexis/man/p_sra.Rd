\name{p_sra}
\alias{p_sra}
\title{p_sra a simple and flexible survivor ratio method function, for modularity.}
\usage{
  p_sra(Pop, Deaths, k = 5, m = 5, l = 5, a = 80, A = 90,
    maxit = 100, reproduce.matlab = FALSE,
    bordercoh = FALSE)
}
\arguments{
  \item{Pop}{The standard internal Population data.frame,
  *after* running p_ecm().}

  \item{Deaths}{after all processing is done. Completed
  triangles.}

  \item{k}{the parameter 'k' from various equations in the
  section on SRM.}

  \item{m}{the parameter 'm' from various equations in the
  section on SRM.}

  \item{l}{the parameter 'l' from various equations in the
  section on SRM.}

  \item{a}{the lower bound for population counts, almost
  always 80 (as in EC), but could be changed}

  \item{A}{controls whether we do SR 90+ or SR 85+ (or
  something else). Default 90.}

  \item{maxit}{maximum number of iterations to optimize the
  improvement coefficient, 'c'.}

  \item{reproduce.matlab}{logical. Do we include the legacy
  matlab kludge? Default \code{FALSE}. This effects results
  very little, and is only used for exact matching of
  output.}

  \item{bordercoh}{logical. Should be keep the border
  cohort between SR and EC zones?}
}
\value{
  Pop standard object, with recent old ages adjusted using
  the survivor ratio method
}
\description{
  This implementation seeks to follow the MP. LDB
  \code{Deaths} must be finalized. Population data is only
  needed to get an estimate for population in an artificial
  open age group, age 90\+ by default, though this may be
  lowered or raised. An argument \code{reproduce.matlab}
  controls whether an unnecessary matlab kludge is
  included. This function calls various other LexisDB
  functions at the moment, such as
  \code{d_addCohortColumn()}, \code{p_addCohortColumn()},
  and \code{p_ecm_findOmega()}. The parameter
  \code{bordercoh} is \code{FALSE} by default, unless
  function called by \code{p_soai()}.
}
\details{
  Years all data present in SR zone could be kept, as long
  as it is contained within the range of years present in
  \code{Pop}. If \code{Pop} is specified as
}
\note{
  This very function will be generalized to account for
  Tadj eventually. Further testing is needed to check
  whether the parameter \code{A} adjusts things properly
  for use with the SR85+ method.
}

