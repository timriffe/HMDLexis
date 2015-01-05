\name{p_srm}
\alias{p_srm}
\title{p_srm a simple survivor ratio method function, for modularity.}
\usage{
  p_srm(Pop, Deaths, k = 5, m = 5, l = 5, MaxAge = 90,
    maxit = 100, reproduce.matlab = FALSE)
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

  \item{MaxAge}{controls whether we do SR 90+ or SR 85+ (or
  something else). Default 90. (85 is untested so far)}

  \item{maxit}{maximum number of iterations to optimize the
  improvement coefficient, 'c'.}

  \item{reproduce.matlab}{logical. Do we include the legacy
  matlab kludge? Default \code{FALSE}. If \code{TRUE} be
  sure to have run \code{p_ecm()} with this argument set to
  \code{TRUE} first.}
}
\value{
  Pop with old ages adjusted using the survivor ratio
  method
}
\description{
  This implementation seeks to follow the MP. LDB
  \code{Deaths} must be finalized. Populations must be in
  single ages and all years, but not necessarily up to age
  130 in all years. An argument \code{reproduce.matlab}
  controls whether an unnecessary matlab kludge is
  included. This function calls various other LexisDB
  functions at the moment, such as
  \code{d_addCohortColumn()}, \code{p_addCohortColumn()},
  \code{p_ecm_findOmega()}, as well as its own auxilliary
  function, \code{p_srm_updatePopGivenCmid()}, which serves
  to speed up a loop for serial dependency.
}
\note{
  This very function will be generalized to account for
  Tadj eventually. Further testing is needed to check
  whether the parameter \code{MaxAge} adjusts things
  properly for use with the SR85+ method.
}

