\name{p_soai}
\alias{p_soai}
\title{redistributes population counts in the open age group}
\usage{
  p_soai(Pop, Deaths, k = 5, l = 5, m = 5, A = 90, a = 80,
    maxit = 100, reproduce.matlab = FALSE, OPENAGE = 130)
}
\arguments{
  \item{Pop}{The standard internal Population data.frame,
  *after* splitting age groups in ages lower than open and
  after redistributing unknown ages.}

  \item{Deaths}{after all processing is done. Completed
  triangles.}

  \item{k}{passed to \code{p_sra()}-- the parameter 'k'
  from various equations in the section on SRM.}

  \item{l}{passed to \code{p_sra()}-- the parameter 'l'
  from various equations in the section on SRM.}

  \item{m}{passed to \code{p_sra()}-- the parameter 'm'
  from various equations in the section on SRM.}

  \item{a}{lower age bound. Default 80, according to MP. If
  population contains some years of data with open ages
  below 80, lower the value of this parameter.}

  \item{A}{passed to \code{p_sra()}-- controls whether we
  do SR 90+ or SR 85+ (or something else). Default 90.}

  \item{maxit}{parameter passed to \code{p_sra()}, which is
  iterative. 100 is more than enough.}

  \item{reproduce.matlab}{logical. Default \code{FALSE}.
  See description.}

  \item{OPENAGE}{pad with 0s out to this age, if
  necessary.}
}
\description{
  There are two methods available in this function, toggled
  by the argument \code{reproduce.matlab}. Users should
  bear in mind a few things. In matlab, the user did not
  directly use this function. Instead it was called under
  strict assumptions within the \code{p_ic()} function
  ecosystem as needed. However, this use was not subject to
  scrutiny or review, and the function design went
  un-noticed for years, until the R translation project.
  Both the old method and a new suggestion are implemented.
  if \code{reproduce.matlab = TRUE}, then we use SR to
  extend all population counts that would be detected as
  the right side of an intercensal estimate. If
  \code{reproduce.matlab = FALSE}, we calculate the entire
  sr-ecm area, and redistribute all open-age population
  counts according to the pattern observated at the open
  age and above in each year of population data with an
  open age group. This differs in a few key ways 1) we rely
  on a mix of EC and SR estimates, depending on the year.
  2) We select all years with open age groups, rather than
  just the right-side of intercensal periods. This is
  innocuous for long stretches of annual data and adds
  little computational overhead, but is also gains us the
  ability to extend open age groups for the first years of
  data, rather than relying on \code{p_ecm_area()}, which
  will likely deprecate. 3) We deal with mid-year
  estimations in a more thorough manner.
}

