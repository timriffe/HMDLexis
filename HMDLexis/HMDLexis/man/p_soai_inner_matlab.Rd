\name{p_soai_inner_matlab}
\alias{p_soai_inner_matlab}
\title{redistributes open age group population counts using SRA method (a single year and sex of data)}
\usage{
  p_soai_inner_matlab(PopYrSex, DeathsSex, k = 5, l = 5,
    m = 5, maxit = 100, reproduce.matlab = FALSE)
}
\arguments{
  \item{PopYrSex}{A single year and sex of the standard
  population object. Single ages only.}

  \item{DeathsSex}{the standard Deaths \code{data.frame},
  for a single sex, all years.}

  \item{k}{the parameter 'k' from various equations in the
  section on SRM.}

  \item{m}{the parameter 'm' from various equations in the
  section on SRM.}

  \item{l}{the parameter 'l' from various equations in the
  section on SRM.}

  \item{maxit}{parameter passed to \code{p_sra()}, which is
  iterative. 100 is more than enough.}

  \item{reproduce.matlab}{logical. Do we include legacy
  matlab kludges? Default \code{FALSE}. This effects
  results very little, and is only used for exact matching
  of output.}
}
\value{
  PopYrSex, with open age group redistributed over higher
  ages, up to and including omega, as detected by
  \code{p_ecm_findOmega()}.
}
\description{
  This function mirrors behavior of the matlab function,
  but not the present live function. The live function was
  found to have a few bugs, which were fixed for purposes
  of testing. Since, this function calls \code{p_sra()}, we
  pass the relevant arguments down the chain. There is one
  shortcoming with this implementation, namely that
  \code{p_sra()} only works with Jan 1 in the current
  implementation. The resultant error is likely very small,
  however. This function is called by \code{p_soai()}.
  Further critiques are given in the code comments and
  description of the parent function, \code{p_soai()}.
}

