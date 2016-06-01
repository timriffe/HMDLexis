\name{p_Date}
\alias{p_Date}
\title{make a Date class column for Population}
\usage{
  p_Date(Pop)
}
\arguments{
  \item{Pop}{the standard Pop object, ideally with open age
  split, e.g., after running \code{p_srecm()} on the data.}
}
\value{
  Pop the same object, but with a date class column added.
}
\description{
  The need for this is evident in cases where we have more
  than one date per year. This can happen, for example, if
  the most recent Input data year is a mid-year estimate.
  In that case, if \code{p_movedata()} is run before
  \code{p_postcensal}, there will be two dates in the most
  recent year, and stuff can misbehave.
}

