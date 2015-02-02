\name{d_sNx1}
\alias{d_sNx1}
\title{split RR deaths in age groups to single ages}
\usage{
  d_sNx1(Deaths)
}
\arguments{
  \item{Deaths}{the standard LexisDB Deaths object.}
}
\value{
  Deaths, with all RR entries in single ages.
}
\description{
  This function iterates over each year and sex of data,
  isolating RR data (grouping to RR in the case of infants)
  and then splitting age intervals greater than one to
  single ages. All data in other Lexis formats are
  unchanged, as are pre-existing single-age data in RR
  format. From an error standpoint, this function may be
  run either before or after \code{d_unk()}. The actual
  work is done by the inner function,
  \code{d_sNx1_inner()}, which itself calls
  \code{d_sNx1_spline()} once data are prepared.
}

