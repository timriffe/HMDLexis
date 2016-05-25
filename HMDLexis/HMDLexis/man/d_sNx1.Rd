\name{d_sNx1}
\alias{d_sNx1}
\title{split RR deaths in age groups to single ages}
\usage{
  d_sNx1(Deaths, MPVERSION = 5)
}
\arguments{
  \item{Deaths}{the standard LexisDB Deaths object.}

  \item{MPVERSION}{which version of the methods protocol
  (experimental)}
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
  \code{d_sNx1_spline()} once data are prepared. If
  \code{MPVERSION < 7} we reproduce the matlab output
  exactly (save for a minor infant triangles bug fix). If
  \code{MPVERSION >= 7} we use a monotonic hyman spline,
  which solves some problems inherent in the MP method (as
  of this writing). Refer to emails from Tim Riffe in
  January, 2015 for more details and case examples.
}

