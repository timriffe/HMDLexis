\name{ypart}
\alias{ypart}
\title{ypart function to determine the proportion of a year passed as of a particular date}
\usage{
  ypart(Year, Month, Day, reproduce.matlab = TRUE)
}
\arguments{
  \item{Year}{4-digit year (string or integer)}

  \item{Month}{month digits (string or integer, 1 or 2
  characters)}

  \item{Day}{Day of month digits (string or integer, 1 or 2
  characters)}

  \item{reproduce.matlab}{logical. Default TRUE. Assume 365
  days in a year.}
}
\description{
  The fraction returned by this is used e.g. for
  intercensal estimates. Function uses 'lubridate' package
  to handle dates elegantly.
}

