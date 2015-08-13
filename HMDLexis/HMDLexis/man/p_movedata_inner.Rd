\name{p_movedata_inner}
\alias{p_movedata_inner}
\title{Internal work function of \code{p_movedata()}}
\usage{
  p_movedata_inner(PopL, PopR, detect.mid.year = TRUE,
    detect.start.end = TRUE, reproduce.matlab = FALSE)
}
\arguments{
  \item{PopL}{The left-side population object. A single
  year of data, including both males and females. Single
  ages only. referred to as C1 in MP.}

  \item{PopR}{The right-side population object. A single
  year of data, including both males and females. Single
  ages only. referred to as C2 in MP. Must be adjacent
  years.}

  \item{detect.mid.year}{logical. if \code{TRUE}, June 30
  or July 1 will always return .5.}

  \item{detect.start.end}{logical. default \code{TRUE}.
  Should Jan 1 always be 0 and Dec 31 always be 1?}

  \item{reproduce.matlab}{logical. Default TRUE. Assume 365
  days in a year.}
}
\description{
  See \code{?p_movedata} for some details on the quirks of
  this function.
}

