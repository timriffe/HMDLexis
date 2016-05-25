\name{p_movedata}
\alias{p_movedata}
\title{Move a series of midyear single-age populations estimates to January 1.}
\usage{
  p_movedata(Pop, detect.mid.year = TRUE,
    detect.start.end = TRUE, reproduce.matlab = FALSE,
    OPENAGE = 130)
}
\arguments{
  \item{Pop}{The standard HMD population \code{data.frame}}

  \item{detect.mid.year}{logical. if \code{TRUE}, June 30
  or July 1 will always return .5.}

  \item{detect.start.end}{logical. default \code{TRUE}.
  Should Jan 1 always be 0 and Dec 31 always be 1?}

  \item{reproduce.matlab}{logical. Default TRUE. Assume 365
  days in a year.}
}
\description{
  This function uses linear interpolation, per the MP, to
  estimate January 1 population. Interpolation takes place
  between like-ages in adjacent years. Time is calculated
  either assuming 365 days per year (\code{reproduce.matlab
  = TRUE}) or with leap-year detection
  (\code{reproduce.matlab = FALSE}). This function covers
  the behavior of the matlab \code{p_movedata()} and
  \code{p_my()}, depending on the settings of
  \code{detect.mid.year} and \code{detect.start.end}. For
  \code{p_my()}, set these two arguments to \code{TRUE},
  otherwise \code{FALSE}. It seems, however, that when a
  statistical office declare June 30 or July 1 that these
  are both intended to be midyear estimates, and so TR
  thinks it wise to just set the defaults to \code{TRUE}
  and not change them. \code{reproduce.matlab} only affects
  the leapyear aspect of estimating year proportions. See
  \code{?ypart} for more info on that. In general, TR
  thinks it would be best to replace this method with a)
  \code{p_ic()} for single-year spans or b) a version of
  \code{p_ic()} for single-year spans that accounts for
  within-cohort variation in the population distribution,
  e.g., from the monthly birth distribution of a given
  cohort. One quirk of this function is that Jan 1 of year
  t is considered adjacent to Dec 31 of year t+1. This
  means that it makes a difference if you run
  \code{p_ey2ny()} before this or not!
}

