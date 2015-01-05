\name{p_ic_prexisting}
\alias{p_ic_prexisting}
\title{p_ic_prexisting connects cohorts between two censuses, distributes error.}
\usage{
  p_ic_prexisting(Pop, Deaths, reproduce.matlab = TRUE)
}
\arguments{
  \item{Pop}{the standard HMD internal population
  data.frame. After unknowns are redistributed.}

  \item{Deaths}{the standard HMD internal deaths
  data.frame, after all processing completed.}
}
\value{
  Pop with new lines for intermediate years.
}
\description{
  Census 1 and census 2 must already be Jan 1st. This
  function does not account for newly born cohorts, nor
  will it extend out to age beyond the open age in census
  2. Territorial adjustments are not accounted for in this
  implementation. Censuses must be in single ages, UNK
  already redistributed.
}

