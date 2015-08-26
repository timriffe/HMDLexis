\name{p_rescale}
\alias{p_rescale}
\title{rescale one population using a second population that has grouped ages, but whose totals are more believable}
\usage{
  p_rescale(Pop, PopN, N = 5)
}
\arguments{
  \item{P1}{a standard HMD Population \code{data.frame} in
  single ages. This one gets rescaled.}

  \item{PN}{a standard HMD Population \code{data.frame},
  that is in N-year age groups. This one is used as the
  standard.}
}
\description{
  Note, any dates that match will get rescaled. Any dates
  that don't match will not get rescaled! For US states, we
  therefore move move midyear official estimates to Jan 1.
  In order to do so they were first split into single ages.
  Here we need them back in 5-year ages. Therefore
  \code{p_agg()} gets run to create \code{P5} just before
  doing this.
}

