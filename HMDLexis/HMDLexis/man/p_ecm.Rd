\name{p_ecm}
\alias{p_ecm}
\title{p_ecm a simple extinct cohort method function, for modularity}
\usage{
  p_ecm(Pop, Deaths, reproduce.matlab = FALSE)
}
\arguments{
  \item{Popualtion}{The standard internal Population
  data.frame, as}

  \item{Deaths}{after all processing is done. Completed
  triangles.}

  \item{reproduce.matlab}{logical, default \code{FALSE}.
  Optional kludge.}
}
\value{
  Population with old ages either filled in or imputed
  using the extinct cohort method
}
\description{
  This implementation seeks to follow the MP. Deaths must
  be finalized. Populations must be in single ages and all
  years. This function covers area B from Figure 6 'Methods
  used for population estimates'.
}

