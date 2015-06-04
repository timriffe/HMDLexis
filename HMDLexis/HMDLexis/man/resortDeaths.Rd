\name{resortDeaths}
\alias{resortDeaths}
\title{resortDeaths a utility function for sorting the standard Deaths data.frame}
\usage{
  resortDeaths(Deaths)
}
\arguments{
  \item{Deaths}{the standard internal Deaths data.frame}
}
\value{
  Deaths the standard internal Deaths data.frame, same
  dimensions but sorted
}
\description{
  The sorting hierarchy is: Lexis within Age within Sex
  within Year.
}

