\name{resortPops}
\alias{resortPops}
\title{resortPops a utility function for sorting the standard Population data.frame}
\usage{
  resortPops(Pop)
}
\arguments{
  \item{Pop}{the standard internal Population data.frame}
}
\value{
  Pop the standard internal Population data.frame, same
  dimensions but sorted
}
\description{
  The sorting hierarchy is: Age within Sex within Month
  within Year.
}

