\name{d_svv}
\alias{d_svv}
\title{a function to split deaths in vertical parallelograms.}
\usage{
  d_svv(Deaths, reproduce.matlab = FALSE)
}
\arguments{
  \item{Deaths}{The standard Deaths data.frame, as
  typically returned by \code{readInputDB}}
}
\value{
  Deaths with VV converted to TL and TU. There will be more
  rows now...
}
\description{
  This method is described in an MP footnote, where we
  state that the method is to just split VV or VH in half.
  That's easy to implement, but the matlab code did much
  more. Since it was never vetted publicly, we do not
  implement the matlab oddities at this time, but perhaps a
  future maintainer will fill in that code chunk. Right now
  the code essentially essentially ignores the
  \code{reproduce.matlab}, argument, following the MP in
  all cases and returning a Warning message if
  \code{reproduce.matlab = TRUE}.
}

