\name{d_unk}
\alias{d_unk}
\title{d_unk distribute deaths of unknown ages}
\usage{
  d_unk(Deaths)
}
\arguments{
  \item{Deaths}{the standard internal Deaths data.frame,
  probably late down the production chain.}
}
\value{
  Deaths rows containing UNK are removed, and these deaths
  have been distributed proportionally over other deaths of
  known age.
}
\description{
  d_unk distribute deaths of unknown ages
}

