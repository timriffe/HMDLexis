\name{d_unk_sex}
\alias{d_unk_sex}
\title{d_unk_sex distribute deaths of unknown sex}
\usage{
  d_unk_sex(Deaths)
}
\arguments{
  \item{Deaths}{the standard internal Deaths data.frame,
  early in the production chain.}
}
\value{
  Deaths rows containing UNK are removed, and these deaths
  have been distributed proportionally over other deaths of
  the same age and known sex.
}
\description{
  d_unk_sex distribute deaths of unknown sex
}

