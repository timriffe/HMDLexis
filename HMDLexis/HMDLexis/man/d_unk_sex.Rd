% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/d_unk_sex.R
\name{d_unk_sex}
\alias{d_unk_sex}
\title{d_unk_sex distribute deaths of unknown sex}
\usage{
d_unk_sex(Deaths)
}
\arguments{
\item{Deaths}{the standard internal Deaths data.frame, early in the production chain.}
}
\value{
Deaths rows containing UNK are removed, and these deaths have been distributed proportionally over other deaths of the same Year, Age, Lexis and AgeInterval and known Sex. If this is not possible for the given ID combination, we split the deaths equally between the sexes.
}

