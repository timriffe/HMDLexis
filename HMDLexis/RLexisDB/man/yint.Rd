\name{yint}
\alias{yint}
\title{yint get interval as fraction of full years}
\usage{
  yint(Day1, Month1, Year1, Day2, Month2, Year2,
    reproduce.matlab = TRUE)
}
\arguments{
  \item{Day1}{Day of first date}

  \item{Month1}{Month of first date}

  \item{Year1}{Year of first date}

  \item{Day1}{Day of second date}

  \item{Month1}{Month of second date}

  \item{Year1}{Year of second date}

  \item{reproduce.matlab}{logical. default \code{TRUE}.
  Assume 365 days in all years?}
}
\value{
  decimal value of year fraction (can be greater than 1)
}
\description{
  Either assume 365 days in the year, or get the precise
  duration.
}

