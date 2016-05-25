\name{computeLDB}
\alias{computeLDB}
\title{compute the standard LDB object for males and females}
\usage{
  computeLDB(Births, Deaths, Pop, Tadj, verbose = FALSE)
}
\arguments{
  \item{Births}{the standard R Births \code{data.frame}.
  All anomalies should be accounted for before this
  function.}

  \item{Deaths}{the standard R Deaths \code{data.frame},
  fully processed and complete.}

  \item{Pop}{the standard R Population \code{data.frame},
  fully processed and complete.}

  \item{Tadj}{the standard R Territorial Adjustment
  \code{data.frame}, as returned by \code{readInputDB()}.}

  \item{verbose.}{logical. do you want lots of info to go
  to the console?}
}
\description{
  the LDB object has columns for Year, Age, Lexis (1,2),
  Population, and Deaths. Lexis = 2 means upper triangle
  deaths and Jan 1 Populations. Lexis = 1 means lower
  triangle deaths and lower age bound birthdays (the bottom
  horizontal line of the Lexis square). This is one way to
  capture births for Age = 0, but Lexis = 1 Population
  counts are otherwise not used in any calculations. TR has
  proposed to instead put DL, DU in columns, and Pop Jan 1,
  Pop Dec 31 in columns, which would obviate the need for
  Tadj from the InputDB downstream. Further, new exposure
  formulas count be accomodated by sd and bbar columns,
  rather than calculating on the fly, then no need to read
  in Monthly births from InputDB for LifeTables calcs. In
  either case, all R \code{data.frame} arguments should be
  fully processed before getting to this stage.
}

