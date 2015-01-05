\name{d_agg}
\alias{d_agg}
\title{a function to aggregate late-registered deaths into the appropriate age and year categories}
\usage{
  d_agg(Deaths.In)
}
\arguments{
  \item{Deaths.In}{This should be the standard Deaths
  data.frame, as typically returned by
  \code{readInputDB()}.}
}
\description{
  InputDB files typically contain a \code{"YearReg"}
  column, but it does not appear to be used either
  implicitly or explicitly in the LexisDB files. Further,
  oftentimes when a particular combination of
  Year-Sex-Age-Lexis is extracted from the data in the
  matlab code it is assumed that this will be a single
  number. In the case of additional late registered events,
  this operation will in fact return a vector, and could
  lead to unexpected results. This function is written to
  explicitly sum-in late registered events prior to further
  matlab operations, and can be considered a pre-processing
  step until such time as the matlab is phased out. In
  order to combine matlab and R in this way, we've decided
  to name the deaths file in the original raw state as
  \code{XYZdeath_pre.txt}, and the output of this function
  as the standard \code{XYZdeath.txt} so that the matlab
  programs will be able to find it. This is a temporary
  kludge for JPN and FIN, hopefully, although there are
  other countries where this procedure may be needed in the
  future, such as Finland or Australia, depending on how
  data come. This function may have other unexpected
  side-effects where redundant deaths have not been
  properly dealt with in the matlab code. For now, this
  function will be called on the CS-side, prior to
  submitting the InputDB, as this function needs to be used
  prior to completion of the LexisDB R programs.
}

