\name{p_addCohortColumn}
\alias{p_addCohortColumn}
\title{p_addCohortColumn adds Cohort column to Population data.frame}
\usage{
  p_addCohortColumn(Pop)
}
\arguments{
  \item{Population}{the standard internal Population
  object, at any point in processing}
}
\value{
  Population same data.frame with an extra column.
}
\description{
  The cohort column will only be added to those entries
  where the age interval is equal to 1. Sometimes it will
  be easier to match populations using cohorts for certain
  LDB operations. This function is an auxiliary function
  available to all others. We assume Populations are
  January 1 for this function.
}

