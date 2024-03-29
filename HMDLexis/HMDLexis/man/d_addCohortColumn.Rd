% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/d_addCohortColumn.R
\name{d_addCohortColumn}
\alias{d_addCohortColumn}
\title{d_addCohortColumn adds Cohort column to Deaths object}
\usage{
d_addCohortColumn(Deaths, TLTUonly = TRUE)
}
\arguments{
\item{Deaths}{the standard Deaths object, at any stage of production}

\item{TLTUonly}{logical, Default = \code{TRUE}. Should Cohorts only be given to triangles?}
}
\description{
This will assign an integer cohort to TL, TU, VV, VH and RV entries. Optionally, it will only do so for triangles. RR and NA Lexis shapes get NA for cohort.
}

