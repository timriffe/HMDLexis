\name{isLeapYear}
\alias{isLeapYear}
\title{a function to determine whether a year is a leap year.}
\usage{
  isLeapYear(Year)
}
\arguments{
  \item{Year}{integer of year to query}
}
\value{
  logical is the Year a leap year or not
}
\description{
  In order to remove lubridate dependency, we self-detect
  leap years and adjust February accordingly.
}

