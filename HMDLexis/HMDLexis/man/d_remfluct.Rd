\name{d_remfluct}
\alias{d_remfluct}
\title{d_remfluct a function called by \code{d_soainew()} to remove outliers prior to fitting the Kannisto survival function}
\usage{
  d_remfluct(deathsRR, N = 25, p1 = 5e-04, p2 = 0.9)
}
\arguments{
  \item{deathsRR}{the universal deaths object AFTER all
  ages below the open age group have been split exclusively
  to triangles.}

  \item{N}{the number of years to take when checking first
  differences for smoothness. The MP states 30, but the
  matlab used 25.}

  \item{p1}{the smoothing parameter for the first cubic
  spline over first differences. In the MP and matlab as
  0.0005.}

  \item{p2}{the smoothing parameter for the second cubic
  spline over all ages of death counts, as seen in matlab
  as 0.9.}
}
\value{
  deaths.out a simple data.frame with 3 columns: Age,
  DeathsForFitting, and Original. These are just for the
  purpose of fitting the Kannisto survival function, so
  everything is in RR format. Ages are from OA up to OA-N
  to OA, Deaths are RR, with detected gaps imputed with
  smooth counts, and Original is an unused indicator of
  which counts were kept vs imputed.
}
\description{
  In order to redistribute deaths in the open age group,
  one first fits the Kannisto survival function to an
  imaginary survival function built backward from death
  counts, similar to the extinct cohort method. Prior to
  backward summing the synthetic survival curve, we look at
  deaths to check for a particular variety of unusual
  fluctuation: that caused by death dearths, as happens
  when small cohorts pass through old ages. If such a
  fluctuation is found, deaths are imputed into the
  affected ages using a cubic spline.
}
\details{
  Some comments: The matlab implementation of this whole
  process is not in sync with MPv5, but attempt here to
  reproduce it. The following quirks have been identified
  at this time: 1) the function only removes death troughs,
  but not spikes, 2) fluctuations of longer than 7 years or
  of only one year are not treated, 3) 25 years are taken
  as reference rather than the 30 stated in the MP, 4) the
  first year in the window cannot be an outlier, 5) only
  smoothed values that differ from the original values by
  more than 10\% are kept. These rules may have been
  reasonable for the three test-countries used, but there
  is no reason to think that they are valid for the whole
  HMD. Further diagnostics are not produced on a regular
  basis, and should be part of the country specialist
  checklist. The particular spline used in matlab was
  \code{csaps()}, and the values of the smoothing parameter
  given in the MP, 0.0005, applies to that particular
  function. One can replicate \code{csaps()} in R by using
  \code{pspline::smooth.Pspline()}, with smoothing
  parameter \code{spar} set to $(1-p)/p$.
}

