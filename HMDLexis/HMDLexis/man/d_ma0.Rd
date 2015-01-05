\name{d_ma0}
\alias{d_ma0}
\title{d_mao a function to add missing ages in deaths with zero counts}
\usage{
  d_ma0(Deaths, sa = 90)
}
\value{
  Deaths the standard Deaths data.frame, with 0s imputed as
  necessary
}
\description{
  The function loops a lot to fill missing age-lexis combos
  with 0s
}

