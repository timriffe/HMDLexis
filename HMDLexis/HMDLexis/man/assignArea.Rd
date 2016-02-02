\name{assignArea}
\alias{assignArea}
\title{assign the correct Area code to each year}
\usage{
  assignArea(Pnew, PopRef)
}
\arguments{
  \item{Pnew}{an interim Population object created by some
  other function}

  \item{PopRef}{the original Population object from which
  Area codes are to be carried over}
}
\description{
  Some newly created Population data, e.g., made by
  \texttt{p_ecm()} needs to have Area assigned to the new
  data. This function just modulates the way of
  accomplishing this.
}

