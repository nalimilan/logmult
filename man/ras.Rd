\name{ras}
\alias{ras}
\title{
  RAS/Deming-Stephan Algorithm for Raking Tables
}
\description{
  Adjust a table by multiplying rows and columns in order to reproduce the provided margins,
  preserving all the odds ratios. This procedure is know as the RAS or Deming-Stephan algorithm,
  as iterative proportional fitting (IPF) or as biproportional fitting.
}
\usage{
ras(tab, row, col, tolerance = .Machine$double.eps)
}
\arguments{
  \item{tab}{a two-way table with only positive or zero entries}
  \item{row}{a vector of strictly positive elements containing the wanted row margins (sums).}
  \item{col}{a vector of strictly positive elements containing the wanted column margins (sums).}
  \item{tolerance}{the convergence criterion to stop iterating.}
}
\details{
  Note that \code{sum(row)} must be equal to \code{sum(col)} for the algorithm to make sense.
}
\value{
  The adjusted table with row sums equal to \code{row} and column sums equal to \code{col}.
}
\author{
  Milan Bouchet-Valat
}
