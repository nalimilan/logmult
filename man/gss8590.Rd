\name{gss8590}
\alias{gss8590}
\docType{data}
\title{
  Education and Occupational Attainment Among Women in the United States, 1985-1990
}
\description{
  Two-way table crossing education and occupational attainment among women
  from the General Social Survey: United States, 1985-1990. This table is used
  by Wong (2001, 2010) to illustrate a log-multiplicative row-column model with two dimensions.
}
\usage{
data(gss8590)
}
\references{
  Wong, R.S-K. (2001). Multidimensional Association Models : A Multilinear Approach
    \emph{Sociol. Methods & Research} 30, 197-240, Table 2.
  Wong, R.S-K. (2010). \emph{Association models}. SAGE. 32, Table 2.3 B.
}
\examples{
  ## see ?rc and ?plot.rc
  # The table reported in Wong (2010) has a cell inconsistent with
  # what was reported in Wong (2001). To fix this:
  data(gss8590)
  tab <- margin.table(gss8590[,,c(2,4)], 1:2)
  tab[2,4] <- 49
}
\keyword{datasets}
