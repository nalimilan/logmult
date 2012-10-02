\name{anoas}
\alias{anoas}
\alias{anoasL}
\title{Analysis of Association Functions}
\description{
  These functions allow performing in a straightforward and efficient way an analysis of association (ANOAS)
  consisting of successive RC(M) or RC(M)-L models from 1 to N dimensions. They fit the models efficiently by
  using scores from the previous model as starting values for the next one.
}
\usage{
anoas(tab, nd = 3, symmetric = FALSE, diagonal = FALSE, ...)

anoasL(tab, nd = 3,
       layer.effect = c("homogeneous.scores", "heterogeneous", "none"),
       symmetric = FALSE,
       diagonal = c("none", "heterogeneous", "homogeneous"), ...)
}
\arguments{
  \item{tab}{ a two-way table, or an object (such as a matrix) that can be coerced into a table;
    if present, dimensions above two will be collapsed as appropriate.}
  \item{nd}{ the number of dimensions to include in the most complex model. Cannot exceed
    \code{min(nrow(tab) - 1, ncol(tab) - 1)} if \code{symmetric} is \code{FALSE} (saturated model),
    and twice this threshold otherwise (quasi-symmetry model).}
  \item{symmetric}{ See \code{\link{rc}} or \code{\link{rcL}}}
  \item{layer.effect}{ See \code{\link{rcL}}.}
  \item{diagonal}{ See \code{\link{rc}} or \code{\link{rcL}}.}
  \item{\dots}{ more arguments to be passed to \code{\link{rc}} or \code{\link{rcL}}.}
}
\details{
  Contrary to most analyses of association in the literature, this function currently does not fit uniform association model
  (\dQuote{U}), nor separate models with only row and column association (\dQuote{R} and \dQuote{C} models), nor log-linear
  row and column association models.

  Currently, no significance test is performed on the models. Please note that it is not correct to test the one-dimension
  association model against the independence model.
}
\value{
  A \code{list} of \code{\link{gnm}} objects. The first element is the independence model, the remaining ones are \code{\link{rc}}
  (for \code{anoas}) or \code{\link{rcL}} (for \code{anoasL}) objects with dimensions from 1 to \code{nd}.
}
\references{
  Wong, R.S-K. (2010). Association models. SAGE: Quantitative Applications in the Social Sciences.
}
\author{
  Milan Bouchet-Valat
}
\seealso{
  \code{\link{rc}}, \code{\link{rcL}}, \code{\link{gnm}}
}
\examples{
  ## Wong (2010), Table 2.6
  data(gss8590)

  # The table used in Wong (2010) is not perfectly consistent
  # with that of Wong (2001)
  tab <- margin.table(gss8590[,,c(2,4)], 1:2)
  tab[2,4] <- 49

  # Results correspond to lines 1, 6 and 11
  results <- anoas(tab, nd=2)
  results

  \dontshow{
      det <- summary(results)
      stopifnot(isTRUE(all.equal(det[[1]], c(12, 6, 2))))
      stopifnot(isTRUE(all.equal(round(det[[2]], 2), c(1373.18, 125.06, 0.60))))
      stopifnot(isTRUE(all.equal(round(det[[4]], 2), c(23.86, 6.44, 0.09))))
      stopifnot(isTRUE(all.equal(round(det[[5]], 2), c(1274.08, 75.51, -15.92))))
      stopifnot(isTRUE(all.equal(round(det[[7]][-1], 2), c(-1248.12, -124.46))))
      stopifnot(isTRUE(all.equal(det[[8]][-1], c(-6, -4))))
  }
}

\keyword{ models }
\keyword{ nonlinear }
