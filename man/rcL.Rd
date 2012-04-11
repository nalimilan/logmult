\name{rcL}
\alias{rcL}
\alias{print.rcL}
\alias{print.rcL.homog}
\title{Fitting Row-Column Association Models With Layer Effect}
\description{
  \code{rcL} fits log-multiplicative row-column association models with layer effect,
  with one or several dimensions. Supported variants include homogeneous or heterogeneous
  scores over the layer variable, and (for square tables, mainly) homogeneous row and column
  scores, possibly combined with separate diagonal parameters.
}
\usage{
  rcL(tab, nd = 1, layer.homogeneous=c("both", "rows", "columns", "none"),
     homogeneous = FALSE, diagonal = FALSE,
     weights = c("marginal", "uniform", "none"),
     std.err = c("none", "jackknife"), family = poisson,
     start=NULL, tolerance = 1e-12, iterMax = 5000, trace = TRUE, ...)
}
\arguments{
  \item{tab}{ a three-way table, or an object (such as a matrix) that can be coerced into a table;
    if present, dimensions above three will be collapsed as appropriate.}
#FIXME
  \item{nd}{ the number of dimensions to include in the model. Cannot exceed
    \code{min(nrow(tab) - 1, ncol(tab) - 1)} if \code{homogeneous} is \code{FALSE} (saturated model),
    and twice this threshold otherwise (quasi-symmetry model).}
  \item{layer.homogeneous} { what scores should be forced to be equal accross layers?}
  \item{homogeneous}{ should row and column scores be constrained to be equal? For rectangular tables,
    applies to common row and column levels.}
  \item{diagonal}{ should the model include parameters specific to each diagonal cell? This amounts to
    taking quasi-independence, rather than independence, as the baseline model.}
  \item{weights}{ what weights should be used when normalizing the scores.}
  \item{std.err}{ whether to compute standard errors for scores. \dQuote{jackknife} is the only currently
     supported method, and is disabled by default because it increases computation time.}
  \item{family}{ a specification of the error distribution and link function
     to be used in the model. This can be a character string naming
     a family function; a family function, or the result of a call
     to a family function. See \code{\link{family}} details of family functions.}
  \item{start}{a vector of starting values for the parameters in the model;
     if a starting value is \code{NA}, the default starting value will
     be used.}
  \item{tolerance}{ a positive numeric value specifying the tolerance level for
     convergence; higher values will speed up the fitting process, but beware of numerical
     instability of estimated scores!}
  \item{iterMax}{ a positive integer specifying the maximum number of main iterations to perform;
     consider raising this value if your model does not converge.}
  \item{trace}{ a logical value indicating whether the deviance
     should be printed after each iteration.}
  \item{\dots}{ more arguments to be passed to \code{\link{gnm}}}
}
\details{
  This function fits log-multiplicative row-column association models with layer effect, usually called (after
  Wong) RC(M)-L models, typically following the equation:
  \deqn{ log F_{ijk} = \lambda + \lambda^I_i + \lambda^J_j + \lambda^K_k + \sum_{m=1}^M { \phi_{km} \mu_{ikm} \nu_{jkm} } }
  where \eqn{F_ijk} is the frequency observed in the cell at the intersection of row i, column j and layer k of
  \code{tab}, and M the number of dimensions. If \code{layer.homogeneous} is set to \sQuote{rows} (resp. \sQuote{columns}),
  then \eqn{\mu_{ikm} = \mu_{im}} (resp. \eqn{\nu_{ikm} = \nu_{im}} for all layers k; both conditions apply if
  \code{layer.homogeneous} is set to \sQuote{both}. See references for detailed information about the
  variants of the model, the degrees of freedom and the identification constraints applied to the scores.

  Actual model fitting is performed using \code{\link{gnm}}, which implements the Newton-Raphson algorithm.
  This function simply ensures correct start values are used, in addition to allowing for identification
  of scores even with several dimensions, computation of their jackknife standard errors, and plotting.
}
\value{
  A \code{rcL} object, with all the components of a \code{\link{gnm}} object, plus an
    \code{assoc} component holding the most relevant association information:
  \item{phi }{The intrisic association parameters, one per dimension and per layer.}
  \item{row }{Row scores, normalized so that their (weighted) mean is 0, their (weighted)
    sum is 1, and their (weighted) cross-dimensional correlation is null.}
  \item{col }{Column scores, normalized so that their (weighted) mean is 0, their (weighted)
    sum is 1, and their (weighted) cross-dimensional correlation is null.}
  \item{weighting }{The name of the weighting method used, reflected by \code{row.weights}
    and \code{col.weights}.}
  \item{row.weights }{The row weights used for the identification of scores, as specified by the
    \code{weights} argument.}
  \item{col.weights }{The column weights used for the identification of scores, as specified by the
    \code{weights} argument.}
  \item{covmat }{The variance-covariance matrix of the parameter's estimators, if \code{std.err}
    was set. It contains first phi coefficients, then normalized row and column scores, then
    sqrt(phi) * row/column scores.}
  \item{covtype }{The method used to compute the variance-covariance matrix. Currently, only
    \dQuote{jackknife} is supported.}
}
\note{
  Default \code{tolerance} value is very low to ensure the model converges to the right point and that
  coefficients are stable even in most complex cases. This comes at a cost, in particular with high-dimensional
  tables and/or when computing jackknife standard errors. If you know what you are doing, raising \code{tolerance}
  to 1e-6 or an intermediate value can lead to substantial computing time reductions; in that case, be sure to run the
  model several times to check for numerical instability of the scores. (Very large jackknife standard errors are
  a good sign that proper convergence was not attained at the chosen tolerance level.)
}
\references{
  Wong, R.S-K. (2010). Association models. SAGE: Quantitative Applications in the Social Sciences.
}
\author{
  Milan Bouchet-Valat
}
\seealso{
  \code{\link{plot.rcL}}, \code{\link{gnm}}
}
\examples{
  ## Wong (2010), Table 4.7 (p. 103), model 9
  data(gender)
  \dontshow{
      # Limit the workload for CRAN checks
      options(cl.cores=2)
  }
  model <- rcL(gender, nd=2, weights="none", std.err="jackknife")
  model$assoc
  se.rcL(model) # Jackknife standard errors are slightly different than their asymptotic counterparts
  # See ?plot.rcL for plotting
  \dontshow{
      stopifnot(all.equal(round(c(model$assoc$phi), d=3),
                          c(3.075, 3.474,  2.949, 2.460,
                            0.539, 1.686, -0.770, 1.891)))
      stopifnot(all.equal(round(c(model$assoc$row), d=3),
                          c(0.640,  0.239, -0.168, -0.711,
                            0.731, -0.217, -0.636,  0.121)))
      stopifnot(all.equal(round(c(model$assoc$col), d=3),
                          c(0.765,  0.250, -0.398, -0.273, -0.344,
                            0.071, -0.480, -0.198, -0.216,  0.824)))
  }
}

\keyword{ models }
\keyword{ nonlinear }
