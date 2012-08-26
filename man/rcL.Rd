\name{rcL}
\alias{rcL}
\alias{print.rcL}
\alias{print.rcL.symm}
\title{Fitting Row-Column Association Models With Layer Effect}
\description{
  \code{rcL} fits log-multiplicative row-column association models with layer effect,
  with one or several dimensions. Supported variants include homogeneous or heterogeneous
  scores over the layer variable, and (for square tables) symmetric (homogeneous)
  row and column scores, possibly combined with separate diagonal parameters.
}
\usage{
rcL(tab, nd = 1,
    layer.effect = c("homogeneous.scores", "heterogeneous", "none"),
    symmetric = FALSE,
    diagonal = c("none", "heterogeneous", "homogeneous"),
    weighting = c("marginal", "uniform", "none"),
    se = c("none", "jackknife"), family = poisson,
    start = NA, etastart = NULL, tolerance = 1e-6,
    iterMax = 5000, trace = TRUE, verbose = TRUE, ...)
}
\arguments{
  \item{tab}{ a three-way table, or an object (such as a matrix) that can be coerced into a table;
    if present, dimensions above three will be collapsed as appropriate.}
  \item{nd}{ the number of dimensions to include in the model. Cannot exceed
    \code{min(nrow(tab) - 1, ncol(tab) - 1)} if \code{symmetric} is \code{FALSE} (saturated model),
    and twice this threshold otherwise (quasi-symmetry model).}
  \item{layer.effect}{ determines the form of the interaction between row-column association and layers.
     See \dQuote{Details} below.}
  \item{symmetric}{ should row and column scores be constrained to be equal? Valid only for square tables.}
  \item{diagonal}{ what type of diagonal-specific parameters to include in the model, if any. This amounts to
    taking quasi-conditional independence, rather than conditional independence, as the baseline model.
    Valid only for square tables.}
  \item{weighting}{ what weights should be used when normalizing the scores.}
  \item{se}{ whether to compute standard errors for scores. \dQuote{jackknife} is the only currently
     supported method, and is disabled by default because it increases computation time.}
  \item{family}{ a specification of the error distribution and link function
     to be used in the model. This can be a character string naming
     a family function; a family function, or the result of a call
     to a family function. See \code{\link{family}} details of family functions.}
  \item{start}{either \code{NA} (the default) to use reasonable starting values, \code{NULL} to use
     random starting values, or a vector of starting values for the parameters in the model;
     if a starting value is \code{NA}, the default starting value will be used.}
  \item{etastart}{ starting values for the linear predictor; set to \code{NULL} to use either default
     starting values (if \code{start = NA}), or random starting values (in all other cases).}
  \item{tolerance}{ a positive numeric value specifying the tolerance level for
     convergence; higher values will speed up the fitting process, but beware of numerical
     instability of estimated scores!}
  \item{iterMax}{ a positive integer specifying the maximum number of main iterations to perform;
     consider raising this value if your model does not converge.}
  \item{trace}{ a logical value indicating whether the deviance
     should be printed after each iteration.}
  \item{verbose}{ a logical value indicating whether progress indicators should be printed,
     including a diagnostic error message if the algorithm restarts.}
  \item{\dots}{ more arguments to be passed to \code{\link{gnm}}}
}
\details{
  This function fits log-multiplicative row-column association models with layer effect, usually called (after
  Wong) RC(M)-L models, typically following the equation:
  \deqn{ log F_{ijk} = \lambda + \lambda^I_i + \lambda^J_j + \lambda^K_k + \sum_{m=1}^M { \phi_{km} \mu_{ikm} \nu_{jkm} } }
  where \eqn{F_ijk} is the frequency observed in the cell at the intersection of row i, column j and layer k of
  \code{tab}, and M the number of dimensions. If \code{layer.effect} is set to \sQuote{heterogeneous}, different scores
  will be computed for each level, which is equivalent to fitting separate RC(M) models on the k two-way tables.
  If it is set to \sQuote{homogeneous.scores}, then \eqn{\mu_{ikm} = \mu_{im}} and \eqn{\nu_{ikm} = \nu_{im}} for all
  layers k: only the \eqn{\phi_{km}} are allowed to vary accross layers. If it is set to \sQuote{none}, then in addition
  to the previous conditions all \eqn{\phi_{km}} are forced to be equal for all layers k, which amounts to a stability
  of the association accross layers. See references for detailed information about the variants of the model, the degrees
  of freedom and the identification constraints applied to the scores.

  Actual model fitting is performed using \code{\link{gnm}}, which implements the Newton-Raphson algorithm.
  This function simply ensures correct start values are used, in addition to allowing for identification
  of scores even with several dimensions, computation of their jackknife standard errors, and plotting.
  The default starting values correspond to the parameters computed for the model without association parameters
  (\dQuote{base model}); association parameters start with random starting values. In some complex cases,
  using \code{start = NA} to get completely random starting values can be more efficient, but it is also
  less stable.
}
\value{
  A \code{rcL} object, with all the components of a \code{\link{gnm}} object, plus an
    \code{assoc.rcL} component holding the most relevant association information:
  \item{phi }{The intrisic association parameters, one per dimension and per layer.}
  \item{row }{Row scores, normalized so that their (weighted) sum is 0, their (weighted)
    sum of squares is 1, and their (weighted) cross-dimensional correlation is null.}
  \item{col }{Column scores, normalized so that their (weighted) sum is 0, their (weighted)
    sum of squares is 1, and their (weighted) cross-dimensional correlation is null.}
  \item{weighting }{The name of the weighting method used, reflected by \code{row.weights}
    and \code{col.weights}.}
  \item{row.weights }{The row weights used for the identification of scores, as specified by the
    \code{weighting} argument.}
  \item{col.weights }{The column weights used for the identification of scores, as specified by the
    \code{weighting} argument.}
  \item{covmat }{The variance-covariance matrix of the parameter's estimators, if \code{se}
    was set. It contains first phi coefficients, then normalized row and column scores, then
    sqrt(phi) * row/column scores.}
  \item{covtype }{The method used to compute the variance-covariance matrix. Currently, only
    \dQuote{jackknife} is supported.}
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
  data(gss7590)
  \dontshow{
      # Limit the workload for CRAN checks
      options(cl.cores=2)
  }
  model <- rcL(gss7590, nd=2, weighting="none", se="jackknife")
  model
  se(model) # Jackknife standard errors are slightly different
                # from their asymptotic counterparts
  # See ?plot.rcL for plotting
  \dontshow{
      stopifnot(all.equal(round(c(model$assoc$phi), d=3),
                          c(3.075, 3.474,  2.949, 2.460,
                            0.539, 1.686, -0.770, 1.891)))
      stopifnot(all.equal(round(c(model$assoc$row[,,1]), d=3),
                          c(0.640,  0.239, -0.168, -0.711,
                            0.731, -0.217, -0.636,  0.121)))
      stopifnot(all.equal(round(c(model$assoc$col[,,1]), d=3),
                          c(0.765,  0.250, -0.398, -0.273, -0.344,
                            0.071, -0.480, -0.198, -0.216,  0.824)))

     # Check scores of heterogeneous model (scores not given by Wong)
     model.heterog <- rcL(gss7590, nd=2, layer.effect="heterogeneous",
                          weighting="none", tolerance=1e-8)
     sep.models <- lapply(1:4, function(i)
         rc(gss7590[,,i], nd=2, weighting="none", tolerance=1e-8))
     # Level 2 does not converge properly because there are too few farmers:
     # values are always slightly different and cannot be checked
     stopifnot(isTRUE(all.equal(model.heterog$assoc$phi[-2,],
                                rbind(sep.models[[1]]$assoc$phi,
                                      sep.models[[3]]$assoc$phi,
                                      sep.models[[4]]$assoc$phi),
                                check.attributes=FALSE)))
     stopifnot(isTRUE(all.equal(model.heterog$assoc$row[,,-2],
                                array(c(sep.models[[1]]$assoc$row,
                                      sep.models[[3]]$assoc$row,
                                      sep.models[[4]]$assoc$row),
                                      dim=c(nrow(gss7590), 2, 3)),
                                check.attributes=FALSE)))
     stopifnot(isTRUE(all.equal(model.heterog$assoc$col[,,-2],
                                array(c(sep.models[[1]]$assoc$col,
                                      sep.models[[3]]$assoc$col,
                                      sep.models[[4]]$assoc$col),
                                      dim=c(ncol(gss7590), 2, 3)),
                                check.attributes=FALSE)))
  }
}

\keyword{ models }
\keyword{ nonlinear }
