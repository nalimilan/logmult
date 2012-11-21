\name{rc}
\alias{rc}
\alias{print.rc}
\alias{print.rc.symm}
\title{Fitting Row-Column Association Models}
\description{
  Fit log-multiplicative row-column association models, also called RC(M) models or Goodman's (1979) Model II,
  with one or several dimensions. Supported variants (for square tables) include symmetric (homogeneous)
  row and column scores, possibly combined with separate diagonal parameters.
}
\usage{
rc(tab, nd = 1, symmetric = FALSE, diagonal = FALSE,
   weighting = c("marginal", "uniform", "none"),
   rowsup = NULL, colsup = NULL,
   se = c("none", "jackknife", "bootstrap"),
   nreplicates = 100, ncpus = getOption("boot.ncpus"),
   family = poisson, weights = NULL,
   start = NA, etastart = NULL, tolerance = 1e-8,
   iterMax = 5000, trace = FALSE, verbose = TRUE, ...)
}
\arguments{
  \item{tab}{ a two-way table, or an object (such as a matrix) that can be coerced into a table;
    if present, dimensions above two will be collapsed.}
  \item{nd}{ the number of dimensions to include in the model. Cannot exceed
    \code{min(nrow(tab) - 1, ncol(tab) - 1)} if \code{symmetric} is \code{FALSE} (saturated model),
    and twice this threshold otherwise (quasi-symmetry model).}
  \item{symmetric}{ should row and column scores be constrained to be equal? Valid only for square tables.}
  \item{diagonal}{ should the model include parameters specific to each diagonal cell? This amounts to
    taking quasi-independence, rather than independence, as the baseline model. Valid only for square tables.}
  \item{weighting}{ what weights should be used when normalizing the scores.}
  \item{rowsup}{ if present, a matrix with the same columns as \code{tab} giving supplementary (passive) rows.
                If \code{symmetric = TRUE}, \code{rowsup} and \code{colsup} must be specified together and
                rows of \code{rowsup} must correspond to columns of \code{colsup}.}
  \item{colsup}{ if present, a matrix with the same rows as \code{tab} giving supplementary (passive) columns.
                See \code{rowsup}.}
  \item{se}{ which method to use to compute standard errors for parameters (see \code{\link{se.assoc}}).}
  \item{nreplicates}{ the number of bootstrap replicates, if enabled.}
  \item{ncpus}{ the number of processes to use for jackknife or bootstrap parallel computing. Defaults to
    the number of cores (see \code{\link[parallel]{detectCores}}), with a maximum of 5, but falls back to 1
    (no parallelization) if packages \code{parallel} or \code{snow} are not available.}
  \item{family}{ a specification of the error distribution and link function
     to be used in the model. This can be a character string naming
     a family function; a family function, or the result of a call
     to a family function. See \code{\link{family}} details of family functions.}
  \item{weights}{ an optional vector of weights to be used in the fitting process.}
  \item{start}{either \code{NA} (the default) to use optimal starting values, \code{NULL} to use
     random starting values, or a vector of starting values for the parameters in the model.}
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
  This function fits log-multiplicative row-column association models, usually called (after
  Goodman) RC(M) models, typically following the equation:
  \deqn{ log F_{ij} = \lambda + \lambda^I_i + \lambda^J_j + \sum_{m=1}^M { \phi_{m} \mu_{im} \nu_{jm} } }
  where \eqn{F_{ij}} is the expected frequency for the cell at the intersection of row i and column j of
  \code{tab}, and M the number of dimensions. See references for detailed information about the
  variants of the model, the degrees of freedom and the identification constraints applied to the scores.

  Actual model fitting is performed using \code{\link{gnm}}, which implements the Newton-Raphson algorithm.
  This function simply ensures correct start values are used, in addition to allowing for identification
  of scores even with several dimensions, computation of their jackknife or bootstrap standard errors, and plotting.
  The default starting values for association parameters are computed using a singular/eigen value decomposition
  from the results of the model without association component (\dQuote{base model}). In some complex cases, using
  \code{start = NULL} to start with random values can be more efficient, but it is also less stable and can converge
  to non-optimal solutions.
}
\value{
  A \code{rc} object, with all the components of a \code{\link{gnm}} object, plus an
    \code{assoc.rc} component holding the most relevant association information:
  \item{phi }{The intrisic association parameters, one per dimension.}
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
  \item{covmat }{The variance-covariance matrix for phi coefficients and normalized row and column
    scores. Only present if \code{se} was not \dQuote{none}.}
  \item{adj.covmats }{An array stacking on its third dimension one variance-covariance matrix for
    the adjusted scores of each layer in the model (used for plotting). Only present if \code{se}
    was not \dQuote{none}.}
  \item{covtype }{The method used to compute the variance-covariance matrix (corresponding to the
    \code{se} argument.}
}
\references{
  Goodman, L.A. (1979). Simple Models for the Analysis of Association in Cross-Classifications
    having Ordered Categories. \emph{J. of the Am. Stat. Association} 74(367), 537-552.

  Becker, M.P., and Clogg, C.C. (1989). Analysis of Sets of Two-Way Contingency Tables Using
    Association Models. \emph{Journal of the American Statistical Association} 84(405), 142-151.

  Goodman, L.A. (1985). The Analysis of Cross-Classified Data Having Ordered and/or Unordered
    Categories: Association Models, Correlation Models, and Asymmetry Models for Contingency
    Tables With or Without Missing Entries. \emph{The Annals of Statistics} 13(1), 10-69.

  Goodman, L.A. (1991). Measures, Models, and Graphical Displays in the Analysis of
    Cross-Classified Data. \emph{J. of the Am. Stat. Association} 86(416), 1085-1111.

  Clogg, C.C., and Shihadeh, E.S. (1994). Statistical Models for Ordinal Variables. Sage: Advanced
      Quantitative Techniques in the Social Sciences (4).

  Wong, R.S-K. (2010). Association models. SAGE: Quantitative Applications in the Social Sciences.
}
\author{
  Milan Bouchet-Valat
}
\seealso{
  \code{\link{plot.rc}}, \code{\link{gnm}}
}
\examples{
  ## Goodman (1991), Table 17.1 (p. 1097)
  data(criminal)
  model <- rc(criminal)

  model$assoc # These are the phi (.07), mu and nu
  model$assoc$row[,1,1] * model$assoc$phi[1,1] # These are the mu'
  model$assoc$col[,1,1] * model$assoc$phi[1,1] # These are the nu'

  ## Becker & Clogg (1989), Table 5 (p. 145)
  # See also ?rcL to run all models in one call
  \dontrun{
  data(color)

  # "Uniform weights" in the authors' terms mean "no weighting" for us
  # See ?rcL for average marginals
  caithness.unweighted <- rc(color[,,1], nd=2, weighting="none",
                             se="jackknife")
  caithness.marginal <- rc(color[,,1], nd=2, weighting="marginal",
                           se="jackknife")
  aberdeen.unweighted <- rc(color[,,2], nd=2, weighting="none",
                            se="jackknife")
  aberdeen.marginal <- rc(color[,,2], nd=2, weighting="marginal",
                          se="jackknife")

  caithness.unweighted
  caithness.marginal
  aberdeen.unweighted
  aberdeen.marginal

  # To see standard errors, either:
  se(caithness.unweighted)

  # and so on...
  # (ours are much smaller for the marginal-weighted case)
  # Or:
  summary(caithness.unweighted)
  }


  ## Clogg & Shihadeh (1994), Tables 3.5a and b (p. 55-61)
  data(gss88)
  model <- rc(gss88)

  # Unweighted scores
  summary(model, weighting="none")
  # Marginally weighted scores
  summary(model, weighting="marginal")
  # Uniformly weighted scores
  summary(model, weighting="uniform")


  ## Wong (2010), Table 2.7 (p. 48-49)
  \dontrun{
  data(gss8590)

  # The table used in Wong (2010) is not perfectly consistent
  # with that of Wong (2001)
  tab <- margin.table(gss8590[,,c(2,4)], 1:2)
  tab[2,4] <- 49

  model <- rc(tab, nd=2, weighting="none", se="jackknife")

  model
  summary(model) # Jackknife standard errors are slightly different
                 # from their asymptotic counterparts

  # Compare with bootstrap standard errors
  model2 <- rc(tab, nd=2, weighting="none", se="bootstrap")
  plot(model, conf.ellipses=0.95)
  summary(model2)
  }
}

\keyword{ models }
\keyword{ nonlinear }
