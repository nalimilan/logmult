\name{rc}
\alias{rc}
\alias{print.rc}
\alias{print.rc.symm}
\title{Fitting Row-Column Association Models}
\description{
  \code{rc} fits log-multiplicative row-column association models first introduced
  by Goodman (1979), with one or several dimensions. Supported variants (for square tables)
  include symmetric (homogeneous) row and column scores, and separate diagonal parameters.
}
\usage{
  rc(tab, nd = 1, symmetric = FALSE, diagonal = FALSE,
     weighting = c("marginal", "uniform", "none"),
     std.err = c("none", "jackknife"), family = poisson,
     start = NA, tolerance = 1e-12, iterMax = 5000, trace = TRUE, ...)
}
\arguments{
  \item{tab}{ a two-way table, or an object (such as a matrix) that can be coerced into a table;
    if present, dimensions above two will be collapsed as appropriate.}
  \item{nd}{ the number of dimensions to include in the model. Cannot exceed
    \code{min(nrow(tab) - 1, ncol(tab) - 1)} if \code{symmetric} is \code{FALSE} (saturated model),
    and twice this threshold otherwise (quasi-symmetry model).}
  \item{symmetric}{ should row and column scores be constrained to be equal? Valid only for square tables.}
  \item{diagonal}{ should the model include parameters specific to each diagonal cell? This amounts to
    taking quasi-independence, rather than independence, as the baseline model. Valid only for square tables.}
  \item{weighting}{ what weights should be used when normalizing the scores.}
  \item{std.err}{ whether to compute standard errors for scores. \dQuote{jackknife} is the only currently
     supported method, and is disabled by default because it increases computation time.}
  \item{family}{ a specification of the error distribution and link function
     to be used in the model. This can be a character string naming
     a family function; a family function, or the result of a call
     to a family function. See \code{\link{family}} details of family functions.}
  \item{start}{either \code{NA} (the default) to use reasonable starting values, \code{NULL} to use
     random starting values, or a vector of starting values for the parameters in the model;
     if a starting value is \code{NA}, the default starting value will be used.}
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
  This function fits log-multiplicative row-column association models, usually called (after
  Goodman) RC(M) models, typically following the equation:
  \deqn{ log F_{ij} = \lambda + \lambda^I_i + \lambda^J_j + \sum_{m=1}^M { \phi_{m} \mu_{im} \nu_{jm} } }
  where \eqn{F_ij} is the frequency observed in the cell at the intersection of row i and column j of
  \code{tab}, and M the number of dimensions. See references for detailed information about the
  variants of the model, the degrees of freedom and the identification constraints applied to the scores.

  Actual model fitting is performed using \code{\link{gnm}}, which implements the Newton-Raphson algorithm.
  This function simply ensures correct start values are used, in addition to allowing for identification
  of scores even with several dimensions, computation of their jackknife standard errors, and plotting.
  The default starting values correspond to the parameters computed for the model without association
  parameters (\dQuote{base model}); association parameters start with values estimated from the base model
  residuals if \code{symmetric = FALSE} (see \code{\link{residSVD}}), and random values otherwise.
  In some complex cases, using \code{start = NA} to get completely random starting values can be more
  efficient, but it is also less stable.
}
\value{
  A \code{rc} object, with all the components of a \code{\link{gnm}} object, plus an
    \code{assoc.rc} component holding the most relevant association information:
  \item{phi }{The intrisic association parameters, one per dimension.}
  \item{row }{Row scores, normalized so that their (weighted) mean is 0, their (weighted)
    sum is 1, and their (weighted) cross-dimensional correlation is null.}
  \item{col }{Column scores, normalized so that their (weighted) mean is 0, their (weighted)
    sum is 1, and their (weighted) cross-dimensional correlation is null.}
  \item{weighting }{The name of the weighting method used, reflected by \code{row.weights}
    and \code{col.weights}.}
  \item{row.weights }{The row weights used for the identification of scores, as specified by the
    \code{weighting} argument.}
  \item{col.weights }{The column weights used for the identification of scores, as specified by the
    \code{weighting} argument.}
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
  \dontshow{
      ## Goodman (1979), scores not reported but present in LEM examples
      # as "GOO91_1G" (where scores reported here come from)
      # No symmetric RC scores seem to be reported in a published paper
      data(occupationalStatus)
      # We use a simpler (collapsed) table
      occupationalStatus[5,]<-colSums(occupationalStatus[5:6,])
      occupationalStatus[,5]<-rowSums(occupationalStatus[,5:6])
      occupationalStatus <- occupationalStatus[-6,-6]

      model <- rc(occupationalStatus, diagonal=TRUE, symmetric=TRUE,
                  weighting="none")
      stopifnot(round(model$assoc$phi[1,1], d=3) == 6.10)
      stopifnot(isTRUE(all.equal(round(c(model$assoc$row), d=3),
                                 c(0.532, 0.438, 0.206, -0.031,
                                   -0.216, -0.426, -0.503))))

      model <- rc(occupationalStatus, diagonal=TRUE, symmetric=TRUE,
                  weighting="uniform")
      stopifnot(round(model$assoc$phi[1,1], d=3) == 0.871)
      stopifnot(isTRUE(all.equal(round(c(model$assoc$row), d=3),
                                 c(1.409, 1.159, 0.544, -0.082,
                                  -0.571, -1.127, -1.332))))

      # Marginal weighted scores from ass2 in LEM are weird: we need
      # to run cor(3) on the fitted table to get correct values, which is not
      # possible for symmetric scores
  }

  ## Goodman (1991), Table 17.1 (p. 1097)
  data(criminal)
  model <- rc(criminal)
  model$assoc # These are the phi (.07), mu and nu
  model$assoc$row[,1,1] * model$assoc$phi[1,1] # These are the mu'
  model$assoc$col[,1,1] * model$assoc$phi[1,1] # These are the nu'
  \dontshow{
      stopifnot(round(model$assoc$phi[1,1], d=2) == 0.07)
      # Scores have reversed signs compared to the orignal article
      stopifnot(isTRUE(all.equal(round(model$assoc$row[,1,1], d=2),
                                 c(1.26, 0.82, -0.56, -1.21),
                                 check.attributes=FALSE)))
      stopifnot(isTRUE(all.equal(round(model$assoc$col[,1,1], d=2),
                                 c(-1.44, -1.30, -0.33, 1.00, 0.89),
                                 check.attributes=FALSE)))
  }

  ## Becker & Clogg (1989), Table 5 (p. 145)
  data(color)
  \dontshow{
      # Limit the workload for CRAN checks
      options(cl.cores=2)
  }
  # "Uniform weights" in the authors' terms mean "no weighting"
  # Average marginals are currently not supported
  caithness.unweighted <- rc(color[,,1], nd=2, weighting="none",
                             std.err="jackknife")
  caithness.marginal <- rc(color[,,1], nd=2, weighting="marginal",
                           std.err="jackknife")
  aberdeen.unweighted <- rc(color[,,2], nd=2, weighting="none",
                            std.err="jackknife")
  aberdeen.marginal <- rc(color[,,2], nd=2, weighting="marginal",
                          std.err="jackknife")
  caithness.unweighted
  caithness.marginal
  aberdeen.unweighted
  aberdeen.marginal
  \dontshow{
      # This is one of the very few articles reporting jackknife standard errors
      # Our results are very close for the unweighted solution;
      # but for marginal weighting their errors are probably not in the same unit,
      # as they are larger than the unweighted ones, while phi are much smaller.
      stopifnot(isTRUE(all.equal(c(round(se.rc(caithness.unweighted)$phi, d=3)),
                                 c(0.282, 0.052))))
      stopifnot(isTRUE(all.equal(c(round(se.rc(aberdeen.unweighted)$phi, d=3)),
                                 c(0.221, 0.040))))
  }


  ## Clogg & Shihadeh (1994), Tables 3.5a and b (p. 55-61)
  data(gss88)
  unweighted <- rc(gss88, weighting="none")$assoc
  marginal <- rc(gss88, weighting="marginal")$assoc
  uniform <- rc(gss88, weighting="uniform")$assoc
  # Row scores
  row.scores <- round(cbind(unweighted$row,
                            unweighted$row * sqrt(unweighted$phi[1,]),
                            marginal$row,
                            marginal$row * sqrt(marginal$phi[1,]),
                            uniform$row,
                            uniform$row * sqrt(uniform$phi[1,])), d=2)
  row.scores
  # Column scores
  col.scores <- round(cbind(unweighted$col,
                            unweighted$col * sqrt(unweighted$phi[1,1]),
                            marginal$col,
                            marginal$col * sqrt(marginal$phi[1,1]),
                            uniform$col,
                            uniform$col * sqrt(uniform$phi[1,1])), d=2)
  col.scores
  \dontshow{
      stopifnot(all.equal(c(row.scores),
                          c(0.57, 0.56, 0.15, -0.01, -0.23, 0.04, -0.42, -0.2, -0.19, -0.16, -0.1,
                            1.61, 1.6, 0.42, -0.04, -0.66, 0.11, -1.19, -0.56, -0.54, -0.46, -0.29,
                            1.78, 1.78, 0.37, -0.17, -0.91, 0, -1.53, -0.78, -0.76, -0.66, -0.47,
                            1.53, 1.52, 0.32, -0.14, -0.78, 0, -1.31, -0.67, -0.65, -0.57, -0.4,
                            1.88, 1.87, 0.48, -0.04,  -0.77, 0.12, -1.38, -0.65, -0.63, -0.53, -0.34, 
                            1.87, 1.87, 0.48, -0.04, -0.77, 0.12, -1.38, -0.65, -0.63, -0.53, -0.34)))
      stopifnot(all.equal(c(col.scores),
                          c(-0.56, -0.41, -0.08, 0.07, 0.41, 0.58,
                            -1.6, -1.16, -0.24, 0.19, 1.16, 1.65,
                            -1.89, -1.39, -0.33, 0.16, 1.29, 1.85,
                            -1.62, -1.19, -0.28, 0.14, 1.1, 1.59,
                            -1.38, -1, -0.21, 0.16, 1, 1.43,
                            -1.37, -1, -0.21, 0.16, 1, 1.42)))
  }

  ## Wong (2010), Table 2.7 (p. 48-49)
  data(gss8590)
  \dontshow{
      # Limit the workload for CRAN checks
      options(cl.cores=2)
  }
  # The table used in Wong (2010) is not perfectly consistent
  # with that of Wong (2001)
  tab <- margin.table(gss8590[,,c(2,4)], 1:2)
  tab[2,4] <- 49
  model <- rc(tab, nd=2, weighting="none", std.err="jackknife")
  model
  \dontshow{
      # A few scores differ from reported results by .001
      stopifnot(isTRUE(all.equal(round(c(model$assoc$phi), d=3),
                                 c(2.601, 1.522))))
      stopifnot(isTRUE(all.equal(round(c(model$assoc$row), d=3),
                                 c(0.743,  0.088, -0.200, -0.632,
                                   0.276, -0.061, -0.777,  0.562))))
      stopifnot(isTRUE(all.equal(round(c(model$assoc$col), d=3),
                                 c(0.765,  0.020, -0.322, -0.55, 0.088,
                                  -0.137, -0.549, -0.120, -0.01, 0.816))))
  }
  se.rc(model) # Jackknife standard errors are slightly different
               # from their asymptotic counterparts
  plot(model, conf.ellipses=TRUE)
}

\keyword{ models }
\keyword{ nonlinear }
