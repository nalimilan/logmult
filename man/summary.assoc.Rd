\name{summary.assoc}
\alias{summary.assocmod}
\alias{print.summary.assocmod}
\title{
  Summarize Association Model Fits
}
\description{
  \code{summary} method for objects of class \code{assocmod}, including \code{\link{rc}},
  \code{\link{rcL}}, \code{\link{rcL.trans}}, \code{\link{hmskew}}, \code{\link{hmskewL}}
  and \code{\link{yrcskew}} models.
}
\usage{
\method{summary}{assocmod}(object, weighting, ...)

\method{print}{summary.assocmod}(x, digits = max(3, getOption("digits") - 4), ...)
}
\arguments{
  \item{object}{an association model of class \code{assocmod}.}
  \item{x}{an object of class \code{summary.gnm}.}
  \item{weighting}{what weights should be used when normalizing the scores.}
  \item{digits}{the number of siginificant digits to use when printing.}
  \item{\dots}{further arguments passed to \code{\link{printCoefmat}} by
      \code{print.summary.assocmod}, and currently ignored by \code{summary.assocmod}.}
}
\details{
  \code{print.summary.assocmod} prints the original call to \code{assoc}; a
  summary of the deviance residuals from the model fit; the coefficients of
  interest of the model; the residual deviance; the residual degrees of freedom; 
  the Schwartz's Bayesian Information Criterion value; the Akaike's An Information
  Criterion value.

  Association coefficients are printed with their standard errors, p-values and significance
  stars. The \dQuote{Estimate} columns contains \emph{normalized scores}, i.e. their (weighted)
  sum is 0, their (weighted) sum of squares is 1, and their (weighted) cross-dimensional
  correlation is null. For models with only one layer (\code{\link{rc}}, \code{\link{hmskew}},
  \code{\link{yrcskew}}), \emph{adjusted} scores are printed in the \dQuote{Adjusted} column:
  these correspond to normalized scores times the square root of the corresponding intrinsic
  association parameter (phi).

  p-values are computed under the assumption that estimators of coefficients are normally distributed,
  even if jackknife of bootstrap are used. This hypothesis can be checked by calling
  \code{link[boot]{plot.boot}} on the \code{boot.results} component of the \code{assoc} object of the models
  (not supported for jackknife). (With a large number of replicates (> 1000), bootstrap confidence intervals
  can be computed using \code{link[boot]{boot.ci}} on the same object.) Checking the values of coefficients in all
  replicates using \code{plot.boot} is a good idea in all cases to ensure there are no dramatic outliers, which
  could indicate some replicates converged to a non optimal solution: in that case, standard errors will be
  artificially large; if it happens, reduce the value of \code{iterMax} when fitting models.

  Note that setting the \code{weighting} argument to a value different from that used at the
  time of the fit discards the computed standard errors, if any.
}
\value{
  An object of class \code{summary.assoc}, with the following components:
  \item{call}{the \code{call} component from object.}
  \item{diagonal}{the \code{diagonal} component from the object's \code{assoc} component.}
  \item{deviance.resid}{the deviance residuals, see \code{\link{residuals.glm}.}}
  \item{coefficients}{a matrix holding the association coefficients estimates, standard
      errors and p-values.}
  \item{diagonal}{a matrix holding the diagonal coefficients, if any.}
  \item{weighting}{the weigthing method used when normalizing the scores.}
  \item{deviance}{the \code{deviance} component from object.}
  \item{chisq}{the Pearson Chi-squared statistic for the model fit.}
  \item{dissim}{the dissimilarity index for the model fit.}
  \item{df.residual}{the \code{df.residual} component from object.}
  \item{bic}{the value of the BIC for the model fit (see \code{\link{extractAIC})}.}
  \item{aic}{the value of the AIC for the model fit (see \code{\link{AIC}},
      \code{\link{extractAIC}}).}
}
\author{
  Milan Bouchet-Valat
}
\seealso{
  \code{\link{assoc}}, \code{\link{plot.assoc}}, \code{\link{rc}}, \code{\link{rcL}},
  \code{\link{rcL.trans}}, \code{\link{hmskew}}, \code{\link{hmskewL}}, \code{\link{yrcskew}}
}

\keyword{ models }
\keyword{ nonlinear }
