\name{summary.unidiff}
\alias{summary.unidiff}
\alias{print.summary.unidiff}
\title{
  Summarize UNIDIFF Model Fits
}
\description{
  \code{summary} method for objects of class \code{unidiff}.
}
\usage{
\method{summary}{unidiff}(object, ...)

\method{print}{summary.unidiff}(x, digits = max(3, getOption("digits") - 4), ...)
}
\arguments{
  \item{object}{an object resulting from a call to \code{\link{unidiff}}}
  \item{x}{an object of class \code{summary.gnm}.}
  \item{digits}{the number of siginificant digits to use when printing.}
  \item{\dots}{further arguments passed to \code{\link{printCoefmat}} by
      \code{print.summary.unidiff}, and currently ignored by \code{summary.unidiff}.}
}
\details{
  \code{print.summary.unidiff} prints the original call to \code{unidiff}; a
  summary of the deviance residuals from the model fit; the coefficients of
  interest of the model; the residual deviance; the residual degrees of freedom; 
  the Schwartz's Bayesian Information Criterion value; the Akaike's An Information
  Criterion value.


  Layer and two-way interaction coefficients are printed with their standard errors,
  quasi-standard errors (see \code{\link[qvcalc]{qvcalc}}), p-values (based on standard
  errors) and significance stars. Constrained coefficients have a value of 0 (by default),
  and 0 standard errors, but still have quasi-standard errors.
}
\value{
  An object of class \code{summary.unidiff}, with the following components:
  \item{call}{the \code{call} component from object.}
  \item{deviance.resid}{the deviance residuals, see \code{\link{residuals.glm}.}}
  \item{layer}{a \code{data.frame} holding the layer coefficients estimates, standard
      errors, quasi-standard errors (see \code{\link[qvcalc]{qvcalc}}) and p-values.}
  \item{phi.layer}{a \code{data.frame} holding the layer coefficients estimates,
      standard errors, and quasi-standard errors (see \code{\link[qvcalc]{qvcalc}})
      multiplied by the intrinsic association coefficient (see \code{\link{iac}})
      for the first layer; p-values are the same as those for the \dQuote{layer} component.}
  \item{interaction}{a \code{data.frame} holding the two-way interaction coefficients
      estimates, standard errors and p-values.}
  \item{deviance}{the \code{deviance} component from object.}
  \item{diagonal}{the \code{diagonal} component from the object's \code{unidiff} component.}
  \item{weighting}{the \code{weighting} component from the object's \code{unidiff} component.}
  \item{chisq}{the Pearson Chi-squared statistic for the model fit.}
  \item{dissim}{the dissimilarity index for the model fit.}
  \item{df.residual}{the \code{df.residual} component from object.}
  \item{bic}{the value of the BIC for the model fit (contrary to the value reported by
      \code{\link{AIC}} and \code{\link{extractAIC}}, the reference is 0 for the
      saturated model).}
  \item{aic}{the value of the AIC for the model fit (contrary to the value reported by
      \code{\link{AIC}} and \code{\link{extractAIC}}, the reference is 0 for the
      saturated model).}
  \item{family}{the \code{family} component from object.}
  \item{dispersion}{the estimated dispersion}
  \item{df}{a 3-vector of the rank of the model; the number of residual
            degrees of freedom; and number of unconstrained coefficients.}
}
\author{
  Milan Bouchet-Valat
}
\seealso{
  \code{\link{unidiff}}, \code{\link{plot.unidiff}}
}

\keyword{ models }
\keyword{ nonlinear }
