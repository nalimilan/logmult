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
  \item{diagonal}{the \code{diagonal} component from the object's \code{unidiff} component.}
  \item{deviance.resid}{the deviance residuals, see \code{\link{residuals.glm}.}}
  \item{layer}{a \code{data.frame} holding the layer coefficients estimates, standard
      errors, quasi-standard errors (see \code{\link[qvcalc]{qvcalc}}) and p-values.}
  \item{interaction}{a \code{data.frame} holding the two-way interaction coefficients
      estimates, standard errors and p-values.}
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
  \code{\link{unidiff}}, \code{\link{plot.unidiff}}
}

\keyword{ models }
\keyword{ nonlinear }
