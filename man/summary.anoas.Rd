\name{summary.anoas}
\alias{summary.anoas}
\alias{print.anoas}
\alias{print.summary.anoas}
\title{Summary and Print Methods for ANOAS objects}
\description{
  These functions print the summary of a list of models fitted using the \code{\link{anoas}} function.
}
\usage{
\method{summary}{anoas}(object, ...)

\method{print}{anoas}(x, ...)

\method{print}{summary.anoas}(x, digits = 1, nsmall = 2, scientific = FALSE, ...)
}
\arguments{
  \item{object}{ an \code{anoas} object.}
  \item{x}{ an \code{anoas} object.}
  \item{digits}{ See ?format.}
  \item{nsmall}{ See ?format.}
  \item{scientific}{ See ?format.}
  \item{\dots}{ more arguments to be passed to further methods (ignored by \code{summary.anoas}).}
}
\details{
  Contrary to most analyses of association in the literature, this function currently does not fit uniform association model
  (\dQuote{U}), nor separate models with only row and column association (\dQuote{R} and \dQuote{C} models), nor log-linear
  row and column association models.

  Currently, no significance test is performed on the models. Please note that it is not correct to test the one-dimension
  association model against the independence model.
}
\value{
  A \code{data.frame} with the following columns:
  \itemize{
      \item{Res. Df}{ the residual number of degrees of freedom of the model.}
      \item{Res. Dev}{ the residual deviance of the model (likelihood ratio Chi-squared statistic, or L-squared).}
      \item{Dev. Indep. (\%)}{ the ratio of the residual deviance of the model over that of the independence model,
          times 100. This measures the share of departure from independence that cannot be explained using the
          number of dimensions of the model.}
      \item{Dissim. (\%)}{ the dissimilarity index of the model's fitted values with regard to the observed data.}
      \item{BIC}{ the Bayesian Information Criterion for the model.}
      \item{AIC}{ Akaike's An Information Criterion for the model.}
      \item{Deviance}{ the reduction in deviance of the model compared to the previous one}
      \item{Df}{ the reduction in the number of degrees of freedom of the model compared to the previous one.}
  }
}
\author{
  Milan Bouchet-Valat
}
\seealso{
  \code{\link{anoas}}, \code{\link{anoasL}}
}

\keyword{ models }
\keyword{ nonlinear }
