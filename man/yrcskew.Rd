\name{yrcskew}
\alias{yrcskew}
\alias{print.yrcskew}
\title{Fitting Yamaguchi RC_SK Skew-Symmetric Association Model}
\description{
  \code{yrcskew} fits a skew-symmetric association model proposed in Yamaguchi (1990)
  to describe asymmetry of square tables. This model can be combined with symmetric
  association models like a quasi-symmetry (the default) or homogeneous RC(M) models.
}
\usage{
yrcskew(tab, nd.symm = NA, nd.skew = 1, diagonal = FALSE,
        weights = c("marginal", "uniform", "none"),
        std.err = c("none", "jackknife"), family = poisson,
        start = NULL, tolerance = 1e-12, iterMax = 15000, trace = TRUE, ...)
}
\arguments{
  \item{tab}{ a two-way table, or an object (such as a matrix) that can be coerced into a table;
    if present, dimensions above two will be collapsed as appropriate.}
  \item{nd.symm}{ the number of dimensions to include in the \emph{symmetric} RC(M) association. Cannot exceed
    \code{2 * min(nrow(tab) - 1, ncol(tab) - 1)}.
    If \code{NA}, a quasi-symmetric association is used instead of a RC(M) model.}
  \item{nd.skew}{ the number of dimensions to include in the \emph{skew-symmetric} RC(M) association.}
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
  The original presented by Yamaguchi (1990), called \dQuote{row-column-effect skew-symmetric
  association (logbilinear) model with full quasi-symmetry (QS+RC_SK)}, combines a skew-symmetric association
  with a quasi-symmetry baseline; it is the variant fitted by default by this function. If \code{nd.symm} is
  set to a positive integer value, though, variants using a RC(M) model to describe the \emph{symmetric association}
  are used, with our without diagonal-specific parameters (depending on the value of the \code{diagonal} argument);
  among them is the HM_RC+RC_SK variant, when \code{nd.symm} is 1.

  These models follow the equation:
  \deqn{ log F_{ij} = q_{ij} + \delta_{i<j} \nu_i (\nu_j - \nu_i) - \delta_{i>j} \nu_j (\nu_i - \nu_j)
                    = q_{ij} + (\delta_{i<j} - \delta_{i>j}) \nu_{min(i,j)} (\nu_{max(i,j)} - \nu_{min(i,j)})}
  where \eqn{F_i} is the frequency observed in the cell at the intersection of row i and column j of
  \code{tab}, and \eqn{q_{ij}} a quasi-symmetric or a RC(M) association. See reference for detailed information
  about the degrees of freedom and the identification constraints applied to the scores.

  Please note that contrary to other association models, this model is sensitive to reorderings of rows and columns.
  You have to take care of passing a table whose categories follow a hierarchical order with a substantive meaning.

  Another model presented in the paper, the homogeneous symmetric and skew-symmetric associations models (HM_(S+SK))
  is not currently supported.

  Actual model fitting is performed using \code{\link{gnm}}, which implements the Newton-Raphson algorithm.
  This function simply ensures correct start values are used, in addition to allowing for identification
  of scores even with several dimensions, computation of their jackknife standard errors, and plotting.
}
\value{
  A \code{yrcskew} object, which is a subclass of an \code{rc.homog} object (see\code{\link{rc}}) unless
    \code{nd.symm} is \code{NA}. In addition to this class, it contains a \code{yrcskew.assoc} component
    holding information about the \emph{skew-symmetric} association:
  \item{phi }{The intrisic association parameters, one per dimension.}
  \item{row }{Row scores, normalized so that their (weighted) mean is 0, their (weighted)
    sum is 1, and their (weighted) cross-dimensional correlation is null.}
  \item{col }{Column scores, normalized so that their (weighted) mean is 0, their (weighted)
    sum is 1, and their (weighted) cross-dimensional correlation is null.}
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
\section{Warning}{
  This family of model sometimes converges to a non-optimal solution, in which case the reported scores
  are wrong. To protect yourself from this problem, you are advised to run the models several times to find out which
  convergence point is the true one. Furthermore, when model converges slowly, restarting the fitting procedure may produce
  much better random starting values.
}
\references{
    Yamaguchi, K. (1990). Some Models for the Analysis of Asymmetric Association in Square Contingency Tables with
      Ordered Categories. \emph{Sociol. Methodology} 20, 181-212.
}
\author{
Milan Bouchet-Valat
}
\seealso{
  \code{\link{plot.yrcskew}}, \code{\link{gnm}}
}
\examples{

  ## Yamaguchi (1990), Table 5, p. 202, and Table 6B, p. 205
  data(ocg1973)

  # Simple homogeneous RC(1) model ("Null skew-symmetry")
  rc.model <- rc(ocg1973, diagonal=TRUE, homogeneous=TRUE, weights="none")
  # Reported phi is slightly different, coefficients agree
  rc.model
  \dontshow{
      stopifnot(isTRUE(all.equal(round(rc.model$assoc$row[,1,1], d=2),
                                 c(0.69, 0.3, -0.1, -0.33, -0.56),
                                 check.attributes=FALSE)))
  }

\dontrun{
  # Warning: model does not always converge, several attempts may be needed
  model <- yrcskew(ocg1973, nd.symm=1, nd.skew=1, diagonal=TRUE, weights="none")
  # We do not get the same results as the author, but the smaller deviance
  # indicates a better fit in our version (!)
  model
}
}
