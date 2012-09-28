\name{unidiff}
\alias{unidiff}
\title{
  Fitting Log-Multiplicative Uniform Difference/Layer Effect Model
}
\description{
  \code{unidiff} fits the log-multiplicative uniform difference model (UNIDIFF,
  see Erikson & Goldthorpe, 1992), also called the log-multiplicative layer effect model (Xie, 1992).
}
\usage{
unidiff(tab, diagonal = c("included", "excluded", "only"),
        constrain = "auto", family = poisson,
        tolerance = 1e-8, iterMax = 5000,
        trace = FALSE, verbose = TRUE, ...)
}
\arguments{
  \item{tab}{a three-way table, or an object (such as a matrix) that can be coerced into a table;
      if present, dimensions above three will be collapsed as appropriate.}
  \item{diagonal}{\code{included} fits the standard model with full two-way interaction;
      \code{excluded} adds to this model diagonal-specific parameters for each years, effectively
      removing the influence of diagonal cells on the layer coefficients; \code{only} fits a model
      without the full two-way interaction, where only diagonal parameters are affected by the layer
      effect (see \dQuote{Details} below).}
  \item{constrain}{(non-eliminated) coefficients to constrain, specified by a regular expression,
      a numeric vector of indices, a logical vector, a character vector of names, or "[?]" to
      select from a Tk dialog. The default constrains to 0 the first layer parameter and
      interaction coefficients for the first row and column of the table.}
  \item{family}{a specification of the error distribution and link function to be used in the model.
      This can be a character string naming a family function; a family function, or the result of
      a call to a family function. See \code{\link{family}} details of family functions.}
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
  The equation of the fitted model is:
  \deqn{ log F_{ijk} = \lambda + \lambda^I_i + \lambda^J_j + \lambda^K_k
                      + \lambda^IK_ik + \lambda^JK_jk
                      + \phi_k \psi^IJ_ij }
  where \eqn{F_ijk} is the expected frequency for the cell at the intersection of row i, column j and
  layer k of \code{tab}. When \code{diagonal = "excluded"}, \eqn{\lambda^IJK_ijk} parameters are added
  but set to 0 when \eqn{i != j} (off-diagonal). When \code{diagonal = "only"}, \eqn{\psi^IJ_ij} is set
  to 0 when \eqn{i != j}.

  Layer coefficients \eqn{\phi_k} are internally exponentiated in the gnm formula, which means the reported
  values are in log scale, with reference 0 for the first year. Interaction coefficients use the
  \dQuote{sum} contrast, also known as \dQuote{effect} coding, except when \code{diagonal} is different from
  \code{included}, in which case \dQuote{treatment} constrast (a.k.a \dQuote{reference} or \dQuote{dummy}
  coding) is used.

  Actual model fitting is performed using \code{\link{gnm}}, which implements the Newton-Raphson algorithm.
  This function simply allows for direct identification of the log-multiplicative parameters by setting the
  appropriate constraints, and improves performance by eliminating less interesting coefficients.
}
\value{
  A \code{unidiff} object, with all the components of a \code{\link{gnm}} object, plus an
      \code{unidiff} component holding the most relevant information:
  \item{diagonal}{the value of the \code{diagonal} argument above.}
  \item{layer}{a \code{\link{qvcalc}} object holding the (log) layer coefficients, their standard
      errors and quasi-standard errors.}
  \item{interaction}{a data frame object holding the two-way interaction coefficients, and their
      standard errors.}
}
\references{
  Erikson, R., and Goldthorpe, J.H. (1992). The Constant Flux: A Study of Class Mobility in Industrial Societies.
      Oxford: Clarendon Press. Ch. 3.

  Xie, Yu (1992). The Log-Multiplicative Layer Effect Model for Comparing Mobility Tables. \emph{Am. Sociol. Rev.}
      57(3):380-395.

  Yaish, M. (1998). Opportunities, Little Change. Class Mobility in Israeli Society, 1974-1991. Ph.D. thesis,
      Nuffield College, University of Oxford.

  Yaish, M. (2004). Class Mobility Trends in Israeli Society, 1974-1991. Lewiston: Edwin Mellen Press.
}
\author{
  Milan Bouchet-Valat
}
\seealso{
  \code{\link{plot.unidiff}}, \code{\link{summary.unidiff}}
}
\examples{
  ## Yaish (1998, 2004)
  data(yaish)

  # Last layer omitted because of low frequencies
  yaish <- yaish[,,-7]

  # Layer (education) must be the third dimension
  yaish <- aperm(yaish, 3:1)

  model <- unidiff(yaish)

  model
  summary(model)
  plot(model)

  \dontshow{
      # The last layer coefficient can sometimes be rounded to -2.249 without this, which breaks the test
      model <- unidiff(yaish, tolerance=1e-8)

      stopifnot(all.equal(round(model$deviance, 3), 200.326))
      stopifnot(all.equal(model$df.residual, 116))
      stopifnot(all.equal(round(model$unidiff$layer$qvframe$estimate, 3),
                          c(0, -0.225, -0.743, -1.039, -2.250)))
      # Interaction coefficients obtained using LEM
      stopifnot(all.equal(round(model$unidiff$interaction[,1], 4),
                          c(0.7084, -0.5925, -0.1608,  0.9028, -0.5433, -0.3146,  1.1473,
                            0.2729, -0.3248, -1.4509,  0.3526,  0.0029, -0.1734,  0.1755,
                            0.5605, -0.2185, -0.1669, -0.1772, -0.4809, -0.9644,  0.0862,
                            1.9604, -0.6909,  0.0896, -0.9153, -0.0035,  0.2453, -0.414,
                            0.7927,  0.2948, -0.3342,  0.3295,  0.0612, -0.5712,  0.2573,
                            0.2574,  0.0482,  0.7824, -0.4677, -0.2085, -0.0015, -0.153)))

  }
}

\keyword{ models }
\keyword{ nonlinear }
