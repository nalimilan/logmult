\name{plot.assoc}
\alias{plot.assoc}
\alias{plot.rc}
\alias{plot.rc.homog}
\alias{plot.hmskew}
\alias{plot.yrcskew}
\alias{plot.rcL}
\alias{plot.rcL.homog}
\title{Plotting Scores from Association Models}
\description{Graphical display of scores from association models in two dimensions}
\usage{
\method{plot}{rc}(x, dim = c(1, 2), what = c("both", "rows", "columns"),
     mass = TRUE, luminosity = length(x$assoc$diagonal > 0),
     conf.ellipses = FALSE, coords = c("cartesian", "polar"),
     rev.axes = c(FALSE, FALSE), cex = par("cex"),
     col = c("blue", "red"), groups = NULL,
     xlim, ylim, xlab, ylab, ...)

\method{plot}{rc.homog}(x, dim = c(1, 2),
     mass = TRUE, luminosity = length(x$assoc$diagonal > 0),
     conf.ellipses = FALSE, coords = c("cartesian", "polar"),
     rev.axes = c(FALSE, FALSE), cex = par("cex"),
     col = "dark grey", groups = NULL,
     xlim, ylim, xlab, ylab, ...)

\method{plot}{hmskew}(x, dim = c(1, 2), what = c("skew-symmetric", "symmetric"),
     mass = TRUE, luminosity = length(x$assoc.hmskew$diagonal > 0),
     arrow = 45, conf.ellipses = FALSE, coords = c("polar", "cartesian"),
     rev.axes = c(FALSE, FALSE), cex = par("cex"),
     col = "dark grey", groups = NULL,
     xlim, ylim, xlab, ylab, ...)

\method{plot}{yrcskew}(x, dim = c(1, 2), what = c("skew-symmetric", "symmetric"),
     mass = TRUE, luminosity = length(x$assoc.yrcskew$diagonal > 0),
     arrow = 45, conf.ellipses = FALSE, coords = c("polar", "cartesian"),
     rev.axes = c(FALSE, FALSE), cex = par("cex"),
     col = "dark grey", groups = NULL,
     xlim, ylim, xlab, ylab, ...)

\method{plot}{rcL}(x, dim = c(1, 2), layer = 1,
     what = c("both", "rows", "columns"),
     mass = TRUE, luminosity = length(x$assoc$diagonal > 0),
     conf.ellipses = FALSE, coords = c("cartesian", "polar"),
     rev.axes = c(FALSE, FALSE), cex = par("cex"),
     col = c("blue", "red"), groups = NULL,
     xlim, ylim, xlab, ylab, ...)

\method{plot}{rcL.homog}(x, dim = c(1, 2), layer = 1,
     mass = TRUE, luminosity = length(x$assoc$diagonal > 0),
     conf.ellipses = FALSE, coords = c("cartesian", "polar"),
     rev.axes = c(FALSE, FALSE), cex = par("cex"),
     col = c("blue", "red"), groups = NULL,
     xlim, ylim, xlab, ylab, ...)

\method{plot}{assoc}(x, dim = c(1, 2), layer = 1,
     what = c("both", "rows", "columns"),
     mass = TRUE, luminosity = length(x$diagonal > 0),
     arrow = NULL, conf.ellipses = FALSE, coords = c("cartesian", "polar"),
     rev.axes = c(FALSE, FALSE), cex = par("cex"),
     col = c("blue", "red"), groups = NULL,
     xlim, ylim, xlab, ylab, ...)
}
\arguments{
  \item{x}{an association model returned by \code{\link{rc}}, \code{\link{hmskew}} or
           \code{\link{yrcskew}}}
  \item{dim}{numerical vector of length 2 indicating the dimensions to plot on horizontal
             and vertical axes respectively; default is first dimension horizontal and second
             dimension vertical.}
  \item{layer}{integer indicating which layer should be represented.}
  \item{what}{for \code{rc} and \code{assoc} objects, whether points corresponding to rows, columns
              or both should be plotted; for \code{hmskew} and \code{yrcskew} objects, what association
              should be plotted.}
  \item{mass}{whether the size of the point symbols should reflect on the mass of the categories;
              this only makes sense when marginal weights were used when fitting the model. The precise
              formula is that the \code{pch} of a symbol is equal to the \code{pch} argument, times the
              weight of the category divided by average weight.}
  \item{luminosity}{whether the luminosity of the symbols should vary according to the diagonal-specific
                    parameters (if \code{diagonal = TRUE} was passed when fitting the model); if \code{TRUE},
                    hue and saturation of \code{col} are taken as a base, and value varies from 0 to 0.8
                    as a linear function of the diagonal parameter values.}
  \item{arrow}{Numeric value indicating the angle at which the polar coordinates system arrow should be plotted;
               use NULL for no arrow.}
  \item{conf.ellipses}{whether to plot 95% confidence ellipses around the points; this is only possible if
                       standard errors were computed when fitting the model.}
  \item{coords}{whether to use a Cartesian or a polar coordinate system; the former makes sense when axes
                offer an interpretation (like in RC(M) models), while the latter are more appropriate when
                only the angle and distance to origin are of interest (like in \code{hmskew} models).}
  \item{rev.axes}{a numeric of length 1 or 2 indicating whether the sign of scores on the axes should be changed;
                  as this sign is arbitrary in RC(M) models, changing it has no incidence on the results and may
                  be more intuitive or consistent with other presentations.}
  \item{cex}{a numeric vector indicating the size of the point symbols, recycled as necessary; the most common
             choice is probably to pass only one value and use \code{mass} to make the size vary.}
  \item{col}{a vector indicating the color of the point symbols, recycled as necessary; as a special case, a
             vector of length 2 can be passed, to indicate the color of row and column points, respectively.
             See also \code{luminosity}.}
  \item{groups}{a vector indicating what point symbols should be used, recycled as necessary; groups will
                use \code{pch} values 21, 24, 22, 23 and 25, in this order, cycling if needed.}
  \item{xlab}{a title for the x axis: see \code{\link{title}}. For RC(M) axes and Cartesian coordinates,
              the default is \dQuote{Dimension N (phi)}; it is empty in other cases where axes have no meaning.}
  \item{ylab}{a title for the y axis: see \code{\link{title}}. For RC(M) axes and Cartesian coordinates,
              the default is \dQuote{Dimension N (phi)}; it is empty in other cases where axes have no meaning.}
  \item{xlim}{numeric vector of length 2, giving the x coordinates range.}
  \item{ylim}{numeric vector of length 2, giving the y coordinates range.}
  \item{...}{Further arguments passed to \code{\link{plot}}.}
}
\details{
  All the functions documented here represent in a two-dimensional space the scores obtained from a
  log-multiplicative association model. They produce \emph{symmetric biplots} in which the coordinates of points
  on both axes are the product of normalized scores and of the square root of the intrinsic association
  coefficient corresponding to each dimension: thus, row and column points share the same \dQuote{unit}
  on all axes (Goodman, 1991, Appendix 2; Wong, 2010, eq. 2.38; Clogg & Shihadeh, 1994, p. 91-92).

  Various convenience options are provided, with reasonable defaults for each model family. In particular,
  you may find it necessary to adapt the \code{cex}, \code{mass}, \code{luminosity} and \code{groups} arguments
  depening on the number of categories to be plotted and to their respective weights. When plotting a RC(2) model,
  a polar coordinate system can be of substantive interest, allowing to interpret at a glance the distance to
  origin as the general strength of the association for a category on both axes (a property that is lost for
  higher-dimensional models).

  The \code{plot.assoc} function is called internally by all others, and may be leveraged for advanced use cases,
  like plotting custom models that do not correspond stricly to the supported types.
}

\references{
  For RC(M) models:

  Goodman, L.A. (1991). Measures, Models, and Graphical Displays in the Analysis of
    Cross-Classified Data. \emph{J. of the Am. Stat. Association} 86(416), 1085-1111.

  Clogg, C.C., and Shihadeh, E.S. (1994). Statistical Models for Ordinal Variables. Sage: Advanced
      Quantitative Techniques in the Social Sciences (4).

  Wong, R.S-K. (2010). Association models. Sage: Quantitative Applications in the Social Sciences (164).


  For van der Heijden & Mooijaart models:

  van der Heijden, P.G.M., and Mooijaart, A. (1995). Some new log bilinear models for the analysis
    of asymmetry in a square contingency table. \emph{Sociol. Methods and Research} 24, 7-29.

}
\seealso{\code{\link{rc}}, \code{\link{hmskew}}, \code{\link{yrcskew}}}
\examples{
  ## Wong (2010), Figures 2.2 and 2.1 (p. 50-51)
  data(gss8590)
  model <- rc(margin.table(gss8590[,,c(2,4)], 1:2),
              nd=2, weights="none", std.err="jackknife")
  plot(model, what="row", rev.axes=c(TRUE, FALSE), conf.ellipses=TRUE)
  plot(model, what="col", rev.axes=c(TRUE, FALSE), conf.ellipses=TRUE)

  ## Wong (2010), Figures 4.1 and 4.2 (p. 108-109)
  data(gender)
  model <- rcL(gender, nd=2, weights="none")
  par(mfrow=c(2, 2))
  for(i in 1:4)
      plot(model, layer=i, what="rows", rev.axes=c(TRUE, FALSE),
           main=rownames(model$assoc$phi)[i], xlim=c(-1.2, 1.2), ylim=c(-1.2, 1.2))
  par(mfrow=c(2, 2))
  for(i in 1:4)
      plot(model, layer=i, what="col", rev.axes=c(TRUE, FALSE),
           main=rownames(model$assoc$phi)[i], xlim=c(-1.4, 1.4), ylim=c(-1, 1))


  ## van der Heijden & Mooijaart (1995), Figure 1c (p. 23)
  data(ocg1973)
  # 5:1 is here to take "Farmers" as reference category (angle 0)
  model <- hmskew(ocg1973[5:1, 5:1], weights="uniform")
  # Reproduce the plot from the original article
  plot(model, coords="cartesian")
  # Use a polar coordinates system, which makes more sense in this setting
  plot(model)
}
