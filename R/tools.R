printModelHeading <- function(x, digits=max(3, getOption("digits") - 4)) {
  cat("Call:\n", deparse(x$call), "\n", sep="", fill=TRUE)

  cat("Deviance Residuals:\n")
  if (x$df.residual > 5) {
      x$deviance.resid <- quantile(x$deviance.resid, na.rm = TRUE)
      names(x$deviance.resid) <- c("Min", "1Q", "Median", "3Q", "Max")
  }
  print(x$deviance.resid, digits=digits, na.print="", print.gap=2)
}

printModelStats <- function(x, digits=max(3, getOption("digits") - 4)) {
  cat("\nDeviance:              ", format(x$deviance, digits),
      "\nPearson chi-squared:   ",
      format(sum(na.omit(c(residuals(x, type="pearson")))^2), digits),
      "\nDissimilarity index:   ",
      format(sum(na.omit(c(abs(residuals(x, "response")))))/sum(na.omit(c(abs(fitted(x)))))/2*100, digits), "%",
      "\nResidual df:           ", x$df.residual,
      "\nBIC:                   ", x$deviance - log(sum(na.omit(c(x$data)))) * x$df.residual,
      "\nAIC:                   ", x$deviance - 2 * x$df.residual, "\n", sep="")
}

get.probs <- function(x) {
  if(inherits(x, "assoc.symm")) {
      # Weight with marginal frequencies, cf. Becker & Clogg (1994), p. 83-84, and Becker & Clogg (1989), p. 144.
      if(x$weighting == "marginal")
          rp <- cp <- prop.table(rowSums(x$row.weights) + rowSums(x$col.weights))
      else if(x$weighting == "uniform")
          rp <- cp <- rep(1/nr, nr)
      else
          rp <- cp <- rep(1, nr)
  }
  else {
      # Weight with marginal frequencies, cf. Becker & Clogg (1994), p. 83-84, and Becker & Clogg (1989), p. 144.
      if(x$weighting == "marginal") {
          rp <- prop.table(rowSums(x$row.weights))
          cp <- prop.table(rowSums(x$col.weights))
      }
      else if(x$weighting == "uniform") {
          rp <- rep(1/nr, nr)
          cp <- rep(1/nc, nc)
      }
      else {
          rp <- rep(1, nr)
          cp <- rep(1, nc)
      }
  }

  list(rp=rp, cp=cp)
}