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
      format(sum(abs(abs(residuals(x, "response"))))/sum(abs(fitted(x)))/2*100, digits), "%",
      "\nResidual df:           ", x$df.residual,
      "\nBIC:                   ", extractAIC(x, k=log(sum(x$data)))[2],
      "\nAIC:                   ", extractAIC(x)[2], "\n", sep="")
}
