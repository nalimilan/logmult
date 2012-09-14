printModelStats <- function(x, digits=max(3, getOption("digits") - 4)) {
  cat("\nDeviance:             ", format(x$deviance, digits),
      "\nPearson chi-squared:  ",
      format(sum(na.omit(c(residuals(x, type="pearson")))^2), digits),
      "\nDissimilarity index:  ",
      format(sum(abs(abs(residuals(x, "response"))))/sum(abs(fitted(x)))/2*100, digits), "%",
      "\nResidual df:          ", x$df.residual,
      "\nBIC:                  ", extractAIC(x, k=log(sum(x$data)))[2],
      "\nAIC:                  ", extractAIC(x)[2], "\n")
}
