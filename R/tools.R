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

get.probs.asymm <- function(weighting, weights, indices) {
  # Weight with marginal frequencies, cf. Becker & Clogg (1994), p. 83-84, and Becker & Clogg (1989), p. 144.
  if(weighting == "marginal")
      prop.table(rowSums(weights[indices,, drop=FALSE]))
  else if(weighting == "uniform")
      rep(1/nrow(weights), nrow(weights))
  else
      rep(1, nrow(weights))
}

get.probs.symm <- function(weighting, row.weights, col.weights, indices) {
  # Weight with marginal frequencies, cf. Becker & Clogg (1994), p. 83-84, and Becker & Clogg (1989), p. 144.
  if(weighting == "marginal")
      prop.table(rowSums(row.weights[indices,, drop=FALSE]) + rowSums(col.weights[indices,, drop=FALSE]))
  else if(weighting == "uniform")
      rep(1/nrow(row.weights), nrow(row.weights))
  else
      rep(1, nrow(row.weights))
}

get.probs <- function(ass) {
  if(inherits(ass, "assoc.symm")) {
      if(length(ass$row.sup) > 0 & length(ass$col.sup) > 0) {
          p <- c(get.probs.symm(ass$weighting, ass$row.weights, ass$col.weights, setdiff(seq(nrow(ass$row)), ass$row.sup)),
                 get.probs.symm(ass$weighting, ass$row.weights, ass$col.weights, ass$row.sup))
          list(rp=p, cp=p)
      }
      else {
          p <- get.probs.symm(ass$weighting, ass$row.weights, ass$col.weights, seq(nrow(ass$row)))
          list(rp=p, cp=p)
      }

  }
  else {
      if(length(ass$row.sup) > 0)
          rp <- c(get.probs.asymm(ass$weighting, ass$row.weights, setdiff(seq(nrow(ass$row)), ass$row.sup)),
                  get.probs.asymm(ass$weighting, ass$row.weights, ass$row.sup))
      else
          rp <- get.probs.asymm(ass$weighting, ass$col.weights)

      if(length(ass$col.sup) > 0)
          cp <- c(get.probs.asymm(ass$weighting, ass$col.weights[setdiff(seq(ncol(ass$col)), ass$col.sup),]),
                  get.probs.asymm(ass$weighting, ass$col.weights[ass$col.sup,]))
      else
          cp <- get.probs.asymm(ass$weighting, ass$col.weights)

      list(rp=rp, cp=cp)
  }
}

goodman.phi <- function(tab, weighting=c("marginal", "uniform", "none"), row.weights=NULL, col.weights=NULL) {
  weighting <- match.arg(weighting)
  stopifnot(is.matrix(tab))

  if(any(is.na(tab)))
      stop("NA cells are currently not supported.")

  if(any(tab == 0))
      stop("Index cannot be computed in the presence of cells with 0 counts; replace them with 1/2 if appropriate.")

  p <- prop.table(tab)

  if(weighting == "marginal") {
      rp <- margin.table(p, 1)
      cp <- margin.table(p, 2)
  }
  else if(weighting == "uniform") {
      rp <- rep(1/nrow(tab), nrow(tab))
      cp <- rep(1/ncol(tab), ncol(tab))
  }
  else {
      rp <- rep(1, nrow(tab))
      cp <- rep(1, ncol(tab))
  }

  if(!is.null(row.weights))
      rp <- row.weights

  if(!is.null(col.weights))
      cp <- col.weights

  rp1 <- prop.table(rp)
  cp1 <- prop.table(cp)
  logp <- log(p)

  lambda <- sweep(logp, 1, rowSums(sweep(logp, 2, cp1, "*")), "-")
  lambda <- sweep(lambda, 2, colSums(sweep(logp, 1, rp1, "*")), "-")
  lambda <- lambda + sum(logp * rp1 %o% cp1)

  sqrt(sum(lambda^2 * rp %o% cp))
}