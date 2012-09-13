## UNIDIFF model (Erikson & Goldthorpe, 1992), or uniform layer effect model (Xie, 1992)

unidiff <- function(tab, diagonal=c("included", "excluded", "only"),
                    constrain="auto", family=poisson,
                    tolerance=1e-6, iterMax=5000,
                    trace=FALSE, verbose=TRUE, ...) {
  diagonal <- match.arg(diagonal)

  tab <- as.table(tab)

  if(length(dim(tab)) < 3)
      stop("tab must have (at least) three dimensions")

  if(diagonal != "included" &&
     (nrow(tab) != ncol(tab) || !all(rownames(tab) == colnames(tab))))
     stop(sprintf("diagonal = %s is only supported for square tables with identical row and column categories",
                  diagonal))

  if(length(dim(tab)) > 3)
      tab <- margin.table(tab, 1:3)

  # When gnm evaluates the formulas, tab will have been converted to a data.frame,
  # with a fallback if both names are empty
  vars <- make.names(names(dimnames(tab)))
  if(length(vars) == 0)
      vars <- c("Var1", "Var2", "Var3")

  f <- sprintf("Freq ~ %s + %s + %s + %s:%s + %s:%s",
               vars[1], vars[2], vars[3], vars[1], vars[3], vars[2], vars[3])

  if(diagonal == "included") {
      f <- sprintf("%s + Mult(Exp(%s), %s:%s)", f, vars[3], vars[1], vars[2])

      if(identical(constrain, "auto"))
          constrain <- sprintf("(Mult\\(Exp\\(.\\), \\Q%s:%s\\E\\)\\Q.%s%s\\E)|(Mult\\(Exp\\(\\Q%s\\E.*\\.(\\Q%s%s\\E:)|(:\\Q%s%s\\E$))",
                               vars[1], vars[2], vars[3], dimnames(tab)[[3]][1],
                               vars[3], vars[1], rownames(tab)[1], vars[2], colnames(tab)[1])
  }
  if(diagonal == "excluded") {
      f <- sprintf("%s + %s:Diag(%s, %s) + Mult(Exp(%s), %s:%s)",
                   f, vars[3], vars[1], vars[2], vars[3], vars[1], vars[2])

      if(identical(constrain, "auto"))
          # Last pattern matches diagonal coefficients
          constrain <- sprintf("(Mult\\(Exp\\(.\\), \\Q%s:%s\\E\\)\\Q.%s%s\\E)|(Mult\\(Exp\\(\\Q%s\\E.*\\.(\\Q%s\\E(\\Q%s\\E):)|(:\\Q%s%s\\E$))|(Mult\\(Exp\\(\\Q%s\\E.*\\.(\\Q%s\\E)$)",
                               vars[1], vars[2], vars[3], dimnames(tab)[[3]][1],
                               vars[3], vars[1], paste(rownames(tab)[1:2], collapse="\\E|\\Q"), vars[2], colnames(tab)[1],
                               vars[3], paste(paste(vars[1], rownames(tab), ":", vars[2],
                                                    rownames(tab), sep=""), collapse="\\E|\\Q"))
  }
  else if(diagonal == "only") {
      f <- sprintf("%s + Mult(Exp(%s), Diag(%s, %s))", f, vars[3], vars[1], vars[2])

      if(identical(constrain, "auto"))
          constrain <- sprintf("(Mult\\(Exp\\(.\\), Diag\\(\\Q%s, %s\\E\\)\\)\\Q.%s%s\\E)|(Mult\\(Exp\\(\\Q%s\\E.*\\.(\\Q%s%s\\E:)|(:\\Q%s%s\\E$))",
                               vars[1], vars[2], vars[3], dimnames(tab)[[3]][1],
                               vars[3], vars[1], rownames(tab)[1], vars[2], colnames(tab)[1])
  }


  # FIXME: we should be able to eliminate 3:1, but this triggers a "numerically singular system" error
#  if(missing(eliminate))
#      eliminate <- eval(parse(text=sprintf("quote(%s:%s)", vars[1], vars[3])))

  # We need to handle ... manually, else they would not be found when modelFormula() evaluates the call
  args <- list(formula=as.formula(f), data=tab, constrain=constrain,
               family=family, tolerance=tolerance, iterMax=iterMax,
               trace=trace, verbose=verbose)
  dots <- as.list(substitute(list(...)))[-1]
  args <- c(args, dots)

  model <- do.call("gnm", args)

  if(is.null(model))
      return(NULL)

  model$unidiff <- list()
  
  model$unidiff$layer <- getContrasts(model, pickCoef(model, sprintf("Mult\\(Exp\\(\\.\\)", vars[3])),
                                      "first", check=FALSE)


  if(diagonal == "included") {
      mat <- tab[,,1]
      nr <- nrow(mat)
      nc <- ncol(mat)
      con <- matrix(0, length(coef(model)), length(mat))
      ind <- pickCoef(model, sprintf("Mult\\(Exp\\(\\Q%s\\E\\)", vars[3]))
      for(i in 1:nr) {
          for(j in 1:nc) {
              mat[] <- 0
              mat[i,] <- -1/nc
              mat[,j] <- -1/nr
              mat[i,j] <- 1 - 1/nc - 1/nr
              mat <- mat + 1/(nr * nc)
              con[ind, (j - 1) * nr + i] <- mat
          }
      }

      colnames(con) <- names(ind)
      model$unidiff$interaction <- gnm:::se(model, con)
  }
  else {
      # Quasi-variances cannot be computed for these coefficients, so hide the warning
      suppressMessages(model$unidiff$interaction <- getContrasts(model, pickCoef(model, sprintf("Mult\\(Exp\\(\\Q%s\\E\\)", vars[3])),
                                                                 ref="first", check=TRUE)$qvframe)
  }

  rownames(model$unidiff$layer$qvframe) <- dimnames(tab)[[3]]
  rownames(model$unidiff$interaction) <- gsub(sprintf("Mult\\(Exp\\(\\Q%s\\E\\), \\.\\)\\.(Diag\\(%s, %s\\))?",
                                                      vars[3], vars[1], vars[2]), "",
                                              rownames(model$unidiff$interaction))

  if(diagonal == "only") {
     # Diag() sorts levels alphabetically, which is not practical
     model$unidiff$interaction$qvframe <- model$unidiff$interaction$qvframe[c(1, 1 + order(rownames(tab))),]
     model$unidiff$interaction$covmat <- model$unidiff$interaction$covmat[c(1, 1 + order(rownames(tab))),]
  }

  model$unidiff$diagonal <- diagonal

  class(model) <- c("unidiff", class(model))

  model$call <- match.call()

  model
}

print.unidiff <- function(x, digits = max(3, getOption("digits") - 4), ...) {
  cat("Call:\n", deparse(x$call), "\n", sep = "", fill = TRUE)

  cat("Layer coefficients:\n")
  layer <- x$unidiff$layer$qvframe[,1]
  names(layer) <- paste(names(dimnames(tab))[3], rownames(x$unidiff$layer$qvframe), sep="")
  print.default(format(layer, digits=digits, ...), quote=FALSE)

  if(x$unidiff$diagonal != "only") {
      cat("\nFull two-way interaction coefficients:\n")
      interaction <- x$data[,,1]
      interaction[] <- x$unidiff$interaction[,1]
  }
  else {
      cat("\nDiagonal interaction coefficients:\n")
      interaction <- x$unidiff$interaction$qvframe[,1]
      names(interaction) <- rownames(x$unidiff$interaction)
  }

  print.default(format(interaction, digits=digits, ...), quote=FALSE)

  cat("\nDeviance:", round(x$deviance, digits=3), "\n")
  cat("Degrees of freedom:", x$df.residual, "\n")
  cat("BIC:", round(extractAIC(x, k=log(sum(x$data)))[2]), " - ",
      "AIC:", round(extractAIC(x)[2]), "\n")

  invisible(x)
}

summary.unidiff <- function(x, ...) {
  layer <- x$unidiff$layer$qvframe[,-4]
  interaction <- x$unidiff$interaction

  layer <- cbind(layer, 2 * pnorm(-abs(layer[,1]/layer[,2])))
  interaction <- cbind(interaction, 2 * pnorm(-abs(interaction[,1]/interaction[,2])))

  colnames(layer) <- c("Estimate", "Std. Error", "Quasi SE", "Pr(>|z|)")
  colnames(interaction) <- c("Estimate", "Std. Error", "Pr(>|z|)")

  res <- list(call=x$call, diagonal=x$unidiff$diagonal,
              deviance.resid=residuals(x, type="deviance"),
              layer=layer, interaction=interaction,
              deviance=x$deviance, df.residual=x$df.residual,
              bic=extractAIC(x, k=log(sum(x$data)))[2],
              aic=extractAIC(x)[2])

  class(res) <- "summary.unidiff"

  res
}

print.summary.unidiff <- function(x, digits = max(3, getOption("digits") - 4), ...) {
  cat("Call:\n", deparse(x$call), "\n", sep = "", fill = TRUE)

  cat("Deviance Residuals:\n")
  if (x$df.residual > 5) {
      x$deviance.resid <- quantile(x$deviance.resid, na.rm = TRUE)
      names(x$deviance.resid) <- c("Min", "1Q", "Median", "3Q", "Max")
  }
  print.default(x$deviance.resid, digits = digits, na = "", print.gap = 2)

  cat("\nLayer coefficients:\n")
  printCoefmat(x$layer, digits, signif.legend=FALSE, ...)

  if(x$diagonal != "only")
      cat("\nFull two-way interaction coefficients:\n")
  else
      cat("\nDiagonal interaction coefficients:\n")

  printCoefmat(x$interaction, digits, has.Pvalue=TRUE, ...)


  cat("\nDeviance:", round(x$deviance, digits=3), "\n")
  cat("Degrees of freedom:", x$df.residual, "\n")
  cat("BIC:", round(x$bic), " - ",
      "AIC:", round(x$aic), "\n")
}

plot.unidiff <- function(x, exponentiate=TRUE, se.type=c("quasi.se", "se"), conf.int=.95,
                         numeric.auto=TRUE, type="o",
                         xlab=names(dimnames(x$data))[3], ylab="Layer coefficient", ...) {
  if(!inherits(x, "unidiff"))
      stop("x must be a unidiff object")

  se.type <- match.arg(se.type)

  qv <- x$unidiff$layer$qvframe

  w <- qnorm((1 - conf.int)/2, lower.tail=FALSE)

  coefs <- qv$estimate

  if(se.type == "quasi.se") {
      tops <- qv$estimate + w * qv$quasiSE
      tails <- qv$estimate - w * qv$quasiSE
   }
   else {
      tops <- qv$estimate + w * qv$SE
      tails <- qv$estimate - w * qvSE
   }

  if(exponentiate) {
      coefs <- exp(coefs)
      tops <- exp(tops)
      tails <- exp(tails)
  }

  range <- max(tops) - min(tails)
  ylim <- c(min(tails) - range/20, max(tops) + range/20)

  # plot() converts x coordinates to numeric if possible, but segments
  # needs a real coordinate, so convert directly
  if(numeric.auto && !any(is.na(as.numeric(rownames(qv)))))
      at <- as.numeric(rownames(qv))
  else
      at <- factor(rownames(qv), levels=rownames(qv))

  plot.default(at, coefs, type=type, xaxt="n",
       ylim=ylim, xlab=xlab, ylab=ylab, ...)
  axis(1, at, labels=at)
  segments(as.numeric(at), tails, as.numeric(at), tops)

  invisible(x)
}
