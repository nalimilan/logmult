## RC-L with constant scores over time, but different scores for rows and columns
## (called homogeneous RC(M) by Wong(2010))

rcL <- function(tab, nd=1, layer.homogeneous=c("both", "rows", "columns", "none"),
               homogeneous=FALSE, diagonal=FALSE,
               weights=c("marginal", "uniform", "none"), std.err=c("none", "jackknife"),
               family=poisson, start=NULL, tolerance=1e-12, iterMax=5000, trace=TRUE, ...) {
  layer.homogeneous <- match.arg(layer.homogeneous)
  weights <- match.arg(weights)
  std.err <- match.arg(std.err)
  tab <- as.table(tab)

  if(length(dim(tab)) < 3)
      stop("tab must have (at least) three dimensions")

  if(is.na(nd) || nd <= 0)
      stop("nd must be strictly positive")

  # FIXME
  stopifnot(layer.homogeneous == "both")

  if(homogeneous && nd/2 > min(nrow(tab), ncol(tab)) - 1)
     stop("Number of dimensions of homogeneous model cannot exceed 2 * (min(nrow(tab), ncol(tab)) - 1)")

  if(!homogeneous && nd > min(nrow(tab), ncol(tab)) - 1)
     stop("Number of dimensions cannot exceed min(nrow(tab), ncol(tab)) - 1")

  # When gnm evaluates the formulas, tab will have been converted to a data.frame,
  # with a fallback if both names are empty
  vars <- make.names(names(dimnames(tab)))
  if(length(vars) == 0)
      vars <- c("Var1", "Var2", "Var3")

  if(diagonal)
      diagstr <- sprintf("+ %s:Diag(%s, %s) ", vars[3], vars[1], vars[2])
  else
      diagstr <- ""

  if(homogeneous) {
      f <- sprintf("Freq ~ %s + %s + %s + %s:%s + %s:%s %s+",
                   vars[1], vars[2], vars[3], vars[1], vars[3], vars[2], vars[3], diagstr)
      for(i in 1:nd)
          f <- paste(f, sprintf("+ Mult(%s, MultHomog(%s, %s), inst = %i)", vars[3], vars[1], vars[2], i))

      if(is.null(start))
          eval(parse(text=sprintf("model <- gnm(%s, data=tab, family=family, tolerance=%e, iterMax=%i, trace=%s, ...)",
                                  f, tolerance, iterMax, if(trace) "TRUE" else "FALSE")))
      else
          eval(parse(text=sprintf("model <- gnm(%s, data=tab, family=family, start=start, tolerance=%e, iterMax=%i, trace=%s, ...)",
                                  f, tolerance, iterMax, if(trace) "TRUE" else "FALSE")))
  }
  else {
      f <- sprintf("Freq ~ %s + %s + %s + %s:%s + %s:%s %s+ instances(Mult(%s, %s, %s), %i)",
                   vars[1], vars[2], vars[3], vars[1], vars[3], vars[2], vars[3], diagstr,
                   vars[3], vars[1], vars[2], nd)

      if(is.null(start)) {
          eval(parse(text=sprintf("model <- gnm(%s, data=tab, family=family, tolerance=%e, iterMax=%i, trace=%s, ...)",
                                   f, tolerance, iterMax, if(trace) "TRUE" else "FALSE")))
      }
      else {
          eval(parse(text=sprintf("model <- gnm(%s, data=tab, family=family, start=start, tolerance=%e, iterMax=%i, trace=%s, ...)",
                                  f, tolerance, iterMax, if(trace) "TRUE" else "FALSE")))
      }
  }

  if(is.null(model))
      return(NULL)

  newclasses <- if(homogeneous) c("rcL.homog", "rcL") else "rcL"
  class(model) <- c(newclasses, class(model))

  model$assoc <- assoc(model, weights=weights)

  if(std.err == "jackknife") {
      cat("Computing jackknife standard errors...\n")
      model$assoc$covmat <- jackknife(1:length(tab), w=tab, theta.assoc, model,
                                      getS3method("assoc", class(model)), NULL,
                                      family, weights)$jack.vcov
      scnames <- t(outer(paste(vars[3], ".", dimnames(tab)[[3]], sep=""),
                         c(t(outer(paste("D", 1:nd, " ", vars[1], ".", sep=""), rownames(tab), paste, sep="")),
                           t(outer(paste("D", 1:nd, " ", vars[2], ".", sep=""), colnames(tab), paste, sep=""))),
                         paste))
      rownames(model$assoc$covmat) <- colnames(model$assoc$covmat) <-
          c(t(outer(paste(vars[3], dimnames(tab)[[3]], sep=""), paste("Dim", 1:nd, sep=""), paste)),
              scnames, paste(scnames, "*", sep=""))

      model$assoc$covtype <- "jackknife"
  }
  else {
      model$assoc$covmat <- numeric(0)
      model$assoc$covtype <- "none"
  }

  model
}


assoc.rcL <- function(model, weights=c("marginal", "uniform", "none"), ...) {
  if(!inherits(model, "gnm"))
      stop("model must be a gnm object")

  # gnm doesn't include coefficients for NA row/columns, so get rid of them too
  tab <- as.table(model$data[!is.na(rownames(model$data)),
                             !is.na(colnames(model$data)),
                             !is.na(dimnames(model$data)[3])])

  # Weight with marginal frequencies, cf. Becker & Clogg (1994), p. 83-84, and Becker & Clogg (1989), p. 144.
  weights <- match.arg(weights)
  if(weights == "marginal") {
      rp <- prop.table(margin.table(tab, 1))
      cp <- prop.table(margin.table(tab, 2))
  }
  else if(weights == "uniform") {
      rp <- rep(1/nrow(tab), nrow(tab))
      cp <- rep(1/ncol(tab), ncol(tab))
  }
  else {
      rp <- rep(1, nrow(tab))
      cp <- rep(1, ncol(tab))
  }

  # When gnm evaluates the formulas, tab will have been converted to a data.frame,
  # with a fallback if both names are empty
  vars <- make.names(names(dimnames(tab)))
  if(length(vars) == 0)
      vars <- c("Var1", "Var2", "Var3")

  # Prepare matrices before filling them
  row <- matrix(NA, nrow(tab), 0)
  col <- matrix(NA, ncol(tab), 0)
  layer <- matrix(NA, dim(tab)[3], 0)

  nd <- 0
  while(TRUE) {
      mu <- coef(model)[pickCoef(model, sprintf("inst = %i.*\\.%s", nd+1, vars[1]))]
      nu <- coef(model)[pickCoef(model, sprintf("inst = %i.*\\.%s", nd+1, vars[2]))]
      phi <- coef(model)[pickCoef(model, sprintf("inst = %i.*\\.%s", nd+1, vars[3]))]

      if(!(length(mu) == nrow(tab) && length(nu) == ncol(tab) && length(phi) == dim(tab)[3]))
        break

      # This is a correct dimension, add it
      nd <- nd + 1

      row <- cbind(row, mu)
      col <- cbind(col, nu)
      layer <- cbind(layer, phi)
  }

  if(nd <= 0) {
    mu <- coef(model)[pickCoef(model, sprintf("\\.%s", vars[1]))]
    nu <- coef(model)[pickCoef(model, sprintf("\\.%s", vars[2]))]
    phi <- coef(model)[pickCoef(model, sprintf("\\.%s", vars[3]))]

      if(length(mu) == nrow(tab) && length(nu) == ncol(tab) && length(phi) == dim(tab)[3]) {
          nd <- 1

          row <- cbind(row, mu)
          col <- cbind(col, nu)
          layer <- cbind(layer, phi)
      }
      else {
        stop("No dimensions found. Are you sure this is a row-column association model with layer effect?")
      }
  }

  if(nd < 1)
      stop("No dimensions found. Are you sure this is a row-column association model with layer effect?")

  dg <- matrix(NA, dim(tab)[3], nrow(tab))
  dg[] <- coef(model)[pickCoef(model, "Diag\\(")]

  # Center
  row <- sweep(row, 2, colSums(sweep(row, 1, rp/sum(rp), "*")), "-")
  col <- sweep(col, 2, colSums(sweep(col, 1, cp/sum(cp), "*")), "-")

  # Scale
  phi.row <- sqrt(colSums(sweep(row^2, 1, rp, "*")))
  phi.col <- sqrt(colSums(sweep(col^2, 1, cp, "*")))
  row <- sweep(row, 2, phi.row, "/")
  col <- sweep(col, 2, phi.col, "/")
  layer <- sweep(layer, 2, phi.row * phi.col, "*")

  # Order dimensions according to phi on first layer category
  ord <- order(abs(layer[1,]), decreasing=TRUE)
  layer <- layer[, ord, drop=FALSE]
  row <- row[, ord, drop=FALSE]
  col <- col[, ord, drop=FALSE]

  # By convention, keep layer coefficients positive for the first layer category
  for(i in 1:nd) {
      if(layer[1,i] < 0) {
          layer[,i] <- -layer[,i]
          row[,i] <- -row[,i]
      }
  }

  # Since the sign of scores is arbitrary, conventionnally choose positive scores
  # for the first row category: this ensures the results are stable when jackknifing.
  for(i in 1:nd) {
      if(row[1,i] < 0) {
          row[,i] <- -row[,i]
          col[,i] <- -col[,i]
      }
  }

  ## Prepare objects
  colnames(row) <- colnames(col) <- colnames(layer) <- paste("Dim", 1:nd, sep="")

  rownames(row) <- rownames(tab)
  rownames(col) <- colnames(tab)
  rownames(layer) <- dimnames(tab)[[3]]

  if(length(dg) > 0) {
      colnames(dg) <- if(all(rownames(tab) == colnames(tab))) rownames(tab)
                      else paste(rownames(tab), colnames(tab), sep=":")
      rownames(dg) <- dimnames(tab)[[3]]
  }

  obj <- list(phi = layer, row = row, col = col, diagonal = dg,
              weighting = weights, row.weights = rp, col.weights = cp)

  class(obj) <- c("rcL.assoc", "assoc")
  obj
}

assoc.rcL.homog <- function(model, weights=c("marginal", "uniform", "none"), ...) {
  if(!inherits(model, "gnm"))
      stop("model must be a gnm object")

  # gnm doesn't include coefficients for NA row/columns, so get rid of them too
  tab <- as.table(model$data[!is.na(rownames(model$data)),
                             !is.na(colnames(model$data)),
                             !is.na(dimnames(model$data)[3])])

  # Weight with marginal frequencies, cf. Becker & Clogg (1994), p. 83-84, and Becker & Clogg (1989), p. 144.
  weights <- match.arg(weights)
  if(weights == "marginal")
      p <- prop.table(margin.table(tab, 1) + margin.table(tab, 2))
  else if(weights == "uniform")
      p <- rep(1/nrow(tab), nrow(tab))
  else
      p <- rep(1, nrow(tab))

  # When gnm evaluates the formulas, tab will have been converted to a data.frame,
  # with a fallback if both names are empty
  vars <- make.names(names(dimnames(tab)))
  if(length(vars) == 0)
      vars <- c("Var1", "Var2", "Var3")

  # Prepare matrices before filling them
  sc <- matrix(NA, nrow(tab), 0)
  layer <- matrix(NA, dim(tab)[3], 0)

  nd <- 0
  while(TRUE) {
      mu <- coef(model)[pickCoef(model, sprintf("inst = %i.*\\.%s\\|%s", nd+1, vars[1], vars[2]))]
      phi <- coef(model)[pickCoef(model, sprintf("inst = %i.*\\.%s", nd+1, vars[3]))]

      if(!(length(mu) == nrow(tab) && length(phi) == dim(tab)[3]))
        break

      # This is a correct dimension, add it
      nd <- nd + 1

      sc <- cbind(sc, mu)
      layer <- cbind(layer, phi)
  }

  if(nd <= 0) {
    mu <- coef(model)[pickCoef(model, sprintf("\\.%s\\|%s", vars[1], vars[2]))]
    phi <- coef(model)[pickCoef(model, sprintf("\\.%s", vars[3]))]

      if(length(mu) == nrow(tab) && length(phi) == dim(tab)[3]) {
          nd <- 1

          sc <- cbind(sc, mu)
          layer <- cbind(layer, phi)
      }
      else {
        stop("No dimensions found. Are you sure this is an homogeneous row-column association model with layer effect?")
      }
  }

  if(nd < 1)
      stop("No dimensions found. Are you sure this is an homogeneous row-column association model with layer effect?")

  dg <- matrix(NA, dim(tab)[3], nrow(tab))
  dg[] <- coef(model)[pickCoef(model, "Diag\\(")]

  # Center
  sc <- sweep(sc, 2, colSums(sweep(sc, 1, p/sum(p), "*")), "-")

  # Scale
  phi <- sqrt(colSums(sweep(sc^2, 1, p/sum(p), "*")))
  sc <- sweep(sc, 2, phi, "/")
  layer <- sweep(layer, 2, phi, "*")

  # Order dimensions according to phi on first layer category
  ord <- order(layer[1,], decreasing=TRUE)
  layer <- layer[, ord, drop=FALSE]
  sc <- sc[, ord, drop=FALSE]

  # By convention, keep layer coefficients positive for the first layer category
  for(i in 1:nd) {
      if(layer[1,i] < 0) {
          layer[,i] <- -layer[,i]
          row[,i] <- -row[,i]
      }
  }

  # Since the sign of scores is arbitrary, conventionnally choose positive scores
  # for the first category: this ensures the results are stable when jackknifing.
  for(i in 1:nd) {
      if(sc[1,i] < 0)
          sc[,i] <- -sc[,i]
  }

  ## Prepare objects
  colnames(sc) <- colnames(layer) <- paste("Dim", 1:nd, sep="")

  rownames(sc) <- rownames(tab)
  rownames(layer) <- dimnames(tab)[[3]]

  if(length(dg) > 0) {
      colnames(dg) <- if(all(rownames(tab) == colnames(tab))) rownames(tab)
                      else paste(rownames(tab), colnames(tab), sep=":")
      rownames(dg) <- dimnames(tab)[[3]]
  }

  obj <- list(phi = layer, row = sc, col = sc, diagonal = dg,
              weighting = weights, row.weights = p, col.weights = p)

  class(obj) <- c("rcL.homog.assoc", "rcL.assoc", "assoc")
  obj
}
