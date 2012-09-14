## RC(M) model

rc <- function(tab, nd=1, symmetric=FALSE, diagonal=FALSE,
               weighting=c("marginal", "uniform", "none"), se=c("none", "jackknife"),
               family=poisson, start=NA, etastart=NULL, tolerance=1e-6, iterMax=5000,
               trace=TRUE, verbose=TRUE, ...) {
  weighting <- match.arg(weighting)
  se <- match.arg(se)
  tab <- as.table(tab)

  if(length(dim(tab)) < 2)
      stop("tab must have (at least) two dimensions")

  if(is.na(nd) || nd <= 0)
      stop("nd must be strictly positive")

  if(symmetric && nrow(tab) != ncol(tab))
      stop("tab must be a square table for symmetric model")

  if(symmetric && !all(rownames(tab) == colnames(tab)))
      stop("tab must have identical row and column names for symmetric model")

  if(symmetric && nd/2 > min(nrow(tab), ncol(tab)) - 1)
      stop("Number of dimensions of symmetric model cannot exceed 2 * (min(nrow(tab), ncol(tab)) - 1)")

  if(!symmetric && nd > min(nrow(tab), ncol(tab)) - 1)
      stop("Number of dimensions cannot exceed min(nrow(tab), ncol(tab)) - 1")

  if(length(dim(tab)) > 2)
      tab <- margin.table(tab, 1:2)

  # When gnm evaluates the formulas, tab will have been converted to a data.frame,
  # with a fallback if both names are empty
  vars <- make.names(names(dimnames(tab)))
  if(length(vars) == 0)
      vars <- c("Var1", "Var2")

  if(diagonal)
      diagstr <- sprintf("+ Diag(%s, %s) ", vars[1], vars[2])
  else
      diagstr <- ""


  nastart <- length(start) == 1 && is.na(start)

  if(symmetric) {
      f <- sprintf("Freq ~ %s + %s %s+ instances(MultHomog(%s, %s), %i)",
                   vars[1], vars[2], diagstr, vars[1], vars[2], nd)


      if(nastart)
          start <- NULL
  }
  else {
      if(nastart) {
          cat("Running base model to find starting values...\n")

          # We need to handle ... manually, else they would not be found when modelFormula() evaluates the call
          args <- list(formula=as.formula(sprintf("Freq ~ %s + %s %s", vars[1], vars[2], diagstr)),
                       data=tab, family=family,
                       tolerance=tolerance, iterMax=iterMax)

          base <- do.call("gnm", c(args, list(...)))

          # residSVD evaluates the variable names in parent.frame(), which uses any object
          # called "vars" in the global environment if not handled like this
          res <- eval(parse(text=sprintf("residSVD(base, %s, %s, %i)", vars[1], vars[2], nd)))
          start <- c(parameters(base), res)

          if(is.null(etastart))
              etastart <- as.numeric(predict(base))

          cat("Running real model...\n")
      }

      f <- sprintf("Freq ~ %s + %s %s+ instances(Mult(%s, %s), %i)",
                   vars[1], vars[2], diagstr, vars[1], vars[2], nd)
  }

  # We need to handle ... manually, else they would not be found when modelFormula() evaluates the call
  args <- list(formula=as.formula(f), data=tab,
               family=family, start=start, etastart=etastart,
               tolerance=tolerance, iterMax=iterMax, trace=trace)
  model <- do.call("gnm", c(args, list(...)))

  if(is.null(model))
      return(NULL)

  newclasses <- if(symmetric) c("rc.symm", "rc") else "rc"
  class(model) <- c(newclasses, class(model))

  model$call <- match.call()

  model$assoc <- assoc(model, weighting=weighting)


  if(se == "jackknife") {
      cat("Computing jackknife standard errors...\n")
      model$assoc$covmat <- jackknife(1:length(tab), w=tab, theta.assoc, model,
                                      getS3method("assoc", class(model)), NULL,
                                      family, weighting, base=base, verbose=verbose)$jack.vcov
      scnames <- c(t(outer(paste("D", 1:nd, " ", vars[1], ".", sep=""), rownames(tab), paste, sep="")),
                   t(outer(paste("D", 1:nd, " ", vars[2], ".", sep=""), colnames(tab), paste, sep="")))
      rownames(model$assoc$covmat) <- colnames(model$assoc$covmat) <-
          c(paste("Dim", 1:nd, sep=""), scnames, paste(scnames, "*", sep=""))

      model$assoc$covtype <- "jackknife"
  }
  else {
      model$assoc$covmat <- numeric(0)
      model$assoc$covtype <- "none"
  }

  model
}

assoc.rc <- function(model, weighting=c("marginal", "uniform", "none"), ...) {
  if(!inherits(model, "gnm"))
      stop("model must be a gnm object")

  # gnm doesn't include coefficients for NA row/columns, so get rid of them too
  if(length(dim(model$data)) == 2)
      tab <- as.table(model$data[!is.na(rownames(model$data)),
                                 !is.na(colnames(model$data))])
  else if(length(dim(model$data)) == 3)
      tab <- as.table(model$data[!is.na(rownames(model$data)),
                                 !is.na(colnames(model$data)),
                                 !is.na(dimnames(model$data)[[3]])])
  else
      stop("Only two and three dimensional tables are supported")

  # Weight with marginal frequencies, cf. Becker & Clogg (1994), p. 83-84, and Becker & Clogg (1989), p. 144.
  weighting <- match.arg(weighting)
  if(weighting == "marginal") {
      rp <- prop.table(margin.table(tab, 1))
      cp <- prop.table(margin.table(tab, 2))
  }
  else if(weighting == "uniform") {
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
      vars <- c("Var1", "Var2")

  # Prepare matrices before filling them
  row <- matrix(NA, nrow(tab), 0)
  col <- matrix(NA, ncol(tab), 0)

  nd <- 0
  while(TRUE) {
      mu <- parameters(model)[pickCoef(model, sprintf("Mult\\(\\., \\Q%s\\E, inst = %i\\)\\.\\Q%s\\E(\\Q%s\\E)$",
                                                vars[2], nd + 1, vars[1],
                                                paste(rownames(tab), collapse="\\E|\\Q")))]
      nu <- parameters(model)[pickCoef(model, sprintf("Mult\\(\\Q%s\\E, \\., inst = %i\\)\\.\\Q%s\\E(\\Q%s\\E)$",
                                                vars[1], nd + 1, vars[2],
                                                paste(colnames(tab), collapse="\\E|\\Q")))]

      if(!(length(mu) == nrow(tab) && length(nu) == ncol(tab)))
          break

      # This is a correct dimension, add it
      nd <- nd + 1

      row <- cbind(row, mu)
      col <- cbind(col, nu)
  }

  if(nd <= 0) {
      mu <- parameters(model)[pickCoef(model, sprintf("Mult\\(\\., \\Q%s\\E)\\.\\Q%s\\E(\\Q%s\\E)$",
                                                vars[2], vars[1],
                                                paste(rownames(tab), collapse="\\E|\\Q")))]
      nu <- parameters(model)[pickCoef(model, sprintf("Mult\\(\\Q%s\\E, \\.)\\.\\Q%s\\E(\\Q%s\\E)$",
                                                vars[1], vars[2],
                                                paste(colnames(tab), collapse="\\E|\\Q")))]

      if(length(mu) == nrow(tab) && length(nu) == ncol(tab)) {
          nd <- 1

          row <- cbind(row, mu)
          col <- cbind(col, nu)
      }
      else {
          stop("No dimensions found. Are you sure this is a row-column association model?")
      }
  }

  if(length(pickCoef(model, "Diag\\(")) > nrow(tab))
      dg <- matrix(parameters(model)[pickCoef(model, "Diag\\(")], ncol=nrow(tab))
  else if(length(pickCoef(model, "Diag\\(")) > 0)
      dg <- matrix(parameters(model)[pickCoef(model, "Diag\\(")], 1, nrow(tab))
  else
      dg <- numeric(0)


  # Center
  row <- sweep(row, 2, colSums(sweep(row, 1, rp/sum(rp), "*")), "-")
  col <- sweep(col, 2, colSums(sweep(col, 1, cp/sum(cp), "*")), "-")

  # Technique proposed in Goodman (1991), Appendix 4
  lambda <- matrix(0, nrow(tab), ncol(tab))
  for(i in 1:nd) {
      lambda <- lambda + row[,i] %o% col[,i]
  }
  lambda0 <- lambda * sqrt(rp %o% cp) # Eq. A.4.3
  sv <- svd(lambda0)
  row[] <- diag(1/sqrt(rp)) %*% sv$u[,1:nd] # Eq. A.4.7
  col[] <- diag(1/sqrt(cp)) %*% sv$v[,1:nd] # Eq. A.4.7
  phi <- sv$d[1:nd]

  # Since the sign of scores is arbitrary, conventionally choose positive scores
  # for the first row category: this ensures the results are stable when jackknifing.
  for(i in 1:nd) {
      if(row[1,i] < 0) {
          row[,i] <- -row[,i]
          col[,i] <- -col[,i]
      }
  }

  ## Prepare objects
  phi <- rbind(c(phi))
  dim(row)[3] <- dim(col)[3] <- 1
  colnames(row) <- colnames(col) <- colnames(phi) <- paste("Dim", 1:nd, sep="")
  rownames(row) <- rownames(tab)
  rownames(col) <- colnames(tab)

  if(length(dg) > 0) {
      # Diag() sorts coefficients alphabetically!
      dg[,order(rownames(tab))] <- dg

      colnames(dg) <- if(all(rownames(tab) == colnames(tab))) rownames(tab)
                      else paste(rownames(tab), colnames(tab), sep=":")

      if(nrow(dg) > 1)
          rownames(dg) <- dimnames(tab)[[3]]
      else
          rownames(dg) <- "All levels"
  }

  obj <- list(phi = phi, row = row, col = col, diagonal = dg,
              weighting = weighting, row.weights = rp, col.weights = cp)

  class(obj) <- c("assoc.rc", "assoc")
  obj
}

## RC(M) model with symmetric row and column scores
assoc.rc.symm <- function(model, weighting=c("marginal", "uniform", "none"), ...) {
  if(!inherits(model, "gnm"))
      stop("model must be a gnm object")

  # gnm doesn't include coefficients for NA row/columns, so get rid of them too
  if(length(dim(model$data)) == 2)
      tab <- as.table(model$data[!is.na(rownames(model$data)),
                                 !is.na(colnames(model$data))])
  else if(length(dim(model$data)) == 3)
      tab <- as.table(model$data[!is.na(rownames(model$data)),
                                 !is.na(colnames(model$data)),
                                 !is.na(dimnames(model$data)[[3]])])
  else
      stop("Only two and three dimensional tables are supported")

  # Weight with marginal frequencies, cf. Becker & Clogg (1994), p. 83-84, and Becker & Clogg (1989), p. 144.
  weighting <- match.arg(weighting)
  if(weighting == "marginal")
      p <- prop.table(margin.table(tab, 1) + margin.table(tab, 2))
  else if(weighting == "uniform")
      p <- rep(1/nrow(tab), nrow(tab))
  else
      p <- rep(1, nrow(tab))

  vars <- make.names(names(dimnames(tab)))
  if(length(vars) == 0)
      vars <- c("Var1", "Var2")

  sc <- matrix(NA, nrow(tab), 0)

  nd <- 0
  while(TRUE) {
      mu <- parameters(model)[pickCoef(model, sprintf("MultHomog\\(\\Q%s\\E\\, \\Q%s\\E\\, inst = %i\\)(\\Q%s\\E)$",
                                                vars[1], vars[2], nd + 1,
                                                paste(c(rownames(tab), colnames(tab)), collapse="\\E|\\Q")))]

      if(!(length(mu) == nrow(tab)))
          break

      # This is a correct dimension, add it
      nd <- nd + 1

      sc <- cbind(sc, mu)
  }

  if(nd <= 0) {
      mu <- parameters(model)[pickCoef(model, sprintf("MultHomog\\(\\Q%s\\E\\, \\Q%s\\E\\)(\\Q%s\\E)$",
                                                vars[1], vars[2],
                                                paste(c(rownames(tab), colnames(tab)), collapse="\\E|\\Q")))]

      if(length(mu) == nrow(tab)) {
          nd <- 1
          sc <- cbind(sc, mu)
      }
      else {
          stop("No dimensions found. Are you sure this is a row-column association model with symmetric row and column scores?")
      }
  }

  if(length(pickCoef(model, "Diag\\(")) > nrow(tab))
      dg <- matrix(parameters(model)[pickCoef(model, "Diag\\(")], ncol=nrow(tab))
  else if(length(pickCoef(model, "Diag\\(")) > 0)
      dg <- matrix(parameters(model)[pickCoef(model, "Diag\\(")], 1, nrow(tab))
  else
      dg <- numeric(0)


  # Center
  sc <- sweep(sc, 2, colSums(sweep(sc, 1, p/sum(p), "*")), "-")

  # Technique proposed in Goodman (1991), Appendix 4, but with eigenvalues decomposition
  lambda <- matrix(0, nrow(tab), ncol(tab))
  for(i in 1:nd)
      lambda <- lambda + sc[,i] %o% sc[,i]
  lambda0 <- lambda * sqrt(p %o% p) # Eq. A.4.3
  eigen <- eigen(lambda0, symmetric=TRUE)
  sc[,1:nd] <- diag(1/sqrt(p)) %*% eigen$vectors[,1:nd] # Eq. A.4.7
  phi <- eigen$values[1:nd]

  # Since the sign of scores is arbitrary, conventionally choose positive scores
  # for the first category: this ensures the results are stable when jackknifing.
  for(i in 1:nd) {
      if(sc[1,i] < 0)
          sc[,i] <- -sc[,i]
  }

  ## Prepare objects
  phi <- rbind(c(phi))
  dim(sc)[3] <- 1
  colnames(sc) <- colnames(phi) <- paste("Dim", 1:nd, sep="")
  rownames(sc) <- rownames(tab)

  if(length(dg) > 0) {
      # Diag() sorts coefficients alphabetically!
      dg[,order(rownames(tab))] <- dg

      colnames(dg) <- if(all(rownames(tab) == colnames(tab))) rownames(tab)
                      else paste(rownames(tab), colnames(tab), sep=":")

      if(nrow(dg) > 1)
          rownames(dg) <- dimnames(tab)[[3]]
      else
          rownames(dg) <- "All levels"
  }

  obj <- list(phi = phi, row = sc, col= sc, diagonal = dg,
              weighting = weighting, row.weights = p, col.weights = p)

  class(obj) <- c("assoc.rc", "assoc.symm", "assoc")
  obj
}

assoc <- function(model, weighting, ...) UseMethod("assoc", model)
