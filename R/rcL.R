## RC(M)-L model Wong(2010))

rcL <- function(tab, nd=1, layer.effect=c("homogeneous.scores", "heterogeneous", "none"),
                symmetric=FALSE, diagonal=c("none", "heterogeneous", "homogeneous"),
                weighting=c("marginal", "uniform", "none"), std.err=c("none", "jackknife"),
                family=poisson, start=NULL, tolerance=1e-12, iterMax=5000, trace=TRUE, ...) {
  layer.effect <- match.arg(layer.effect)
  diagonal <- match.arg(diagonal)
  weighting <- match.arg(weighting)
  std.err <- match.arg(std.err)
  tab <- as.table(tab)

  if(length(dim(tab)) < 3)
      stop("tab must have (at least) three dimensions")

  if(is.na(nd) || nd <= 0)
      stop("nd must be strictly positive")

  if(symmetric && nd/2 > min(nrow(tab), ncol(tab)) - 1)
      stop("Number of dimensions of symmetric model cannot exceed 2 * (min(nrow(tab), ncol(tab)) - 1)")

  if(!symmetric && nd > min(nrow(tab), ncol(tab)) - 1)
      stop("Number of dimensions cannot exceed min(nrow(tab), ncol(tab)) - 1")

  if(length(dim(tab)) > 3)
      tab <- margin.table(tab, 1:3)

  # When gnm evaluates the formulas, tab will have been converted to a data.frame,
  # with a fallback if both names are empty
  vars <- make.names(names(dimnames(tab)))
  if(length(vars) == 0)
      vars <- c("Var1", "Var2", "Var3")

  if(diagonal != "none")
      diagstr <- sprintf("+ Diag(%s, %s) ", vars[1], vars[2])
  else
      diagstr <- ""


  f1 <- sprintf("Freq ~ %s + %s + %s + %s:%s + %s:%s",
                vars[1], vars[2], vars[3], vars[1], vars[3], vars[2], vars[3])

  base <- NULL

  if(is.null(start)) {
      cat("Running base model to find starting values...\n")

      # We need to handle ... manually, else they would not be found when modelFormula() evaluates the call
      args <- list(formula=eval(as.formula(paste(f1, diagstr))), data=tab,
                   family=family,
                   tolerance=1e-3, iterMax=iterMax)
      dots <- as.list(substitute(list(...)))[-1]
      args <- c(args, dots)

      base <- do.call("gnm", args)
  }

  if(symmetric) {
      if(layer.effect == "homogeneous.scores") {
          f2 <- ""

          for(i in 1:nd)
              f2 <- paste(f2, sprintf("+ Mult(%s, MultHomog(%s, %s), inst = %i)",
                                      vars[3], vars[1], vars[2], i))

          if(is.null(start))
              start <- c(coef(base), rep(NA, nd * (dim(tab)[3] + nrow(tab))))
      }
      else if(layer.effect == "heterogeneous") {
          stop("Symmetric association with heterogeneous layer effect is currently not supported")

          f2 <- ""

          for(i in 1:nd)
              f2 <- paste(f2, sprintf("+ MultHomog(%s:%s, %s:%s, inst = %i)", 
                                      vars[3], vars[1], vars[3], vars[2], i))

          if(is.null(start))
              start <- c(coef(base), rep(NA, nd * nrow(tab)))
      }
      else {
          f2 <- sprintf("+ instances(MultHomog(%s, %s), %i)", vars[1], vars[2], nd)

          if(is.null(start))
              start <- c(coef(base), rep(NA, nd * nrow(tab)))
      }
  }
  else {
      if(layer.effect == "homogeneous.scores") {
          f2 <- sprintf("+ instances(Mult(%s, %s, %s), %i)",
                        vars[3], vars[1], vars[2], nd)

          if(is.null(start))
              start <- c(coef(base), rep(NA, nd * (nrow(tab) + ncol(tab) + dim(tab)[3])))
      }
      else if(layer.effect == "heterogeneous") {
          f2 <- sprintf("+ instances(Mult(%s:%s, %s:%s), %i)",
                        vars[3], vars[1], vars[3], vars[2], nd)

          if(is.null(start))
              start <- c(coef(base), rep(NA, nd * dim(tab)[3] * (nrow(tab) + ncol(tab))))
      }
      else {
          f2 <- sprintf("+ instances(Mult(%s, %s), %i)",
                        vars[1], vars[2], nd)

          if(is.null(start))
              start <- c(coef(base), rep(NA, nd * (nrow(tab) + ncol(tab))))
      }
  }

  # Heterogeneous diagonal parameters can make the convergence really slow unless
  # correct starting values are used
  if(!is.null(base) && diagonal == "heterogeneous") {
      cat("Running model with homogeneous diagonal parameters...\n")

      # We need to handle ... manually, else they would not be found when modelFormula() evaluates the call
      args <- list(formula=as.formula(paste(f1, diagstr, f2)),
                   data=tab, family=family, start=start,
                   tolerance=1e-3, iterMax=iterMax, trace=trace)
      dots <- as.list(substitute(list(...)))[-1]
      args <- c(args, dots)

      base2 <- do.call("gnm", args)

      start <- c(head(coef(base2), length(coef(base)) - nrow(tab)),
                 rep(NA, nrow(tab) * dim(tab)[3]), tail(coef(base2),
                 length(start) - length(coef(base))))
  }

  if(diagonal == "heterogeneous")
      diagstr <- sprintf("+ %s:Diag(%s, %s) ", vars[3], vars[1], vars[2])

  if(!is.null(base))
      cat("Running real model...\n")

  # We need to handle ... manually, else they would not be found when modelFormula() evaluates the call
  args <- list(formula=eval(as.formula(paste(f1, diagstr, f2))), data=substitute(tab),
               family=substitute(family), start=start,
               tolerance=tolerance, iterMax=iterMax, trace=trace)
  dots <- as.list(substitute(list(...)))[-1]
  args <- c(args, dots)

  model <- do.call("gnm", args)

  if(is.null(model))
      return(NULL)

  newclasses <- if(symmetric) c("rcL.symm", "rcL") else "rcL"
  class(model) <- c(newclasses, class(model))

  if(layer.effect == "none") {
      if(symmetric)
          model$assoc <- assoc.rc.symm(model, weighting=weighting)
      else
          model$assoc <- assoc.rc(model, weighting=weighting)
  }
  else {
      model$assoc <- assoc(model, weighting=weighting)
  }

  class(model$assoc) <- if(symmetric) c("assoc.rcL", "assoc.symm", "assoc")
                        else c("assoc.rcL", "assoc")


  if(std.err == "jackknife") {
      cat("Computing jackknife standard errors...\n")
      model$assoc$covmat <- jackknife(1:length(tab), w=tab, theta.assoc, model,
                                      getS3method("assoc", class(model)), NULL,
                                      family, weighting)$jack.vcov
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


assoc.rcL <- function(model, weighting=c("marginal", "uniform", "none"), ...) {
  if(!inherits(model, "gnm"))
      stop("model must be a gnm object")

  # gnm doesn't include coefficients for NA row/columns, so get rid of them too
  tab <- as.table(model$data[!is.na(rownames(model$data)),
                             !is.na(colnames(model$data)),
                             !is.na(dimnames(model$data)[3])])

  nr <- nrow(tab)
  nc <- ncol(tab)
  nl <- dim(tab)[3]

  # Weight with marginal frequencies, cf. Becker & Clogg (1994), p. 83-84, and Becker & Clogg (1989), p. 144.
  weighting <- match.arg(weighting)
  if(weighting == "marginal") {
      rp <- prop.table(margin.table(tab, 1))
      cp <- prop.table(margin.table(tab, 2))
  }
  else if(weighting == "uniform") {
      rp <- rep(1/nr, nr)
      cp <- rep(1/nc, nc)
  }
  else {
      rp <- rep(1, nr)
      cp <- rep(1, nc)
  }

  # When gnm evaluates the formulas, tab will have been converted to a data.frame,
  # with a fallback if both names are empty
  vars <- make.names(names(dimnames(tab)))
  if(length(vars) == 0)
      vars <- c("Var1", "Var2", "Var3")

  # Find out the number of dimensions
  nd <- 0
  while(length(pickCoef(model, sprintf("Mult\\(.*\\Q%s\\E.*inst = %i\\)", vars[3], nd + 1))) > 0)
      nd <- nd + 1

  homogeneous <- TRUE

  # Only one dimension, or none
  if(nd <= 0) {
      mu <- coef(model)[pickCoef(model, sprintf("Mult\\(.*\\Q%s\\E.*\\).*[.:]\\Q%s\\E(\\Q%s\\E)$",
                                                vars[3], vars[1],
                                                paste(rownames(tab), collapse="\\E|\\Q")))]
      nu <- coef(model)[pickCoef(model, sprintf("Mult\\(.*\\Q%s\\E.*\\).*[.:]\\Q%s\\E(\\Q%s\\E)$",
                                                vars[3], vars[2],
                                                paste(colnames(tab), collapse="\\E|\\Q")))]
      phi <- coef(model)[pickCoef(model, sprintf("Mult\\(.*[.:]\\Q%s\\E(\\Q%s\\E)$",
                                                 vars[3],
                                                 paste(dimnames(tab)[[3]], collapse="\\E|\\Q")))]

      # Homogeneous scores for rows and/or columns
      if(length(phi) == nl &&
         (length(mu) == nr || length(mu) == nr * nl) &&
         (length(nu) == nc || length(nu) == nc * nl)) {
          nd <- 1

          if(length(mu) == nr) {
              row <- array(mu, dim=c(nr, 1, 1))
          }
          else {
              row <- array(mu, dim=c(nr, nl, 1))
              homogeneous <- FALSE
          }

          if(length(nu) == nc) {
              row <- array(nu, dim=c(nc, 1, 1))
          }
          else {
              row <- array(nu, dim=c(nc, nl, 1))
              homogeneous <- FALSE
          }

          layer <- matrix(phi, nl, 1)
      }
      else {
          stop("No dimensions found. Are you sure this is a row-column association model with layer effect?")
      }
  }

  # Several dimensions: prepare arrays before filling them
  row <- array(NA, dim=c(nr, nd, nl))
  col <- array(NA, dim=c(nc, nd, nl))
  layer <- matrix(NA, nl, nd)

  for(i in 1:nd) {
      mu <- coef(model)[pickCoef(model, sprintf("Mult\\(.*\\Q%s\\E.*inst = %i.*[.:]\\Q%s\\E(\\Q%s\\E)$",
                                                vars[3], i, vars[1],
                                                paste(rownames(tab), collapse="\\E|\\Q")))]
      nu <- coef(model)[pickCoef(model, sprintf("Mult\\(.*\\Q%s\\E.*inst = %i.*[.:]\\Q%s\\E(\\Q%s\\E)$",
                                                vars[3], i, vars[2],
                                                paste(colnames(tab), collapse="\\E|\\Q")))]
      phi <- coef(model)[pickCoef(model, sprintf("Mult\\(.*inst = %i.*\\.\\Q%s\\E(\\Q%s\\E)$",
                                                 i, vars[3],
                                                 paste(dimnames(tab)[[3]], collapse="\\E|\\Q")))]

      # Homogeneous scores for rows and/or columns
      if(length(phi) == nl &&
         (length(mu) == nr || length(mu) == nr * nl) &&
         (length(nu) == nc || length(nu) == nc * nl)) {
          if(length(mu) == nr)
              row[,i,1] <- mu
          else {
              row[,i,] <- t(matrix(mu, nl, nr))
              homogeneous <- FALSE
          }

          if(length(nu) == nc)
              col[,i,1] <- nu
          else {
              col[,i,] <- t(matrix(nu, nl, nc))
              homogeneous <- FALSE
          }

          layer[,i] <- phi
      }
      # Fully heterogeneous scores
      else if(length(phi) == 0 &&
              length(mu) == nr * nl &&
              length(nu) == nc * nl) {
          homogeneous <- FALSE
          row[,i,] <- t(matrix(mu, nl, nr))
          col[,i,] <- t(matrix(nu, nl, nc))
          layer[,i] <- 1
      }
      else {
          stop("Invalid dimensions found. Are you sure this is a row-column association model with layer effect?")
      }
  }

  if(length(pickCoef(model, "Diag\\(")) > nr)
      dg <- matrix(coef(model)[pickCoef(model, "Diag\\(")], nl, nr)
  else if(length(pickCoef(model, "Diag\\(")) > 0)
      dg <- matrix(coef(model)[pickCoef(model, "Diag\\(")], 1, nr)
  else
      dg <- numeric(0)


  # Center
  row <- sweep(row, 2:3, margin.table(sweep(row, 1, rp/sum(rp), "*"), 2:3), "-")
  col <- sweep(col, 2:3, margin.table(sweep(col, 1, cp/sum(cp), "*"), 2:3), "-")

  if(homogeneous) {
      # Scale
      phi.row <- sqrt(margin.table(sweep(row[,,1, drop=FALSE]^2, 1, rp, "*"), 2))
      phi.col <- sqrt(margin.table(sweep(col[,,1, drop=FALSE]^2, 1, cp, "*"), 2))
      row <- sweep(row[,,1, drop=FALSE], 2, phi.row, "/")
      col <- sweep(col[,,1, drop=FALSE], 2, phi.col, "/")
      layer <- sweep(layer, 2, phi.row * phi.col, "*")

      # Order dimensions according to phi on first layer category
      ord <- order(abs(layer[1,]), decreasing=TRUE)
      layer <- layer[,ord, drop=FALSE]
      row <- row[,ord,, drop=FALSE]
      col <- col[,ord,, drop=FALSE]
  }
  else {
      for(l in 1:nl) {
          # Technique proposed in Goodman (1991), Appendix 4
          lambda <- matrix(0, nr, nc)
          for(i in 1:nd) {
              lambda <- lambda + layer[l,i] * row[,i,l] %o% col[,i,l]
          }
          lambda0 <- lambda * sqrt(rp %o% cp) # Eq. A.4.3
          sv <- svd(lambda0)
          row[,,l] <- diag(1/sqrt(rp)) %*% sv$u[,1:nd] # Eq. A.4.7
          col[,,l] <- diag(1/sqrt(cp)) %*% sv$v[,1:nd] # Eq. A.4.7
          layer[l,] <- sv$d[1:nd]
      }
  }

  # By convention, keep layer coefficients positive for the first layer category
  if(homogeneous) {
      for(i in 1:nd) {
          if(layer[1,i] < 0) {
              layer[,i] <- -layer[,i]
              row[,i,] <- -row[,i,]
          }
      }
  }

  # Since the sign of scores is arbitrary, conventionally choose positive scores
  # for the first row category: this ensures the results are stable when jackknifing.
  for(i in 1:dim(row)[3]) {
      for(j in 1:nd) {
          if(row[1,j,i] < 0) {
              row[,j,i] <- -row[,j,i]
              col[,j,i] <- -col[,j,i]
          }
      }
  }

  ## Prepare objects
  rownames(row) <- rownames(tab)
  rownames(col) <- colnames(tab)
  colnames(row) <- colnames(col) <- colnames(layer) <- paste("Dim", 1:nd, sep="")
  rownames(layer) <- dimnames(tab)[[3]]

  if(!homogeneous)
      dimnames(row)[[3]] <- dimnames(col)[[3]] <- dimnames(tab)[[3]]
  else
      dimnames(row)[[3]] <- dimnames(col)[[3]] <- "All levels"

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

  obj <- list(phi = layer, row = row, col = col, diagonal = dg,
              weighting = weighting, row.weights = rp, col.weights = cp)

  class(obj) <- c("assoc.rcL", "assoc")
  obj
}

assoc.rcL.symm <- function(model, weighting=c("marginal", "uniform", "none"), ...) {
  if(!inherits(model, "gnm"))
      stop("model must be a gnm object")

  # gnm doesn't include coefficients for NA row/columns, so get rid of them too
  tab <- as.table(model$data[!is.na(rownames(model$data)),
                             !is.na(colnames(model$data)),
                             !is.na(dimnames(model$data)[3])])

  nr <- nrow(tab)
  nc <- ncol(tab)
  nl <- dim(tab)[3]

  # Weight with marginal frequencies, cf. Becker & Clogg (1994), p. 83-84, and Becker & Clogg (1989), p. 144.
  weighting <- match.arg(weighting)
  if(weighting == "marginal")
      p <- prop.table(margin.table(tab, 1) + margin.table(tab, 2))
  else if(weighting == "uniform")
      p <- rep(1/nr, nr)
  else
      p <- rep(1, nr)

  # When gnm evaluates the formulas, tab will have been converted to a data.frame,
  # with a fallback if both names are empty
  vars <- make.names(names(dimnames(tab)))
  if(length(vars) == 0)
      vars <- c("Var1", "Var2", "Var3")


  # Find out the number of dimensions
  nd <- 0
  while(length(pickCoef(model, paste("Mult\\(.*MultHomog\\(.*inst =", nd + 1))) > 0)
      nd <- nd + 1

  homogeneous <- TRUE

  # Only one dimension, or none
  if(nd <= 0) {
      mu <- coef(model)[pickCoef(model, sprintf("Mult\\(.*MultHomog\\(.*\\.\\Q%s\\E\\|\\Q%s\\E(\\Q%s\\E)$",
                                                vars[1], vars[2],
                                                paste(c(rownames(tab), colnames(tab)), collapse="\\E|\\Q")))]
      phi <- coef(model)[pickCoef(model, sprintf("Mult\\(.*MultHomog\\(.*\\.\\Q%s\\E(\\Q%s\\E)$", vars[3],
                                                 paste(dimnames(tab)[[3]], collapse="\\E|\\Q")))]

      # Homogeneous scores
      if(length(phi) == nl &&
         (length(mu) == nr || length(mu) == nr * nl)) {
          nd <- 1

          if(length(mu) == nr) {
              sc <- array(mu, dim=c(nr, 1, 1))
          }
          else {
              sc <- array(mu, dim=c(nr, nl, 1))
              homogeneous <- FALSE
          }

          layer <- matrix(phi, nl, 1)
      }
      else {
          stop("No dimensions found. Are you sure this is a symmetric row-column association model with layer effect?")
      }
  }

  # Several dimensions: prepare arrays before filling them
  sc <- array(NA, dim=c(nr, nd, nl))
  layer <- matrix(NA, nl, nd)

  for(i in 1:nd) {
      mu <- coef(model)[pickCoef(model, sprintf("Mult\\(.*inst = %i.*MultHomog\\(.*\\.\\Q%s\\E\\|\\Q%s\\E(\\Q%s\\E)$",
                                                i, vars[1], vars[2],
                                                paste(c(rownames(tab), colnames(tab)), collapse="\\E|\\Q")))]
      phi <- coef(model)[pickCoef(model, sprintf("Mult\\(.*MultHomog\\(.*inst = %i.*\\.\\Q%s\\E(\\Q%s\\E)$",
                                                 i, vars[3],
                                                 paste(dimnames(tab)[[3]], collapse="\\E|\\Q")))]

      # Homogeneous scores
      if(length(phi) == nl &&
         (length(mu) == nr || length(mu) == nr * nl)) {
          if(length(mu) == nr)
              sc[,i,1] <- mu
          else {
              sc[,i,] <- t(matrix(mu, nl, nr))
              homogeneous <- FALSE
          }

          layer[,i] <- phi
      }
      # Fully heterogeneous scores
      else if(length(phi) == 0 &&
              length(mu) == nr * nl) {
          homogeneous <- FALSE
          sc[,i,] <- t(matrix(mu, nl, nr))
          layer[,i] <- 1
      }
      else {
          stop("Invalid dimensions found. Are you sure this is a symmetric row-column association model with layer effect?")
      }
  }

  if(length(pickCoef(model, "Diag\\(")) > nr)
      dg <- matrix(coef(model)[pickCoef(model, "Diag\\(")], nl, nr)
  else if(length(pickCoef(model, "Diag\\(")) > 0)
      dg <- matrix(coef(model)[pickCoef(model, "Diag\\(")], 1, nr)
  else
      dg <- numeric(0)

  # Center
  sc <- sweep(sc, 2:3, margin.table(sweep(sc, 1, p/sum(p), "*"), 2:3), "-")

  if(homogeneous) {
      # Scale
      phi <- sqrt(margin.table(sweep(sc[,,1, drop=FALSE]^2, 1, p/sum(p), "*"), 2))
      sc <- sweep(sc[,,1, drop=FALSE], 2, phi, "/")
      layer <- sweep(layer, 2, phi, "*")

      # Order dimensions according to phi on first layer category
      ord <- order(layer[1,], decreasing=TRUE)
      layer <- layer[,ord, drop=FALSE]
      sc <- sc[,ord,, drop=FALSE]
  }
  else {
      for(l in 1:nl) {
          # Technique proposed in Goodman (1991), Appendix 4, but with eigenvalues decomposition
          lambda <- matrix(0, nr, nc)
          for(i in 1:nd)
              lambda <- lambda + sc[,i,l] %o% sc[,i,l]
          lambda0 <- lambda * sqrt(p %o% p) # Eq. A.4.3
          eigen <- eigen(lambda0, symmetric=TRUE)
          sc[,,l] <- diag(1/sqrt(p)) %*% eigen$vectors[,1:nd] # Eq. A.4.7
          layer[l,] <- eigen$values[1:nd]
      }
  }

  # Since the sign of scores is arbitrary, conventionally choose positive scores
  # for the first category: this ensures the results are stable when jackknifing.
  for(i in 1:dim(sc)[3]) {
      for(j in 1:nd) {
          if(sc[1,j,i] < 0)
              sc[,j,i] <- -sc[,j,i]
      }
  }

  ## Prepare objects
  rownames(sc) <- rownames(tab)
  colnames(sc) <- colnames(layer) <- paste("Dim", 1:nd, sep="")
  rownames(layer) <- dimnames(tab)[[3]]

  if(!homogeneous)
      dimnames(sc)[[3]] <- dimnames(tab)[[3]]
  else
      dimnames(sc)[[3]] <- "All levels"


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

  obj <- list(phi = layer, row = sc, col = sc, diagonal = dg,
              weighting = weighting, row.weights = p, col.weights = p)

  class(obj) <- c("assoc.rcL", "assoc.symm", "assoc")
  obj
}
