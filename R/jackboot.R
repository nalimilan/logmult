# Run jackknife or bootstrap replicates of the model
jackboot <- function(se, ncpus, nreplicates, tab, model, assoc1, assoc2,
                     weighting, family, weights, base, ...) {
  cat("Computing", se, "standard errors...\n")

  if(is.null(ncpus))
      ncpus <- if(require(parallel)) min(parallel::detectCores(), 5) else 1

  if(se == "jackknife") {
      jack <- jackknife((1:length(tab))[!is.na(tab)], jackknife.assoc,
                        w=tab[!is.na(tab)], ncpus=ncpus,
                        model=model, assoc1=assoc1, assoc2=assoc2,
                        weighting=weighting, family=family, weights=weights, ...,
                        base=base, verbose=FALSE)

      covmat <- jack$vcov
      jack.results <- list(bias=jack$bias, values=jack$values)
      boot.results <- numeric(0)
  }
  else {
      boot.results <- boot::boot(1:sum(tab, na.rm=TRUE), boot.assoc,
                                 R=nreplicates, ncpus=ncpus, parallel="snow",
                                 args=list(model=model, assoc1=assoc1, assoc2=assoc2,
                                           weighting=weighting, family=family,
                                           weights=weights, ..., base=base))

      covmat <- cov(boot.results$t, use="na.or.complete")
      jack.results <- numeric(0)
  }


  # When gnm evaluates the formulas, tab will have been converted to a data.frame,
  # with a fallback if both names are empty
  vars <- make.names(names(dimnames(tab)))
  if(length(vars) == 0)
      vars <- c("Var1", "Var2")

  mat.names <- function(ass) {
      if(length(rownames(ass$phi)) > 0)
          lnames <- paste(vars[3], rownames(ass$phi), ":", sep="")
      else
          lnames <- ""

      nd <- ncol(ass$phi)

      scnames <- t(outer(lnames,
                         c(t(outer(paste("Dim", 1:nd, sep=""),
                                   paste(vars[1], rownames(ass$row), sep=""),
                                   paste, sep=":")),
                           t(outer(paste("Dim", 1:nd, sep=""),
                                   paste(vars[2], rownames(ass$col), sep=""),
                                   paste, sep=":"))),
                         paste, sep=""))

      c(t(outer(lnames, paste("Dim", 1:nd, sep=""), paste, sep="")),
                scnames,
                paste(scnames, "*", sep=""))
  }

  ass1 <- assoc1(model, weighting=weighting)

  nl <- nrow(ass1$phi)
  nd <- ncol(ass1$phi)
  nr <- nrow(ass1$row)
  nc <- nrow(ass1$col)

  int <- seq.int(1, nl * nd + 2 * nl * nd * (nr + nc))
  covmat1 <- covmat[int, int]

  rownames(covmat1) <- colnames(covmat1) <- mat.names(ass1)

  if(length(boot.results) > 0) {
      boot.results1 <- boot.results
      boot.results1$t0 <- boot.results$t0[int]
      boot.results1$t <- boot.results$t[,int]
      names(boot.results1$t0) <- colnames(boot.results1$t) <- rownames(covmat1)
  }
  else {
      boot.results1 <- numeric(0)
  }

  if(length(jack.results) > 0) {
      jack.results1 <- jack.results
      jack.results1$bias <- jack.results$bias[int]
      jack.results1$values <- jack.results$values[,int]
      names(jack.results1$bias) <- colnames(jack.results1$values) <- rownames(covmat1)
  }
  else {
      jack.results1 <- numeric(0)
  }

  if(!is.null(assoc2)) {
      covmat2 <- covmat[-int, -int]
      rownames(covmat2) <- colnames(covmat2) <- mat.names(assoc2(model, weighting=weighting))

      if(length(boot.results) > 0) {
          boot.results2 <- boot.results
          boot.results2$t0 <- boot.results$t0[-int]
          boot.results2$t <- boot.results$t[,-int]
          names(boot.results2$t0) <- colnames(boot.results2$t) <- rownames(covmat2)
      }
      else {
          boot.results2 <- numeric(0)
      }

      if(length(jack.results) > 0) {
          jack.results2 <- jack.results
          jack.results2$bias <- jack.results$bias[-int]
          jack.results2$values <- jack.results$values[,-int]
          names(jack.results2$bias) <- colnames(jack.results2$values) <- rownames(covmat2)
      }
      else {
          jack.results2 <- numeric(0)
      }
  }
  else {
      covmat2 <- numeric(0)
      boot.results2 <- numeric(0)
      jack.results2 <- numeric(0)
  }


  if(is.null(assoc2))
      list(covmat=covmat1, boot.results=boot.results1, jack.results=jack.results1)
  else
      list(covmat1=covmat1, boot.results1=boot.results1, jack.results1=jack.results1,
           covmat2=covmat2, boot.results2=boot.results2, jack.results2=jack.results2)
}

# Additional arguments are needed so that update() finds them even when using parLapply
jackknife.assoc <- function(x, model, repl.verbose=FALSE, ...) {
  tab <- model$data

  if(repl.verbose) {
      iter <- which(!1:length(tab) %in% x)
      if(length(iter) == 1)
          cat(sprintf("Iteration for cell %i of %i\n", iter, length(tab)))
      else
          cat("Initial iteration\n")
  }

  if(sum(tab[-x], na.rm=TRUE) > 0) {
      mat <- tab
      mat[] <- -1
      mat[x] <- 0
      tab <- tab + mat
  }

  replicate.assoc(model, tab, repl.verbose=repl.verbose, ...)
}

boot.assoc <- function(data, indices, args) {
  tab <- args$model$data

  # Create a table from the indices - one index identifies an observation in the original table,
  # following the cumulative sum, from 1 to sum(tab)
  tab[!is.na(tab) & tab > 0] <- tapply(tabulate(indices, nbins=sum(tab, na.rm=TRUE)),
                                      rep.int(1:length(tab), ifelse(is.na(tab), 0, tab)),
                                      sum)

  # Basic sanity check
  stopifnot(sum(tab, na.rm=TRUE) == sum(args$model$data, na.rm=TRUE))

  # We need to pass all arguments through "args" to prevent them
  # from being catched by boot(), especially "weights"
  args$tab <- tab
  do.call(replicate.assoc, args)
}

# Replicate model with new data, and combine assoc components into a vector
replicate.assoc <- function(model.orig, tab, assoc1, assoc2, weighting, ...,
                            base=NULL, repl.verbose=FALSE) {
  library(assoc)

  # Models can generate an error if they fail repeatedly
  # Remove warnings because we handle them below
  model <- tryCatch(suppressWarnings(update(model.orig, tab=tab,
                                            start=parameters(model.orig),
                                            etastart=as.numeric(predict(model.orig)),
                                            verbose=repl.verbose, trace=repl.verbose, se="none")),
                    error=function(e) NULL)

  if(is.null(model) || !model$converged) {
      if(is.null(model)) {
          cat("Model replicate failed.\nData was:\n")
          model <- model.orig
      }
      else {
          cat("Model replicate did not converge.\nData was:\n")
      }


      print(tab)
      cat(sprintf("Trying again with different starting values...\n"))

      model <- tryCatch(suppressWarnings(update(model, tab=tab, start=NA, etastart=NULL,
                                                verbose=TRUE, trace=TRUE, se="none")),
                        error=function(e) NULL)
  }

  if(is.null(model) || !model$converged) {
      cat(sprintf("Model failed again. Trying one last time with random starting values...\n"))

      # Without the quote(NULL), update.gnm() does call$start <- NULL, which removes it,
      # and eventually restores the default value (NA)
      model <- tryCatch(suppressWarnings(update(model, tab=tab, start=quote(NULL), etastart=NULL,
                                                verbose=TRUE, trace=TRUE, se="none")),
                        error=function(e) NULL)
  }

  if(is.null(model) || !model$converged) {
      warning("Model failed to converge three times: ignoring the results of this replicate. Standard errors may not be completely accurate. Consider raising the value of iterMax.", immediate.=TRUE)

      # This NA value is skipped when computing variance-covariance matrix
      return(NA)
  }

  ass1 <- assoc1(model, weighting=weighting)
  ass1.orig <- assoc1(model.orig, weighting=weighting)

  ret <- if(inherits(ass1, c("assoc.hmskew", "assoc.hmskewL"))) find.stable.scores.hmskew(ass1, ass1.orig)
         else find.stable.scores(ass1, ass1.orig)

  # For double association models like some hmskew and yrcskew variants
  if(!is.null(assoc2)) {
      ass2 <- assoc2(model, weighting=weighting)
      ass2.orig <- assoc2(model.orig, weighting=weighting)

      ret <- c(ret, if(inherits(ass1, c("assoc.hmskew", "assoc.hmskewL"))) find.stable.scores.hmskew(ass1, ass1.orig)
                    else find.stable.scores(ass1, ass1.orig))
  }

  ret
}

# Originally based on procrustes() from package vegan 2.0-4, by
# Jari Oksanen, F. Guillaume Blanchet, Roeland Kindt, Pierre Legendre,
# Peter R. Minchin, R. B. O'Hara, Gavin L. Simpson, Peter Solymos, M. Henry H. Stevens, Helene Wagner.
# License GPL-2.
# In this very simplified version, we perform no centering nor scaling (handled manually with weighting)
procrustes <- function (X, Y) {
  XY <- crossprod(X, Y)
  sol <- svd(XY)

  A <- sol$v %*% t(sol$u)
  Yrot <- Y %*% A

  list(Yrot = Yrot, rotation = A, svd=sol)
}

# Compute distance between adjusted scores for this replicate to that of the original model
# We choose the permutation and sign of dimensions that minimizes the sum of squares,
# weighted by the inverse of row frequencies
find.stable.scores <- function(ass, ass.orig) {
  weights <- 1

  nd <- ncol(ass$phi)
  nr <- nrow(ass$row)
  nc <- nrow(ass$col)
  nl <- nrow(ass$phi)

  sc <- adj.orig <- array(NA, dim=c(nr + nc, nd, nl))
  phi <- ass$phi

  sc[1:nr,,] <- ass$row
  sc[-(1:nr),,] <- ass$col

  adj.orig[1:nr,,] <- ass.orig$row
  adj.orig[-(1:nr),,] <- ass.orig$col

  adj <- sc
  adj <- sweep(adj, 3:2, sqrt(abs(ass$phi)), "*")
  adj.orig <- sweep(adj.orig, 3:2, sqrt(abs(ass.orig$phi)), "*")

  # If phi is negative, change sign of columns for adjusted scores so that the interpretation
  # is consistent with positive phi for plotting (plot.assoc() does the same)
  # For symmetric models, we change the sign of adjusted column scores since this makes computations easier
  adj[-(1:nr),,] <- sweep(adj[-(1:nr),, , drop=FALSE], 3:2, sign(ass$phi), "*")
  adj.orig[-(1:nr),,] <- sweep(adj.orig[-(1:nr),, , drop=FALSE], 3:2, sign(ass$phi), "*")

  # Transform arrays to matrices with one column per dimension
  adj <- aperm(adj, c(1, 3, 2))
  adj.orig <- aperm(adj.orig, c(1, 3, 2))
  dim(adj) <- dim(adj.orig) <- c((nr + nc) * nl, nd)
  procr <- procrustes(adj.orig, adj)

  adj <- procr$Yrot
  dim(adj) <- c(nr + nc, nl, nd)
  adj <- aperm(adj, c(1, 3, 2))

  for(l in 1:nl) {
      # phi for rows and column are identical since rotation is the same for both,
      # so take an average in case there are rounding errors
      phi[l,] <- margin.table(sweep(adj[,,l]^2, 1,
                                    c(ass.orig$row.weights,
                                      ass.orig$col.weights), "*"), 2)/2 * sign(ass$phi[l,])

      sc[1:nr,,l] <- sweep(adj[1:nr,,l, drop=FALSE], 2, sqrt(abs(phi[l,])), "/")
      sc[-(1:nr),,l] <- sweep(adj[-(1:nr),,l, drop=FALSE], 2, sqrt(abs(phi[l,])) * sign(ass$phi[l,]), "/")
  }

  # Sanity check 1: rebuild the association matrix from normalized scores and compare with the original
  for(l in 1:nl) {
      lambda <- lambda.sav <- matrix(0, nr, nc)

      for(i in 1:nd) {
          # Heterogeneous scores
          if(dim(ass$row)[3] > 1) {
              lambda <- lambda + phi[l, i] * sc[1:nr, i, l] %o% sc[-(1:nr), i, l]
              lambda.sav <- lambda.sav + ass$phi[l, i] * ass$row[, i, l] %o% ass$col[, i, l]
          }
          # Homogeneous scores
          else {
              lambda <- lambda + phi[l, i] * sc[1:nr, i, l] %o% sc[-(1:nr), i, l]
              lambda.sav <- lambda.sav + ass$phi[l, i] * ass$row[, i, 1] %o% ass$col[, i, 1]
          }
      }

      stopifnot(isTRUE(all.equal(lambda, lambda.sav, check.attr=FALSE, tolerance=1e-8)))
  }

  # Sanity check 2: rebuild the association matrix from adjusted scores and compare with the original
  for(l in 1:nl) {
     lambda <- lambda.sav <- matrix(0, nr, nc)

      for(i in 1:nd) {
          # Heterogeneous scores
          if(dim(ass$row)[3] > 1) {
              lambda <- lambda + adj[1:nr, i, l] %o% adj[-(1:nr), i, l]
              lambda.sav <- lambda.sav + ass$phi[l, i] * ass$row[, i, l] %o% ass$col[, i, l]
          }
          # Homogeneous scores
          else {
              lambda <- lambda + adj[1:nr, i, l] %o% adj[-(1:nr), i, l]
              lambda.sav <- lambda.sav + ass$phi[l, i] * ass$row[, i, 1] %o% ass$col[, i, 1]
          }
      }

      stopifnot(isTRUE(all.equal(lambda, lambda.sav, check.attr=FALSE, tolerance=1e-8)))
  }

  row <- sc[1:nr,,, drop=FALSE]
  col <- sc[-(1:nr),,, drop=FALSE]
  adjrow <- adj[1:nr,,, drop=FALSE]
  adjcol <- adj[-(1:nr),,, drop=FALSE]

  ret <- numeric(nl * nd + 2 * nl * nd * (nr + nc))
  ret[seq(1, nd * nl)] <- t(phi)

  # We replicate normalized scores for all layers even if they are homogeneous for simplicity
  for(l in 1:nl) {
      int <- nl * nd + seq((l - 1) * nd * (nr + nc) + 1,
                           l * nd * (nr + nc))
      ret[int] <- c(row[,, l], col[,, l])
  }

  for(l in 1:nl) {
      int <- nl * nd + nl * nd * (nr + nc) + seq((l - 1) * nd * (nr + nc) + 1,
                                                 l * nd * (nr + nc))
      ret[int] <- c(adjrow[,, l], adjcol[,, l])
  }

  ret
}

# Simplified version of permutations() from the gtools 2.7.0 package,
# by Gregory R. Warnes.
# Original version by Bill Venables and cited by Matthew
# Wiener (mcw@ln.nimh.nih.gov) in an email to R-help dated
# Tue, 14 Dec 1999 09:11:32 -0500 (EST) in response to
# Alex Ahgarin <datamanagement@email.com>
perms <- function(n) {
  sub <- function(n, v) {
      if(n == 1) return(matrix(v, n, 1))

      X <- NULL

      for(i in 1:n)
        X <- rbind(X, cbind(v[i], Recall(n - 1, v[-i])))

      X
  }

  sub(n, 1:n)
}

# This class of models is special because dimensions are paired,
# and a rotation is needed rather than changing signs
find.stable.scores.hmskew <- function(ass, ass.orig) {
  weights <- 1

  nd <- ncol(ass$phi)
  nr <- nrow(ass$row)
  nc <- nrow(ass$col)
  nl <- nrow(ass$phi)

  phi <- ass$phi

  # Repeat scores once for each layer
  sc <- array(ass$row, dim=c(nr, nd, nl))
  adj <- array(ass$row, dim=c(nr, nd, nl))
  adj.orig <- array(ass.orig$row, dim=c(nr, nd, nl))

  # Compute adjusted scores
  adj <- sweep(adj, 3:2, sqrt(abs(ass$phi)), "*")
  # Where phi is negative, change signs on second dimension of each pair to get the same effect
  adj[,seq(1, nd, by=2),] <- sweep(adj[,seq(1, nd, by=2),, drop=FALSE], 3:2,
                                   sign(ass$phi)[,seq(1, nd, by=2)], "*")
  adj.orig[,seq(1, nd, by=2),] <- sweep(adj.orig[,seq(1, nd, by=2),, drop=FALSE], 3:2,
                                   sign(ass.orig$phi)[,seq(1, nd, by=2)], "*")

  # Rotate scores and return the sum of squares to the old scores
  rot <- function(angle, adj.tmp, dim) {
      rotmat <- matrix(c(cos(angle), sin(angle), -sin(angle), cos(angle)), 2, 2)
      adj.rot <- adj.tmp[, c(dim, dim + 1), , drop=FALSE]
      adj.rot[] <- apply(adj.rot, 3, "%*%", rotmat)
      sum(sweep((adj.rot - adj.orig[, c(dim,  dim + 1), , drop=FALSE])^2, 1, weights, "*"))
  }

  perms <- perms(nd/2)
  nperms <- nrow(perms)
  angles <- matrix(NA, nd/2, nperms)
  sq <- numeric(nperms)

  for(j in 1:nperms) {
      order <- rep(2 * perms[j,] - 1, each=2) + c(0, 1)
      adj.tmp <- adj[, order, , drop=FALSE]

      # Find optimal rotation for each pair of dimensions and apply it
      for(i in seq.int(1, nd, by=2)) {
          angle <- optim(0, rot, gr=NULL, adj.tmp, i, method="BFGS")$par
          angles[i, j] <- angle
          rotmat <- matrix(c(cos(angle), sin(angle), -sin(angle), cos(angle)), 2, 2)
          adj.tmp[, c(i, i + 1),] <- apply(adj.tmp[, c(i, i + 1), , drop=FALSE], 3, "%*%", rotmat)       
      }

      sq[i] <- sum(sweep((adj.orig - adj.tmp)^2, 1, weights, "*"))
  }

  # Find out which permutation allows for the smaller sum of squares,
  # choosing the best rotation for each pair of dimensions
  best.perm <- which.min(sq)
  order <- rep(2 * perms[best.perm,] - 1, 2) + c(0, 1)
  adj <- adj[, order, , drop=FALSE]
  # Phi does not need to be recomputed from adjusted scores since it does not change after a rotation
  phi <- phi[, order, drop=FALSE]

  for(i in seq.int(1, nd, by=2)) {
      best.angle <- angles[i, best.perm]
      rotmat <- matrix(c(cos(best.angle), sin(best.angle), -sin(best.angle), cos(best.angle)), 2, 2)

      sc[, c(i, i + 1),] <- apply(sc[, c(i, i + 1), , drop=FALSE], 3, "%*%", rotmat)
      adj[, c(i, i + 1),] <- apply(adj[, c(i, i + 1), , drop=FALSE], 3, "%*%",  rotmat)
  }

  # Sanity check 1: rebuild the association matrix from normalized scores and compare with the original
  for(l in 1:nl) {
     lambda <- lambda.sav <- matrix(0, nr, nc)

      for(i in seq.int(1, nd, by=2)) {
          # Heterogeneous scores
          if(dim(ass$row)[3] > 1) {
              # We introduce sign(phi) separately because adjusted scores cannot take it into account
              lambda <- lambda + phi[l, i] * (sc[, i + 1, l] %o% sc[, i, l] -
                                              sc[, i, l] %o% sc[, i + 1, l])
              lambda.sav <- lambda.sav + ass.sav$phi[l, i] * (ass.sav$row[, i + 1, l] %o% ass$row[, i, l] -
                                                              ass$row[, i, l] %o% ass$row[, i + 1, l])
          }
          # Homogeneous scores
          else {
              # We introduce sign(phi) separately because adjusted scores cannot take it into account
              lambda <- lambda + phi[l, i] * (sc[, i + 1, l] %o% sc[, i, l] -
                                              sc[, i, l] %o% sc[, i + 1, l])
              lambda.sav <- lambda.sav + ass$phi[l, i] * (ass$row[, i + 1, 1] %o% ass$row[, i, 1] -
                                                          ass$row[, i, 1] %o% ass$row[, i + 1, 1])
          }
      }

      stopifnot(isTRUE(all.equal(lambda, lambda.sav, check.attr=FALSE, tolerance=1e-8)))
  }

  # Sanity check 2: rebuild the association matrix from adjusted scores and compare with the original
  for(l in 1:nl) {
     lambda <- lambda.sav <- matrix(0, nr, nc)

      for(i in seq.int(1, nd, by=2)) {
          # Heterogeneous scores
          if(dim(ass$row)[3] > 1) {
              lambda <- lambda + adj[, i + 1, l] %o% adj[, i, l] -
                                 adj[, i, l] %o% adj[, i + 1, l]
              lambda.sav <- lambda.sav + ass$phi[l, i] * (ass$row[, i + 1, l] %o% ass$row[, i, l] -
                                                          ass$row[, i, l] %o% ass$row[, i + 1, l])
          }
          # Homogeneous scores
          else {
              lambda <- lambda + adj[, i + 1, l] %o% adj[, i, l] -
                                 adj[, i, l] %o% adj[, i + 1, l]
              lambda.sav <- lambda.sav + ass$phi[l, i] * (ass$row[, i + 1, 1] %o% ass$row[, i, 1] -
                                                          ass$row[, i, 1] %o% ass$row[, i + 1, 1])
          }
      }

      stopifnot(isTRUE(all.equal(lambda, lambda.sav, check.attr=FALSE, tolerance=1e-8)))
  }

  ret <- numeric(nl * nd + 2 * nl * nd * (nr + nc))
  ret[seq(1, nd * nl)] <- t(phi)

  # We replicate normalized scores for all layers even if they are homogeneous for simplicity
  for(l in 1:nl) {
      int <- nl * nd + seq((l - 1) * nd * (nr + nc) + 1,
                           l * nd * (nr + nc))
      ret[int] <- rep(sc[,, l], 2)
  }

  for(l in 1:nl) {
      int <- nl * nd + nl * nd * (nr + nc) + seq((l - 1) * nd * (nr + nc) + 1,
                                                 l * nd * (nr + nc))
      ret[int] <- rep(adj[,, l], 2)
  }

  ret
}

