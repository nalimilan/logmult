lambda <- function(tab, rp=rep(1/nrow(tab), nrow(tab)), cp=rep(1/ncol(tab), ncol(tab))) {
  logp <- log(prop.table(tab))

  lambda <- sweep(logp, 1, rowSums(sweep(logp, 2, cp, "*")), "-")
  lambda <- sweep(lambda, 2, colSums(sweep(logp, 1, rp, "*")), "-")
  lambda + sum(logp * rp %o% cp)
}

iac <- function(tab, cell=FALSE,
                weighting=c("marginal", "uniform", "none"),
                component=c("total", "symmetric", "antisymmetric"),
                shrink=FALSE,
                row.weights=NULL, col.weights=NULL) {
  weighting <- match.arg(weighting)
  component <- match.arg(component)

  if(!length(dim(tab)) %in% 2:3) {
      stop("Only two- and three-way tables are supported.")
  }

  if(!(is.null(row.weights) || !is.null(col.weights)) && !missing(weighting))
      warning("Argument 'weighting' is ignored when custom row/column weights are specified.")

  if(any(is.na(tab))) {
      warning("NA cells are currently not supported, returning NA.")
      return(NA)
  }

  if(any(tab == 0)) {
      tab <- tab + 0.5
      warning("Cells with zero counts found: adding 0.5 to each cell of the table.")
  }

  if(weighting == "marginal") {
      p <- prop.table(tab)
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
      rp <- prop.table(row.weights)

  if(!is.null(col.weights))
      cp <- prop.table(col.weights)

  if(shrink) {
      if(length(dim(tab)) != 3)
          stop("shrink=TRUE only makes sense for three-way tables")

      if(component != "total")
          stop("shrink=TRUE is currently only supported with component=\"total\"")

      return(iac_shrunk(tab, rp, cp))
  }
  else if(length(dim(tab)) == 3) {
      rp <- margin.table(tab, 1)
      cp <- margin.table(tab, 2)

      if(!is.null(row.weights))
          rp <- row.weights

      if(!is.null(col.weights))
          cp <- col.weights

      if(weighting == "marginal")
          res <- apply(tab, 3, iac, cell=cell, component=component,
                       row.weights=rp, col.weights=cp)
      else
          res <- apply(tab, 3, iac, cell=cell, weighting=weighting, component=component,
                       row.weights=row.weights, col.weights=col.weights)

      if(cell) {
          dim(res) <- dim(tab)
          dimnames(res) <- dimnames(tab)
      }

      return(res)
  }

  if(any(tab == 0)) {
      if(all(tab == 0)) {
          warning("Table contains only empty cells, returning NA.")
          return(NA)
      }
      else if(sum(tab == 0)/length(tab) > .25) {
          warning("More than 25% of cells are empty, the value of the index may not be reliable.")
      }
  }

  if(component %in% c("symmetric", "antisymmetric"))
      rp <- cp <- (rp + cp)/2

  rp1 <- prop.table(rp)
  cp1 <- prop.table(cp)

  l <- lambda(tab, rp1, cp1)

  if(component == "symmetric")
      l <- (l + t(l))/2
  else if(component == "antisymmetric")
      l <- (l - t(l))/2

  lambdasq <- l^2 * rp %o% cp

  if(cell)
      lambdasq
  else
      sqrt(sum(lambdasq))
}

# Data frame of all log-odds ratios that can be computed from tab,
# and of their variances and weights
lor <- function(tab, rp, cp) {
    ltab <- log(tab)
    res <- data.frame(LOR=rep(NA, length(tab)^2), V=NA, W=NA)
    r <- 0
    for(i in 1:nrow(tab)) {
        for(j in 1:ncol(tab)) {
            for(k in 1:nrow(tab)) {
                for(l in 1:ncol(tab)) {
                    LOR <- ltab[i, j] + ltab[k, l] - ltab[i, l] - ltab[k, j]
                    V <- 1/tab[i, j] + 1/tab[k, l] + 1/tab[i, l] + 1/tab[k, j] # Zhou (2015, p. 324)
                    W <- rp[i]*cp[j]*rp[k]*cp[l] # extension, to weight LORs
                    r <- r + 1
                    res[r,] <- c(LOR, V, W)
                }
            }
        }
    }
    res
}

# Function adapted from code provided by Zhou (2015, section 3)
# y contains the observed log-odds ratios, and sigma.sq their variances
shrink <- function(y, sigma.sq) {
  k <- length(y)
  # solve tau.sq, weight, and mu
  tau.sq <- 0
  weight <- 1/(sigma.sq+tau.sq)
  mu <- weighted.mean(y, weight)
  for(i in 1:100) {
      sum.sq <- sum(weight*(y-mu)^2)
      sum.sq.deriv <- -sum((weight^2)*(y-mu)^2)
      tau.sq.update <- max(0, tau.sq + (k-1-sum.sq)/sum.sq.deriv)
      weight.update <- 1/(sigma.sq+tau.sq.update)
      mu.update <- weighted.mean(y, weight.update)

      if(is.finite(abs(tau.sq.update-tau.sq)) && abs(tau.sq.update-tau.sq) < 1.0e-8)
          break()

      tau.sq <- tau.sq.update
      weight <- weight.update
      mu <- mu.update
  }
  if(i == 100)
      stop("algorithm did not converge")

  # point estimates for theta
  shrink <- (k-3)/(k-1)*sigma.sq*weight
  theta <- mu+(1-shrink)*(y-mu)
  theta
}

iac_shrunk <- function(tab, rp, cp) {
  lors <- array(NA, dim=c(nrow(tab)^2*ncol(tab)^2, 3, dim(tab)[[3]]))
  colnames(lors) <- c("LOR", "V", "W")
  dimnames(lors)[[3]] <- dimnames(tab)[[3]]

  for(i in seq(dim(tab)[3])) {
      x <- lor(tab[,,i], rp, cp)
      lors[, "LOR", i] <- x$LOR
      lors[, "V", i] <- x$V
      lors[, "W", i] <- x$W
  }

  lors_shrunk <- apply(lors, 1, function(x) shrink(x["LOR",], x["V",]))
  res <- sqrt(rowSums(lors_shrunk^2 * t(lors[,"W",])))
  names(res) <- dimnames(lors)[[3]]
  res/sqrt(4 * sum(rp %o% cp))
}
