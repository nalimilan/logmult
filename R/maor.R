
maor <- function(tab, phi=FALSE,
                 weighting=c("marginal", "uniform", "none"), norm=2,
                 row.weights=NULL, col.weights=NULL) {
  weighting <- match.arg(weighting)

  if(!length(dim(tab)) %in% 2:3) {
      stop("Only two- and three-way tables are supported.")
  }

  if(!norm %in% 1:2)
     stop("'norm' must be 1 or 2")

  if(norm != 2 && !phi)
      stop("Only norm=2 is currently supported when phi=FALSE")

  if(!(is.null(row.weights) || !is.null(col.weights)) && !missing(weighting))
      warning("Argument 'weighting' is ignored when custom row/column weights are specified.")

  if(any(is.na(tab)))
      stop("NA cells are currently not supported.")

  if(length(dim(tab)) == 3) {
      rp <- margin.table(tab, 1)
      cp <- margin.table(tab, 2)

      if(any(tab == 0)) {
          tab <- tab + 0.5
          warning("Cells with zero counts found: adding 0.5 to each cell of the table.")
      }

      if(weighting == "marginal")
          return(apply(tab, 3, maor,
                       phi=phi, norm=norm,
                       row.weights=rp, col.weights=cp))
      else
          return(apply(tab, 3, maor,
                       phi=phi, weighting=weighting, norm=norm,
                       row.weights=row.weights, col.weights=col.weights))
  }

  if(any(tab == 0)) {
      tab <- tab + 0.5
      warning("Cells with zero counts found: adding 0.5 to each cell of the table.")
  }

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
      rp <- prop.table(row.weights)

  if(!is.null(col.weights))
      cp <- prop.table(col.weights)



  rp1 <- prop.table(rp)
  cp1 <- prop.table(cp)
  logp <- log(p)

  lambda <- sweep(logp, 1, rowSums(sweep(logp, 2, cp1, "*")), "-")
  lambda <- sweep(lambda, 2, colSums(sweep(logp, 1, rp1, "*")), "-")
  lambda <- lambda + sum(logp * rp1 %o% cp1)

  phi.norm <- sum(abs(lambda)^norm * rp %o% cp)

  if(phi)
      phi.norm^(1/norm)
  else
      exp((4/sum((rp1 * (1 - rp1)) %o% (cp1 * (1 - cp1))) * phi.norm)^(1/norm))
}