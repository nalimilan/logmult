# See Becker-Clogg (1989) test

library(logmult)
data(color)

rcm <- rc(color[,,1], 2, weighting="marginal", start=NA)
rcu <- rc(color[,,1], 2, weighting="uniform", start=NA)
rcn <- rc(color[,,1], 2, weighting="none", start=NA)

phim <- maor(fitted(rcm), TRUE, "marginal", norm=2)
phiu <- maor(fitted(rcu), TRUE, "uniform", norm=2)
phin <- maor(fitted(rcn), TRUE, "none", norm=2)

stopifnot(all.equal(phim, sqrt(sum((rcm$assoc$phi)^2))))
stopifnot(all.equal(phiu, sqrt(sum(abs(rcu$assoc$phi)^2))))
stopifnot(all.equal(phin, sqrt(sum(abs(rcn$assoc$phi)^2))))



# Test for phi computed from UNIDIFF two-way interaction coefficients
data(yaish)
tab <- aperm(yaish[,,-7], 3:1)

# 1-norm
# Currently disabled until supported correctly
# u1m <- unidiff(tab, weighting="marginal", norm=1)
# rp <- prop.table(margin.table(tab, 1))
# cp <- prop.table(margin.table(tab, 2))
# stopifnot(all.equal(u1m$unidiff$phi,
#                     maor(fitted(u1m)[,,1], TRUE, "marginal", norm=1, rp, cp)))
# stopifnot(all.equal(u1m$unidiff$phi,
#                     maor(fitted(u1m)[,,1], TRUE, "marginal", norm=1, rp, cp)))
# 
# u1u <- unidiff(tab, weighting="uniform", norm=1)
# stopifnot(all.equal(u1u$unidiff$phi,
#                     maor(fitted(u1u)[,,1], TRUE, "uniform", norm=1)))
# stopifnot(all.equal(u1u$unidiff$phi * exp(u1u$unidiff$layer$qvframe$estimate[2]),
#                     maor(fitted(u1u)[,,2], TRUE, "uniform", norm=1)))
# 
# u1n <- unidiff(tab, weighting="none", norm=1)
# stopifnot(all.equal(u1n$unidiff$phi,
#                     maor(fitted(u1n)[,,1], TRUE, "none", norm=1)))
# stopifnot(all.equal(u1n$unidiff$phi * exp(u1u$unidiff$layer$qvframe$estimate[2]),
#                     maor(fitted(u1n)[,,2], TRUE, "none", norm=1)))


# 2-norm
u2m <- unidiff(tab, weighting="marginal", norm=2)
rp <- prop.table(margin.table(tab, 1))
cp <- prop.table(margin.table(tab, 2))
stopifnot(all.equal(u2m$unidiff$phi,
                    maor(fitted(u2m)[,,1], TRUE, "marginal", norm=2, rp, cp)))
stopifnot(all.equal(u2m$unidiff$phi * exp(u2m$unidiff$layer$qvframe$estimate[2]),
                    maor(fitted(u2m)[,,2], TRUE, "marginal", norm=2, rp, cp)))

u2u <- unidiff(tab, weighting="uniform", norm=2)
stopifnot(all.equal(u2u$unidiff$phi,
                    maor(fitted(u2u)[,,1], TRUE, "uniform", norm=2)))
stopifnot(all.equal(u2u$unidiff$phi * exp(u2u$unidiff$layer$qvframe$estimate[2]),
                    maor(fitted(u2u)[,,2], TRUE, "uniform", norm=2)))

u2n <- unidiff(tab, weighting="none", norm=2)
stopifnot(all.equal(u2n$unidiff$phi,
                    maor(fitted(u2n)[,,1], TRUE, "none", norm=2)))
stopifnot(all.equal(u2n$unidiff$phi * exp(u2n$unidiff$layer$qvframe$estimate[2]),
                    maor(fitted(u2n)[,,2], TRUE, "none", norm=2)))


###
## Comparison with mean/sum of all spanning odds ratios
###
# Can be set to arbitrary values
nr <- 4
nc <- 5

norm <- 2

or <- function(tab) {
    or1 <- function(i, j, i2, j2) (tab[i, j] * tab[i2, j2]) / (tab[i, j2] * tab[i2, j])
    or <- array(NA, c(nrow(tab), ncol(tab), nrow(tab), ncol(tab)))

    for(i in 1:nrow(tab))
      for(j in 1:ncol(tab))
        for(i2 in 1:nrow(tab))
          for(j2 in 1:ncol(tab))
            if(i2 != i && j2 != j)
              or[i, j, i2, j2] <- or1(i, j, i2, j2)
    or
}

wlor2 <- function(tab) {
    rp <- prop.table(margin.table(tab, 1)) * nrow(tab)
    cp <- prop.table(margin.table(tab, 2)) * ncol(tab)
    wlor2 <- w <- array(NA, c(nrow(tab), ncol(tab), nrow(tab), ncol(tab)))

    for(i in 1:nrow(tab))
      for(j in 1:ncol(tab))
        for(i2 in 1:nrow(tab))
          for(j2 in 1:ncol(tab))
            if(i2 != i && j2 != j) {
              wlor2[i, j, i2, j2] <- log((tab[i, j] * tab[i2, j2]) / (tab[i, j2] * tab[i2, j]))^norm
              w[i, j, i2, j2] <- rp[i] * cp[j] * rp[i2] * cp[j2]
            }

    wlor2 * w / sum(w, na.rm=TRUE)
}

## Unweighted mean
# General case
res <- replicate(10, {
    tab <- matrix(rpois(nr*nc, 1000), nr, nc) + .5
    rp <- rep(1/nr, nr)
    cp <- rep(1/nc, nc)
    c(maor(tab, weighting="uniform", norm=norm),
      maor(tab, norm=norm, row.weights=rp, col.weights=cp),
      exp(mean(abs(log(or(tab)))^norm, na.rm=TRUE)^(1/norm)))
})

stopifnot(all.equal(res[1,], res[2,]))
stopifnot(all.equal(res[1,], res[3,]))

# 2x2 table (equality with single odds ratio)
res <- replicate(10, {
    tab <- matrix(rpois(2*2, 1000), 2, 2) + .5
    rp <- cp <- rep(1/2, 2)
    c(maor(tab, weighting="uniform", norm=norm),
      maor(tab, norm=norm, row.weights=rp, col.weights=cp),
      exp(sqrt(log(tab[1,1] * tab[2,2] / (tab[1,2] * tab[2,1]))^2)))
})

stopifnot(all.equal(res[1,], res[2,]))
stopifnot(all.equal(res[1,], res[3,]))


## Unweighted sum
# General case
res <- replicate(10, {
    tab <- matrix(rpois(nr*nc, 1000), nr, nc) + .5
    rp <- rep(1, nr)
    cp <- rep(1, nc)
    c(maor(tab, weighting="none", norm=norm),
      exp((sum(abs(log(or(tab)))^norm, na.rm=TRUE)/((nr-1) * (nc-1)))^(1/norm)))
})

stopifnot(all.equal(res[1,], res[2,]))

# 2x2 table (equality with single odds ratio)
res <- replicate(10, {
    tab <- matrix(rpois(2*2, 1000), 2, 2) + .5
    rp <- cp <- rep(1, 2)
    c(maor(tab, weighting="none", norm=norm),
      exp(sqrt(4 * log(tab[1,1] * tab[2,2] / (tab[1,2] * tab[2,1]))^2)))
})

stopifnot(all.equal(res[1,], res[2,]))


## Marginal-weighted mean
# General case
res <- replicate(10, {
    tab <- matrix(rpois(nr*nc, 100), nr, nc) +.5
    rp <- prop.table(margin.table(tab, 1))
    cp <- prop.table(margin.table(tab, 2))
    c(maor(tab, weighting="marginal", norm=norm),
      maor(tab, norm=norm, row.weights=rp, col.weights=cp),
      exp(abs(sum(wlor2(tab), na.rm=TRUE))^(1/norm)))
})

stopifnot(all.equal(res[1,], res[2,]))
stopifnot(all.equal(res[1,], res[3,]))

# 2x2 table (equality with single odds ratio)
res <- replicate(10, {
    tab <- matrix(rpois(2*2, 100), 2, 2) + .5
    rp <- prop.table(margin.table(tab, 1))
    cp <- prop.table(margin.table(tab, 2))
    c(maor(tab, weighting="marginal", norm=norm),
      maor(tab, norm=norm, row.weights=rp, col.weights=cp),
      exp(abs(log(tab[1,1] * tab[2,2] / (tab[1,2] * tab[2,1])))))
})

stopifnot(all.equal(res[1,], res[2,]))
stopifnot(all.equal(res[1,], res[3,]))