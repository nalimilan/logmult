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
# For this, use two tabs with identical association but different margins,
# and stack them to fit the UNIDIFF model
a <- matrix(1:9, 3, 3)
b <- a
b[,1] <- b[,1] * 5
b[1,] <- b[1,] * 2

tab <- array(c(a, b), c(3, 3, 2))

# 1-norm
u1m <- unidiff(tab, weighting="marginal", norm=1)
stopifnot(all.equal(u1m$unidiff$phi,
                    maor(a, TRUE, "marginal", norm=1,
                         prop.table(margin.table(tab, 1)),
                         prop.table(margin.table(tab, 2)))))
stopifnot(all.equal(u1m$unidiff$phi,
                    maor(b, TRUE, "marginal", norm=1,
                         prop.table(margin.table(tab, 1)),
                         prop.table(margin.table(tab, 2)))))

u1u <- unidiff(tab, weighting="uniform", norm=1)
stopifnot(all.equal(u1u$unidiff$phi, maor(a, TRUE, "uniform", norm=1)))
stopifnot(all.equal(u1u$unidiff$phi, maor(b, TRUE, "uniform", norm=1)))

u1n <- unidiff(tab, weighting="none", norm=1)
stopifnot(all.equal(u1n$unidiff$phi, maor(a, TRUE, "none", norm=1)))
stopifnot(all.equal(u1n$unidiff$phi, maor(b, TRUE, "none", norm=1)))


# 2-norm
u2m <- unidiff(tab, weighting="marginal", norm=2)
stopifnot(all.equal(u2m$unidiff$phi,
                    maor(a, TRUE, "marginal", norm=2,
                         prop.table(margin.table(tab, 1)),
                         prop.table(margin.table(tab, 2)))))
stopifnot(all.equal(u2m$unidiff$phi,
                    maor(b, TRUE, "marginal", norm=2,
                         prop.table(margin.table(tab, 1)),
                         prop.table(margin.table(tab, 2)))))

u2u <- unidiff(tab, weighting="uniform", norm=2)
stopifnot(all.equal(u2u$unidiff$phi, maor(a, TRUE, "uniform", norm=2)))
stopifnot(all.equal(u2u$unidiff$phi, maor(b, TRUE, "uniform", norm=2)))

u2n <- unidiff(tab, weighting="none", norm=2)
stopifnot(all.equal(u2n$unidiff$phi, maor(a, TRUE, "none", norm=2)))
stopifnot(all.equal(u2n$unidiff$phi, maor(b, TRUE, "none", norm=2)))