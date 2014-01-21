# See Becker-Clogg (1989) test

library(logmult)
data(color)

rcm <- rc(color[,,1], 2, weighting="marginal")
rcu <- rc(color[,,1], 2, weighting="uniform")
rcn <- rc(color[,,1], 2, weighting="none")

phim <- goodman.phi(fitted(rcm), "marginal")
phiu <- goodman.phi(fitted(rcu), "uniform")
phin <- goodman.phi(fitted(rcn), "none")

stopifnot(all.equal(phim, sqrt(sum(rcm$assoc$phi^2))))
stopifnot(all.equal(phiu, sqrt(sum(rcu$assoc$phi^2))))
stopifnot(all.equal(phin, sqrt(sum(rcn$assoc$phi^2))))


# Test for phi computed from UNIDIFF two-way interaction coefficients
# For this, use two tabs with identical association but different margins,
# and stack them to fit the UNIDIFF model
a <- matrix(1:9, 3, 3)
b <- a
b[,1] <- b[,1] * 5
b[1,] <- b[1,] * 2

tab <- array(c(a, b), c(3, 3, 2))

um <- unidiff(tab, weighting="marginal")
stopifnot(all.equal(um$unidiff$phi,
                    goodman.phi(a, "marginal",
                                prop.table(margin.table(tab, 1)),
                                prop.table(margin.table(tab, 2)))))
stopifnot(all.equal(um$unidiff$phi,
                    goodman.phi(b, "marginal",
                                prop.table(margin.table(tab, 1)),
                                prop.table(margin.table(tab, 2)))))

uu <- unidiff(tab, weighting="uniform")
stopifnot(all.equal(uu$unidiff$phi, goodman.phi(a, "uniform")))
stopifnot(all.equal(uu$unidiff$phi, goodman.phi(b, "uniform")))

un <- unidiff(tab, weighting="none")
stopifnot(all.equal(un$unidiff$phi, goodman.phi(a, "none")))
stopifnot(all.equal(un$unidiff$phi, goodman.phi(b, "none")))