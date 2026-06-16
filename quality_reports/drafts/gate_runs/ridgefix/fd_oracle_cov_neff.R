#!/usr/bin/env Rscript
# FD-ORACLE (cov-path lift): the ridge total differential C:dOmega + (tr(C)/n_eff) tr(dOmega)
# with an EXPLICIT n_eff != n (i.e. testing that swapping n -> n_eff in the ridge denominator
# keeps the analytic EE channel EXACT, since lambda(Omega) = tr(Omega)/n_eff and n_eff is a
# fixed scalar). The package wires .ridge_neff = n_eff_edid(unit_weights, all-TRUE, n); on the
# cov path unit_weights is NULL today so n_eff == n -- here we FORCE a non-trivial n_eff to
# confirm the differential is exact at the value n_eff actually takes, at tol 1e-8.
suppressMessages(pkgload::load_all("/Users/pcostag/Documents/GitHub/did", quiet=TRUE))
options(digits = 12)
set.seed(7)

# exact smooth inverse-variance weight functional and its coupling C = dtheta/dM
th <- function(M, mbar) { Minv <- solve(0.5*(M+t(M))); w <- drop(Minv %*% rep(1,nrow(M))); w <- w/sum(w); sum(w*mbar) }
Csmooth <- function(M, mbar) { Minv <- solve(0.5*(M+t(M))); w <- drop(Minv %*% rep(1,nrow(M)))
  w <- w/sum(w); theta <- sum(w*mbar); q <- drop(Minv %*% (mbar - theta)); -0.5*(outer(q,w)+outer(w,q)) }
mk_pd <- function(H) { S <- 0.3*crossprod(matrix(rnorm(H*H),H))/H + diag(seq(1.5,1.5+H-1)); 0.5*(S+t(S)) }

cat("=== FD-ORACLE cov-path ridge total differential at n_eff (tol 1e-8) ===\n\n")
worst <- 0
for (n_full in c(40L, 120L, 300L)) {
 for (neff_factor in c(1.0, 0.42, 0.08)) {          # n_eff = neff_factor * n_full (1.0 = unweighted)
  n_eff <- neff_factor * n_full
  ridge <- function(M) M + (sum(diag(M))/n_eff) * diag(nrow(M))    # lambda = tr(M)/n_eff
  for (H in c(3L, 5L)) {
    M0 <- mk_pd(H); mbar <- rnorm(H); Mr <- ridge(M0)
    C <- Csmooth(Mr, mbar); trC <- sum(diag(C))
    dE <- matrix(rnorm(H*H),H); dE <- 0.5*(dE+t(dE))
    ana <- sum(C*dE) + (trC/n_eff) * sum(diag(dE))                 # C:dOmega + (trC/n_eff) tr dOmega
    h <- 1e-6
    fd <- (th(ridge(M0 + h*dE), mbar) - th(ridge(M0 - h*dE), mbar)) / (2*h)
    rel <- abs(ana - fd) / max(abs(fd), 1e-300)
    worst <- max(worst, rel)
    cat(sprintf("  n=%3d  n_eff=%6.1f (n/n_eff=%5.2f)  H=%d | ana=%+.6e fd=%+.6e  rel=%.2e %s\n",
                n_full, n_eff, n_full/n_eff, H, ana, fd, rel, if (rel < 1e-8) "OK" else "FAIL"))
  }
 }
}
cat(sprintf("\nworst rel = %.3e  -> %s\n", worst, if (worst < 1e-8) "PASS" else "FAIL"))
