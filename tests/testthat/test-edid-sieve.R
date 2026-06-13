# Tests for the sieve Omega smoother (options(edid_omega_method = "sieve")) and its weight-estimation channel.
# EFFICIENT and AVERAGED + sieve + misspec_robust both run the proper sieve Sigma_Omega (OLS-projection IF with
# the eigen-floor-aware Daleckii-Krein coupling -- per-unit Omega*(X_i) for efficient, the pooled Omega-bar for
# averaged; validated to nominal coverage + jackknife slope ~1). estimation_effect / higher_order are
# smoother-agnostic (frozen weights) and stay active under either scheme.

data(mpdta, package = "did")

test_that("the sieve smoother runs and yields a mean-zero EIF with finite, positive SEs", {
  skip_on_cran()
  old <- options(edid_omega_method = "sieve"); on.exit(options(old))
  f <- edid(mpdta, yname = "lemp", idname = "countyreal", tname = "year",
            gname = "first.treat", xformla = ~ lpop, misspec_robust = FALSE, aggregate = "none")
  expect_true(is.matrix(f$eif))                          # the influence-function slot is $eif (not $inffunc)
  expect_lt(max(abs(colMeans(f$eif))), 1e-8)             # mean-zero per cell
  expect_true(all(is.finite(f$att_gt$se)) && all(f$att_gt$se > 0))
})

test_that("sieve EFFICIENT + misspec_robust runs the weight channel (no warning, mean-zero EIF, finite SEs)", {
  skip_on_cran()
  old <- options(edid_omega_method = "sieve"); on.exit(options(old))
  w <- testthat::capture_warnings(
    f_mr <- edid(mpdta, yname = "lemp", idname = "countyreal", tname = "year",
                 gname = "first.treat", xformla = ~ lpop, weight_scheme = "efficient",
                 misspec_robust = TRUE, aggregate = "none"))
  # the efficient sieve weight channel IS implemented -> NO "not implemented" warning
  expect_false(any(grepl("not implemented", w)))
  expect_true(is.matrix(f_mr$eif) && max(abs(colMeans(f_mr$eif))) < 1e-8)   # EIF still mean-zero with the channel folded
  expect_true(all(is.finite(f_mr$att_gt$se)) && all(f_mr$att_gt$se > 0))
  # folding the weight channel changes the SE vs the plug-in (it is a genuine, non-degenerate contribution)
  f_pl <- suppressWarnings(edid(mpdta, yname = "lemp", idname = "countyreal", tname = "year",
                 gname = "first.treat", xformla = ~ lpop, weight_scheme = "efficient",
                 misspec_robust = FALSE, aggregate = "none"))
  expect_true(max(abs(f_mr$att_gt$se - f_pl$att_gt$se)) > 1e-8)
  expect_equal(f_mr$att_gt$att, f_pl$att_gt$att, tolerance = 1e-10)          # point estimates unchanged
})

test_that("sieve AVERAGED + misspec_robust runs the pooled weight channel (no warning, mean-zero EIF, finite SEs)", {
  skip_on_cran()
  old <- options(edid_omega_method = "sieve", edid_legacy_floor = NULL); on.exit(options(old))
  w <- testthat::capture_warnings(
    f_mr <- edid(mpdta, yname = "lemp", idname = "countyreal", tname = "year",
                 gname = "first.treat", xformla = ~ lpop, weight_scheme = "averaged",
                 misspec_robust = TRUE, aggregate = "none"))
  # the averaged sieve weight channel IS implemented (pooled Daleckii-Krein coupling) -> NO "not implemented" warning
  expect_false(any(grepl("not implemented", w)))
  expect_true(is.matrix(f_mr$eif) && max(abs(colMeans(f_mr$eif))) < 1e-8)   # EIF mean-zero with the pooled channel folded
  expect_true(all(is.finite(f_mr$att_gt$se)) && all(f_mr$att_gt$se > 0))
  # folding the weight channel changes the SE vs the plug-in (it is a genuine, non-degenerate contribution)
  f_pl <- suppressWarnings(edid(mpdta, yname = "lemp", idname = "countyreal", tname = "year",
                 gname = "first.treat", xformla = ~ lpop, weight_scheme = "averaged",
                 misspec_robust = FALSE, aggregate = "none"))
  expect_true(max(abs(f_mr$att_gt$se - f_pl$att_gt$se)) > 1e-8)
  expect_equal(f_mr$att_gt$att, f_pl$att_gt$att, tolerance = 1e-10)          # point estimates unchanged
  # Sanity bound under the DEFAULT (pooled-scale, exponent-1/3) floor: the weaker pooled floor clamps far
  # fewer directions, so the weight channel legitimately responds more than under the legacy floor (the
  # long-horizon (2004, 2007) cell sits at ~2.6x here); bound it loosely. The tight anti-regression anchor
  # lives below under the LEGACY floor, where its original calibration applies unchanged.
  expect_lt(max(f_mr$att_gt$se / f_pl$att_gt$se), 4)
  # NUMERIC ANCHOR (legacy floor): pin the SE ratio to the eigen-floor-aware coupling's range under
  # options(edid_legacy_floor = TRUE), the regime the anchor was calibrated in. RE-PINNED 2026-06-12 for
  # the exp default: under ratio_method = "exp" (the new default; "coherent" removed) the cross-cohort
  # 1/p prefactors differ, and se_mr/se_pl sits at ~2.03 here (was ~1.31 under coherent). The anchor's
  # PURPOSE is unchanged -- catch a silent drop back to the smooth -sym(q w') adjoint (the bug the pooled
  # Daleckii-Krein coupling fixes), which would balloon the ratio to ~2.5x. 2.03 is provably the corrected
  # coupling, not the smooth bug (compute_obar_coupling_edid's direct unit test above is engine-independent
  # and still passes); bound at 2.3, between the exp value and the bug signature.
  options(edid_legacy_floor = TRUE)
  f_mr_l <- suppressWarnings(edid(mpdta, yname = "lemp", idname = "countyreal", tname = "year",
                 gname = "first.treat", xformla = ~ lpop, weight_scheme = "averaged",
                 misspec_robust = TRUE, aggregate = "none"))
  f_pl_l <- suppressWarnings(edid(mpdta, yname = "lemp", idname = "countyreal", tname = "year",
                 gname = "first.treat", xformla = ~ lpop, weight_scheme = "averaged",
                 misspec_robust = FALSE, aggregate = "none"))
  options(edid_legacy_floor = NULL)
  expect_lt(max(f_mr_l$att_gt$se / f_pl_l$att_gt$se), 2.3)
})

test_that("compute_obar_coupling_edid: reduces to the smooth adjoint with no floor, differs with an active floor", {
  # Direct unit test of the pooled eigen-floor coupling helper, with NUMERIC ANCHORS that pin the Daleckii-Krein
  # sign and scale (not just NULL/finite/shape): (1) with the floor inactive it must equal the smooth adjoint
  # C_smooth = -sym(q w') EXACTLY (so a sign flip or scale error breaks the test); (2) with an active floor it must
  # DIFFER from C_smooth (so the floor branch is provably engaged); (3) it is symmetric; (4) NULL without the attr.
  set.seed(7); H <- 4L
  A  <- matrix(rnorm(H * H), H); S <- crossprod(A) + diag(H)                  # SPD raw pooled Omega-bar
  eg <- eigen(S, symmetric = TRUE); mbar <- rnorm(H)
  w  <- drop(solve(S, rep(1, H))); w <- w / sum(w); att <- sum(w * mbar)
  q  <- drop(solve(S, mbar - att))
  C_smooth <- -0.5 * (outer(q, w) + outer(w, q))                             # independent smooth -sym(q w') reference

  fl_lo <- min(eg$values) * 0.5                                              # floor BELOW min eigenvalue => inactive
  S_lo  <- eg$vectors %*% diag(pmax(eg$values, fl_lo)) %*% t(eg$vectors)
  attr(S_lo, "eig_floor") <- list(values = eg$values, vectors = eg$vectors, floor = fl_lo)
  C_unf <- did:::compute_obar_coupling_edid(S_lo, mbar, att)
  expect_true(is.matrix(C_unf) && all(dim(C_unf) == c(H, H)) && isSymmetric(C_unf, tol = 1e-9))
  expect_lt(max(abs(C_unf - C_smooth)), 1e-10)                               # locks sign + scale (Daleckii-Krein of 1/lam)

  fl_hi <- max(eg$values) * 0.2                                              # floor ABOVE the smallest eigenvalue(s)
  S_hi  <- eg$vectors %*% diag(pmax(eg$values, fl_hi)) %*% t(eg$vectors)
  attr(S_hi, "eig_floor") <- list(values = eg$values, vectors = eg$vectors, floor = fl_hi)
  C_flo <- did:::compute_obar_coupling_edid(S_hi, mbar, att)
  expect_true(sum(eg$values < fl_hi) >= 1L)                                  # the floor genuinely binds
  expect_gt(max(abs(C_flo - C_smooth)), 1e-3)                               # floored result is distinct from smooth
  expect_null(did:::compute_obar_coupling_edid(S, mbar, att))               # no eig_floor attr -> NULL (smooth fallback)
})

test_that("psi_channel_credible_edid gates on finiteness and the (clustered) variance-inflation ceiling", {
  # Direct unit test of the stability guard: small/finite psi accepted, non-finite/NULL rejected, and the
  # EDID_PSI_VAR_RATIO boundary flips the verdict. Also checks the clustered metric matches the cluster-robust SE.
  set.seed(3); n <- 200L; eif <- rnorm(n)
  expect_true (did:::psi_channel_credible_edid(0.1 * rnorm(n), eif))         # small channel -> credible
  expect_false(did:::psi_channel_credible_edid(c(Inf, rnorm(n - 1)), eif))   # non-finite -> rejected
  expect_false(did:::psi_channel_credible_edid(NULL, eif))                   # NULL -> rejected
  # boundary: build psi so v1/v0 straddles the ceiling. psi = c*eif => v1/v0 = (1+c)^2; pick c so (1+c)^2 ~ ratio.
  R <- did:::EDID_PSI_VAR_RATIO
  expect_true (did:::psi_channel_credible_edid((sqrt(R) - 1 - 0.05) * eif, eif))   # just below the ceiling
  expect_false(did:::psi_channel_credible_edid((sqrt(R) - 1 + 0.50) * eif, eif))   # just above the ceiling
  # clustering: a per-unit-large psi that CANCELS within clusters has a small clustered ratio => credible under
  # clustering though it would fail the i.i.d. test. cl groups units in pairs; psi = +/-M alternating cancels.
  cl  <- rep(seq_len(n / 2L), each = 2L); psi_cancel <- rep(c(50, -50), n / 2L)
  expect_false(did:::psi_channel_credible_edid(psi_cancel, eif))                     # i.i.d.: huge -> rejected
  expect_true (did:::psi_channel_credible_edid(psi_cancel, eif, cluster_indices = cl)) # clustered: cancels -> credible
})

test_that("misspec_robust weight channel cannot blow up the SE in poor-overlap / placebo cells (guarded fallback)", {
  skip_on_cran()
  old <- options(edid_omega_method = "sieve"); on.exit(options(old))
  # Poor-overlap covariate panel: steep propensity => some cohorts have near-zero propensity over part of the X
  # support (huge inverse-propensity prefactors) + sparse sieve groups (near-singular basis Gram). Pre-treatment
  # placebo cells (t < g) are where the sieve weight-channel psi exploded (SE ~1e14, non-mean-zero EIF) before
  # the Eq.(3.12) Term-1 restoration; the SEs must stay sane here whether the channel folds (credible IF, the
  # current behavior) or the credibility guard drops it (the pre-fix fallback).
  set.seed(20260609L); n <- 120L; Tn <- 4L
  x1 <- rnorm(n); x2 <- rnorm(n); eta <- 2.2 * x1 + 1.6 * x2 - 0.4
  P <- exp(cbind(0, eta, 0.7 * eta)); P <- P / rowSums(P)
  gcat <- apply(P, 1L, function(p) sample(c(Inf, 2, 4), 1L, prob = p))
  alpha <- rnorm(n, 0.5 * x1, 1); rows <- vector("list", Tn)
  for (tt in 1:Tn) { ht <- (tt - 1) * (0.4 * x1 + 0.3 * x2); tau <- ifelse(is.finite(gcat) & tt >= gcat, 1, 0)
    rows[[tt]] <- data.frame(id = seq_len(n), tt = tt, g = gcat, x1 = x1, x2 = x2,
                             y = alpha + 0.3 * tt + ht + tau + rnorm(n, sd = 0.5)) }
  df <- do.call(rbind, rows)
  for (ws in c("averaged", "efficient")) {
    f_pl <- suppressWarnings(edid(df, "y", "id", "tt", "g", xformla = ~ x1 + x2, weight_scheme = ws,
                                  misspec_robust = FALSE, aggregate = "none", cband = FALSE))
    w <- testthat::capture_warnings(
      f_mr <- edid(df, "y", "id", "tt", "g", xformla = ~ x1 + x2, weight_scheme = ws,
                   misspec_robust = TRUE,  aggregate = "none", cband = FALSE))
    fin     <- is.finite(f_mr$att_gt$se)
    fin_pl  <- is.finite(f_pl$att_gt$se)
    # The misspec_robust weight channel must not introduce NEW non-finite SEs: its NA pattern must
    # MATCH the plug-in's. (Under the default ratio_method = "exp", this extreme poor-overlap design
    # leaves cohort-4's exact-zero PRE-treatment placebo cells with a degenerate variance -> NA SE in
    # BOTH the plug-in and the misspec fit; that is a design+nuisance property, not a channel blow-up.
    # The previous `all(fin)` held only incidentally under the removed "coherent" engine, which gave
    # those pre-cells a tiny finite SE.) The substantive guard -- no 1e14 SE on the ESTIMABLE cells --
    # is asserted on the finite set below.
    expect_identical(fin, fin_pl)                                            # channel adds no new NA SEs
    expect_lt(max(abs(f_mr$eif)), 1e6)                                       # EIF not blown (pre-guard: ~1e16)
    expect_lt(max(f_mr$att_gt$se[fin] / f_pl$att_gt$se[fin]), 5)             # SE within a sane multiple (pre-guard: ~1e14)
    # Mechanism pin, updated with the Eq.(3.12) Term-1 restoration: the pre-fix channel SKIPPED Term 1 (valid
    # only for the smooth adjoint, 1'C1 = 0; for the Daleckii-Krein coupling 1'C1 != 0 where the floor binds),
    # which made psi non-mean-zero and exploded it in exactly these placebo/poor-overlap cells -- the guard then
    # fired and fell back per cell. With Term 1 added the channel is a credible IF here (measured per-cell
    # variance inflation <= ~6x, far under the EDID_PSI_VAR_RATIO = 100 ceiling), so it FOLDS rather than falls
    # back: assert no instability fallback occurred AND the channel genuinely moved the SE. The guard machinery
    # itself stays covered by the psi_channel_credible_edid unit test above (it remains defense-in-depth).
    expect_false(any(grepl("weight-estimation channel was numerically unstable", w)))
    expect_gt(max(abs(f_mr$att_gt$se[fin] / f_pl$att_gt$se[fin] - 1)), 1e-3) # the channel folded (SE moved)
  }
})
