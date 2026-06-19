# edid-hausman.R
# Hausman-type specification test of PT-All vs PT-Post for edid fits
# (Section 5.1, Theorem 5.1 of Chen, Sant'Anna & Xie 2025), plus the shared
# internals used by the Section-5 toolkit (edid_sargan, edid_frontier,
# edid_adaptive).

# ---------------------------------------------------------------------------
# Shared toolkit internals
# ---------------------------------------------------------------------------

# Validate that two edid fits are comparable: same sample (n, unit order, cohort
# assignment = the data fingerprint available on the fit), same design (periods,
# cohorts, anticipation), and same clustering. `require_pt = TRUE` additionally
# warns unless (unrestricted, restricted) = (PT-Post, PT-All), the pairing the
# Section-5 results are stated for.
.edid_toolkit_check_fits <- function(fit_unrestricted, fit_restricted, require_pt = TRUE) {
  if (!inherits(fit_unrestricted, "edid_fit") || !inherits(fit_restricted, "edid_fit")) {
    stop("`fit_unrestricted` and `fit_restricted` must be `edid_fit` objects returned by edid().",
         call. = FALSE)
  }
  fu <- fit_unrestricted; fr <- fit_restricted
  if (!identical(fu$n, fr$n)) {
    stop("The two fits have different sample sizes (n = ", fu$n, " vs ", fr$n,
         "); they must be estimated on the same data.", call. = FALSE)
  }
  if (!is.null(fu$all_units) && !is.null(fr$all_units) && !identical(fu$all_units, fr$all_units)) {
    stop("The two fits carry different unit identifiers / unit order; the per-unit influence ",
         "functions cannot be differenced. Fit both estimators on the same data.", call. = FALSE)
  }
  if (!is.null(fu$unit_cohorts) && !is.null(fr$unit_cohorts) &&
      !isTRUE(all.equal(fu$unit_cohorts, fr$unit_cohorts))) {
    stop("The two fits assign units to different cohorts (data fingerprint mismatch); ",
         "they must be estimated on the same data.", call. = FALSE)
  }
  if (!isTRUE(all.equal(fu$time_periods, fr$time_periods)) ||
      !isTRUE(all.equal(fu$treatment_groups, fr$treatment_groups))) {
    stop("The two fits have different time periods or treatment cohorts.", call. = FALSE)
  }
  if (!identical(fu$anticipation, fr$anticipation)) {
    stop("The two fits use different `anticipation` values.", call. = FALSE)
  }
  cu <- fu$cluster_indices; cr <- fr$cluster_indices
  if (is.null(cu) != is.null(cr) || (!is.null(cu) && !identical(cu, cr))) {
    stop("The two fits use different cluster assignments (`clustervars`); the estimator-",
         "difference covariance must be computed under a common clustering.", call. = FALSE)
  }
  if (require_pt &&
      !(identical(fu$pt_assumption, "post") && identical(fr$pt_assumption, "all"))) {
    warning("Expected `fit_unrestricted` from edid(pt_assumption = \"post\") and ",
            "`fit_restricted` from edid(pt_assumption = \"all\") (the conservative vs efficient ",
            "pairing of Section 5 of Chen, Sant'Anna & Xie 2025). Got pt_assumption = \"",
            fu$pt_assumption, "\" vs \"", fr$pt_assumption, "\"; results are only meaningful if ",
            "the restricted fit imposes strictly more moment restrictions.", call. = FALSE)
  }
  invisible(TRUE)
}

# ---------------------------------------------------------------------------
# Leg-health guards for the Section-5 toolkit (broken-leg + few-cluster)
# ---------------------------------------------------------------------------

# Number of clusters feeding a fit's cluster-robust covariance (G; n when unclustered,
# i.e. unit-level / iid). The few-cluster guard compares this to EDID_FEWCLUSTER_MIN.
.edid_n_clusters <- function(fit) {
  ci <- fit$cluster_indices
  if (is.null(ci)) (fit$n %||% length(unique(fit$all_units))) else length(unique(ci))
}

# Effective (Kish ESS) sample size of the units feeding a fit's aggregation IFs,
# for the weight-dispersion noise floor of .edid_if_diff_quadform(). Equals the
# full unit count `n` exactly when the fit is unweighted (unit_weights NULL),
# making the noise floor a no-op there (byte-identical legacy behaviour). Reuses
# the same n_eff_edid() machinery the cov-ridge fix added; every unit feeds the
# event-study / overall aggregation IF, so the active mask is all-TRUE.
.edid_overid_n_eff <- function(fit) {
  n <- fit$n %||% length(fit$all_units)
  n_eff_edid(fit$unit_weights, rep(TRUE, n), n_full = n)
}

# Is a fit a numerically degenerate leg? A non-rejection (or a point estimate) built on
# such a leg is HOLLOW -- the restricted/efficient fit carried extreme propensity ratios,
# a non-credible weight channel, or cross-cohort hedges that carry the estimand, or its
# reported SEs are non-finite. Reads the fit's own $diagnostics read-out (set by edid())
# plus the aggregation SEs the test will difference. Returns a list with a logical
# `broken` and a character vector of `reasons` (empty when healthy). Older fits without
# $diagnostics degrade gracefully to the SE-finiteness check only.
.edid_leg_health <- function(fit, label) {
  reasons <- character(0L)
  d <- fit$diagnostics
  if (!is.null(d)) {
    if (isTRUE(d$unstable)) {
      if ((d$n_extreme_ratio %||% 0L) > 0L)
        reasons <- c(reasons, sprintf("extreme propensity ratios in %d estimation step(s)", d$n_extreme_ratio))
      if ((d$n_psi_unstable %||% 0L) > 0L)
        reasons <- c(reasons, sprintf("weight-estimation channel not a credible IF in %d cell(s)", d$n_psi_unstable))
      if (isTRUE(d$net_hedge_flag))
        reasons <- c(reasons, sprintf("cross-cohort hedges carry the estimand (net hedge mass %.2f ~ gross %.2f)",
                                      d$net_hedge_mass %||% NA_real_, d$gross_hedge_mass %||% NA_real_))
    }
  }
  # Non-finite reported SEs on the headline (overall) aggregation: a hard breakage signal
  # available even without $diagnostics.
  ov <- fit$overall
  se_ov <- if (!is.null(ov)) suppressWarnings(ov$overall.se) else NULL
  if (!is.null(se_ov) && length(se_ov) && any(!is.finite(se_ov))) {
    reasons <- c(reasons, "non-finite reported standard error on the overall aggregation")
  }
  list(broken = length(reasons) > 0L, reasons = unique(reasons), label = label)
}

# Dynamic (event-study) AGGTEobj of a fit: reuse the stored one, else aggregate
# on the fly with the same machinery edid() uses (aggte_edid -> did::aggte).
.edid_dynamic_aggte <- function(fit) {
  a <- fit$event_study
  if (is.null(a)) a <- aggte_edid(fit, type = "dynamic", na.rm = TRUE)
  if (is.null(a) || !inherits(a, "AGGTEobj")) {
    stop("Could not construct the event-study aggregation for an edid fit.", call. = FALSE)
  }
  a
}

# Per-unit influence functions + estimates of the requested parameter vector.
# parameter = "event_study": the ES(e) vector over e_set (default: all finite
# post-treatment e's); parameter = "overall": the scalar ES_avg (the average of
# ES(e) over e >= 0, i.e. the dynamic AGGTEobj's overall). The IFs are the same
# aggregation IFs that vcov.edid_fit() reads (did::aggte's per-element
# influence functions, including the cohort-share weight-estimation correction).
.edid_param_ifs <- function(fit, parameter, e_set = NULL) {
  a <- .edid_dynamic_aggte(fit)
  g <- .edid_agg_if(a)
  if (identical(parameter, "overall")) {
    if (is.null(a$overall.att) || is.null(g$overall)) {
      stop("The dynamic aggregation does not expose the overall ES_avg influence function.",
           call. = FALSE)
    }
    return(list(e = NA_real_, est = a$overall.att,
                IF = matrix(as.numeric(g$overall), ncol = 1L)))
  }
  if (is.null(g$egt) || is.null(a$egt)) {
    stop("The dynamic aggregation does not expose per-e influence functions.", call. = FALSE)
  }
  e_all <- a$egt
  post  <- e_all[e_all >= 0 & is.finite(a$att.egt)]
  if (is.null(e_set)) {
    e_set <- post
  } else {
    bad <- setdiff(e_set, post)
    if (length(bad) > 0L) {
      stop("`e_set` contains event times not available as finite post-treatment ES(e) ",
           "coordinates in this fit: ", paste(bad, collapse = ", "), call. = FALSE)
    }
  }
  e_set <- sort(unique(e_set))
  if (length(e_set) == 0L) stop("No post-treatment event times available.", call. = FALSE)
  idx <- match(e_set, e_all)
  list(e = e_set, est = as.numeric(a$att.egt[idx]),
       IF = as.matrix(g$egt)[, idx, drop = FALSE])
}

# Default e_set for a two-fit comparison: the intersection of both fits'
# finite post-treatment event times.
.edid_shared_e_set <- function(fit_unrestricted, fit_restricted, e_set = NULL) {
  if (!is.null(e_set)) return(sort(unique(e_set)))
  eu <- .edid_param_ifs(fit_unrestricted, "event_study")$e
  er <- .edid_param_ifs(fit_restricted,  "event_study")$e
  shared <- intersect(eu, er)
  if (length(shared) == 0L) {
    stop("The two fits share no finite post-treatment event times.", call. = FALSE)
  }
  sort(shared)
}

# Rank-aware quadratic form for the IF-difference Hausman statistic:
#   H = n * d' D^+ d,   D = n * Var-hat(xi)  (cluster-robust when clustered),
# with df = rank(D) by eigenvalue thresholding and the Moore-Penrose
# pseudoinverse on the rank-deficient branch. The full-rank branch uses the
# exact inverse (numerically identical to solve()). Estimating D from the
# per-unit IF differences xi_i = psi_U,i - psi_R,i makes it positive
# semi-definite in finite samples (the footnote to eqn (5.3) of Chen,
# Sant'Anna & Xie 2025). Caveat (Andrews 1987): the chi^2(rank) limit on the
# rank-deficient branch additionally requires the estimated rank to be
# consistent for the true rank; the eigenvalue threshold is the standard
# practical device but is not a formal guarantee.
#
# `v_scale` is the absolute variance scale of the constituent estimators (the
# largest per-coordinate asymptotic variance of either fit; 1 if unknown). It
# feeds the degenerate-contrast guard, the joint-path mirror of the scalar
# .edid_scalar_hausman() guard: when the two estimators coincide (e.g. both
# fits pinned to the same just-identified moments by the thin-cohort guard),
# xi is pure float dust (differences ~1e-17, D entries ~1e-33) and the
# RELATIVE eigenvalue threshold below would still "find" rank in that noise,
# returning an arbitrary large H with a tiny p-value. A D that is negligible
# on the ABSOLUTE scale of the estimators' own variances carries no testable
# contrast: report H = 0, df = 0, p = 1 with degenerate = TRUE instead.
#
# `n_eff` is the effective (Kish ESS) sample size of the units feeding xi
# (`n_eff_edid(unit_weights, ...)`); it equals `n` exactly on the unweighted
# path. UNDER DISPERSED WEIGHTS + THIN COHORTS the cluster-robust D-hat has
# downward-biased small eigenvalues (it is an average over n_eff << n
# independent pieces, but each direction's sampling floor is set by n_eff, not
# n). The bare pseudoinverse `/ ev$values[pos]` then OVER-AMPLIFIES those
# statistically-unreliable directions and inflates H spuriously (the weighted
# Bailey-Goodman-Bacon over-id: H = 290.6, p < 2e-16, even though the faithful
# weighted pre-trend is clean, p ~ 0.43). We therefore RIDGE the spectrum at a
# weight-dispersion-aware sampling-noise floor before inverting: each eigenvalue
# is lifted to `max(lambda_k, floor)` with
#   floor = mx * max(sqrt(eps), c * sqrt(disp / n_eff)),  disp = max(0, n/n_eff - 1),
# i.e. the over-amplified directions are damped rather than discarded (the rank,
# hence the chi^2 df, is preserved -- a smooth Tikhonov lift, not a lumpy hard
# rank cut on a smoothly-decaying spectrum). The scale `sqrt(disp / n_eff)` is
# the relative sampling SD of an estimated-covariance eigenvalue built from
# n_eff effective pieces, inflated by the excess weight dispersion `disp` (the
# part of the small-eigenvalue bias that the raw 1/n averaging does not see); the
# constant c = 1 is the bare random-matrix noise-edge coefficient (NOT a size-tuned
# knob): validated to control over-id size (~0.01, slightly conservative, -> nominal as
# n_eff grows) at power equal to any larger c on a clean-PT, thin-cohort, dispersed-weight
# DGP (see quality_reports/drafts/gate_runs/overid/{overid_fix,val_power-size-c,
# c1_lock_and_Klarge_scope}.md). INVARIANTS, all by
# construction: (i) UNWEIGHTED / uniform weights => n_eff = n => disp = 0 =>
# floor = mx*sqrt(eps), so NO eigenvalue is lifted and the code takes the
# ORIGINAL solve()/pseudoinverse branch byte-for-byte (the regularizer is a
# no-op); (ii) for any fixed weight distribution disp -> const and n_eff -> inf,
# so floor -> mx*sqrt(eps): the lift VANISHES asymptotically and the chi^2(rank)
# limit on well-identified designs is untouched (the project-wide invariant that
# every Omega/D regularizer is asymptotically negligible); (iii) the degenerate-
# contrast guard (H = 0, df = 0, p = 1) still fires first, before any lift.
EDID_OVERID_DISP_C <- 1.0  # noise-floor constant = the bare MP noise-edge coefficient (uniform, not size-tuned)
# Bell-McCaffrey (2002) / Pustejovsky-Tipton (2018, AHT) Satterthwaite EFFECTIVE DEGREES OF FREEDOM of
# the cluster-robust IF-difference covariance D-hat. D-hat is a sandwich (meat = sum over the G
# independent units/clusters of u_g u_g', u_g = the cluster's IF-difference sum); its effective df is
# NOT n but ~ the effective number of independent clusters feeding it, which is below G when a few
# clusters dominate. With the per-cluster "leverages" w_g = s * u_g' D-hat^+ u_g (s the cluster-robust
# scale; sum_g w_g = rank(D-hat) by construction since s*sum_g u_g u_g' = D-hat), the Satterthwaite
# match gives
#   m_hat = (sum_g w_g)^2 / sum_g w_g^2 = rk^2 / sum_g w_g^2.
# Properties: (i) clustered with G clusters => m_hat <= G, and m_hat < G-1 under cluster imbalance
# (the few-cluster over-rejection regime); (ii) iid with many balanced units (G = n) => each w_i ~ rk/n
# so m_hat ~ n -> the F reference below -> chi^2 (asymptotically negligible, NO-OP on the already-sized
# unweighted path); (iii) dispersed unit weights concentrate the leverages => m_hat ~ Kish ESS << n,
# supplying exactly the finite-sample correction the weighted over-id needs. References: Bell &
# McCaffrey (2002) Survey Methodology; Pustejovsky & Tipton (2018) JBES (AHT test, eqs 12-13);
# Imbens & Kolesar (2016) REStat.
.edid_overid_satdf <- function(xi, V, lam, cluster_indices, n, rk) {
  if (is.null(cluster_indices)) { U <- xi; G <- nrow(xi); cf <- 1 }
  else { U <- rowsum(xi, cluster_indices); G <- nrow(U); cf <- if (G > 1L) G / (G - 1) else 1 }
  s <- cf / n
  proj <- U %*% V                                        # G x rk : column k = u_g' v_k
  w    <- s * drop((proj * proj) %*% (1 / lam))          # w_g = s * sum_k (u_g'v_k)^2 / lam_k
  sw2  <- sum(w * w)
  if (!is.finite(sw2) || sw2 <= 0) return(Inf)
  (sum(w)^2) / sw2
}

# Rank-aware IF-difference quadratic form H = n d' D^+ d, df = rank(D-hat), referred to the AHT
# (approximate Hotelling T^2) F distribution with the Satterthwaite effective df m_hat above:
#   H * (m_hat - rk + 1) / (m_hat * rk)  ~  F(rk, m_hat - rk + 1).
# This is the exact Hotelling rescaling of a quadratic form in an ESTIMATED covariance (Bell-McCaffrey /
# Pustejovsky-Tipton), and -> chi^2(rk) as m_hat -> infinity, so it is asymptotically negligible (the
# project-wide invariant). It corrects the finite-sample over-rejection of the chi^2 reference that the
# cluster-robust / dispersed-weight D-hat (effective df = #clusters / Kish ESS, NOT n) suffers; on the
# many-balanced-iid-units path it is a numerical no-op (m_hat ~ n). The dispersed-weight eigen-ridge is
# retained underneath as a numerical conditioning safeguard.
.edid_if_diff_quadform <- function(d, xi, n, cluster_indices, v_scale = 1, n_eff = n,
                                   rel_tol = 0) {
  xi <- as.matrix(xi)
  D  <- n * cluster_cov_edid(xi, cluster_indices, n)     # = E_n[xi xi'] when iid
  if (any(!is.finite(D))) {
    return(list(statistic = NA_real_, df = NA_integer_, p_value = NA_real_, D = D,
                degenerate = NA, m_eff = NA_real_, df2 = NA_real_))
  }
  auto   <- identical(rel_tol, "auto")
  rt_num <- if (auto) 0 else as.numeric(rel_tol)
  if (!is.finite(n_eff) || n_eff <= 0) n_eff <- n
  eps_D  <- .Machine$double.eps^0.5
  # Absolute degenerate guard (estimators coincide; D ~ 0 relative to the parameter variance scale)
  # ONLY on the bare convention (rel_tol == 0) -> edid_hausman / edid_sargan are byte-identical. For
  # rel_tol > 0 / "auto" (edid_overid) this v_scale-relative guard is SKIPPED: it fires spuriously when
  # D has genuine rank but small magnitude relative to a large efficient-fit v_scale (the Bailey-GB
  # panel), masking rank-deficiency as a false p = 1. There the eigen-based logic below decides
  # (lambda_max ~ 0 -> genuine p = 1; relative floor drops all real directions -> NA).
  if (rt_num == 0 && !auto && max(abs(D)) <= eps_D * max(v_scale, 1)) {
    return(list(statistic = 0, df = 0L, p_value = 1, D = D, degenerate = TRUE,
                m_eff = NA_real_, df2 = NA_real_))
  }
  ev  <- eigen(D, symmetric = TRUE)
  mx  <- max(ev$values, 0)
  rk_bare <- if (mx > 0) sum(ev$values > mx * sqrt(.Machine$double.eps)) else 0L
  if (rk_bare == 0L) {                                   # lambda_max ~ 0: D genuinely null (coincide)
    return(list(statistic = 0, df = 0L, p_value = 1, D = D, degenerate = TRUE,
                m_eff = NA_real_, df2 = NA_real_))
  }
  G_eff <- if (is.null(cluster_indices)) n_eff else length(unique(cluster_indices))
  # SATURATION guard (over-id path only) -- a DEFENSIVE guard for genuinely FEW-CLUSTER designs. When the
  # bare numerical rank reaches the cluster-robust rank ceiling G_eff - 1 (a centered sandwich of G_eff
  # cluster scores has rank <= G_eff - 1), D-hat is full-rank-for-its-cluster-budget => there is NO null
  # space to anchor a noise floor => the over-identifying dimension meets/exceeds the cluster budget and
  # the JOINT over-id is not reliably estimable. Return NA (rank_deficient), never a misleading p; read the
  # per-cell breakdown / edid_sargan. This fires ONLY when the user clusters coarsely AND the over-id
  # dimension is large relative to the cluster count. It does NOT fire under the DEFAULT unit-level
  # clustering (G_eff = n units, thousands of pieces, comfortably supports the fixed low-rank over-id) --
  # e.g. Bailey-Goodman-Bacon, clustered at the county=unit level per the original paper (G_eff ~ 3059 >>
  # structural rank 83), computes a normal joint J (df 15, p ~ 5.6e-5, REJECTS). Cannot fire on the
  # validated designs (clustered no-cov r_bare = 6 << G_eff - 1; unclustered cov/mpdta G_eff = n_eff).
  # Gated on `auto`: rel_tol = 0 (edid_hausman/edid_sargan) is byte-identical; a numeric rel_tol is respected.
  if (auto && G_eff >= 2L && rk_bare >= G_eff - 1L) {
    return(list(statistic = NA_real_, df = 0L, p_value = NA_real_, D = D, degenerate = TRUE,
                rank_deficient = TRUE, m_eff = NA_real_, df2 = NA_real_))
  }
  # Rank threshold. rel_tol = 0 -> the byte-identical numerical cut mx*sqrt(eps) (edid_hausman /
  # edid_sargan). rel_tol = "auto" (edid_overid default) -> the EFFECTIVE-RANK relative floor
  # r_bare / n_eff, n_eff = Kish effective sample size. This is a DIVISION OF LABOR with the AHT F below:
  # the FLOOR determines the RANK (which spectral directions are genuine over-id content vs the covariate
  # decaying noise tail), and the AHT F handles the few-cluster sampling RELIABILITY of the surviving
  # directions via m = G_eff - 1. The denominator for the floor is n_eff (the spectral noise scale of the
  # covariate smear), NOT G_eff. [Validated empirically, round-3b MC on GENUINELY clustered data (ICC 0.44):
  # n_eff recovers the true rank (df ~ 6) with nominal SIZE and strong POWER (.95) on clustered cov AND
  # no-cov; the G_eff denominator instead COLLAPSES the rank (df ~ 1.3), under-rejects (.008), and loses
  # ~40% power (.54) -- because r_bare over-counts the numerator, so r_bare/G_eff overshoots the noise edge.
  # The cluster-robust over-resolution that motivated G_eff (Bailey) is handled NOT by the denominator but
  # by the SATURATION guard above (r_bare = G_eff - 1 -> NA), which fires before this floor.] NO-OP on the
  # clean low-rank no-covariate spectrum (r_bare small => tiny floor below the genuine eigenvalues, even at
  # small n_eff); lands in the spectral gap on the covariate decaying spectrum. A fixed numeric rel_tol
  # overrides.
  rel_eff <- if (auto) rk_bare / n_eff else rt_num
  tol <- mx * max(sqrt(.Machine$double.eps), rel_eff)
  pos <- ev$values > tol
  rk  <- sum(pos)
  if (rk == 0L) {
    # The relative floor dropped EVERY real direction (rel_eff >= 1: extreme rank-deficiency).
    # Uncomputable -> NA (rank_deficient = TRUE), never a misleading p = 1. (For rel_tol = 0 this is
    # unreachable -- rk_bare >= 1 guarantees rk >= 1 at the bare cut -- so the bare convention is safe.)
    if (rel_eff > sqrt(.Machine$double.eps)) {
      return(list(statistic = NA_real_, df = 0L, p_value = NA_real_, D = D, degenerate = TRUE,
                  rank_deficient = TRUE, m_eff = NA_real_, df2 = NA_real_))
    }
    return(list(statistic = 0, df = 0L, p_value = 1, D = D, degenerate = TRUE,
                m_eff = NA_real_, df2 = NA_real_))
  }
  V   <- ev$vectors[, pos, drop = FALSE]
  # NO eigen-ridge. The finite-sample SIZING is carried entirely by the AHT effective-df F below, which
  # is the principled, asymptotically-negligible correction and is NOMINAL + power-preserving for the
  # dispersed-weight / few-cluster regimes. The previous dispersed-weight eigen-ridge (a fixed eigenvalue
  # lift) is REMOVED because it DOUBLE-corrected with the F and drove the weighted size to ~0 (over-
  # conservative, no power: MC dispersed size collapsed to 0.000-0.012 vs the F-only 0.036-0.058 nominal);
  # and any fixed lift also over-corrects healthy dispersion (a no-op lift is impossible at a fixed
  # floor). A genuinely near-singular D-hat (thin-cohort / collinear-moment blow-up, the Bailey regime)
  # is NOT ridged into a spurious "clean" verdict; it is FLAGGED by m_sat << G_eff (and the few-cluster /
  # leg-unstable guards), where the protocol reads the localized edid_sargan rather than the diffuse
  # joint -- the honest treatment. Only the numerical rank threshold (mx*sqrt(eps), above) is applied.
  lam <- ev$values[pos]
  H   <- as.numeric(n * crossprod(d, V %*% (crossprod(V, d) / lam)))
  # AHT (approximate Hotelling T^2) effective-df F reference. The cluster-robust D-hat is a sandwich
  # built from G_eff INDEPENDENT pieces -- the number of CLUSTERS when clustered, else the Kish effective
  # sample size n_eff of the (possibly weighted) units -- so its reliability, hence the finite-sample
  # reference, is governed by G_eff, NOT n:
  #   H * (m - rk + 1) / (m * rk)  ~  F(rk, m - rk + 1),   m = G_eff - 1.
  # This is the exact Hotelling rescaling of a quadratic form in an ESTIMATED covariance; it -> chi^2(rk)
  # as G_eff -> inf (a numerical no-op for many balanced iid units; asymptotically negligible, the
  # project invariant) and removes the chi^2 over-rejection of the few-cluster / dispersed-weight D-hat.
  # Validated nominal on correct spec (clustered 0.32 -> 0.05; dispersed 0.07 -> 0.05; unweighted no-op),
  # power-preserving. Refs: Bell-McCaffrey (2002); Pustejovsky-Tipton (2018, JBES); Imbens-Kolesar (2016).
  # m_sat is the Bell-McCaffrey/Satterthwaite LEVERAGE df (rk^2 / sum_g w_g^2, w_g the per-cluster/unit
  # leverage of D-hat); m_sat << G_eff FLAGS a D-hat dominated by a few high-leverage clusters/units
  # (weak overlap / severe imbalance), where even the F reference is fragile and trimming is the remedy.
  # (G_eff defined above, where it also sets the relative rank floor's denominator.)
  m_sat <- .edid_overid_satdf(xi, V, ev$values[pos], cluster_indices, n, rk)
  m     <- G_eff - 1
  if (is.finite(m) && m > rk) {
    df2 <- m - rk + 1
    p   <- stats::pf(H * df2 / (m * rk), df1 = rk, df2 = df2, lower.tail = FALSE)
  } else {
    # too few effective pieces (G_eff - 1 <= rk): F denominator df <= 1, uninformative; fall back to
    # the chi^2 reference and flag via df2 = NA (callers warn, mirroring the few-cluster guard).
    df2 <- NA_real_
    p   <- stats::pchisq(H, df = rk, lower.tail = FALSE)
  }
  list(statistic = H, df = rk, p_value = p, D = D, degenerate = FALSE,
       m_eff = m, m_sat = m_sat, df2 = df2)
}

# Scalar Hausman component H = n d^2 / D with the degenerate-D guard of
# Theorem 5.2 (D > 0 is required; xi ~ 0 makes the statistic 0/0, so report
# H = 0 / p = 1 instead of NaN). The guard is relative to the parameter's
# asymptotic variance scale.
#
# The scalar statistic divides by a single, directly-estimated variance D (a
# weighted mean of squares), NOT by an inverted small eigenvalue, so it does NOT
# suffer the joint test's pseudoinverse over-amplification: the weighted
# Bailey-Goodman-Bacon scalar statistics are already sane (ES_avg H = 3.7,
# p = 0.05) at the exact same dispersed weights that send the joint statistic to
# H = 290.6. The weight-dispersion noise floor is therefore an over-id-JOINT
# device, not needed in 1-D; we accept `n_eff` for signature parity with
# .edid_if_diff_quadform (and so callers pass it uniformly) but the scalar H is
# left BYTE-IDENTICAL -- adding a 1-D floor would only ever shrink an already-
# sane statistic and could mask genuine per-coordinate evidence.
.edid_scalar_hausman <- function(d, xi_vec, n, cluster_indices, v_scale = 1, n_eff = n) {
  D <- as.numeric(n * cluster_cov_edid(matrix(xi_vec, ncol = 1L), cluster_indices, n))
  eps_D <- .Machine$double.eps^0.5
  if (!is.finite(D) || D <= eps_D * max(v_scale, 1)) {
    return(list(D = D, H = 0, p_value = 1, degenerate = TRUE))
  }
  H <- n * d^2 / D
  list(D = D, H = H, p_value = stats::pchisq(H, df = 1, lower.tail = FALSE),
       degenerate = FALSE)
}

# ---------------------------------------------------------------------------
# edid_hausman
# ---------------------------------------------------------------------------

#' Hausman test of PT-All against PT-Post for edid fits
#'
#' Implements the Hausman-type specification test of Theorem 5.1 in Chen,
#' Sant'Anna & Xie (2025): it compares the efficient event-study estimator
#' \eqn{\widehat{ES}} (consistent and semiparametrically efficient under
#' PT-All) with the conservative just-identified estimator
#' \eqn{\widecheck{ES}} of eqns (5.1)-(5.2) (consistent under PT-Post alone),
#' via the statistic of eqn (5.3),
#' \deqn{\widehat{H} = n\,(\widehat{ES} - \widecheck{ES})'\,\widehat{D}^{-1}\,
#'   (\widehat{ES} - \widecheck{ES}),}
#' where \eqn{\widehat{D}} is estimated from the per-unit difference of the two
#' estimators' influence functions, \eqn{\xi_i = \psi_{U,i} - \psi_{R,i}} --- the
#' positive semi-definite rendering noted in the footnote to eqn (5.3). Under
#' PT-All, \eqn{\widehat{H} \overset{d}{\to} \chi^2(|\mathcal{E}|)}; rejection
#' is evidence against the additional moment restrictions that PT-All imposes
#' beyond PT-Post.
#'
#' @param fit_unrestricted An \code{edid_fit} from
#'   \code{edid(..., pt_assumption = "post")}: the conservative just-identified
#'   estimator (the paper's staggered just-identification corollary),
#'   consistent under PT-Post alone. Its event-study aggregation (via
#'   \code{did::aggte}) includes the cohort-share weight-estimation
#'   influence-function correction, matching the conservative estimator used in
#'   the paper's empirical application.
#' @param fit_restricted An \code{edid_fit} from
#'   \code{edid(..., pt_assumption = "all")}: the efficient estimator under
#'   PT-All. Both fits must be estimated on the same data with the same
#'   clustering.
#' @param parameter \code{"event_study"} (default) for the joint test over the
#'   post-treatment event-study coefficients \eqn{ES(e), e \in \mathcal{E}}, or
#'   \code{"overall"} for the scalar test on \eqn{ES_{\mathrm{avg}}} (the
#'   average of \eqn{ES(e)} over \eqn{e \ge 0}).
#' @param e_set Numeric vector of post-treatment event times defining
#'   \eqn{\mathcal{E}}, or \code{NULL} (default: the intersection of the two
#'   fits' finite post-treatment event times). Ignored for
#'   \code{parameter = "overall"}.
#' @param data The panel data used to fit the two legs, or \code{NULL} (default),
#'   in which case the data expression stored in \code{fit_restricted$call} is
#'   re-evaluated in the caller's environment. Both legs are refit in the
#'   \strong{efficient plug-in configuration} (all estimation-effect channels
#'   off) before the contrast is formed, so the over-identification statistic
#'   uses the efficient inverse-variance covariance (Andrews, Chen and Tecchio
#'   2025) rather than any misspecification-robust SE the fits may report; the
#'   point estimates, hence the contrast \eqn{d}, are unchanged. Supply
#'   \code{data} explicitly when the original object is no longer reachable.
#'
#' @details
#' The joint statistic uses \eqn{df = \mathrm{rank}(\widehat{D})} by eigenvalue
#' thresholding with a Moore-Penrose pseudoinverse on the rank-deficient branch
#' (the generically full-rank case reproduces the exact-inverse statistic with
#' \eqn{df = |\mathcal{E}|}). Following Andrews (1987), the
#' \eqn{\chi^2(\mathrm{rank})} limit under rank deficiency additionally
#' requires the estimated rank to be consistent; the threshold is the standard
#' practical device, not a formal guarantee. The covariance \eqn{\widehat{D}}
#' is cluster-robust when the fits carry cluster assignments.
#'
#' \strong{Finite-sample reference (AHT effective-df F).} \eqn{\widehat{D}} is a
#' sandwich estimate built from \eqn{G_{\mathrm{eff}}} independent pieces --- the
#' number of clusters when clustered, else the Kish effective sample size
#' \eqn{n_{\mathrm{eff}}} of the (possibly weighted) units --- so its reliability,
#' and hence the reference distribution, is governed by \eqn{G_{\mathrm{eff}}},
#' \emph{not} \eqn{n}. The \eqn{\chi^2} reference therefore over-rejects with few
#' clusters or dispersed weights. \eqn{\widehat{H}} is instead referred to the
#' approximate Hotelling \eqn{T^2} (AHT) F distribution,
#' \deqn{\widehat{H}\,\frac{m - df + 1}{m\,df} \;\sim\; F(df,\; m - df + 1),
#'   \qquad m = G_{\mathrm{eff}} - 1,}
#' the exact Hotelling rescaling of a quadratic form in an estimated covariance
#' (Bell & McCaffrey 2002; Pustejovsky & Tipton 2018; Imbens & Kolesar 2016). It
#' converges to \eqn{\chi^2(df)} as \eqn{G_{\mathrm{eff}} \to \infty} (a numerical
#' no-op for many balanced i.i.d. units; asymptotically negligible), and removes
#' the finite-sample over-rejection of the few-cluster / dispersed-weight
#' \eqn{\widehat{D}}. When \eqn{G_{\mathrm{eff}} - 1 \le df} the F denominator df
#' is \eqn{\le 1} and the test is uninformative; the p-value then falls back to
#' \eqn{\chi^2} and \code{df2} is \code{NA} (flagged like the few-cluster guard).
#' The Bell-McCaffrey/Satterthwaite \emph{leverage} effective df
#' \eqn{\widehat m_{\mathrm{sat}} = df^2 / \sum_g w_g^2} (\eqn{w_g} the per-cluster
#' leverage of \eqn{\widehat{D}}) is reported as \code{m_sat}: when
#' \eqn{\widehat m_{\mathrm{sat}} \ll G_{\mathrm{eff}}} the covariance is dominated
#' by a few high-leverage units/clusters (weak overlap / severe imbalance), where
#' even the F reference is fragile --- trim overlap (\code{trim_level}) and read
#' the localized \code{\link{edid_sargan}} rather than the diffuse joint statistic.
#'
#' The returned object also reports the scalar per-coordinate statistics
#' \eqn{H_{\theta,n} = n(\widehat\theta_U - \widehat\theta_R)^2/\widehat{D}}
#' of eqn (5.5) for each \eqn{ES(e)} and for \eqn{ES_{\mathrm{avg}}}, with a
#' degenerate-\eqn{\widehat{D}} guard (coordinates where the two estimators
#' coincide report \eqn{H = 0}, \eqn{p = 1}). The joint statistic carries the
#' same guard: when the two estimators coincide on every coordinate --- e.g.
#' both fits pinned to the same just-identified moments by the thin-cohort
#' guard, so \eqn{\widehat{D}} is numerical noise relative to the estimators'
#' own variances --- the joint contrast is degenerate and is reported as
#' \eqn{H = 0}, \eqn{df = 0}, \eqn{p = 1} (with a message and
#' \code{degenerate = TRUE}) rather than ranking the noise.
#'
#' @return An object of class \code{edid_hausman}: a list with elements
#'   \code{statistic}, \code{df}, \code{p_value} (the joint test, AHT effective-df
#'   F p-value), \code{m_eff} (the AHT effective df used,
#'   \eqn{G_{\mathrm{eff}} - 1}: clusters minus one, or Kish \eqn{n_{\mathrm{eff}}}
#'   minus one when unclustered), \code{m_sat} (the Bell-McCaffrey/Satterthwaite
#'   leverage effective df, a fragility diagnostic; \code{m_sat << m_eff} flags
#'   weak overlap / severe imbalance), \code{df2} (the F denominator df
#'   \eqn{m_{\mathrm{eff}} - df + 1}; \code{NA} when it fell back to \eqn{\chi^2}),
#'   \code{degenerate} (\code{TRUE} when the joint contrast was degenerate and
#'   the \eqn{H = 0}, \eqn{df = 0}, \eqn{p = 1} guard applied), \code{d}
#'   (the estimate difference vector, unrestricted minus restricted), \code{D}
#'   (the estimated asymptotic covariance of \eqn{\sqrt{n}\,d}), \code{scalar}
#'   (data.frame of per-coordinate eqn (5.5) statistics, including an
#'   \code{ES_avg} row), \code{parameter}, \code{e_set}, \code{n},
#'   \code{clustered}, plus the round-3 sanity guards:
#'   \code{leg_unstable} (\code{TRUE} when a constituent fit is numerically
#'   degenerate -- extreme propensity ratios, a non-credible weight channel,
#'   cross-cohort hedges carrying the estimand, or non-finite reported SEs --
#'   so a non-rejection is \emph{hollow}; a loud warning is also emitted),
#'   \code{leg_reasons} (the per-leg breakage descriptions; empty when healthy),
#'   \code{few_clusters} (\code{TRUE} when the fits carry fewer than 5 clusters,
#'   so the cluster-robust statistic is unreliable -- a few-cluster artifact,
#'   not PT evidence), and \code{n_clusters}.
#'
#' @references Chen, X., Sant'Anna, P. H. C., & Xie, H. (2025). Efficient
#'   Difference-in-Differences and Event Study Estimators. Section 5.1,
#'   Theorem 5.1. \cr
#'   Hausman, J. A. (1978). Specification Tests in Econometrics.
#'   \emph{Econometrica}, 46(6), 1251-1271. \cr
#'   Andrews, D. W. K. (1987). Asymptotic Results for Generalized Wald Tests.
#'   \emph{Econometric Theory}, 3(3), 348-358. \cr
#'   Bell, R. M., & McCaffrey, D. F. (2002). Bias Reduction in Standard Errors
#'   for Linear Regression with Multi-Stage Samples. \emph{Survey Methodology},
#'   28(2), 169-181. \cr
#'   Pustejovsky, J. E., & Tipton, E. (2018). Small-Sample Methods for
#'   Cluster-Robust Variance Estimation and Hypothesis Testing in Fixed Effects
#'   Models. \emph{Journal of Business & Economic Statistics}, 36(4), 672-683. \cr
#'   Imbens, G. W., & Kolesar, M. (2016). Robust Standard Errors in Small Samples:
#'   Some Practical Advice. \emph{Review of Economics and Statistics}, 98(4), 701-712.
#'
#' @seealso \code{\link{edid}}, \code{\link{edid_sargan}},
#'   \code{\link{edid_frontier}}, \code{\link{edid_adaptive}}
#'
#' @examples
#' \donttest{
#' df <- data.frame(
#'   id   = rep(1:120, each = 6),
#'   time = rep(1:6, 120),
#'   g    = rep(sample(c(3, 5, Inf), 120, replace = TRUE), each = 6)
#' )
#' df$y <- rnorm(120)[df$id] + 0.2 * df$time + 1 * (df$time >= df$g) +
#'   rnorm(nrow(df), 0, 0.5)
#' fit_R <- edid(df, "y", "id", "time", "g", pt_assumption = "all",
#'               aggregate = "event_study", cband = FALSE)
#' fit_U <- edid(df, "y", "id", "time", "g", pt_assumption = "post",
#'               aggregate = "event_study", cband = FALSE)
#' edid_hausman(fit_U, fit_R)
#' }
#'
#' @export
edid_hausman <- function(fit_unrestricted, fit_restricted,
                         parameter = c("event_study", "overall"),
                         e_set = NULL, data = NULL) {
  parameter <- match.arg(parameter)
  .edid_toolkit_check_fits(fit_unrestricted, fit_restricted)

  # The over-identification contrast uses the EFFICIENT plug-in influence function, not the
  # misspecification-robust one (Andrews, Chen & Tecchio 2025, Sec 5): refit BOTH legs in the bare plug-in
  # configuration (all estimation-effect channels off). Point estimates are unchanged (the channels are
  # variance-only), so the contrast d is identical; only its covariance reverts to the efficient one. This
  # makes the test invariant to how the legs were fit (e.g. the covariate default folds psi_Omega). `data`
  # is recovered from the fits' call when not supplied.
  data <- .edid_recover_data(fit_restricted, data, parent.frame())
  fit_unrestricted <- .edid_plugin_refit(fit_unrestricted, data, parent.frame())
  fit_restricted   <- .edid_plugin_refit(fit_restricted,   data, parent.frame())

  n  <- fit_restricted$n
  ci <- fit_restricted$cluster_indices
  n_eff <- .edid_overid_n_eff(fit_restricted)   # Kish ESS for the over-id noise floor (== n unweighted)

  # Broken-leg sanity guard (a hollow non-rejection footgun, Nguyen/Bailey-GB/ACA gate
  # evidence): if EITHER leg is a numerically degenerate fit -- extreme propensity ratios,
  # a non-credible weight channel, cross-cohort hedges carrying the estimand, or non-finite
  # SEs -- the Hausman contrast inherits the breakage and a non-rejection means nothing.
  # Warn loudly, naming the unhealthy leg(s) and reason(s), and flag the result object
  # ($leg_unstable / $leg_reasons) so a hollow p ~ 0.3-0.95 is not mistaken for a pass.
  health_U <- .edid_leg_health(fit_unrestricted, "fit_unrestricted (PT-Post)")
  health_R <- .edid_leg_health(fit_restricted,  "fit_restricted (PT-All)")
  leg_unstable <- health_U$broken || health_R$broken
  leg_reasons  <- character(0L)
  if (leg_unstable) {
    msgs <- c(if (health_R$broken) sprintf("%s: %s", health_R$label, paste(health_R$reasons, collapse = "; ")),
              if (health_U$broken) sprintf("%s: %s", health_U$label, paste(health_U$reasons, collapse = "; ")))
    leg_reasons <- msgs
    warning(paste0(
      "edid_hausman: a constituent fit is numerically degenerate, so this test is HOLLOW -- ",
      "a non-rejection here carries no evidence (the contrast inherits the broken leg's ",
      "instability). ", paste(msgs, collapse = " | "),
      ". Repair the fit (e.g. moment_set = \"own\" to drop cross-cohort pairs, weight_scheme = ",
      "\"averaged\", or a low-dimensional covariate index) before reading the p-value; the result ",
      "is returned with $leg_unstable = TRUE."), call. = FALSE)
  }

  if (parameter == "event_study") {
    e_set <- .edid_shared_e_set(fit_unrestricted, fit_restricted, e_set)
    pU <- .edid_param_ifs(fit_unrestricted, "event_study", e_set)
    pR <- .edid_param_ifs(fit_restricted,  "event_study", e_set)
  } else {
    pU <- .edid_param_ifs(fit_unrestricted, "overall")
    pR <- .edid_param_ifs(fit_restricted,  "overall")
    e_set <- NULL
  }

  d  <- pU$est - pR$est                  # unrestricted minus restricted
  xi <- pU$IF - pR$IF                    # per-unit IF difference (n x |E|)

  # Absolute variance scale of the two estimators (largest per-coordinate
  # asymptotic variance of either fit), for the degenerate-contrast guard.
  v_scale <- max(diag(as.matrix(n * cluster_cov_edid(pU$IF, ci, n))),
                 diag(as.matrix(n * cluster_cov_edid(pR$IF, ci, n))))
  joint <- .edid_if_diff_quadform(d, xi, n, ci, v_scale = v_scale, n_eff = n_eff)
  if (isTRUE(joint$degenerate)) {
    message("edid_hausman: the two estimators coincide (the IF-difference covariance is at ",
            "numerical-noise scale relative to the estimators' own variances), so the joint ",
            "contrast is degenerate and is reported as H = 0, df = 0, p = 1. This is expected ",
            "when the fits carry no over-identifying content to test -- e.g. when the ",
            "thin-cohort guard has pinned every cell to its just-identified moment.")
  }

  # Few-cluster guard: with G < EDID_FEWCLUSTER_MIN clusters the cluster-robust D is too
  # noisy / degenerate for the chi-square reference (the ACA gate's 2-3-state cohorts give
  # H ~ 350, p ~ 1e-73 -- a few-cluster artifact, not PT evidence). Flag the statistic as
  # unreliable (message + field); it is still returned.
  n_clusters <- .edid_n_clusters(fit_restricted)
  few_clusters <- !is.null(ci) && n_clusters < EDID_FEWCLUSTER_MIN
  if (few_clusters) {
    warning(sprintf(paste0(
      "edid_hausman: only %d cluster(s) at the clustering level (clustervars). The cluster-robust ",
      "statistic is UNRELIABLE below %d clusters -- the cluster-level covariance is too noisy / ",
      "degenerate for the chi-square reference, so the p-value is not trustworthy (flagged ",
      "$few_clusters = TRUE). Report unit-level (unclustered) toolkit statistics, or aggregate to ",
      "a coarser level with enough clusters."), n_clusters, EDID_FEWCLUSTER_MIN), call. = FALSE)
  }

  # Scalar eqn (5.5) statistics: each ES(e) plus ES_avg (always included).
  oU <- .edid_param_ifs(fit_unrestricted, "overall")
  oR <- .edid_param_ifs(fit_restricted,  "overall")
  lab_e <- if (parameter == "event_study") pU$e else numeric(0L)
  rows  <- vector("list", length(lab_e) + 1L)
  for (j in seq_along(lab_e)) {
    vR <- as.numeric(n * cluster_cov_edid(pR$IF[, j, drop = FALSE], ci, n))
    sc <- .edid_scalar_hausman(d[j], xi[, j], n, ci, v_scale = vR, n_eff = n_eff)
    rows[[j]] <- data.frame(
      parameter = sprintf("ES(%g)", lab_e[j]), e = lab_e[j],
      theta_U = pU$est[j], theta_R = pR$est[j], difference = d[j],
      D = sc$D, H = sc$H, p_value = sc$p_value, stringsAsFactors = FALSE)
  }
  d_ov  <- oU$est - oR$est
  xi_ov <- oU$IF[, 1L] - oR$IF[, 1L]
  vR_ov <- as.numeric(n * cluster_cov_edid(oR$IF, ci, n))
  sc_ov <- .edid_scalar_hausman(d_ov, xi_ov, n, ci, v_scale = vR_ov, n_eff = n_eff)
  rows[[length(rows)]] <- data.frame(
    parameter = "ES_avg", e = NA_real_,
    theta_U = oU$est, theta_R = oR$est, difference = d_ov,
    D = sc_ov$D, H = sc_ov$H, p_value = sc_ov$p_value, stringsAsFactors = FALSE)
  scalar <- do.call(rbind, rows)
  rownames(scalar) <- NULL

  out <- list(
    statistic  = joint$statistic,
    df         = joint$df,
    p_value    = joint$p_value,
    m_eff      = joint$m_eff,                 # AHT effective df used (G_eff - 1: #clusters or Kish n_eff)
    m_sat      = joint$m_sat,                 # Bell-McCaffrey/Satterthwaite leverage df (fragility diagnostic)
    df2        = joint$df2,                   # F denominator df = m_eff - df + 1 (NA => fell back to chi^2)
    degenerate = isTRUE(joint$degenerate),
    d          = stats::setNames(d, if (parameter == "event_study") sprintf("e=%g", pU$e) else "overall"),
    D          = joint$D,
    scalar     = scalar,
    parameter  = parameter,
    e_set      = e_set,
    n          = n,
    clustered  = !is.null(ci),
    leg_unstable = isTRUE(leg_unstable),   # a constituent fit is numerically degenerate -> hollow test
    leg_reasons  = leg_reasons,            # per-leg breakage descriptions (empty when healthy)
    few_clusters = isTRUE(few_clusters),   # G < EDID_FEWCLUSTER_MIN -> cluster-robust stat unreliable
    n_clusters   = n_clusters,
    alpha      = fit_restricted$alpha %||% 0.05
  )
  class(out) <- c("edid_hausman", "list")
  out
}

#' @describeIn edid_hausman Print method.
#' @param x an \code{edid_hausman} object
#' @param digits number of significant digits to print
#' @param ... ignored
#' @export
print.edid_hausman <- function(x, digits = 4, ...) {
  cat("\nHausman test of PT-All vs PT-Post (Chen, Sant'Anna & Xie 2025, Theorem 5.1)\n")
  param_lab <- if (identical(x$parameter, "event_study")) {
    sprintf("event study, E = {%s}", paste(x$e_set, collapse = ", "))
  } else "overall ES_avg"
  cat(sprintf("  Parameter: %s%s\n", param_lab,
              if (isTRUE(x$clustered)) " (cluster-robust)" else ""))
  ref <- if (!is.null(x$df2) && is.finite(x$df2))
    sprintf("F(%s, %s) [AHT, effective df m = %s]", format(x$df), format(round(x$df2, 1)),
            format(round(x$m_eff, 1)))
  else sprintf("chi^2(%s)", format(x$df))
  cat(sprintf("  H = %s on %s df (rank of D-hat), p-value = %s  [ref: %s]\n",
              format(x$statistic, digits = digits), format(x$df),
              format.pval(x$p_value, digits = digits), ref))
  if (!is.null(x$df2) && !is.finite(x$df2) && !isTRUE(x$degenerate)) {
    cat("  NOTE: the cluster-robust D-hat has effective df <= the test dimension (a few clusters\n")
    cat("  dominate); the AHT F reference is uninformative, so the p-value falls back to chi^2 and\n")
    cat("  is NOT trustworthy (treat like the few-cluster guard).\n")
  } else if (!is.null(x$m_sat) && is.finite(x$m_sat) && !is.null(x$m_eff) && is.finite(x$m_eff) &&
             x$m_sat < 0.5 * max(x$m_eff, x$df + 1)) {
    cat(sprintf("  NOTE: the covariance D-hat is dominated by a few high-leverage units/clusters\n"))
    cat(sprintf("  (Bell-McCaffrey leverage df m_sat = %.1f << effective df %.1f): weak overlap / severe\n",
                x$m_sat, x$m_eff))
    cat("  imbalance, so even the F reference is fragile. Consider overlap trimming (trim_level) and\n")
    cat("  read the localized Sargan rather than the diffuse joint statistic.\n")
  }
  if (isTRUE(x$degenerate)) {
    cat("  NOTE: degenerate contrast -- the two estimators coincide, so there is no\n")
    cat("  over-identifying content to test (H = 0, df = 0, p = 1 by construction).\n")
  }
  if (isTRUE(x$leg_unstable)) {
    cat("  WARNING: HOLLOW TEST -- a constituent fit is numerically degenerate, so a\n")
    cat("  non-rejection carries no evidence. Reason(s):\n")
    for (r in x$leg_reasons) cat("    - ", r, "\n", sep = "")
  }
  if (isTRUE(x$few_clusters)) {
    cat(sprintf("  WARNING: only %d cluster(s) -- the cluster-robust statistic is unreliable\n", x$n_clusters))
    cat("  below ", EDID_FEWCLUSTER_MIN, " clusters (few-cluster artifact, not PT evidence).\n", sep = "")
  }
  cat("\nPer-parameter scalar statistics (eqn 5.5):\n")
  tab <- x$scalar
  num <- vapply(tab, is.numeric, logical(1L))
  tab[num] <- lapply(tab[num], function(z) signif(z, digits))
  print(tab, row.names = FALSE)
  cat("\nH0: parallel trends holds across all groups and pre-treatment periods (PT-All).\n")
  invisible(x)
}
