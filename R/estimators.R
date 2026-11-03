## =============================================================================
## estimators.R
## Hand-coded ATE estimators that SHARE one Super Learner nuisance library, so that
## differences between them come from the combination rule and the fitting scheme (and,
## for split-based estimators, from the split), not from separate implementations.
##
## Estimators implemented (all for the average treatment effect, ATE = E[Y^1 - Y^0]):
##   - g-computation (outcome-regression plug-in)      : singly robust
##   - IPTW (Horvitz-Thompson / Hajek)                 : singly robust
##   - AIPW  (augmented IPW = one-step estimator)      : doubly robust
##   - TMLE  (targeted maximum likelihood)             : doubly robust
## Each doubly-robust estimator is run with three nuisance strategies:
##   (F)  full-sample  : nuisances fit on all data, predicted in-sample  -> "classic" AIPW/TMLE
##   (X)  cross-fit    : out-of-fold nuisances (V folds)                 -> DML (AIPW) / CV-TMLE
##   (D)  double-cross-fit : g and Q fit on DISJOINT folds               -> DC-AIPW / DC-TMLE
##
## THEORY SPINE:
## AIPW, TMLE and DML are three routes to the SAME efficient influence function (EIF)
## for the ATE under the nonparametric model:
##
##   D(O) = [ A/g(X) - (1-A)/(1-g(X)) ] * ( Y - Q_A(X) )  +  ( Q_1(X) - Q_0(X) )  - psi
##            \___________________________/   \_________/     \_________________/
##             clever covariate H(A,X)         residual        outcome contrast
##
##   * AIPW/one-step  : psi_hat = psi_plugin + mean( EIF evaluated at nuisances )  (solve E[D]=0 by adding the mean)
##   * TMLE           : fluctuate Q so that the empirical mean of the H*(Y-Q) part is 0, then plug in
##   * DML            : identical AIPW score, but nuisances are estimated OUT-OF-FOLD (Neyman-orthogonal
##                      moment + cross-fitting => remainder is second-order, valid CIs with ML nuisances)
## All three are asymptotically equivalent (Kennedy 2022/2023; Hines et al. 2022; Diaz 2020).
## =============================================================================

suppressMessages({
  library(SuperLearner)
})

## ---- 0. small helpers -------------------------------------------------------
.bound <- function(x, lo, hi) pmin(pmax(x, lo), hi)

## how much propensity truncation actually did: range before clipping, number clipped
.clip_report <- function(p, gbound) c(raw_min = min(p), raw_max = max(p),
                                      n_clipped = sum(p < gbound[1] | p > gbound[2]))

## default Super Learner library (shared by g and Q for a FAIR comparison)
default_SL <- c("SL.mean", "SL.glm", "SL.glmnet", "SL.ranger", "SL.earth")

make_folds <- function(n, V = 5, seed = 1L) {
  set.seed(seed)
  sample(rep(seq_len(V), length.out = n))
}

## Fit one Super Learner and return a predictor closure.
.fit_sl <- function(Y, X, family, SL.library) {
  fit <- SuperLearner(Y = Y, X = X, family = family, SL.library = SL.library,
                      cvControl = list(V = 5L))
  function(newX) as.numeric(predict(fit, newX, onlySL = TRUE)$pred)
}

## ---- 1. nuisance estimation strategies --------------------------------------
## Every strategy returns the SAME object: out-of-(the-relevant)-sample predictions
##   pi = g(X) = P(A=1|X),  Q1 = E[Y|A=1,X],  Q0 = E[Y|A=0,X],  QA = Q at observed A
## so the downstream estimators are identical regardless of nuisance strategy.

## (F) full-sample: fit on all data, predict in-sample  -> classic AIPW / TMLE
nuisance_full <- function(Y, A, X, SL.library = default_SL,
                          family_y = gaussian(), gbound = c(.025, .975)) {
  Xdf <- as.data.frame(X)
  g   <- .fit_sl(A, Xdf, binomial(), SL.library)
  q   <- .fit_sl(Y, cbind(A = A, Xdf), family_y, SL.library)
  praw <- g(Xdf)
  pi  <- .bound(praw, gbound[1], gbound[2])
  Q1  <- q(cbind(A = 1, Xdf)); Q0 <- q(cbind(A = 0, Xdf))
  list(pi = pi, Q1 = Q1, Q0 = Q0, QA = ifelse(A == 1, Q1, Q0), clip = .clip_report(praw, gbound))
}

## (X) cross-fit (DML / CV-TMLE): V folds; g and Q trained on the SAME V-1 folds,
## predicted on the held-out fold => out-of-fold nuisances for every observation.
nuisance_crossfit <- function(Y, A, X, V = 5L, SL.library = default_SL,
                              family_y = gaussian(), gbound = c(.025, .975), seed = 1L) {
  Xdf <- as.data.frame(X); n <- length(Y); fold <- make_folds(n, V, seed)
  pi <- Q1 <- Q0 <- rep(NA_real_, n)
  for (v in seq_len(V)) {
    tr <- fold != v; te <- fold == v
    g <- .fit_sl(A[tr], Xdf[tr, , drop = FALSE], binomial(), SL.library)
    q <- .fit_sl(Y[tr], cbind(A = A[tr], Xdf[tr, , drop = FALSE]), family_y, SL.library)
    pi[te] <- g(Xdf[te, , drop = FALSE])
    Q1[te] <- q(cbind(A = 1, Xdf[te, , drop = FALSE]))
    Q0[te] <- q(cbind(A = 0, Xdf[te, , drop = FALSE]))
  }
  clip <- .clip_report(pi, gbound)
  pi <- .bound(pi, gbound[1], gbound[2])
  list(pi = pi, Q1 = Q1, Q0 = Q0, QA = ifelse(A == 1, Q1, Q0), fold = fold, clip = clip)
}

## (D) double cross-fit (DC-AIPW / DC-TMLE; Newey & Robins 2018; Zivich & Breskin 2021):
## split into V = 3 folds and use ONE fixed rotation of roles, as in Zivich & Breskin's
## reference code: observations in fold v are scored with g fit on fold v+1 alone and Q
## fit on fold v+2 alone (indices mod 3). Every scored observation is excluded from both
## training sets, and the g and Q used for it never share training data. Averaging the
## predictions of several role assignments would break this: each averaged g and Q would
## then use both training folds.
nuisance_double <- function(Y, A, X, V = 3L, SL.library = default_SL,
                            family_y = gaussian(), gbound = c(.025, .975), seed = 1L) {
  stopifnot(V == 3L)
  Xdf <- as.data.frame(X); n <- length(Y); fold <- make_folds(n, V, seed)
  pi <- Q1 <- Q0 <- rep(NA_real_, n)
  for (v in seq_len(V)) {
    te   <- fold == v
    tr_g <- fold == (v %% V) + 1L            # the next fold trains g
    tr_q <- fold == ((v + 1L) %% V) + 1L     # the fold after that trains Q
    g <- .fit_sl(A[tr_g], Xdf[tr_g, , drop = FALSE], binomial(), SL.library)
    q <- .fit_sl(Y[tr_q], cbind(A = A[tr_q], Xdf[tr_q, , drop = FALSE]), family_y, SL.library)
    pi[te] <- g(Xdf[te, , drop = FALSE])
    Q1[te] <- q(cbind(A = 1, Xdf[te, , drop = FALSE]))
    Q0[te] <- q(cbind(A = 0, Xdf[te, , drop = FALSE]))
  }
  clip <- .clip_report(pi, gbound)
  pi <- .bound(pi, gbound[1], gbound[2])
  list(pi = pi, Q1 = Q1, Q0 = Q0, QA = ifelse(A == 1, Q1, Q0), fold = fold, clip = clip)
}

## ---- 2. estimators (agnostic to which nuisance strategy produced `nu`) -------

## g-computation: average of the individual outcome contrasts (plug-in).
est_gcomp <- function(Y, A, nu) {
  psi <- mean(nu$Q1 - nu$Q0)
  list(estimate = psi, se = NA_real_, ci = c(NA, NA), method = "g-computation")
}

## IPTW (Hajek / stabilized): reweight observed outcomes by inverse propensity.
est_iptw <- function(Y, A, nu) {
  pi <- nu$pi
  m1 <- sum(A * Y / pi) / sum(A / pi)
  m0 <- sum((1 - A) * Y / (1 - pi)) / sum((1 - A) / (1 - pi))
  list(estimate = m1 - m0, se = NA_real_, ci = c(NA, NA), method = "IPTW (Hajek)")
}

## AIPW / one-step = plug-in + mean of the EIF. With cross-fit nu this IS DML.
## SE from the influence curve (no bootstrap needed).
est_aipw <- function(Y, A, nu) {
  pi <- nu$pi; Q1 <- nu$Q1; Q0 <- nu$Q0
  psi_i <- (Q1 - Q0) + (A / pi) * (Y - Q1) - ((1 - A) / (1 - pi)) * (Y - Q0)  # EIF + psi
  est <- mean(psi_i)
  ic  <- psi_i - est                      # influence curve
  se  <- sqrt(var(ic) / length(Y))
  list(estimate = est, se = se, ci = est + c(-1, 1) * 1.96 * se,
       ic = ic, method = "AIPW / one-step")
}

## TMLE: targeted logistic fluctuation of Q via the clever covariate H(A,X).
## Continuous Y is scaled to [0,1] (Gruber & van der Laan) so the logistic update is valid;
## binary Y uses y_range = c(0,1) (the scaling is then the identity).
est_tmle <- function(Y, A, nu, y_range = NULL) {
  pi <- nu$pi
  if (is.null(y_range)) y_range <- range(Y)
  a <- y_range[1]; b <- y_range[2]; span <- b - a
  if (!(span > 0)) stop("est_tmle: the outcome range is zero, so it cannot be rescaled to [0,1]")
  Ys  <- (Y      - a) / span
  Q1s <- .bound((nu$Q1 - a) / span, 1e-5, 1 - 1e-5)
  Q0s <- .bound((nu$Q0 - a) / span, 1e-5, 1 - 1e-5)
  QAs <- ifelse(A == 1, Q1s, Q0s)
  ## clever covariate and its values under A=1 / A=0
  H  <- A / pi - (1 - A) / (1 - pi)
  H1 <- 1 / pi; H0 <- -1 / (1 - pi)
  ## fluctuation parameter eps: logistic reg of Ys on H with offset logit(QAs), no intercept
  ## (warnings about non-integer "successes" are expected for a [0,1]-scaled outcome;
  ##  convergence and the empirical score residual are returned below and checked)
  fit <- suppressWarnings(glm(Ys ~ -1 + H, offset = qlogis(QAs), family = binomial()))
  eps <- coef(fit)
  ## update the initial fits
  Q1s_star <- plogis(qlogis(Q1s) + eps * H1)
  Q0s_star <- plogis(qlogis(Q0s) + eps * H0)
  QAs_star <- plogis(qlogis(QAs) + eps * H)
  ## back to original scale, take the mean contrast
  Q1_star <- Q1s_star * span + a; Q0_star <- Q0s_star * span + a
  est <- mean(Q1_star - Q0_star)
  ## influence curve on the original scale (note the span factor from the residual term)
  ic <- span * H * (Ys - QAs_star) + (Q1_star - Q0_star) - est
  se <- sqrt(var(ic) / length(Y))
  list(estimate = est, se = se, ci = est + c(-1, 1) * 1.96 * se,
       ic = ic, eps = eps, converged = isTRUE(fit$converged),
       score_resid = mean(H * (Ys - QAs_star)), method = "TMLE")
}

## ---- 3. convenience: run the whole panel on one dataset ---------------------
## Returns a tidy data.frame of estimator x nuisance-strategy with estimate, SE, CI.
run_panel <- function(Y, A, X, family_y = gaussian(), y_range = NULL,
                      V = 5L, Vd = 3L, SL.library = default_SL, seed = 1L,
                      truth = NA_real_) {
  binary_y <- isTRUE(all.equal(sort(unique(Y)), c(0, 1)))
  if (binary_y) y_range <- c(0, 1)
  nu_F <- nuisance_full(Y, A, X, SL.library, family_y)
  nu_X <- nuisance_crossfit(Y, A, X, V,  SL.library, family_y, seed = seed)
  nu_D <- nuisance_double(Y, A, X, Vd, SL.library, family_y, seed = seed)

  mk <- function(res, label, strat) data.frame(
    estimator = label, nuisance = strat,
    estimate = res$estimate, se = res$se,
    ci_lo = res$ci[1], ci_hi = res$ci[2],
    bias = res$estimate - truth, stringsAsFactors = FALSE)

  out <- rbind(
    mk(est_gcomp(Y, A, nu_F), "g-computation", "full"),
    mk(est_iptw (Y, A, nu_F), "IPTW",          "full"),
    mk(est_aipw (Y, A, nu_F), "AIPW",          "full"),
    mk(est_tmle (Y, A, nu_F, y_range), "TMLE",  "full"),
    mk(est_aipw (Y, A, nu_X), "AIPW = DML",    "cross-fit"),
    mk(est_tmle (Y, A, nu_X, y_range), "TMLE (CV-TMLE)", "cross-fit"),
    mk(est_aipw (Y, A, nu_D), "AIPW (DC-AIPW)", "double-cross-fit"),
    mk(est_tmle (Y, A, nu_D, y_range), "TMLE (DC-TMLE)", "double-cross-fit")
  )
  attr(out, "truth") <- truth
  out
}
