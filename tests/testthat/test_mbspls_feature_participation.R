feature_toy = function(n = 80L, seed = 1L, signal = 1) {
  set.seed(seed)
  f = rnorm(n)
  X = list(
    A = cbind(a1 = signal * f + rnorm(n), a2 = rnorm(n)),
    B = cbind(b1 = signal * f + rnorm(n), b2 = rnorm(n)),
    C = cbind(c1 = signal * f + rnorm(n))
  )
  W = list(A = c(1, 0), B = c(0.8, 0.6), C = 1)
  T_mat = vapply(names(X), function(b) drop(X[[b]] %*% W[[b]]), numeric(n))
  list(X = X, T = T_mat)
}

test_that("mbspls_feature_heldout_tests computes the stated statistic", {
  tr = feature_toy(seed = 1L)
  te = feature_toy(seed = 2L)
  res = mbspls_feature_heldout_tests(te$X, tr$X, te$T, tr$T, n_perm = 199L)
  expect_identical(names(res), c("block", "feature", "r", "sign", "train_r", "p"))
  expect_identical(res$feature, c("a1", "a2", "b1", "b2", "c1"))
  expect_true(all(res$p >= 1 / 200 & res$p <= 1))

  # By hand for feature a1: composite of the held-out scores of B and C, signed by training correlations.
  sg = sign(cor(tr$T))
  zs = function(v) (v - mean(v)) / sd(v)
  comp = zs(sg["A", "B"] * zs(te$T[, "B"]) + sg["A", "C"] * zs(te$T[, "C"]))
  comp_tr = drop(scale(tr$T)[, c("B", "C")] %*% sg["A", c("B", "C")])
  s1 = sign(cor(tr$X$A[, "a1"], comp_tr))
  expect_equal(res[feature == "a1", r], s1 * cor(te$X$A[, "a1"], comp))
  expect_equal(res[feature == "a1", train_r], cor(tr$X$A[, "a1"], comp_tr))
  # Signal features are detected, noise features are not.
  expect_true(all(res[feature %in% c("a1", "b1", "c1"), p] <= 0.01))

  # Cross-modality composite: A and B share a modality, so the composite of A is C alone.
  cm = mbspls_feature_heldout_tests(te$X, tr$X, te$T, tr$T, n_perm = 99L,
    modality = c(A = "m1", B = "m1", C = "m2"))
  s1c = sign(cor(tr$X$A[, "a1"], scale(tr$T)[, "C"] * sg["A", "C"]))
  expect_equal(cm[feature == "a1", r], s1c * cor(te$X$A[, "a1"], zs(te$T[, "C"])))

  # Strata: the function runs and returns valid p-values.
  st = mbspls_feature_heldout_tests(te$X, tr$X, te$T, tr$T, n_perm = 99L, strata = rep(1:2, 40))
  expect_true(all(st$p > 0 & st$p <= 1))

  # A constant held-out feature gets r = NA and p = 1.
  te2 = te
  te2$X$B[, "b2"] = 1
  cst = mbspls_feature_heldout_tests(te2$X, tr$X, te$T, tr$T, n_perm = 99L)
  expect_true(is.na(cst[feature == "b2", r]))
  expect_identical(cst[feature == "b2", p], 1)

  expect_error(mbspls_feature_heldout_tests(te$X, tr$X, te$T[, 1:2], tr$T, n_perm = 9L))
  bad = te$T
  colnames(bad) = c("A", "B", "Z")
  expect_error(mbspls_feature_heldout_tests(te$X, tr$X, bad, bad, n_perm = 9L), "same named columns")
  expect_error(mbspls_feature_heldout_tests(te$X, tr$X, te$T, tr$T, n_perm = 9L, modality = c(A = "m")), "every block")
})

test_that("mbspls_feature_heldout_tests checks that blocks and columns match", {
  tr = feature_toy(seed = 1L)
  te = feature_toy(seed = 2L)
  # Swapped columns would give each feature the orientation of the other.
  swapped = te$X
  swapped$A = swapped$A[, c("a2", "a1")]
  expect_error(mbspls_feature_heldout_tests(swapped, tr$X, te$T, tr$T, n_perm = 9L), "block 'A'")
  renamed = te$X
  colnames(renamed$B) = c("b1", "other")
  expect_error(mbspls_feature_heldout_tests(renamed, tr$X, te$T, tr$T, n_perm = 9L), "block 'B'")
  # Blocks in a different order.
  expect_error(
    mbspls_feature_heldout_tests(te$X, tr$X[c("B", "A", "C")], te$T, tr$T, n_perm = 9L),
    "same blocks in the same order"
  )
  expect_error(
    mbspls_feature_heldout_tests(te$X[c("A", "B")], tr$X, te$T, tr$T, n_perm = 9L),
    "same blocks in the same order"
  )
  # Without column names the features are numbered and the check is skipped.
  strip = function(X) lapply(X, unname)
  un = mbspls_feature_heldout_tests(strip(te$X), strip(tr$X), te$T, tr$T, n_perm = 9L)
  expect_identical(un$feature, c("A_1", "A_2", "B_1", "B_2", "C_1"))
  named = mbspls_feature_heldout_tests(te$X, tr$X, te$T, tr$T, n_perm = 9L)
  expect_equal(un$r, named$r)
})

test_that("mbspls_feature_heldout_tests is valid under the null", {
  skip_on_cran()
  p = vapply(seq_len(150), function(i) {
    tr = feature_toy(n = 60L, seed = 1000L + i)
    te = feature_toy(n = 40L, seed = 5000L + i)
    te$X$A[, "a2"] = rnorm(40) # independent of the other blocks
    mbspls_feature_heldout_tests(te$X, tr$X, te$T, tr$T, n_perm = 199L)[feature == "a2", p]
  }, numeric(1L))
  expect_gt(mean(p), 0.4)
  expect_lt(mean(p), 0.6)
  expect_lt(mean(p <= 0.05), 0.1)
})

test_that("as_mbspls_feature_participation counts replication and controls the FDR", {
  centres = data.frame(
    site = rep(c("a", "b", "c", "d"), each = 3), block = rep(c("X", "X", "Y"), 4),
    feature = rep(c("x1", "x2", "y1"), 4), n = rep(c(50, 40, 30, 20), each = 3),
    r = c(0.40, 0.02, 0.35, 0.38, -0.05, 0.30, 0.33, 0.04, 0.28, 0.31, 0.01, 0.2),
    p = c(0.0001, 0.45, 0.002, 0.0004, 0.62, 0.03, 0.004, 0.4, 0.06, 0.001, 0.7, 0.2)
  )
  ref = data.frame(block = c("X", "X", "Y"), feature = c("x1", "x2", "y1"), weight = c(-0.9, 0, 1))
  res = as_mbspls_feature_participation(centres, reference = ref)
  expect_s3_class(res, "mbspls_feature_participation")
  expect_identical(res$r, 2L)
  f = res$features
  # x1: p sorted 0.0001, 0.0004, 0.001, 0.004 -> partial conjunction at r = 2: max(4e-4, 3 * 4e-4) = 0.0012.
  expect_equal(f[feature == "x1", p_r], 0.0012)
  expect_identical(f[feature == "x1", r_hat], 4L)
  expect_equal(f$q_r, p.adjust(f$p_r, method = "BH"))
  expect_identical(f$replicated, f$q_r <= 0.05)
  expect_identical(f[feature == "x2", replicated], FALSE)
  expect_equal(f[feature == "x1", heldout_r], sum(c(50, 40, 30, 20) * c(0.40, 0.38, 0.33, 0.31)) / 140)
  expect_identical(f$direction, c(-1, NA, 1))
  expect_identical(f$weight, c(-0.9, 0, 1))
  expect_output(print(res), "Replicated at 2 or more centres")

  by = as_mbspls_feature_participation(centres, r = 3L, fdr = "BY")
  expect_equal(by$features$q_r, p.adjust(by$features$p_r, method = "BY"))
  none = as_mbspls_feature_participation(centres[centres$site != "d" | centres$feature != "y1", ], fdr = "none")
  expect_identical(nrow(none$centres), 12L)
  expect_identical(none$centres[site == "d" & feature == "y1", p], 1)
  expect_error(as_mbspls_feature_participation(centres[, names(centres) != "p"]), "'p'")
  expect_error(as_mbspls_feature_participation(centres, r = 5L))
})

test_that("autoplot.mbspls_feature_participation draws a signature", {
  centres = data.frame(
    site = rep(c("a", "b", "c"), each = 3), block = rep(c("X", "X", "Y"), 3),
    feature = rep(c("x1", "x2", "y1"), 3), n = rep(c(50, 40, 30), each = 3),
    r = c(0.40, 0.02, 0.35, 0.38, -0.05, 0.30, 0.33, 0.04, 0.28),
    p = c(0.0001, 0.45, 0.002, 0.0004, 0.62, 0.03, 0.04, 0.4, 0.06)
  )
  ref = data.frame(block = c("X", "X", "Y"), feature = c("x1", "x2", "y1"), weight = c(-0.9, 0, 1))
  res = as_mbspls_feature_participation(centres, reference = ref, r = 2L)
  p = ggplot2::autoplot(res, feature_labels = c(x1 = "First"), unit = c("site", "sites"))
  expect_s3_class(p, "ggplot")
  pts = ggplot2::layer_data(p, 3L)
  expect_identical(nrow(pts), 3L)
  # The signed held-out correlation uses the reference direction.
  expect_true(any(abs(pts$x + sum(c(50, 40, 30) * c(0.40, 0.38, 0.33)) / 120) < 1e-12))
  expect_true(any(grepl("sites", levels(p$data$status))))
  expect_identical(nrow(ggplot2::layer_data(ggplot2::autoplot(res, top = 1L), 3L)), 2L)
  p2 = ggplot2::autoplot(res, compare = list(Second = res), label = "First")
  expect_identical(levels(p2$data$analysis), c("First", "Second"))
  expect_identical(nrow(ggplot2::layer_data(p2, 3L)), 6L)
  expect_error(ggplot2::autoplot(res, compare = list(res)), "named list")
})
