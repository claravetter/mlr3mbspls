participation_toy = function(n_site = 40L, sites = c("s1", "s2", "s3", "s4"), seed = 1L) {
  set.seed(seed)
  n = n_site * length(sites)
  f = rnorm(n)
  e = function() 0.5 * rnorm(n)
  df = data.frame(
    x1 = f + e(), x2 = f + e(), x3 = rnorm(n),
    y1 = f + e(), y2 = f + e(),
    u1 = rnorm(n), u2 = rnorm(n), u3 = rnorm(n),
    site = rep(sites, each = n_site)
  )
  # Centre effects that shift a feature's mean within one block.
  df$x3 = df$x3 + rep(c(0, 3, -2, 5)[seq_along(sites)], each = n_site)
  task = mlr3::as_task_regr(cbind(df, target = rnorm(n)), target = "target", id = "participation_toy")
  task$col_roles$feature = setdiff(task$col_roles$feature, "site")
  list(task = task, blocks = list(X = c("x1", "x2", "x3"), Y = c("y1", "y2"), U = c("u1", "u2", "u3")))
}

participation_learner = function(blocks) {
  mbspls_graph_learner(
    learner = mlr3::lrn("regr.featureless"),
    blocks = blocks,
    ncomp = 1L,
    c_matrix = matrix(c(sqrt(3), sqrt(2), sqrt(3)), ncol = 1L, dimnames = list(names(blocks), NULL)),
    impute_within_blocks = TRUE,
    bootstrap = FALSE,
    bootstrap_selection = FALSE
  )
}

test_that("mbspls_partial_conjunction matches the Bonferroni partial conjunction", {
  p = cbind(A = c(0.02, 0.001, 0.3, 0.004), B = c(0.2, 0.5, 0.04, 0.9))
  res = mbspls_partial_conjunction(p)
  expect_equal(unname(res$p_pc[, "A"]), c(0.004, 0.012, 0.04, 0.3))
  expect_equal(unname(res$p_pc[, "B"]), c(0.16, 0.6, 1, 1))
  expect_identical(res$r_hat, c(A = 3L, B = 0L))
  expect_equal(res$level, 0.05)

  sim = mbspls_partial_conjunction(p, simultaneous = TRUE)
  expect_equal(sim$level, 0.025)
  expect_identical(sim$r_hat, c(A = 2L, B = 0L))

  # A vector is one block; the curve is monotone.
  v = mbspls_partial_conjunction(c(0.03, 0.01, 0.02), alpha = 0.05)
  expect_equal(unname(v$p_pc[, 1L]), c(0.03, 0.04, 0.04))
  expect_identical(unname(v$r_hat), 3L)
  expect_identical(unname(mbspls_partial_conjunction(c(0.03, 0.01, 0.02), alpha = 0.035)$r_hat), 1L)
  expect_true(all(diff(v$p_pc[, 1L]) >= 0))
  expect_error(mbspls_partial_conjunction(c(0.1, NA)))
  expect_error(mbspls_partial_conjunction(c(0.1, 1.2)), "p-values")
})

test_that("mbspls_partial_conjunction counts p-values below the fixed threshold", {
  p = cbind(A = c(0.02, 0.001, 0.3, 0.004), B = c(0.2, 0.5, 0.04, 0.9))
  res = mbspls_partial_conjunction(p)
  # The existing elements come first and are unchanged.
  expect_identical(names(res), c("p_pc", "r_hat", "level", "r_count"))
  # Threshold 0.05 / 4 = 0.0125: 0.001 and 0.004 of A.
  expect_identical(res$r_count, c(A = 2L, B = 0L))
  # Simultaneous: 0.025 / 4 = 0.00625.
  expect_identical(mbspls_partial_conjunction(p, simultaneous = TRUE)$r_count, c(A = 2L, B = 0L))
  expect_identical(unname(mbspls_partial_conjunction(c(0.03, 0.01, 0.02))$r_count), 1L)

  # The count is the number of p-values at most level / S, never above r_hat,
  # and agrees with r_hat on the any-centre test.
  set.seed(3)
  pm = matrix(stats::runif(10L * 200L)^4, 10L)
  for (sim in c(FALSE, TRUE)) {
    rr = mbspls_partial_conjunction(pm, simultaneous = sim)
    expect_identical(unname(rr$r_count), as.integer(colSums(pm <= rr$level / 10)))
    expect_true(all(rr$r_count <= rr$r_hat))
    expect_identical(rr$r_count >= 1L, rr$r_hat >= 1L)
  }
})

test_that("the count stays valid when the hypotheses depend on the data", {
  # Two centres. The rule tested at one centre is trained on the other and
  # chooses a strong signal when the other centre's noise p-value is below
  # 0.05, and the noise otherwise. A centre is non-null if it tests the signal.
  p_two = cbind(c(0.03, 1e-9))
  expect_identical(unname(mbspls_partial_conjunction(p_two)$r_hat), 2L)
  expect_identical(unname(mbspls_partial_conjunction(p_two)$r_count), 1L)

  set.seed(11)
  R = 20000L
  q1 = stats::runif(R)
  q2 = stats::runif(R)
  null2 = q1 >= 0.05
  null1 = q2 >= 0.05
  p1 = ifelse(null1, q1, 1e-12)
  p2 = ifelse(null2, q2, 1e-12)
  N = (!null1) + (!null2)
  res = mbspls_partial_conjunction(rbind(p1, p2))
  # P(r_hat > N) is 2 * 0.05 * 0.95 = 0.095; P(r_count > N) is at most 0.05.
  expect_gt(mean(res$r_hat > N), 0.085)
  expect_lt(mean(res$r_count > N), 0.055)
})

test_that("mbspls_participation detects participating blocks and returns valid p-values", {
  testthat::skip_if_not("regr.featureless" %in% mlr3::mlr_learners$keys())
  toy = participation_toy()
  gl = participation_learner(toy$blocks)
  res = mbspls_participation(toy$task, gl, site = "site", n_perm = 199L, seed = 3L)

  expect_s3_class(res, "mbspls_participation")
  expect_identical(names(res$centres), c("site", "n", "block", "r", "train_r", "p"))
  expect_identical(nrow(res$centres), 12L)
  expect_identical(res$blocks$block, c("X", "Y", "U"))
  expect_true(all(res$centres$p >= 1 / 200 & res$centres$p <= 1))
  expect_identical(dim(res$partial_conjunction), c(4L, 3L))

  # X and Y share a strong factor at every centre; U is independent noise.
  expect_true(all(res$centres[block %in% c("X", "Y"), p] <= 0.05))
  expect_identical(res$blocks[block == "X", r_hat], 4L)
  expect_identical(res$blocks[block == "Y", r_hat], 4L)
  expect_identical(res$blocks[block == "U", r_hat_simultaneous], 0L)
  expect_true(all(res$blocks$r_hat_simultaneous <= res$blocks$r_hat))
  expect_output(print(res), "lower 95% confidence bound")

  # Fixed-threshold counts: S = 4 centres, K = 3 blocks.
  expect_identical(
    names(res$blocks),
    c("block", "S", "median_r", "n_significant", "r_hat", "r_hat_simultaneous", "r_count", "r_count_simultaneous")
  )
  for (b in c("X", "Y", "U")) {
    pb = res$centres[block == b, p]
    expect_identical(res$blocks[block == b, r_count], sum(pb <= 0.05 / 4))
    expect_identical(res$blocks[block == b, r_count_simultaneous], sum(pb <= 0.05 / 12))
  }
  expect_true(all(res$blocks$r_count <= res$blocks$r_hat))
  expect_true(all(res$blocks$r_count_simultaneous <= res$blocks$r_hat_simultaneous))
  expect_output(print(res), "r_count: centres with p <= 0.0125")
})

test_that("mbspls_participation computes the stated statistic", {
  testthat::skip_if_not("regr.featureless" %in% mlr3::mlr_learners$keys())
  toy = participation_toy()
  gl = participation_learner(toy$blocks)
  res = mbspls_participation(toy$task, gl, site = "site", n_perm = 99L, seed = 2L)

  # Recompute centre s2 by hand.
  site = toy$task$backend$data(rows = toy$task$row_ids, cols = "site")$site
  tr = toy$task$row_ids[site != "s2"]
  te = toy$task$row_ids[site == "s2"]
  fit = gl$clone(deep = TRUE)$train(toy$task, row_ids = tr)
  scores = function(rows) {
    ev = mbspls_eval_new_data(fit, toy$task$clone(deep = TRUE)$filter(rows))
    m = ev$T_mat[, paste0("LV1_", ev$blocks)]
    colnames(m) = ev$blocks
    m
  }
  t_tr = scores(tr)
  t_te = scores(te)
  sg = sign(cor(t_tr))
  z = scale(t_te)
  for (k in colnames(z)) {
    others = setdiff(colnames(z), k)
    r = cor(z[, k], z[, others] %*% sg[k, others])
    expect_equal(res$centres[site == "s2" & block == k, r], drop(r), tolerance = 1e-10)
  }
})

test_that("mbspls_participation accepts a fitting function and is reproducible", {
  testthat::skip_if_not("regr.featureless" %in% mlr3::mlr_learners$keys())
  toy = participation_toy(n_site = 30L, sites = c("a", "b", "c"))
  gl = participation_learner(toy$blocks)
  fit_fun = function(task, row_ids) gl$clone(deep = TRUE)$train(task, row_ids = row_ids)
  site = toy$task$backend$data(rows = toy$task$row_ids, cols = "site")$site

  r1 = mbspls_participation(toy$task, gl, site = "site", n_perm = 99L, seed = 11L)
  r2 = mbspls_participation(toy$task, fit_fun, site = site, n_perm = 99L, seed = 11L)
  expect_equal(r1$centres, r2$centres)
  expect_equal(r1$blocks, r2$blocks)

  # Stratified permutations run and stay valid.
  strata = rep(c("g1", "g2"), length.out = length(site))
  r3 = mbspls_participation(toy$task, gl, site = "site", strata = strata, n_perm = 99L, seed = 11L)
  expect_true(r3$strata)
  expect_true(all(r3$centres$p >= 1 / 100 & r3$centres$p <= 1))
  expect_equal(r3$centres$r, r1$centres$r)

  expect_error(mbspls_participation(toy$task, gl, site = "nope", n_perm = 9L), "not in the task")
  expect_error(mbspls_participation(toy$task, gl, site = rep("a", 90L), n_perm = 9L), "two centres")
  expect_error(mbspls_participation(toy$task, mlr3::lrn("regr.featureless"), site = "site"), "GraphLearner")
})

test_that("mbspls_participation is valid under the null", {
  testthat::skip_on_cran()
  testthat::skip_if_not("regr.featureless" %in% mlr3::mlr_learners$keys())
  # Three independent blocks: every p-value is exact; their mean is about 0.5.
  p = unlist(lapply(1:6, function(rep) {
    set.seed(100 + rep)
    n = 90
    df = data.frame(matrix(rnorm(n * 6), n), site = rep(c("a", "b", "c"), each = 30), target = rnorm(n))
    task = mlr3::as_task_regr(df, target = "target")
    task$col_roles$feature = setdiff(task$col_roles$feature, "site")
    blocks = list(A = c("X1", "X2"), B = c("X3", "X4"), C = c("X5", "X6"))
    gl = mbspls_graph_learner(
      learner = mlr3::lrn("regr.featureless"), blocks = blocks, ncomp = 1L,
      c_matrix = matrix(sqrt(2), 3L, 1L), impute_within_blocks = TRUE,
      bootstrap = FALSE, bootstrap_selection = FALSE
    )
    mbspls_participation(task, gl, site = "site", n_perm = 99L, seed = rep)$centres$p
  }))
  expect_length(p, 54L)
  expect_gt(mean(p), 0.35)
  expect_lt(mean(p), 0.65)
})

test_that("block-wise k-NN imputation uses only the block's own features", {
  set.seed(5)
  n = 60
  df = data.frame(a1 = rnorm(n), a2 = rnorm(n), a3 = rnorm(n), b1 = rnorm(n), b2 = rnorm(n), target = rnorm(n))
  df$a2[c(3, 7, 11)] = NA
  df[c(20, 21), c("a1", "a2", "a3")] = NA # whole block missing
  task = mlr3::as_task_regr(df, target = "target")
  blocks = list(A = c("a1", "a2", "a3"), B = c("b1", "b2"))

  impute = function(d, graph) {
    graph = graph$clone(deep = TRUE)
    graph$train(mlr3::as_task_regr(d, target = "target"))[[1L]]$data()
  }
  g_block = impute_knn_graph(k = 5, blocks = blocks)
  expect_true(all(c("impute_num_knn_block_A", "impute_num_knn_block_B", "impute_num_knn_unblocked") %in% g_block$ids()))

  d2 = df
  d2[c(3, 7, 11, 20, 21), c("b1", "b2")] = d2[c(3, 7, 11, 20, 21), c("b1", "b2")] + 10
  out1 = impute(df, g_block)
  out2 = impute(d2, g_block)
  expect_false(anyNA(out1))
  expect_equal(out1$a2, out2$a2)
  expect_equal(out1$a1[20:21], rep(mean(df$a1, na.rm = TRUE), 2L))

  # Imputing from all features lets block B change block A's imputed values.
  g_all = impute_knn_graph(k = 5)
  df_all = df[-c(20, 21), ]
  d2_all = d2[-c(20, 21), ]
  expect_false(isTRUE(all.equal(impute(df_all, g_all)$a2, impute(d2_all, g_all)$a2)))

  # The preprocessing graph passes the option through.
  pre = mbspls_preproc_graph(blocks = blocks, impute_within_blocks = TRUE)
  expect_true("impute_num_knn_block_A" %in% pre$ids())
  expect_true("impute_num_knn" %in% mbspls_preproc_graph(blocks = blocks)$ids())
})

test_that("Gower k-NN learners can fall back to the training marginal", {
  set.seed(9)
  n = 30
  df = data.frame(x1 = rnorm(n), x2 = rnorm(n), y = rnorm(n))
  task = mlr3::as_task_regr(df, target = "y")
  new = data.frame(x1 = c(NA, 0.1), x2 = c(NA, -0.2), y = c(0, 0))
  new_task = mlr3::as_task_regr(new, target = "y")

  lr = mlr3::lrn("regr.knngower", k = 3L)
  lr$train(task)
  expect_error(lr$predict(new_task), "no eligible neighbors")
  lr$param_set$values$no_neighbors = "marginal"
  pr = lr$predict(new_task)
  expect_equal(pr$response[1L], mean(df$y))
  lr$param_set$values$no_neighbors = "error"
  expect_equal(pr$response[2L], lr$predict(mlr3::as_task_regr(new[2L, ], target = "y"))$response)

  dfc = data.frame(x1 = rnorm(n), x2 = rnorm(n), y = factor(rep(c("a", "b", "b"), 10L)))
  lc = mlr3::lrn("classif.knngower", k = 3L, predict_type = "prob", no_neighbors = "marginal")
  lc$train(mlr3::as_task_classif(dfc, target = "y"))
  pc = lc$predict(mlr3::as_task_classif(data.frame(x1 = NA_real_, x2 = NA_real_, y = factor("a", levels = c("a", "b"))), target = "y"))
  expect_equal(as.character(pc$response), "b")
  expect_equal(unname(pc$prob[1L, ]), c(1, 2) / 3)
})

test_that("as_mbspls_participation reproduces the summary of mbspls_participation", {
  testthat::skip_if_not("regr.featureless" %in% mlr3::mlr_learners$keys())
  toy = participation_toy()
  gl = participation_learner(toy$blocks)
  res = mbspls_participation(toy$task, gl, site = "site", n_perm = 99L, seed = 5L)
  rebuilt = as_mbspls_participation(as.data.frame(res$centres), alpha = res$alpha, n_perm = res$n_perm)
  expect_equal(rebuilt$centres, res$centres)
  expect_equal(rebuilt$blocks, res$blocks)
  expect_equal(rebuilt$partial_conjunction, res$partial_conjunction)
  expect_s3_class(ggplot2::autoplot(res), "ggplot")
  expect_s3_class(ggplot2::autoplot(res, type = "heldout"), "ggplot")
})

test_that("as_mbspls_participation fills missing cells and checks its input", {
  centres = data.frame(
    site = c("a", "a", "b", "c", "c"), block = c("X", "Y", "X", "X", "Y"),
    n = c(50, 50, 40, 30, 30), r = c(0.4, 0.1, 0.35, 0.3, 0),
    p = c(0.001, 0.2, 0.004, 0.01, 0.6)
  )
  res = as_mbspls_participation(centres)
  expect_s3_class(res, "mbspls_participation")
  expect_identical(nrow(res$centres), 6L)
  filled = res$centres[site == "b" & block == "Y"]
  expect_identical(filled$p, 1)
  expect_identical(filled$n, 40L)
  expect_true(all(is.na(res$centres$train_r)))
  expect_identical(res$blocks$S, c(3L, 3L))
  # X: 3 * 0.001, 2 * 0.004 and 0.01 are all below 0.05 / 2.
  expect_identical(res$blocks$r_hat, c(3L, 0L))
  expect_identical(res$blocks$r_hat_simultaneous, c(3L, 0L))
  # X: all three p-values are below 0.05 / 3; 0.01 is above 0.05 / 6.
  expect_identical(res$blocks$r_count, c(3L, 0L))
  expect_identical(res$blocks$r_count_simultaneous, c(2L, 0L))
  expect_output(print(res), "3 centres\\.")

  # One external validation cohort (S = 1).
  one = as_mbspls_participation(centres[centres$site == "a", ])
  expect_identical(one$blocks$r_hat, c(1L, 0L))
  expect_identical(one$blocks$r_count, c(1L, 0L))
  expect_identical(one$blocks$r_count_simultaneous, c(1L, 0L))

  expect_error(as_mbspls_participation(centres[, names(centres) != "p"]), "'p'")
  expect_error(as_mbspls_participation(rbind(centres, centres[1L, ])), "more than one row")
  bad = centres
  bad$p[1L] = 2
  expect_error(as_mbspls_participation(bad), "p-values")
})

test_that("autoplot.mbspls_participation draws centres, bounds and the held-out drop", {
  centres = data.frame(
    site = rep(c("a", "b", "c", "d"), each = 2), block = rep(c("X_block", "Y"), 4),
    n = rep(c(120, 80, 60, 30), each = 2),
    r = c(0.42, 0.05, 0.38, -0.02, 0.29, 0.11, 0.35, 0.08),
    train_r = c(0.45, 0.30, 0.44, 0.31, 0.47, 0.28, 0.46, 0.30),
    p = c(0.0001, 0.31, 0.0013, 0.55, 0.012, 0.20, 0.03, 0.34)
  )
  res = as_mbspls_participation(centres)
  p = ggplot2::autoplot(res)
  expect_s3_class(p, "ggplot")
  pts = ggplot2::layer_data(p, 2L)
  expect_identical(nrow(pts), 8L)
  # The block with the larger held-out correlations is drawn on top.
  expect_true(all(as.numeric(pts$y[pts$x > 0.25]) == 2))
  expect_true(all(as.numeric(pts$y[pts$x < 0.25]) == 1))
  # The bound column. X_block has p = 0.0001, 0.0013, 0.012, 0.03 at S = 4
  # centres, K = 2 blocks. Count: 2 below 0.05 / 8, 3 below 0.05 / 4.
  # Partial conjunction: 3 (0.024 <= 0.05 / 2 < 0.03) and 4. Y: 0.
  lab_of = function(pl) sort(ggplot2::layer_data(pl, 3L)$label)
  expect_identical(lab_of(p), c("0", "2"))
  expect_identical(lab_of(ggplot2::autoplot(res, bound = "unadjusted")), c("0", "3"))
  expect_identical(lab_of(ggplot2::autoplot(res, rule = "partial_conjunction")), c("0", "3"))
  expect_identical(lab_of(ggplot2::autoplot(res, bound = "unadjusted", rule = "partial_conjunction")), c("0", "4"))
  expect_error(ggplot2::autoplot(res, rule = "fdr"))
  # Objects without the count columns are summarised again from their centres.
  old = res
  old$blocks = old$blocks[, !c("r_count", "r_count_simultaneous")]
  expect_identical(lab_of(ggplot2::autoplot(old)), c("0", "2"))
  expect_identical(lab_of(ggplot2::autoplot(old, rule = "partial_conjunction")), c("0", "3"))
  expect_length(ggplot2::autoplot(res, bound = "none")$layers, 2L)
  p2 = ggplot2::autoplot(res, compare = list(Other = res), label = "Main")
  expect_length(p2$layers, 5L)
  expect_identical(ggplot2::layer_data(p2, 4L)$label, c("Main", "Other"))

  h = ggplot2::autoplot(res, type = "heldout", labels = c(Y = "Second block"))
  hd = ggplot2::layer_data(h, 2L)
  expect_identical(nrow(hd), 4L)
  w = c(120, 80, 60, 30)
  expect_true(any(abs(hd$x - sum(w * c(0.42, 0.38, 0.29, 0.35)) / sum(w)) < 1e-12))
  expect_true(any(abs(hd$x - mean(c(0.45, 0.44, 0.47, 0.46))) < 1e-12))
  expect_true(any(abs(hd$x - mean(c(0.30, 0.31, 0.28, 0.30))) < 1e-12))

  # The unit word appears in every label.
  us = ggplot2::autoplot(res, unit = c("site", "sites"))
  expect_identical(
    ggplot2::get_labs(us)$x,
    "Correlation with the other blocks\u2019 composite at the held-out site"
  )
  expect_identical(ggplot2::get_labs(us)$size, "Site n")
  expect_identical(levels(us$data$sig), c("p \u2264 0.05 at that site", "p > 0.05"))
  uh = ggplot2::autoplot(res, type = "heldout", unit = c("site", "sites"))
  expect_identical(levels(uh$layers[[2L]]$data$kind), c("Training sites", "Held-out site"))
  expect_error(ggplot2::autoplot(res, unit = "site"))

  # A given order (top to bottom) overrides the held-out ordering.
  po = ggplot2::autoplot(res, order = c("Y", "X_block"))
  expect_identical(levels(po$data$label), c("X block", "Y"))
  expect_true(all(as.numeric(ggplot2::layer_data(po, 2L)$y[ggplot2::layer_data(po, 2L)$x > 0.25]) == 1))
  expect_error(ggplot2::autoplot(res, order = "Y"), "every block")

  # Common limits widen the axis; the bound column moves to the right of them.
  pl = ggplot2::autoplot(res, limits = c(-0.5, 0.9))
  expect_true(min(ggplot2::layer_data(pl, 3L)$x) > 0.9)
  expect_error(ggplot2::autoplot(res, limits = c(1, 0)))

  # Two analyses in one held-out plot: one colour each, offset within the rows.
  res2 = as_mbspls_participation(transform(centres, r = r + 0.05, train_r = train_r + 0.02))
  hc = ggplot2::autoplot(res, type = "heldout", compare = list(Second = res2), label = "First",
    colours = c("#2a78d6", "#eb6834"))
  hl = ggplot2::layer_data(hc, 2L)
  expect_identical(nrow(hl), 8L)
  expect_setequal(unique(hl$colour), c("#2a78d6", "#eb6834"))
  expect_equal(sort(unique(hl$y)), c(1 - 0.18, 1 + 0.18, 2 - 0.18, 2 + 0.18))
  expect_identical(levels(hc$layers[[2L]]$data$analysis), c("First", "Second"))
  expect_error(ggplot2::autoplot(res, type = "heldout", compare = list(B = res2, C = res2)), "one colour per analysis")

  no_train = as_mbspls_participation(centres[, names(centres) != "train_r"])
  expect_error(ggplot2::autoplot(no_train, type = "heldout"), "train_r")
  expect_error(ggplot2::autoplot(res, compare = list(res)), "named list")
  expect_error(ggplot2::autoplot(res, labels = "a"), "named by block")
})
