pair_toy = function(n = 80L, seed = 1L, signal = 1) {
  set.seed(seed)
  f = rnorm(n)
  cbind(A = signal * f + rnorm(n), B = signal * f + rnorm(n), C = rnorm(n))
}

test_that("mbspls_pair_heldout_tests computes the stated statistic", {
  tr = pair_toy(seed = 1L)
  te = pair_toy(seed = 2L)
  res = mbspls_pair_heldout_tests(te, tr, n_perm = 199L)
  expect_identical(names(res), c("a", "b", "n", "train_r", "r", "p"))
  expect_identical(res$a, c("A", "A", "B"))
  expect_identical(res$b, c("B", "C", "C"))
  expect_true(all(res$n == 80L))
  expect_true(all(res$p >= 1 / 200 & res$p <= 1))
  # By hand for A-C: held-out correlation signed by the training correlation.
  s = sign(cor(tr[, "A"], tr[, "C"]))
  expect_equal(res[a == "A" & b == "C", r], s * cor(te[, "A"], te[, "C"]))
  expect_equal(res[a == "A" & b == "C", train_r], cor(tr[, "A"], tr[, "C"]))
  # The associated pair is detected, the noise pairs are not.
  expect_lte(res[a == "A" & b == "B", p], 0.01)
  expect_gt(min(res[b == "C", p]), 0.01)

  # Strata: the function runs and returns valid p-values.
  st = mbspls_pair_heldout_tests(te, tr, n_perm = 99L, strata = rep(1:2, 40))
  expect_true(all(st$p > 0 & st$p <= 1))

  # A constant held-out score gives r = NA and p = 1 for its pairs.
  te2 = te
  te2[, "C"] = 1
  cst = mbspls_pair_heldout_tests(te2, tr, n_perm = 99L)
  expect_true(all(is.na(cst[b == "C", r])))
  expect_identical(cst[b == "C", p], c(1, 1))

  bad = te
  colnames(bad) = c("A", "C", "B")
  expect_error(mbspls_pair_heldout_tests(bad, tr, n_perm = 9L), "same unique column names")
  expect_error(mbspls_pair_heldout_tests(unname(te), unname(tr), n_perm = 9L), "same unique column names")
  expect_error(mbspls_pair_heldout_tests(te, tr[, 1:2], n_perm = 9L))
  expect_error(mbspls_pair_heldout_tests(te, tr, n_perm = 9L, strata = 1:3), "strata")
})

test_that("mbspls_pair_heldout_tests is valid under the null", {
  skip_on_cran()
  p = vapply(seq_len(200), function(i) {
    tr = pair_toy(n = 60L, seed = 1000L + i)
    te = pair_toy(n = 40L, seed = 5000L + i)
    mbspls_pair_heldout_tests(te, tr, n_perm = 199L)[a == "A" & b == "C", p]
  }, numeric(1L))
  expect_gt(mean(p), 0.4)
  expect_lt(mean(p), 0.6)
  expect_lt(mean(p <= 0.05), 0.1)
})

test_that("as_mbspls_pairs decides simultaneously over the pairs", {
  pairs = data.frame(a = c("A", "A", "B"), b = c("B", "C", "C"), r = c(0.45, 0.02, 0.30),
    p = c(0.0002, 0.44, 0.02), n = 60, train_r = c(0.5, -0.1, 0.25))
  x = as_mbspls_pairs(pairs, modality = c(A = "MRI", B = "MRI", C = "Genotype"))
  expect_s3_class(x, "mbspls_pairs")
  expect_equal(x$level, 0.05 / 3)
  s = x$summary
  expect_identical(s$significant, c(TRUE, FALSE, FALSE)) # 0.02 > 0.05 / 3
  expect_identical(s$same_modality, c(TRUE, FALSE, FALSE))
  expect_equal(s$p, pairs$p)
  expect_equal(s[1, lo], tanh(atanh(0.45) - qnorm(0.975) / sqrt(57)))
  expect_equal(s[1, hi], tanh(atanh(0.45) + qnorm(0.975) / sqrt(57)))
  expect_output(print(x), "3 pairs of 3 blocks, 1 site")

  # Several sites: the fixed-threshold count over sites; a missing row counts as p = 1.
  multi = data.frame(site = rep(c("s1", "s2", "s3"), each = 2), a = "A", b = rep(c("B", "C"), 3),
    r = c(0.5, 0.1, 0.4, 0.0, 0.45, 0.05), p = c(1e-4, 0.3, 1e-3, 0.6, 0.01, 0.4), n = rep(c(50, 40, 30), each = 2))
  multi = multi[-6, ]
  y = as_mbspls_pairs(multi, min_sites = 2L)
  expect_equal(y$level, 0.05 / (2 * 3))
  expect_identical(y$S, 3L)
  ab = y$summary[b == "B"]
  expect_identical(ab$r_count, 2L) # 1e-4 and 1e-3 <= 0.0083; 0.01 is not
  expect_true(ab$significant)
  expect_equal(ab$heldout_r, sum(c(50, 40, 30) * c(0.5, 0.4, 0.45)) / 120)
  expect_true(is.na(ab$p))
  expect_identical(y$pairs[site == "s3" & b == "C", p], 1)
  expect_false(y$summary[b == "C", significant])
  expect_equal(y$summary[b == "C", heldout_r], sum(c(50, 40) * c(0.1, 0)) / 90)

  expect_error(as_mbspls_pairs(pairs[, c("a", "b", "r")]), "'p'")
  expect_error(as_mbspls_pairs(rbind(pairs, pairs)), "once per site")
  expect_error(as_mbspls_pairs(pairs, modality = c(A = "MRI")), "every block")
})

test_that("autoplot.mbspls_pairs draws a network and a dumbbell", {
  pairs = data.frame(a = c("A", "A", "A", "B", "B", "C"), b = c("B", "C", "D", "C", "D", "D"),
    r = c(0.4, 0.05, 0.6, 0.0, 0.1, 0.7), p = c(0.001, 0.4, 1e-4, 0.5, 0.3, 1e-4),
    n = 50, train_r = c(0.5, 0.1, 0.55, 0.05, 0.2, 0.6))
  x = as_mbspls_pairs(pairs, modality = c(A = "MRI", B = "MRI", C = "Genotype", D = "Clinical"))
  p1 = ggplot2::autoplot(x)
  expect_s3_class(p1, "ggplot")
  edges = ggplot2::layer_data(p1, 1L)
  expect_identical(nrow(edges), 6L)
  expect_identical(nrow(ggplot2::layer_data(p1, 3L)), 4L)
  p2 = ggplot2::autoplot(x, type = "pairs", title = "B")
  expect_s3_class(p2, "ggplot")
  expect_true(all(c("Same modality", "MRI with Genotype", "MRI with Clinical", "Genotype with Clinical") %in% levels(p2$data$group)))
  # Custom layout, labels and groups; no modality.
  lay = data.frame(block = c("A", "B", "C", "D"), x = c(0, 0, 2, 2), y = c(1, 0, 1, 0))
  p3 = ggplot2::autoplot(as_mbspls_pairs(pairs), layout = lay, labels = c(C = "Block\nC"))
  pts = ggplot2::layer_data(p3, 3L)
  expect_equal(sort(pts$x), c(0, 0, 2, 2))
  gr = data.frame(a = pairs$a, b = pairs$b, group = rep(c("first", "second"), each = 3))
  p4 = ggplot2::autoplot(x, type = "pairs", groups = gr, labels = c(C = "Block\nC"))
  expect_identical(levels(p4$data$group), c("first", "second"))
  expect_true(any(grepl("Block C", levels(p4$data$pair))))
  expect_error(ggplot2::autoplot(x, layout = lay[1:2, ]), "layout")
  expect_error(ggplot2::autoplot(x, legend_labels = "only one"), "two entries")
})

test_that("edge labels that would collide are moved apart along their edges", {
  # Two crossing edges share their midpoint (1, 0.5).
  e = data.table::data.table(x = c(0, 0), y = c(1, 0), x_b = c(2, 2), y_b = c(0, 1))
  out = mlr3mbspls:::.mbspls_pairs_edge_labels(data.table::copy(e), thr = 0.3)
  d = sqrt(diff(out$mx)^2 + diff(out$my)^2)
  expect_gt(d, 0.3)
  # Each label stays on its own edge.
  on_edge = function(i) {
    t = (out$mx[i] - e$x[i]) / (e$x_b[i] - e$x[i])
    abs(e$y[i] + t * (e$y_b[i] - e$y[i]) - out$my[i]) < 1e-12
  }
  expect_true(on_edge(1L) && on_edge(2L))
  far = mlr3mbspls:::.mbspls_pairs_edge_labels(data.table::data.table(x = c(0, 0), y = c(0, 2), x_b = c(1, 1), y_b = c(0, 2)), thr = 0.3)
  expect_equal(far$mx, c(0.5, 0.5))
})
