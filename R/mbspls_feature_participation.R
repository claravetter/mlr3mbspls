#' Held-out feature participation at one centre
#'
#' @description
#' Feature-level counterpart of the block test in [mbspls_participation()].
#' For every feature of every block it tests whether the feature's held-out
#' values are associated with the composite of the other blocks' held-out
#' scores, with weights and preprocessing learned at other centres.
#'
#' @details
#' For block `k`, the composite is the sum of the other blocks' held-out
#' scores, standardised within the centre and signed by the training score
#' correlations, as in [mbspls_participation()]. The statistic of feature `f`
#' of block `k` is `s_f * cor(x_f, z_k)`, where `x_f` are the feature's
#' held-out values and `s_f` is the sign of the feature's correlation with the
#' training composite. Its one-sided p-value compares the statistic with
#' `n_perm` permutations of the composite within the centre (within `strata`,
#' if given). This is equivalent to permuting the feature and lets all features
#' of a block share one set of permutations.
#'
#' The p-value is exact, conditionally on the training data, when the rows of
#' the centre are exchangeable and the feature is independent of the other
#' blocks' scores there. As for the block test, the held-out values of a
#' feature must depend only on its own block's data (impute within blocks).
#'
#' In MB-sPLS the weight of a feature can be non-zero only if its covariance
#' with the composite of the other blocks is non-zero, because the sparse
#' update soft-thresholds these covariances. A replicated held-out association
#' is therefore the counterpart of a non-zero weight that does not depend on
#' the sparsity budget.
#'
#' @param X_test,X_train (named `list` of `matrix`) Preprocessed held-out and
#'   training data, one matrix per block, with the same block names in the same
#'   order and, within every block, the same columns in the same order. When
#'   both matrices of a block have column names, they must be identical, so
#'   that permuted or renamed columns cannot receive the orientation of another
#'   feature.
#' @param T_test,T_train (`matrix`) Held-out and training block scores of one
#'   component, one column per block, named as `X_test`.
#' @param n_perm (`integer(1)`) Permutations. Default `9999`.
#' @param strata (`NULL` | `vector`) Optional strata of the held-out rows;
#'   permutations then stay within strata.
#' @param modality (`NULL` | named `character()`) Optional modality of every
#'   block. The composite of block `k` then uses only the blocks of other
#'   modalities (cross-modality participation).
#'
#' @return A `data.table` with one row per feature: `block`, `feature`, the
#'   oriented held-out correlation `r`, its orientation `sign`, the signed
#'   training correlation `train_r` and the one-sided p-value `p`. Constant
#'   features get `r = NA` and `p = 1`.
#'
#' @seealso [as_mbspls_feature_participation()], [mbspls_participation()]
#' @examples
#' set.seed(1)
#' n = 60
#' f = rnorm(n)
#' X = list(A = cbind(a1 = f + rnorm(n), a2 = rnorm(n)), B = cbind(b1 = f + rnorm(n)))
#' W = list(A = c(1, 0), B = 1)
#' T_mat = sapply(names(X), function(b) drop(X[[b]] %*% W[[b]]))
#' mbspls_feature_heldout_tests(X, X, T_mat, T_mat, n_perm = 199L)
#' @export
mbspls_feature_heldout_tests = function(
  X_test,
  X_train,
  T_test,
  T_train,
  n_perm = 9999L,
  strata = NULL,
  modality = NULL
) {
  checkmate::assert_list(X_test, types = "matrix", min.len = 2L, names = "unique")
  checkmate::assert_list(X_train, types = "matrix", min.len = 2L, names = "unique")
  if (!identical(names(X_test), names(X_train))) {
    stop("`X_test` and `X_train` need the same blocks in the same order.", call. = FALSE)
  }
  checkmate::assert_matrix(T_test, mode = "numeric", min.cols = 2L)
  checkmate::assert_matrix(T_train, mode = "numeric", ncols = ncol(T_test))
  blocks = colnames(T_test)
  if (is.null(blocks) || !identical(blocks, colnames(T_train)) ||
    !all(blocks %in% names(X_test)) || !all(blocks %in% names(X_train))) {
    stop("`T_test` and `T_train` need the same named columns, one per block of `X_test` and `X_train`.", call. = FALSE)
  }
  n = nrow(T_test)
  for (b in blocks) {
    if (nrow(X_test[[b]]) != n || nrow(X_train[[b]]) != nrow(T_train) || ncol(X_test[[b]]) != ncol(X_train[[b]])) {
      stop(sprintf("The data of block '%s' do not match the scores.", b), call. = FALSE)
    }
    cn_test = colnames(X_test[[b]])
    cn_train = colnames(X_train[[b]])
    if (!is.null(cn_test) && !is.null(cn_train) && !identical(cn_test, cn_train)) {
      stop(sprintf(
        "The columns of block '%s' differ between `X_test` and `X_train` (names or order).", b
      ), call. = FALSE)
    }
  }
  n_perm = checkmate::assert_int(n_perm, lower = 1L, coerce = TRUE)
  if (!is.null(strata) && length(strata) != n) {
    stop("`strata` needs one value per held-out row.", call. = FALSE)
  }
  if (!is.null(modality)) {
    checkmate::assert_character(modality, any.missing = FALSE, names = "unique")
    if (!all(blocks %in% names(modality))) {
      stop("`modality` must name every block.", call. = FALSE)
    }
  }

  K = length(blocks)
  sg = suppressWarnings(sign(stats::cor(T_train)))
  sg[!is.finite(sg)] = 0
  z = lapply(seq_len(K), function(k) .mbspls_participation_std(T_test[, k]))
  z_train = scale(T_train)
  z_train[!is.finite(z_train)] = 0
  perms = if (n >= 3L) .mbspls_permutations(n, n_perm, strata)

  out = lapply(seq_len(K), function(k) {
    b = blocks[k]
    others = setdiff(seq_len(K), k)
    if (!is.null(modality)) others = others[modality[blocks[others]] != modality[[b]]]
    xte = X_test[[b]]
    xtr = X_train[[b]]
    feats = colnames(xte)
    if (is.null(feats)) feats = paste0(b, "_", seq_len(ncol(xte)))
    res = data.table::data.table(block = b, feature = feats, r = NA_real_, sign = 1, train_r = NA_real_, p = 1)
    if (!length(others)) {
      return(res)
    }
    comp_train = drop(z_train[, others, drop = FALSE] %*% sg[k, others])
    tr_r = suppressWarnings(as.numeric(stats::cor(xtr, comp_train)))
    sgn = sign(tr_r)
    sgn[!is.finite(sgn) | sgn == 0] = 1
    data.table::set(res, j = "sign", value = sgn)
    data.table::set(res, j = "train_r", value = tr_r)
    if (n < 3L) {
      return(res)
    }
    zo = vapply(others, function(j) if (is.null(z[[j]])) numeric(n) else z[[j]], numeric(n))
    comp = .mbspls_participation_std(drop(matrix(zo, n) %*% sg[k, others]))
    if (is.null(comp)) {
      return(res)
    }
    sds = apply(xte, 2L, stats::sd)
    ok = is.finite(sds) & sds > 1e-12
    if (!any(ok)) {
      return(res)
    }
    xs = scale(xte[, ok, drop = FALSE])
    obs = drop(crossprod(xs, comp)) / (n - 1) * sgn[ok]
    null = crossprod(xs, matrix(comp[perms], n, n_perm)) / (n - 1) * sgn[ok]
    pv = (1 + rowSums(null >= obs - 1e-10)) / (n_perm + 1)
    data.table::set(res, i = which(ok), j = "r", value = obs)
    data.table::set(res, i = which(ok), j = "p", value = pv)
    res
  })
  data.table::rbindlist(out)
}

# Permutation matrix (n x n_perm), within strata if given.
.mbspls_permutations = function(n, n_perm, strata = NULL) {
  if (is.null(strata)) {
    return(vapply(seq_len(n_perm), function(i) sample.int(n), integer(n)))
  }
  idx = split(seq_len(n), strata)
  vapply(seq_len(n_perm), function(i) {
    o = seq_len(n)
    for (g in idx) o[g] = g[sample.int(length(g))]
    o
  }, integer(n))
}


#' Feature participation across centres
#'
#' @description
#' Combines centre-level feature tests, such as those of
#' [mbspls_feature_heldout_tests()] at every held-out centre, into a
#' replication count per feature and a false discovery rate decision across
#' features: which features of a signature participate at `r` or more centres.
#'
#' @details
#' For every feature, the centre-level p-values give Bonferroni partial
#' conjunction p-values (see [mbspls_partial_conjunction()]), valid under
#' arbitrary dependence between centres. `r_hat` is the lower `1 - alpha`
#' bound for the number of centres at which the feature participates. The
#' partial conjunction p-values for "at least `r` centres" are then adjusted
#' across all features with `fdr` (Benjamini-Hochberg by default, or
#' Benjamini-Yekutieli for arbitrary dependence between features), and
#' `replicated` marks the features with an adjusted value of at most `alpha`.
#'
#' @param centres (`data.frame`) One row per centre and feature, with columns
#'   `site`, `block`, `feature`, `n`, the oriented held-out correlation `r`
#'   and the p-value `p`, and optionally `train_r`. A feature without a row at
#'   some centre is recorded there with `p = 1`.
#' @param reference (`NULL` | `data.frame`) Optional reference signature, for
#'   example the model fitted to all centres, with columns `block`, `feature`
#'   and `weight`, and optionally `direction` (the sign of the feature in the
#'   signature, used to sign the held-out correlations in plots; by default the
#'   sign of a non-zero weight).
#' @param r (`NULL` | `integer(1)`) Number of centres for the replication
#'   decision. Default `NULL`: a majority, `ceiling(S / 2)`.
#' @param alpha (`numeric(1)`) Level of the bounds and of the false discovery
#'   rate. Default `0.05`.
#' @param fdr (`character(1)`) `"BH"` (default), `"BY"` or `"none"`.
#'
#' @return An object of class `mbspls_feature_participation`, a list with
#'   `features` (one row per feature: `block`, `feature`, `S`, `heldout_r`
#'   weighted by centre size, mean `train_r`, `n_significant`, `r_hat`,
#'   `p_r`, `q_r`, `replicated`, and `weight` and `direction` from
#'   `reference`), `centres`, and the settings `r`, `alpha`, `fdr` and `S`.
#'
#' @seealso [mbspls_feature_heldout_tests()],
#'   [autoplot.mbspls_feature_participation()]
#' @examples
#' centres = data.frame(
#'   site = rep(c("a", "b", "c"), each = 3), block = rep(c("X", "X", "Y"), 3),
#'   feature = rep(c("x1", "x2", "y1"), 3), n = rep(c(50, 40, 30), each = 3),
#'   r = c(0.40, 0.02, 0.35, 0.38, -0.05, 0.30, 0.33, 0.04, 0.28),
#'   p = c(0.0001, 0.45, 0.002, 0.0004, 0.62, 0.03, 0.04, 0.4, 0.06)
#' )
#' as_mbspls_feature_participation(centres)
#' @export
as_mbspls_feature_participation = function(
  centres,
  reference = NULL,
  r = NULL,
  alpha = 0.05,
  fdr = c("BH", "BY", "none")
) {
  checkmate::assert_data_frame(centres, min.rows = 1L)
  missing_cols = setdiff(c("site", "block", "feature", "n", "r", "p"), names(centres))
  if (length(missing_cols)) {
    stop(sprintf("`centres` lacks the column(s) %s.", paste0("'", missing_cols, "'", collapse = ", ")), call. = FALSE)
  }
  checkmate::assert_number(alpha, lower = 0, upper = 1)
  fdr = match.arg(fdr)
  d = data.table::data.table(
    site = as.character(centres$site), block = as.character(centres$block),
    feature = as.character(centres$feature), n = as.integer(centres$n),
    r = as.numeric(centres$r), p = as.numeric(centres$p),
    train_r = if ("train_r" %in% names(centres)) as.numeric(centres$train_r) else NA_real_
  )
  if (anyDuplicated(d, by = c("site", "block", "feature"))) {
    stop("`centres` has more than one row for some centre and feature.", call. = FALSE)
  }
  if (any(d$p < 0 | d$p > 1, na.rm = TRUE)) {
    stop("`p` must contain p-values in [0, 1].", call. = FALSE)
  }
  sites = unique(d$site)
  S = length(sites)
  r_target = if (is.null(r)) ceiling(S / 2) else r
  r_target = checkmate::assert_int(r_target, lower = 1L, upper = S, coerce = TRUE)
  keys = unique(d[, c("block", "feature")])
  grid = data.table::CJ(site = sites, key = seq_len(nrow(keys)), sorted = FALSE)
  grid[, `:=`(block = keys$block[key], feature = keys$feature[key])][, key := NULL]
  site_n = d[!is.na(n), list(n_site = max(n)), by = "site"]
  d = d[grid, on = c("site", "block", "feature")]
  d[is.na(n), n := site_n$n_site[match(site, site_n$site)]]
  d[is.na(p), p := 1]

  feats = d[, {
    pc = cummax(pmin(1, (S - seq_len(S) + 1) * sort(p)))
    ok = !is.na(r)
    list(
      S = .N,
      heldout_r = if (any(ok)) sum((n * r)[ok]) / sum(n[ok]) else NA_real_,
      train_r = if (all(is.na(train_r))) NA_real_ else mean(train_r, na.rm = TRUE),
      n_significant = sum(p <= alpha),
      r_hat = sum(pc <= alpha),
      p_r = pc[r_target]
    )
  }, by = c("block", "feature")]
  feats[, q_r := if (fdr == "none") p_r else stats::p.adjust(p_r, method = fdr)]
  feats[, replicated := q_r <= alpha]
  feats[, `:=`(weight = NA_real_, direction = NA_real_)]
  if (!is.null(reference)) {
    checkmate::assert_data_frame(reference)
    if (!all(c("block", "feature", "weight") %in% names(reference))) {
      stop("`reference` needs the columns 'block', 'feature' and 'weight'.", call. = FALSE)
    }
    ref = data.table::as.data.table(reference)
    m = match(paste(feats$block, feats$feature, sep = "\r"), paste(ref$block, ref$feature, sep = "\r"))
    feats[, weight := as.numeric(ref$weight[m])]
    dir = if ("direction" %in% names(ref)) as.numeric(ref$direction[m]) else rep(NA_real_, nrow(feats))
    from_weight = is.na(dir) & !is.na(feats$weight) & feats$weight != 0
    dir[from_weight] = sign(feats$weight[from_weight])
    feats[, direction := dir]
  }
  structure(list(
    features = feats[], centres = d[], r = r_target, alpha = alpha, fdr = fdr, S = S
  ), class = "mbspls_feature_participation")
}


#' @export
print.mbspls_feature_participation = function(x, ...) {
  f = x$features
  cat(sprintf(
    "Held-out feature participation: %d features in %d blocks, %d centres.\n",
    nrow(f), length(unique(f$block)), x$S
  ))
  cat(sprintf(
    "Replicated at %d or more centres (%s, level %g): %d features.\n",
    x$r, if (x$fdr == "none") "unadjusted" else paste(x$fdr, "FDR"), x$alpha, sum(f$replicated)
  ))
  print(f[, list(features = .N, replicated = sum(replicated)), by = "block"])
  invisible(x)
}


#' Plot a feature-level signature
#'
#' @description
#' Draws the result of [as_mbspls_feature_participation()]: one row per
#' feature, grouped by block, with its held-out correlation with the composite
#' of the other blocks (signed by the reference direction, if known). Filled
#' points mark features that participate at `r` or more centres after the
#' false discovery rate adjustment. If the object carries reference weights,
#' circles mark features with a non-zero weight in the reference signature and
#' triangles features with a zero weight. With `compare`, further analyses are
#' drawn in further columns with the same rows.
#'
#' @param object (`mbspls_feature_participation`) The result to plot.
#' @param compare (`NULL` | named `list()`) Further
#'   `mbspls_feature_participation` objects, named by their labels.
#' @param label (`NULL` | `character(1)`) Column label of `object`. Default
#'   `"Main"` when `compare` is given.
#' @param blocks (`NULL` | `character()`) Blocks to show, in this order.
#'   Default: all blocks of `object`.
#' @param top (`NULL` | `integer(1)`) Show at most this many features per
#'   block, those with the largest absolute held-out correlation in any
#'   analysis. Default `NULL` shows all.
#' @param labels (`NULL` | named `character()`) Display names of the blocks.
#' @param feature_labels (`NULL` | named `character()`) Display names of the
#'   features, named by feature.
#' @param colours (`character(2)`) Colours of replicated and of the other
#'   features.
#' @param unit (`character(2)`) Singular and plural of the held-out unit as
#'   used in the labels. Default `c("centre", "centres")`.
#' @param base_size (`numeric(1)`) Base font size in points. Default `11`.
#' @param ... Not used.
#'
#' @return A [ggplot2::ggplot()] object.
#' @seealso [as_mbspls_feature_participation()]
#' @examples
#' centres = data.frame(
#'   site = rep(c("a", "b", "c"), each = 3), block = rep(c("X", "X", "Y"), 3),
#'   feature = rep(c("x1", "x2", "y1"), 3), n = rep(c(50, 40, 30), each = 3),
#'   r = c(0.40, 0.02, 0.35, 0.38, -0.05, 0.30, 0.33, 0.04, 0.28),
#'   p = c(0.0001, 0.45, 0.002, 0.0004, 0.62, 0.03, 0.04, 0.4, 0.06)
#' )
#' ref = data.frame(block = c("X", "X", "Y"), feature = c("x1", "x2", "y1"), weight = c(0.9, 0, 1))
#' ggplot2::autoplot(as_mbspls_feature_participation(centres, reference = ref, r = 2))
#' @export
autoplot.mbspls_feature_participation = function(
  object,
  compare = NULL,
  label = NULL,
  blocks = NULL,
  top = NULL,
  labels = NULL,
  feature_labels = NULL,
  colours = c("#2a78d6", "#898781"),
  unit = c("centre", "centres"),
  base_size = 11,
  ...
) {
  checkmate::assert_class(object, "mbspls_feature_participation")
  checkmate::assert_character(unit, len = 2L, any.missing = FALSE, min.chars = 1L)
  unit_word = if (object$r == 1L) unit[1L] else unit[2L]
  if (!is.null(compare)) {
    checkmate::assert_list(compare, min.len = 1L)
    if (is.null(names(compare)) || any(!nzchar(names(compare))) ||
      !all(vapply(compare, inherits, logical(1L), what = "mbspls_feature_participation"))) {
      stop("`compare` must be a named list of `mbspls_feature_participation` objects.", call. = FALSE)
    }
  }
  checkmate::assert_string(label, null.ok = TRUE)
  checkmate::assert_character(blocks, any.missing = FALSE, null.ok = TRUE)
  checkmate::assert_int(top, lower = 1L, null.ok = TRUE)
  checkmate::assert_character(colours, len = 2L, any.missing = FALSE)
  checkmate::assert_number(base_size, lower = 1)
  objs = c(list(object), compare)
  names(objs) = c(if (is.null(label)) "Main" else label, names(compare))
  inks = .mbspls_participation_inks()

  d = data.table::rbindlist(lapply(names(objs), function(a) {
    f = data.table::copy(objs[[a]]$features)
    f[, analysis := a]
    f
  }), use.names = TRUE)
  if (is.null(blocks)) blocks = unique(object$features$block)
  d = d[block %in% blocks]
  d[, x := ifelse(is.na(direction), heldout_r, direction * heldout_r)]
  # Rows: features of the main analysis, ordered within blocks by its signed held-out correlation.
  ord = d[analysis == names(objs)[1L]][order(match(block, blocks), -x)]
  if (!is.null(top)) {
    strength = d[, list(s = max(abs(heldout_r), na.rm = TRUE)), by = c("block", "feature")]
    keep = strength[order(-s), utils::head(.SD, top), by = "block"]
    ord = ord[paste(block, feature) %in% paste(keep$block, keep$feature)]
  }
  d = d[paste(block, feature) %in% paste(ord$block, ord$feature)]
  flab = .nice_label(ord$feature)
  if (!is.null(feature_labels)) {
    hit = ord$feature %in% names(feature_labels)
    flab[hit] = feature_labels[ord$feature[hit]]
  }
  d[, row := factor(paste(block, feature, sep = "\r"), levels = rev(paste(ord$block, ord$feature, sep = "\r")))]
  blab = .nice_label(blocks)
  if (!is.null(labels)) {
    hit = blocks %in% names(labels)
    blab[hit] = labels[blocks[hit]]
  }
  d[, block_lab := factor(blab[match(block, blocks)], levels = blab)]
  d[, analysis := factor(analysis, levels = names(objs))]
  status_levels = c(
    sprintf("Participates at %d+ %s (%s)", object$r, unit_word, if (object$fdr == "none") "unadjusted" else paste(object$fdr, "FDR")),
    "Not established"
  )
  d[, status := factor(ifelse(replicated %in% TRUE, status_levels[1L], status_levels[2L]), levels = status_levels)]
  has_weight = any(!is.na(d$weight))
  wlev = c("Weight \u2260 0 in the reference fit", "Weight = 0")
  d[, wt := factor(ifelse(!is.na(weight) & weight != 0, wlev[1L], wlev[2L]), levels = wlev)]
  row_labels = stats::setNames(rev(flab), rev(paste(ord$block, ord$feature, sep = "\r")))

  p = ggplot2::ggplot(d, ggplot2::aes(x = x, y = row)) +
    ggplot2::geom_vline(xintercept = 0, colour = inks[["rule"]], linewidth = 0.4) +
    ggplot2::geom_segment(ggplot2::aes(x = 0, xend = x, yend = row), colour = inks[["rule"]], linewidth = 0.4, na.rm = TRUE)
  p = if (has_weight) {
    p + ggplot2::geom_point(ggplot2::aes(fill = status, shape = wt), colour = inks[["ink2"]], size = 1.9, stroke = 0.45, na.rm = TRUE) +
      ggplot2::scale_shape_manual(values = stats::setNames(c(21, 24), wlev), drop = FALSE, name = NULL)
  } else {
    p + ggplot2::geom_point(ggplot2::aes(fill = status), shape = 21, colour = inks[["ink2"]], size = 1.9, stroke = 0.45, na.rm = TRUE)
  }
  p + ggplot2::scale_fill_manual(values = stats::setNames(c(colours[1L], "white"), status_levels), drop = FALSE, name = NULL) +
    ggplot2::scale_y_discrete(labels = row_labels) +
    ggplot2::facet_grid(
      block_lab ~ analysis, scales = "free_y", space = "free_y",
      switch = "y"
    ) +
    ggplot2::guides(fill = ggplot2::guide_legend(order = 1L, override.aes = list(shape = 21))) +
    ggplot2::labs(x = "Held-out correlation with the other blocks\u2019 composite", y = NULL) +
    .mbspls_participation_theme(base_size) +
    ggplot2::theme(
      strip.placement = "outside", strip.text.y.left = ggplot2::element_text(angle = 0, hjust = 1, face = "bold"),
      strip.text.x = ggplot2::element_text(face = "bold"), panel.spacing.y = grid::unit(0.35, "lines"),
      axis.text.y = ggplot2::element_text(size = ggplot2::rel(0.8))
    )
}
