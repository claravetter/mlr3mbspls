#' Held-out association of every pair of blocks at one centre
#'
#' @description
#' Pairwise counterpart of [mbspls_participation()]. For every pair of blocks
#' it tests whether their held-out scores, computed with weights and
#' preprocessing learned elsewhere, are associated at the held-out centre or
#' cohort. Pairwise results show which blocks carry the participation of a
#' block, for example whether two blocks derived from the same measurement
#' participate only through each other.
#'
#' @details
#' For blocks `a` and `b` (in column order, `a` before `b`), the statistic is
#' `s * cor(t_a, t_b)` on the held-out scores, where `s` is the sign of the
#' training correlation of the two scores (1 if it is zero or undefined). Its
#' one-sided p-value compares the statistic with `n_perm` permutations of the
#' held-out score of `b` within the centre (within `strata`, if given),
#' `(1 + #{perm >= obs}) / (n_perm + 1)`. The p-value is exact, conditionally
#' on the training data, when the rows of the centre are exchangeable and the
#' two held-out scores are independent there.
#'
#' @param T_test,T_train (`matrix`) Held-out and training block scores of one
#'   component, one column per block, with the same column names in the same
#'   order.
#' @param n_perm (`integer(1)`) Permutations. Default `9999`.
#' @param strata (`NULL` | `vector`) Optional strata of the held-out rows;
#'   permutations then stay within strata.
#'
#' @return A `data.table` with one row per pair: `a`, `b`, the number of
#'   held-out rows `n`, the training correlation `train_r`, the oriented
#'   held-out correlation `r` and the one-sided p-value `p`. Pairs with a
#'   constant held-out score get `r = NA` and `p = 1`.
#'
#' @seealso [as_mbspls_pairs()], [mbspls_participation()]
#' @examples
#' set.seed(1)
#' f = rnorm(80)
#' T_mat = cbind(A = f + rnorm(80), B = f + rnorm(80), C = rnorm(80))
#' mbspls_pair_heldout_tests(T_mat[41:80, ], T_mat[1:40, ], n_perm = 199L)
#' @export
mbspls_pair_heldout_tests = function(T_test, T_train, n_perm = 9999L, strata = NULL) {
  checkmate::assert_matrix(T_test, mode = "numeric", min.cols = 2L)
  checkmate::assert_matrix(T_train, mode = "numeric", ncols = ncol(T_test))
  blocks = colnames(T_test)
  if (is.null(blocks) || anyDuplicated(blocks) || !identical(blocks, colnames(T_train))) {
    stop("`T_test` and `T_train` need the same unique column names, one per block, in the same order.", call. = FALSE)
  }
  n_perm = checkmate::assert_int(n_perm, lower = 1L, coerce = TRUE)
  n = nrow(T_test)
  if (!is.null(strata) && length(strata) != n) {
    stop("`strata` needs one value per held-out row.", call. = FALSE)
  }
  K = length(blocks)
  idx = utils::combn(K, 2L)
  tr = suppressWarnings(stats::cor(T_train))
  z = lapply(seq_len(K), function(k) .mbspls_participation_std(T_test[, k]))
  perms = if (n >= 3L) .mbspls_permutations(n, n_perm, strata)
  res = data.table::data.table(a = blocks[idx[1L, ]], b = blocks[idx[2L, ]], n = n,
    train_r = tr[cbind(idx[1L, ], idx[2L, ])], r = NA_real_, p = 1)
  for (j in seq_len(ncol(idx))) {
    za = z[[idx[1L, j]]]
    zb = z[[idx[2L, j]]]
    if (is.null(za) || is.null(zb) || n < 3L) next
    s = sign(res$train_r[j])
    if (!is.finite(s) || s == 0) s = 1
    obs = s * sum(za * zb) / (n - 1)
    null = s * drop(crossprod(za, matrix(zb[perms], n, n_perm))) / (n - 1)
    data.table::set(res, i = j, j = "r", value = obs)
    data.table::set(res, i = j, j = "p", value = (1 + sum(null >= obs - 1e-10)) / (n_perm + 1))
  }
  res[]
}


#' Block pairs across centres
#'
#' @description
#' Summarises pairwise held-out tests, such as those of
#' [mbspls_pair_heldout_tests()] at one external cohort or at every held-out
#' centre, with a decision for each pair that is simultaneous over all pairs.
#'
#' @details
#' With `P` pairs and `S` sites, a pair counts at a site if its p-value is at
#' most `alpha / (P * S)`. `r_count` is the number of sites at which the pair
#' counts, and a pair is `significant` if `r_count >= min_sites`. With one
#' site this is the Bonferroni test over the pairs at level `alpha / P`. With
#' several sites, `r_count` is the fixed-threshold replication count, which
#' remains a valid lower confidence bound when the rule tested at each centre
#' is learned from the other centres. `heldout_r` is the size-weighted mean of
#' the held-out correlations, and `lo` and `hi` an approximate 95% Fisher-z
#' interval based on the total `n`. `train_r` is the mean training correlation
#' over the sites.
#'
#' @param pairs (`data.frame`) One row per pair (and site), with columns `a`,
#'   `b`, the oriented held-out correlation `r` and the p-value `p`, and
#'   optionally `n`, `train_r` and `site`. Without `site`, all rows belong to
#'   one cohort. A pair without a row at some site is recorded there with
#'   `p = 1`.
#' @param alpha (`numeric(1)`) Family-wise level over all pairs. Default `0.05`.
#' @param modality (`NULL` | named `character()`) Optional modality of every
#'   block, used to mark pairs within one modality and to group pairs in plots.
#' @param min_sites (`integer(1)`) Number of sites a pair must count at to be
#'   significant. Default `1`.
#'
#' @return An object of class `mbspls_pairs`: a list with `pairs` (one row per
#'   site and pair), `summary` (one row per pair: `a`, `b`, `S`, `n`,
#'   `train_r`, `heldout_r`, `p` for a single site, `lo`, `hi`, `r_count`,
#'   `significant` and, with `modality`, `same_modality`), `alpha`, `level`,
#'   `min_sites`, `modality`, `S` and `P`.
#'
#' @seealso [mbspls_pair_heldout_tests()], [autoplot.mbspls_pairs()]
#' @examples
#' pairs = data.frame(a = c("A", "A", "B"), b = c("B", "C", "C"),
#'   r = c(0.45, 0.02, 0.05), p = c(0.0002, 0.44, 0.38), n = 60, train_r = c(0.5, 0.1, 0.08))
#' as_mbspls_pairs(pairs, modality = c(A = "MRI", B = "MRI", C = "Genotype"))
#' @export
as_mbspls_pairs = function(pairs, alpha = 0.05, modality = NULL, min_sites = 1L) {
  checkmate::assert_data_frame(pairs, min.rows = 1L)
  d = data.table::as.data.table(pairs)
  miss = setdiff(c("a", "b", "r", "p"), names(d))
  if (length(miss)) {
    stop(sprintf("`pairs` needs the column(s) %s.", paste0("'", miss, "'", collapse = ", ")), call. = FALSE)
  }
  alpha = checkmate::assert_number(alpha, lower = 0, upper = 1)
  if (!"site" %in% names(d)) d[, site := "all"]
  if (!"n" %in% names(d)) d[, n := NA_real_]
  if (!"train_r" %in% names(d)) d[, train_r := NA_real_]
  d = d[, list(site = as.character(site), a = as.character(a), b = as.character(b), n = as.numeric(n),
    train_r = as.numeric(train_r), r = as.numeric(r), p = as.numeric(p))]
  if (anyDuplicated(d[, c("site", "a", "b")])) {
    stop("Each pair may appear only once per site.", call. = FALSE)
  }
  keys = unique(d[, c("a", "b")])
  sites = unique(d$site)
  P = nrow(keys)
  S = length(sites)
  min_sites = checkmate::assert_int(min_sites, lower = 1L, upper = S, coerce = TRUE)
  full = data.table::CJ(site = sites, k = seq_len(P))
  full[, `:=`(a = keys$a[k], b = keys$b[k])][, k := NULL]
  d = merge(full, d, by = c("site", "a", "b"), all.x = TRUE, sort = FALSE)
  d[is.na(p), p := 1]
  level = alpha / (P * S)
  summ = d[, {
    ok = is.finite(r)
    w = if (all(is.na(n[ok]))) rep(1, sum(ok)) else n[ok]
    nn = if (all(is.na(n))) NA_real_ else sum(n, na.rm = TRUE)
    hr = if (any(ok)) sum(w * r[ok]) / sum(w) else NA_real_
    list(S = .N, n = nn, train_r = if (all(is.na(train_r))) NA_real_ else mean(train_r, na.rm = TRUE),
      heldout_r = hr, p = if (.N == 1L) p else NA_real_, r_count = sum(p <= level))
  }, by = c("a", "b")]
  ci = function(r, n) if (is.finite(r) && is.finite(n) && n > 3 && abs(r) < 1) tanh(atanh(r) + c(-1, 1) * stats::qnorm(0.975) / sqrt(n - 3)) else c(NA_real_, NA_real_)
  b95 = t(mapply(ci, summ$heldout_r, summ$n))
  summ[, `:=`(lo = b95[, 1L], hi = b95[, 2L], significant = r_count >= min_sites)]
  if (!is.null(modality)) {
    checkmate::assert_character(modality, any.missing = FALSE, names = "unique")
    blocks = unique(c(keys$a, keys$b))
    if (!all(blocks %in% names(modality))) {
      stop("`modality` must name every block.", call. = FALSE)
    }
    summ[, same_modality := unname(modality[a] == modality[b])]
  }
  data.table::setcolorder(summ, intersect(c("a", "b", "S", "n", "train_r", "heldout_r", "p", "lo", "hi", "r_count", "significant", "same_modality"), names(summ)))
  structure(list(pairs = d[], summary = summ[], alpha = alpha, level = level, min_sites = min_sites,
    modality = modality, S = S, P = P), class = "mbspls_pairs")
}


#' @export
print.mbspls_pairs = function(x, ...) {
  blocks = unique(c(x$summary$a, x$summary$b))
  cat(sprintf("Held-out block pairs: %d pairs of %d blocks, %d site%s.\n", x$P, length(blocks), x$S, if (x$S == 1L) "" else "s"))
  rule = if (x$S == 1L) {
    sprintf("p <= %s (alpha / %d pairs)", format(signif(x$level, 3)), x$P)
  } else {
    sprintf("p <= %s at %d or more of %d sites", format(signif(x$level, 3)), x$min_sites, x$S)
  }
  cat(sprintf("Significant, simultaneously over the pairs (%s): %d.\n", rule, sum(x$summary$significant)))
  cols = intersect(c("a", "b", "heldout_r", "p", "r_count", "significant", "same_modality"), names(x$summary))
  print(x$summary[, cols, with = FALSE])
  invisible(x)
}


#' Plot block pairs
#'
#' @description
#' Draws the result of [as_mbspls_pairs()]. `type = "network"` shows the
#' blocks as nodes and the held-out correlation of every pair as an edge:
#' solid if the pair is significant simultaneously over all pairs, dashed
#' otherwise, with the width following the correlation. `type = "pairs"`
#' shows, for every pair, the training correlation (open) and the held-out
#' correlation (filled, with its 95% interval), grouped by the kind of pair
#' when the blocks have modalities.
#'
#' @param object (`mbspls_pairs`) The result to plot.
#' @param type (`character(1)`) `"network"` (default) or `"pairs"`.
#' @param labels (`NULL` | named `character()`) Display names of the blocks.
#'   Line breaks are allowed; in `type = "pairs"` they are replaced by spaces.
#' @param layout (`NULL` | `data.frame`) Node positions for `type = "network"`,
#'   with columns `block`, `x` and `y`. Default: the blocks on a circle,
#'   grouped by modality.
#' @param groups (`NULL` | `data.frame`) Groups of pairs for `type = "pairs"`,
#'   with columns `a`, `b` and `group`; the order of first appearance sets the
#'   order of the panels. Default: "Same modality" and "<modality> with
#'   <modality>" when the object has modalities, otherwise no grouping.
#' @param legend_labels (`NULL` | `character()`) Legend entries: for
#'   `type = "network"` the significant and the other pairs (two entries), for
#'   `type = "pairs"` the training, the significant held-out and the other
#'   held-out values (three entries).
#' @param colours (`NULL` | `character()`) Node colours, one per modality in the
#'   order of their first appearance in `modality` (or one colour without
#'   modalities).
#' @param base_size (`numeric(1)`) Base font size in points. Default `9.5`.
#' @param title,subtitle (`NULL` | `character(1)`) Plot title and subtitle.
#' @param ... Not used.
#'
#' @return A [ggplot2::ggplot()] object.
#' @seealso [as_mbspls_pairs()]
#' @examples
#' pairs = data.frame(a = c("A", "A", "B"), b = c("B", "C", "C"),
#'   r = c(0.45, 0.02, 0.05), p = c(0.0002, 0.44, 0.38), n = 60, train_r = c(0.5, 0.1, 0.08))
#' x = as_mbspls_pairs(pairs, modality = c(A = "MRI", B = "MRI", C = "Genotype"))
#' ggplot2::autoplot(x)
#' ggplot2::autoplot(x, type = "pairs")
#' @export
autoplot.mbspls_pairs = function(
  object,
  type = c("network", "pairs"),
  labels = NULL,
  layout = NULL,
  groups = NULL,
  legend_labels = NULL,
  colours = NULL,
  base_size = 9.5,
  title = NULL,
  subtitle = NULL,
  ...
) {
  type = match.arg(type)
  s = data.table::copy(object$summary)
  blocks = unique(c(s$a, s$b))
  lab = stats::setNames(blocks, blocks)
  if (!is.null(labels)) {
    checkmate::assert_character(labels, names = "unique")
    lab[intersect(names(labels), blocks)] = labels[intersect(names(labels), blocks)]
  }
  pal = .mbspls_pairs_palette()
  lvl = format(signif(object$level, 3))
  if (type == "network") {
    return(.mbspls_pairs_network(object, s, blocks, lab, layout, legend_labels, colours, base_size, title, subtitle, pal, lvl))
  }
  .mbspls_pairs_dumbbell(object, s, lab, groups, legend_labels, base_size, title, subtitle, pal, lvl)
}

.mbspls_pairs_palette = function() {
  c(ink = "#0b0b0b", ink2 = "#52514e", muted = "#898781", grid = "#e1e0d9", blue = "#2a78d6", orange = "#eb6834")
}

.mbspls_pairs_theme = function(base_size, pal) {
  ggplot2::theme_minimal(base_size = base_size) +
    ggplot2::theme(
      text = ggplot2::element_text(colour = pal[["ink"]]), axis.text = ggplot2::element_text(colour = pal[["ink2"]]),
      panel.grid.minor = ggplot2::element_blank(), panel.grid.major = ggplot2::element_line(colour = pal[["grid"]], linewidth = 0.3),
      plot.title.position = "plot", plot.title = ggplot2::element_text(face = "bold", size = ggplot2::rel(1.1)),
      plot.subtitle = ggplot2::element_text(colour = pal[["ink2"]], size = ggplot2::rel(0.95)),
      legend.position = "bottom", legend.title = ggplot2::element_blank(), legend.text = ggplot2::element_text(colour = pal[["ink2"]])
    )
}

# Default node layout: a circle starting at the top, blocks grouped by modality.
.mbspls_pairs_layout = function(blocks, modality) {
  ord = if (is.null(modality)) blocks else blocks[order(match(modality[blocks], unique(modality[blocks])), match(blocks, blocks))]
  ang = pi / 2 - 2 * pi * (seq_along(ord) - 1) / length(ord)
  data.table::data.table(block = ord, x = cos(ang), y = sin(ang))
}

# Edge-label positions: midpoints; labels closer than `thr` are moved apart along their own edges.
.mbspls_pairs_edge_labels = function(e, thr) {
  e[, `:=`(mx = (x + x_b) / 2, my = (y + y_b) / 2)]
  n = nrow(e)
  if (n < 2L) {
    return(e)
  }
  for (i in seq_len(n - 1L)) {
    for (j in (i + 1L):n) {
      if (sqrt((e$mx[i] - e$mx[j])^2 + (e$my[i] - e$my[j])^2) < thr) {
        for (k in c(i, j)) {
          tt = if (k == i) 0.35 else 0.65
          data.table::set(e, i = k, j = "mx", value = e$x[k] + tt * (e$x_b[k] - e$x[k]))
          data.table::set(e, i = k, j = "my", value = e$y[k] + tt * (e$y_b[k] - e$y[k]))
        }
      }
    }
  }
  e
}

.mbspls_pairs_network = function(object, s, blocks, lab, layout, legend_labels, colours, base_size, title, subtitle, pal, lvl) {
  modality = object$modality
  pos = if (is.null(layout)) .mbspls_pairs_layout(blocks, modality) else data.table::as.data.table(layout)
  if (!all(c("block", "x", "y") %in% names(pos)) || !all(blocks %in% pos$block)) {
    stop("`layout` needs the columns 'block', 'x' and 'y' and a row for every block.", call. = FALSE)
  }
  pos = pos[match(blocks, block), list(block, x = as.numeric(x), y = as.numeric(y))]
  rx = diff(range(pos$x))
  ry = diff(range(pos$y))
  span = max(rx, ry, 1e-8)
  cx = mean(range(pos$x))
  pos[, `:=`(label = unname(lab[block]),
    hj = ifelse(x < cx - 0.05 * span, 1, ifelse(x > cx + 0.05 * span, 0, 0.5)))]
  pos[, `:=`(nx = ifelse(hj == 1, -0.059 * span, ifelse(hj == 0, 0.059 * span, 0)),
    ny = ifelse(hj == 0.5, ifelse(y <= mean(range(y)), -0.12 * span, 0.12 * span), 0))]
  if (is.null(modality)) {
    pos[, modality := "Block"]
    fills = stats::setNames(if (is.null(colours)) pal[["blue"]] else colours[1L], "Block")
  } else {
    mods = unique(unname(modality[blocks]))
    pos[, modality := factor(unname(modality[block]), levels = mods)]
    cols = if (is.null(colours)) c(pal[["blue"]], pal[["orange"]], pal[["ink2"]], "#3a9a5b", "#8a5fbf", "#c9a227") else colours
    if (length(cols) < length(mods)) stop("`colours` needs one colour per modality.", call. = FALSE)
    fills = stats::setNames(cols[seq_along(mods)], mods)
  }
  leg = if (is.null(legend_labels)) {
    c(if (object$S == 1L) {
      sprintf("Held-out p \u2264 %s (simultaneous)", lvl)
    } else {
      sprintf("p \u2264 %s at \u2265 %d of %d sites", lvl, object$min_sites, object$S)
    }, "Not significant")
  } else {
    legend_labels
  }
  if (length(leg) != 2L) stop("`legend_labels` needs two entries for type = \"network\".", call. = FALSE)
  e = merge(s, pos[, list(a = block, x, y)], by = "a")
  e = merge(e, pos[, list(b = block, x_b = x, y_b = y)], by = "b")
  e = .mbspls_pairs_edge_labels(e, thr = 0.15 * span)
  e[, `:=`(edge = factor(ifelse(significant, leg[1L], leg[2L]), levels = leg), w = pmax(heldout_r, 0.02))]
  e[!is.finite(w), w := 0.02]
  lab_args = list(data = e, mapping = ggplot2::aes(x = mx, y = my, label = sprintf("%.2f", heldout_r), colour = edge), size = base_size / 3.3,
    label.padding = grid::unit(1.2, "pt"), fill = "white", show.legend = FALSE)
  if ("border.colour" %in% names(formals(ggplot2::geom_label))) lab_args$border.colour = NA else lab_args$label.size = 0
  p = ggplot2::ggplot() +
    ggplot2::geom_segment(data = e, ggplot2::aes(x = x, y = y, xend = x_b, yend = y_b, linewidth = w, colour = edge, linetype = edge), lineend = "round") +
    do.call(ggplot2::geom_label, lab_args) +
    ggplot2::geom_point(data = pos, ggplot2::aes(x = x, y = y, fill = modality), shape = 21, size = 5.2, colour = "white", stroke = 1.2) +
    ggplot2::geom_text(data = pos, ggplot2::aes(x = x + nx, y = y + ny, label = label, hjust = hj), size = base_size / 2.97, colour = pal[["ink"]], lineheight = 0.9) +
    ggplot2::scale_linewidth_continuous(range = c(0.3, 2.6), guide = "none") +
    ggplot2::scale_colour_manual(values = stats::setNames(c(pal[["ink"]], pal[["muted"]]), leg), drop = FALSE) +
    ggplot2::scale_linetype_manual(values = stats::setNames(c("solid", "22"), leg), drop = FALSE) +
    ggplot2::scale_fill_manual(values = fills, breaks = names(fills)) +
    ggplot2::coord_cartesian(xlim = range(pos$x) + c(-0.43, 0.43) * span, ylim = range(pos$y) + c(-0.18, 0.1) * span, clip = "off") +
    ggplot2::labs(title = title, subtitle = subtitle) +
    ggplot2::guides(fill = if (is.null(modality)) "none" else ggplot2::guide_legend(order = 1, override.aes = list(size = 4)),
      colour = ggplot2::guide_legend(order = 2), linetype = ggplot2::guide_legend(order = 2)) +
    .mbspls_pairs_theme(base_size, pal) +
    ggplot2::theme(axis.text = ggplot2::element_blank(), axis.title = ggplot2::element_blank(), panel.grid.major = ggplot2::element_blank(),
      legend.box = "vertical", legend.spacing.y = grid::unit(0, "pt"), legend.margin = ggplot2::margin(0, 0, 0, 0))
  p
}

.mbspls_pairs_dumbbell = function(object, s, lab, groups, legend_labels, base_size, title, subtitle, pal, lvl) {
  flat = gsub("\\s*\n\\s*", " ", lab)
  s[, pair := paste(flat[a], "\u2013", flat[b])]
  modality = object$modality
  if (!is.null(groups)) {
    g = data.table::as.data.table(groups)
    if (!all(c("a", "b", "group") %in% names(g))) stop("`groups` needs the columns 'a', 'b' and 'group'.", call. = FALSE)
    g = g[, list(a = as.character(a), b = as.character(b), group = as.character(group))]
    s = merge(s, g, by = c("a", "b"), all.x = TRUE, sort = FALSE)
    if (anyNA(s$group)) stop("`groups` needs a row for every pair.", call. = FALSE)
    glev = unique(g$group)
  } else if (!is.null(modality)) {
    mods = unique(unname(modality))
    s[, group := ifelse(modality[a] == modality[b], "Same modality", {
      ma = match(modality[a], mods)
      mb = match(modality[b], mods)
      paste(mods[pmin(ma, mb)], "with", mods[pmax(ma, mb)])
    })]
    cross = unique(s$group[s$group != "Same modality"])
    glev = c(if ("Same modality" %in% s$group) "Same modality", cross[order(match(sub(" with .*", "", cross), mods), match(sub(".* with ", "", cross), mods))])
  } else {
    s[, group := ""]
    glev = ""
  }
  s[, group := factor(group, levels = glev)]
  ord = s[order(group, -heldout_r), pair]
  s[, pair := factor(pair, levels = rev(ord))]
  leg = if (is.null(legend_labels)) {
    c("Training (in-sample)", sprintf("Held-out, p \u2264 %s", lvl), "Held-out, not significant")
  } else {
    legend_labels
  }
  if (length(leg) != 3L) stop("`legend_labels` needs three entries for type = \"pairs\".", call. = FALSE)
  long = rbind(s[is.finite(train_r), list(pair, group, r = abs(train_r), what = leg[1L])],
    s[significant == TRUE, list(pair, group, r = heldout_r, what = leg[2L])],
    s[significant == FALSE, list(pair, group, r = heldout_r, what = leg[3L])])
  long[, what := factor(what, levels = leg)]
  xr = range(c(0, s$lo, s$hi, abs(s$train_r), s$heldout_r), na.rm = TRUE)
  xb = pretty(xr, n = 4)
  xr = range(xb) # the axis ends at labelled breaks and includes every interval
  p = ggplot2::ggplot(s, ggplot2::aes(y = pair)) +
    ggplot2::geom_vline(xintercept = 0, colour = pal[["muted"]], linewidth = 0.3) +
    ggplot2::geom_segment(data = s[is.finite(train_r)], ggplot2::aes(x = abs(train_r), xend = heldout_r, yend = pair), colour = pal[["grid"]], linewidth = 1.1) +
    ggplot2::geom_errorbar(data = s[is.finite(lo)], ggplot2::aes(xmin = lo, xmax = hi), width = 0, orientation = "y", colour = pal[["muted"]], linewidth = 0.45) +
    ggplot2::geom_point(data = long, ggplot2::aes(x = r, shape = what, colour = what), size = 2.3, stroke = 0.8) +
    ggplot2::scale_shape_manual(values = stats::setNames(c(1, 16, 16), leg), drop = FALSE) +
    ggplot2::scale_colour_manual(values = stats::setNames(c(pal[["ink2"]], pal[["ink"]], pal[["muted"]]), leg), drop = FALSE) +
    ggplot2::scale_x_continuous(limits = xr, breaks = xb) +
    ggplot2::labs(title = title, subtitle = subtitle, x = "Correlation", y = NULL) +
    ggplot2::guides(shape = ggplot2::guide_legend(ncol = 1), colour = ggplot2::guide_legend(ncol = 1)) +
    .mbspls_pairs_theme(base_size, pal) +
    ggplot2::theme(panel.grid.major.y = ggplot2::element_blank(), legend.margin = ggplot2::margin(0, 0, 0, 0))
  if (length(glev) > 1L || nzchar(glev[1L])) {
    p = p + ggplot2::facet_grid(group ~ ., scales = "free_y", space = "free_y") +
      ggplot2::theme(strip.text.y = ggplot2::element_text(angle = 0, hjust = 0, face = "bold", colour = pal[["ink2"]], size = ggplot2::rel(0.84)),
        panel.spacing.y = grid::unit(6, "pt"))
  }
  p
}
