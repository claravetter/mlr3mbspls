#' Plot held-out block participation
#'
#' @description
#' Draws the result of [mbspls_participation()] or [as_mbspls_participation()]
#' with \pkg{ggplot2}.
#' * `type = "centres"` shows one row per block. Each point is the block's
#'   held-out correlation at one centre. It is filled if that centre's p-value
#'   is at most `alpha`, and its area is proportional to the centre's size.
#'   Columns on the right give the lower bounds for the number of centres at
#'   which the block participates, for `object` and for every analysis in
#'   `compare`: by default the fixed-threshold counts (`rule = "count"`), or
#'   the partial conjunction bounds (`rule = "partial_conjunction"`); see
#'   [mbspls_participation()] for when each is valid.
#' * `type = "heldout"` contrasts, for every block, the correlation with the
#'   composite of the other blocks in the training centres (the mean over
#'   folds) with the correlation at the held-out centres (the mean weighted by
#'   centre size). A large drop marks an association that holds only
#'   in-sample. With `compare`, every analysis is drawn in its own colour,
#'   next to the others in each block's row, with training correlations as
#'   open and held-out correlations as filled circles.
#'
#' Both types order the blocks by their size-weighted mean held-out
#' correlation, with the strongest at the top, unless `order` is given (for
#' example to draw two analyses with the same order).
#'
#' @param object (`mbspls_participation`) The result to plot.
#' @param type (`character(1)`) `"centres"` (default) or `"heldout"`.
#' @param bound (`character(1)`) The bound printed by `type = "centres"`:
#'   `"simultaneous"` (default; adjusted for the number of blocks),
#'   `"unadjusted"` or `"none"`. Which bound is printed is set by `rule`.
#' @param compare (`NULL` | named `list()`) Further `mbspls_participation`
#'   objects, for example sensitivity analyses or another cohort, named by
#'   their labels. With `type = "centres"` their bounds are printed as further
#'   columns; with `type = "heldout"` they are drawn next to `object`. Blocks
#'   are matched by name, and blocks that `object` lacks are left out.
#' @param label (`NULL` | `character(1)`) Name of `object`: the header of its
#'   bound column (`type = "centres"`, default `"Lower\nbound"`), or its legend
#'   entry with `compare` and `type = "heldout"` (default `"Main"`).
#' @param labels (`NULL` | named `character()`) Display names of the blocks,
#'   named by block. Other blocks are shown with underscores replaced by
#'   spaces.
#' @param order (`NULL` | `character()`) All blocks of `object`, from top to
#'   bottom. Default `NULL` orders them by held-out correlation.
#' @param limits (`NULL` | `numeric(2)`) Range of correlations on the axis of
#'   `type = "centres"`, for example to align several plots. It is widened to
#'   include zero and all points. Default `NULL` uses the range of the data.
#' @param colours (`character()`) Colours of the centres with `p <= alpha`
#'   and of the others (`type = "centres"`), or of the held-out and of the
#'   training correlations (`type = "heldout"`). With `compare` and
#'   `type = "heldout"`, one colour per analysis, `object` first.
#' @param unit (`character(2)`) Singular and plural of the held-out unit as
#'   used in the labels, for example `c("site", "sites")` or
#'   `c("center", "centers")`. Default `c("centre", "centres")`.
#' @param base_size (`numeric(1)`) Base font size in points. Default `11`.
#' @param rule (`character(1)`) The bound printed by `type = "centres"`:
#'   `"count"` (default), the number of centres with `p <= alpha / S`
#'   (`r_count`) or `p <= alpha / (K S)` (`r_count_simultaneous`), which is
#'   valid for leave-one-centre-out results; or `"partial_conjunction"`,
#'   `r_hat` or `r_hat_simultaneous`, which is valid when the null hypotheses
#'   do not depend on the training data. Objects created before the counts
#'   existed are summarised again from their `centres`.
#' @param ... Not used.
#'
#' @return A [ggplot2::ggplot()] object, which can be modified further.
#'
#' @seealso [mbspls_participation()], [as_mbspls_participation()]
#' @examples
#' centres = data.frame(
#'   site = rep(c("a", "b", "c", "d"), each = 2),
#'   block = rep(c("Clinical_scales", "MRI"), 4),
#'   n = rep(c(120, 80, 60, 30), each = 2),
#'   r = c(0.42, 0.05, 0.38, -0.02, 0.29, 0.11, 0.35, 0.08),
#'   train_r = c(0.45, 0.30, 0.44, 0.31, 0.47, 0.28, 0.46, 0.30),
#'   p = c(0.0001, 0.31, 0.0013, 0.55, 0.012, 0.20, 0.03, 0.34)
#' )
#' res = as_mbspls_participation(centres)
#' ggplot2::autoplot(res)
#' ggplot2::autoplot(res, type = "heldout", labels = c(MRI = "Structural MRI"))
#' @export
autoplot.mbspls_participation = function(
  object,
  type = c("centres", "heldout"),
  bound = c("simultaneous", "unadjusted", "none"),
  compare = NULL,
  label = NULL,
  labels = NULL,
  order = NULL,
  limits = NULL,
  colours = c("#2a78d6", "#898781"),
  unit = c("centre", "centres"),
  base_size = 11,
  rule = c("count", "partial_conjunction"),
  ...
) {
  checkmate::assert_class(object, "mbspls_participation")
  type = match.arg(type)
  bound = match.arg(bound)
  rule = match.arg(rule)
  if (!is.null(compare)) {
    checkmate::assert_list(compare, min.len = 1L)
    if (is.null(names(compare)) || any(!nzchar(names(compare)))) {
      stop("`compare` must be a named list.", call. = FALSE)
    }
    if (!all(vapply(compare, inherits, logical(1L), what = "mbspls_participation"))) {
      stop("Every element of `compare` must be an `mbspls_participation` object.", call. = FALSE)
    }
  }
  checkmate::assert_string(label, null.ok = TRUE)
  checkmate::assert_character(labels, any.missing = FALSE, null.ok = TRUE)
  if (!is.null(labels) && is.null(names(labels))) {
    stop("`labels` must be named by block.", call. = FALSE)
  }
  checkmate::assert_numeric(limits, len = 2L, any.missing = FALSE, finite = TRUE, sorted = TRUE, null.ok = TRUE)
  checkmate::assert_character(colours, min.len = 2L, any.missing = FALSE)
  checkmate::assert_character(unit, len = 2L, any.missing = FALSE, min.chars = 1L)
  checkmate::assert_number(base_size, lower = 1)

  centres = data.table::as.data.table(object$centres)
  blocks = as.character(object$blocks$block)
  lab = .nice_label(blocks)
  if (!is.null(labels)) {
    hit = blocks %in% names(labels)
    lab[hit] = labels[blocks[hit]]
  }
  if (anyDuplicated(lab)) {
    stop("The block labels must be unique.", call. = FALSE)
  }
  if (is.null(order)) {
    levels_y = lab[.mbspls_participation_block_order(centres, blocks)]
  } else {
    checkmate::assert_character(order, any.missing = FALSE, unique = TRUE)
    if (!setequal(order, blocks)) {
      stop("`order` must contain every block of `object` exactly once.", call. = FALSE)
    }
    levels_y = lab[match(rev(order), blocks)]
  }
  label_of = function(b) factor(lab[match(b, blocks)], levels = levels_y)

  if (type == "heldout") {
    return(.mbspls_participation_plot_heldout(
      centres, label_of, levels_y, colours, unit, base_size, compare, label
    ))
  }

  cols = c(list(object), compare)
  heads = c(if (is.null(label)) "Lower\nbound" else label, names(compare))
  if (bound == "none") {
    cols = list()
    heads = character()
  }
  .mbspls_participation_plot_centres(
    centres, object$alpha, label_of, levels_y, cols, heads, bound, colours, unit, base_size, limits,
    rule = rule
  )
}

# The printed bound of every block of one analysis. Objects created before the
# counts existed are summarised again from their centre-level results.
.mbspls_participation_bound_values = function(obj, bound, rule) {
  b = data.table::as.data.table(obj$blocks)
  col = paste0(if (rule == "count") "r_count" else "r_hat", if (bound == "simultaneous") "_simultaneous" else "")
  if (!col %in% names(b)) {
    b = data.table::as.data.table(as_mbspls_participation(obj$centres, alpha = obj$alpha)$blocks)
  }
  data.table::data.table(block = as.character(b$block), value = as.character(b[[col]]))
}

.mbspls_participation_capitalise = function(x) {
  paste0(toupper(substr(x, 1L, 1L)), substr(x, 2L, nchar(x)))
}

# Neutral inks shared by both plot types.
.mbspls_participation_inks = function() {
  c(ink = "#0b0b0b", ink2 = "#52514e", rule = "#c3c2b7")
}

.mbspls_participation_theme = function(base_size, top = 6) {
  ggplot2::theme_minimal(base_size = base_size) +
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(),
      legend.position = "bottom",
      legend.title = ggplot2::element_text(size = ggplot2::rel(0.9)),
      plot.title.position = "plot",
      plot.margin = ggplot2::margin(top, 8, 4, 4)
    )
}

# Indices of `blocks` in increasing order of the size-weighted mean held-out
# correlation; ggplot2 draws the first level of a discrete y axis at the bottom.
.mbspls_participation_block_order = function(centres, blocks) {
  d = centres[!is.na(r), list(h = sum(n * r) / sum(n)), by = "block"]
  order(d$h[match(blocks, d$block)], na.last = FALSE)
}

# About three legend sizes, chosen from 1, 2, 3 and 5 times a power of ten.
.mbspls_participation_size_breaks = function(n) {
  n = n[is.finite(n)]
  if (!length(n)) {
    return(ggplot2::waiver())
  }
  nice = sort(as.vector(outer(c(1, 2, 3, 5), 10^(0:7))))
  cand = nice[nice >= min(n) & nice <= max(n)]
  if (length(cand) <= 3L) {
    return(if (length(cand)) cand else unique(signif(range(n), 2L)))
  }
  ends = cand[c(1L, length(cand))]
  mid = cand[which.min(abs(log(cand) - mean(log(ends))))]
  unique(c(ends[1L], mid, ends[2L]))
}

.mbspls_participation_plot_centres = function(
  centres, alpha, label_of, levels_y, cols, heads, bound, colours, unit, base_size, limits = NULL,
  rule = "count"
) {
  inks = .mbspls_participation_inks()
  d = centres[!is.na(r)]
  a = format(alpha)
  sig_levels = c(sprintf("p \u2264 %s at that %s", a, unit[1L]), sprintf("p > %s", a))
  d[, label := label_of(block)]
  d[, sig := factor(ifelse(p <= alpha, sig_levels[1L], sig_levels[2L]), levels = sig_levels)]

  lo = min(0, d$r, limits)
  hi = max(0, d$r, limits)
  span = hi - lo
  if (!is.finite(span) || span <= 0) span = 1
  breaks = pretty(c(lo, hi))
  breaks = breaks[breaks >= lo - 1e-8 * span & breaks <= hi + 1e-8 * span]
  n_cols = length(cols)
  x_col = hi + span * (0.11 + 0.17 * (seq_len(n_cols) - 1L))
  x_max = if (n_cols) x_col[n_cols] + 0.05 * span else hi + 0.05 * span
  fs = base_size / ggplot2::.pt

  p = ggplot2::ggplot(d, ggplot2::aes(x = r, y = label)) +
    ggplot2::geom_vline(xintercept = 0, colour = inks[["rule"]], linewidth = 0.4) +
    ggplot2::geom_point(
      ggplot2::aes(size = n, colour = sig, shape = sig),
      alpha = 0.9, stroke = 0.6
    ) +
    ggplot2::scale_colour_manual(values = stats::setNames(colours, sig_levels), drop = FALSE) +
    ggplot2::scale_shape_manual(values = stats::setNames(c(16, 1), sig_levels), drop = FALSE) +
    ggplot2::scale_size_area(
      max_size = 4.2, breaks = .mbspls_participation_size_breaks(d$n),
      name = paste(.mbspls_participation_capitalise(unit[1L]), "n")
    ) +
    ggplot2::scale_x_continuous(limits = c(lo - 0.11 * span, x_max), breaks = breaks) +
    ggplot2::scale_y_discrete(drop = FALSE) +
    ggplot2::guides(
      colour = ggplot2::guide_legend(order = 1L, title = NULL),
      shape = ggplot2::guide_legend(order = 1L, title = NULL),
      size = ggplot2::guide_legend(order = 2L)
    ) +
    ggplot2::labs(
      x = sprintf("Correlation with the other blocks\u2019 composite at the held-out %s", unit[1L]),
      y = NULL
    )

  top = 6
  if (n_cols) {
    bd = data.table::rbindlist(lapply(seq_len(n_cols), function(j) {
      v = .mbspls_participation_bound_values(cols[[j]], bound, rule)
      v[, col := j]
      v
    }))
    bd = bd[!is.na(label_of(block))]
    bd[, `:=`(label = label_of(block), x = x_col[col])]
    p = p +
      ggplot2::geom_text(
        data = bd[col == 1L], ggplot2::aes(x = x, y = label, label = value),
        inherit.aes = FALSE, size = 0.95 * fs, colour = inks[["ink"]], fontface = "bold"
      ) +
      ggplot2::annotate(
        "text", x = x_col, y = length(levels_y) + 0.8, label = heads, vjust = 0,
        size = 0.8 * fs, colour = inks[["ink2"]], lineheight = 0.85
      ) +
      ggplot2::coord_cartesian(clip = "off")
    if (n_cols > 1L) {
      p = p + ggplot2::geom_text(
        data = bd[col > 1L], ggplot2::aes(x = x, y = label, label = value),
        inherit.aes = FALSE, size = 0.95 * fs, colour = inks[["ink2"]]
      )
    }
    top = 6 + 10 * max(lengths(strsplit(heads, "\n", fixed = TRUE)))
  }
  p + .mbspls_participation_theme(base_size, top = top)
}

# Training correlation (mean over folds) and held-out correlation (weighted by
# centre size) of every block.
.mbspls_participation_heldout_summary = function(centres) {
  centres[, list(
    train = if (all(is.na(train_r))) NA_real_ else mean(train_r, na.rm = TRUE),
    heldout = if (all(is.na(r))) NA_real_ else sum((n * r)[!is.na(r)]) / sum(n[!is.na(r)])
  ), by = "block"]
}

.mbspls_participation_plot_heldout = function(
  centres, label_of, levels_y, colours, unit, base_size, compare = NULL, label = NULL
) {
  objs = c(list(centres), lapply(compare, function(o) data.table::as.data.table(o$centres)))
  names(objs) = c(if (is.null(label)) "Main" else label, names(compare))
  if (any(vapply(objs, function(x) all(is.na(x$train_r)), logical(1L)))) {
    stop("`type = \"heldout\"` needs the training correlations `train_r`.", call. = FALSE)
  }
  inks = .mbspls_participation_inks()
  kinds = c(paste("Training", unit[2L]), paste("Held-out", unit[1L]))
  theme = .mbspls_participation_theme(base_size)
  x_lab = "Correlation with the other blocks\u2019 composite"

  if (length(objs) == 1L) {
    d = .mbspls_participation_heldout_summary(centres)
    d[, label := label_of(block)]
    long = data.table::rbindlist(list(
      d[, list(label, value = train, kind = kinds[1L])],
      d[, list(label, value = heldout, kind = kinds[2L])]
    ))
    long[, kind := factor(kind, levels = kinds)]
    return(ggplot2::ggplot(d, ggplot2::aes(y = label)) +
      ggplot2::geom_segment(
        ggplot2::aes(x = heldout, xend = train, yend = label),
        colour = inks[["rule"]], linewidth = 0.6, na.rm = TRUE
      ) +
      ggplot2::geom_point(
        data = long, ggplot2::aes(x = value, colour = kind, fill = kind),
        shape = 21, size = 2.3, stroke = 0.7, na.rm = TRUE
      ) +
      ggplot2::scale_colour_manual(values = stats::setNames(colours[c(2L, 1L)], kinds)) +
      ggplot2::scale_fill_manual(values = stats::setNames(c("white", colours[1L]), kinds)) +
      ggplot2::scale_y_discrete(drop = FALSE) +
      ggplot2::expand_limits(x = 0) +
      ggplot2::labs(x = x_lab, y = NULL, colour = NULL, fill = NULL) +
      theme)
  }

  # Several analyses: one colour each, side by side within each block's row.
  m = length(objs)
  if (length(colours) < m) {
    stop("`colours` needs one colour per analysis when `type = \"heldout\"` is used with `compare`.", call. = FALSE)
  }
  off = seq(0.18, -0.18, length.out = m)
  d = data.table::rbindlist(lapply(seq_len(m), function(j) {
    sj = .mbspls_participation_heldout_summary(objs[[j]])
    sj[, label := label_of(block)]
    sj = sj[!is.na(label)]
    sj[, `:=`(analysis = names(objs)[j], y = as.numeric(label) + off[j])]
    sj
  }))
  d[, analysis := factor(analysis, levels = names(objs))]
  long = data.table::rbindlist(list(
    d[, list(y, analysis, value = train, kind = kinds[1L])],
    d[, list(y, analysis, value = heldout, kind = kinds[2L])]
  ))
  long[, kind := factor(kind, levels = kinds)]
  ggplot2::ggplot(d, ggplot2::aes(y = y)) +
    ggplot2::geom_segment(
      ggplot2::aes(x = heldout, xend = train, yend = y),
      colour = inks[["rule"]], linewidth = 0.5, na.rm = TRUE
    ) +
    ggplot2::geom_point(
      data = long, ggplot2::aes(x = value, colour = analysis, shape = kind),
      size = 2, stroke = 0.8, na.rm = TRUE
    ) +
    ggplot2::scale_colour_manual(values = stats::setNames(colours[seq_len(m)], names(objs))) +
    ggplot2::scale_shape_manual(values = stats::setNames(c(1, 16), kinds)) +
    ggplot2::scale_y_continuous(
      breaks = seq_along(levels_y), labels = levels_y,
      expand = ggplot2::expansion(add = 0.5)
    ) +
    ggplot2::guides(
      shape = ggplot2::guide_legend(order = 1L, override.aes = list(colour = inks[["ink2"]])),
      colour = ggplot2::guide_legend(order = 2L)
    ) +
    ggplot2::expand_limits(x = 0) +
    ggplot2::labs(x = x_lab, y = NULL, colour = NULL, shape = NULL) +
    theme
}
