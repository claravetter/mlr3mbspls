#' Held-out block participation test across centres
#'
#' @description
#' Tests, for every block of an MB-sPLS component, whether the block
#' *participates* at each centre: whether its held-out score is associated
#' with the composite of the other blocks' held-out scores when the weights are
#' learned at the other centres (leave one centre out). Counting the centres
#' whose p-value passes a Bonferroni threshold then gives a lower confidence
#' bound for the number of centres at which the block participates. A partial
#' conjunction bound, which is valid when the null hypotheses do not depend on
#' the training data, is reported as well.
#'
#' @details
#' For held-out centre `s`, `learner` is trained on all other centres. Let
#' `z_b` be the held-out scores of block `b`, standardised within centre `s`,
#' and `r_kj` the correlations between the training scores of blocks `k` and
#' `j`. The statistic for block `k` is the correlation of `z_k` with the
#' composite `sum_{j != k} sign(r_kj) z_j`. Its one-sided p-value compares the
#' statistic with `n_perm` permutations of `z_k` within the centre (within
#' `strata`, if given), as `(1 + #{permuted >= observed}) / (n_perm + 1)`
#' (see [mb_permutation_pvalue()]).
#'
#' Conditional on the training centres, this p-value is exact when the rows of
#' the held-out centre are exchangeable (within strata) and block `k` is
#' independent of the other blocks there, whatever the training data, the
#' tuning or the sparsity. Two conditions on the preprocessing make this hold:
#' * The held-out scores of block `k` must depend only on block `k`'s data.
#'   Impute within blocks (`impute_within_blocks = TRUE` in
#'   [mbspls_graph_learner()], or [impute_knn_graph()] with `blocks`) and keep
#'   covariate corrections within blocks.
#' * Any centring of the held-out data must be constant within the centre. The
#'   statistic is invariant to it, so centre effects that shift a block's mean
#'   do not create false participation.
#'
#' The p-values of different centres are dependent, because the training sets
#' overlap. Moreover, whether block `k` participates at centre `s` depends on
#' the weights learned at the other centres, so the set of centres at which
#' the null hypothesis holds is itself random. Two bounds are reported
#' ([mbspls_partial_conjunction()]):
#' * `r_count` is the number of centres with `p <= alpha / S`. It is a lower
#'   `1 - alpha` confidence bound for the number of centres at which the block
#'   participates in the rule fitted without that centre, whatever the
#'   dependence between centres and although the hypotheses depend on the
#'   data: counting a centre wrongly requires a true-null p-value at most the
#'   fixed threshold `alpha / S`, which has conditional probability at most
#'   `alpha / S` at each centre. `r_count_simultaneous` uses `alpha / (K S)`
#'   for `K` blocks, so that the counts of all blocks hold jointly.
#' * `r_hat` is the Bonferroni partial conjunction bound, the largest `r`
#'   whose partial conjunction p-value
#'   `max_{r' <= r} min(1, (S - r' + 1) p_(r'))` is at most `alpha`;
#'   `r_hat_simultaneous` uses `alpha / K`. It is never smaller than
#'   `r_count` and is valid under arbitrary dependence when the null
#'   hypotheses are fixed in advance: for example when one frozen signature
#'   is tested at independent centres, or for the number of centres at which
#'   block `k` is not independent of the other blocks at all (which implies
#'   non-participation for every fitted rule). For the participation of rules
#'   learned at the other centres it can exceed its target with probability
#'   above `alpha`, because two centres can make each other's null hypotheses
#'   false.
#'
#' `r_count >= 1` and `r_hat >= 1` coincide: both are the any-centre test
#' `S min_s p_s <= alpha`. With a single external validation cohort (`S = 1`)
#' all bounds reduce to the test in that cohort.
#'
#' A centre whose held-out score of block `k` (or composite) is constant gets
#' `p = 1` for that block. Centres are processed in decreasing order of size.
#'
#' @param task ([mlr3::Task]) Data of all centres.
#' @param learner The model fitted at each set of training centres: an
#'   untrained [mlr3pipelines::GraphLearner] containing one [PipeOpMBsPLS]
#'   (for example from [mbspls_graph_learner()] with fixed sparsity budgets);
#'   an [mlr3tuning::AutoTuner] wrapping such a learner, which then tunes
#'   within each training set; or a function `function(task, row_ids)` that
#'   returns a trained GraphLearner.
#' @param site (`character(1)` | `vector`) Name of the column of the task's
#'   backend that identifies centres, or one centre label per task row in the
#'   order of `task$row_ids`.
#' @param component (`integer(1)`) Component to test. Default `1`.
#' @param n_perm (`integer(1)`) Permutations per centre. Default `9999`.
#' @param alpha (`numeric(1)`) Level of the lower confidence bounds. Default
#'   `0.05`.
#' @param strata (`NULL` | `character(1)` | `vector`) Optional strata, given
#'   like `site`. Permutations then exchange rows only within the same stratum
#'   of the held-out centre (e.g. diagnostic group), which tests participation
#'   beyond differences between strata.
#' @param seed (`NULL` | `integer(1)`) Optional non-negative seed. It seeds the
#'   training via [set.seed()] and the permutations of each centre through an
#'   independent stream ([mb_rng_streams()]), so that the permutations do not
#'   depend on the random numbers the training consumes.
#' @param mbspls_id (`NULL` | `character(1)`) Id of the MB-sPLS node, if the
#'   graph contains more than one.
#'
#' @return An object of class `mbspls_participation`, a list with
#' \describe{
#'   \item{`centres`}{`data.table` with one row per centre and block: `site`,
#'     `n`, `block`, held-out correlation `r`, training correlation `train_r`
#'     of the block with its training composite, and one-sided p-value `p`.}
#'   \item{`blocks`}{`data.table` with one row per block: number of centres
#'     `S`, the median held-out correlation `median_r`, the number of centres
#'     with `p <= alpha` (`n_significant`), the partial conjunction bounds
#'     `r_hat` and `r_hat_simultaneous`, and the counts `r_count` (centres with
#'     `p <= alpha / S`) and `r_count_simultaneous` (`p <= alpha / (K S)`).}
#'   \item{`partial_conjunction`}{Matrix of partial conjunction p-values,
#'     `r = 1, ..., S` (rows) by block (columns).}
#'   \item{`alpha`, `n_perm`, `component`, `strata`}{The settings used.}
#' }
#'
#' @references
#' Benjamini, Y. and Heller, R. (2008). Screening for partial conjunction
#' hypotheses. *Biometrics*, 64(4), 1215-1222.
#'
#' @seealso [mbspls_partial_conjunction()], [as_mbspls_participation()],
#'   [autoplot.mbspls_participation()], [impute_knn_graph()],
#'   [mbspls_eval_new_data()]
#'
#' @examples
#' \donttest{
#' library(mlr3)
#' library(mlr3pipelines)
#' set.seed(1)
#' n = 120
#' site = rep(c("a", "b", "c", "d"), each = n / 4)
#' f = rnorm(n)
#' df = data.frame(
#'   x1 = f + rnorm(n), x2 = f + rnorm(n),
#'   y1 = f + rnorm(n), y2 = rnorm(n),
#'   u1 = rnorm(n), u2 = rnorm(n),
#'   site = site
#' )
#' task = mlr3cluster::as_task_clust(df, id = "toy")
#' task$col_roles$feature = setdiff(task$col_roles$feature, "site")
#' blocks = list(X = c("x1", "x2"), Y = c("y1", "y2"), U = c("u1", "u2"))
#' gl = mbspls_graph_learner(
#'   blocks = blocks, ncomp = 1L, c_matrix = matrix(sqrt(2), 3, 1),
#'   impute_within_blocks = TRUE, bootstrap = FALSE,
#'   bootstrap_selection = FALSE
#' )
#' res = mbspls_participation(task, gl, site = "site", n_perm = 999L, seed = 1L)
#' res
#' }
#' @export
mbspls_participation = function(
  task,
  learner,
  site,
  component = 1L,
  n_perm = 9999L,
  alpha = 0.05,
  strata = NULL,
  seed = NULL,
  mbspls_id = NULL
) {
  checkmate::assert_class(task, "Task")
  component = checkmate::assert_int(component, lower = 1L, coerce = TRUE)
  n_perm = checkmate::assert_int(n_perm, lower = 1L, coerce = TRUE)
  checkmate::assert_number(alpha, lower = 0, upper = 1)
  if (!is.null(seed)) {
    seed = checkmate::assert_int(seed, lower = 0L, coerce = TRUE)
  }
  if (!is.null(mbspls_id)) checkmate::assert_string(mbspls_id)
  fit = .mbspls_participation_fitter(learner, mbspls_id)

  row_ids = task$row_ids
  site = .mbspls_participation_labels(task, site, "site")
  if (anyNA(site)) {
    stop("`site` must not contain missing values.", call. = FALSE)
  }
  if (!is.null(strata)) {
    strata = .mbspls_participation_labels(task, strata, "strata")
    if (anyNA(strata)) {
      stop("`strata` must not contain missing values.", call. = FALSE)
    }
  }
  sizes = table(site)
  sites = names(sizes)[order(-as.integer(sizes), names(sizes))]
  S = length(sites)
  if (S < 2L) {
    stop("At least two centres are required.", call. = FALSE)
  }

  if (!is.null(seed)) set.seed(seed)
  stream_seed = if (is.null(seed)) sample.int(.Machine$integer.max, 1L) else seed
  streams = mb_rng_streams(S, stream_seed)

  per_centre = vector("list", S)
  for (i in seq_len(S)) {
    held_out = site == sites[i]
    model = fit(task, row_ids[!held_out])
    t_train = .mbspls_participation_scores(model, task, row_ids[!held_out], component, mbspls_id)
    t_test = .mbspls_participation_scores(model, task, row_ids[held_out], component, mbspls_id)
    st = if (is.null(strata)) NULL else strata[held_out]
    tests = with_rng_stream_local(streams[[i]], function() {
      .mbspls_participation_centre(t_train, t_test, n_perm, st)
    })
    tests[, `:=`(site = sites[i], n = sum(held_out))]
    per_centre[[i]] = tests
  }

  centres = data.table::rbindlist(per_centre, use.names = TRUE)
  as_mbspls_participation(
    centres,
    alpha = alpha,
    component = component,
    n_perm = n_perm,
    strata = !is.null(strata)
  )
}


#' Participation summary from centre-level results
#'
#' @description
#' Builds an `mbspls_participation` object from held-out results computed
#' elsewhere, for example by a leave-one-centre-out pipeline that does not fit
#' into [mbspls_participation()], or for one external validation cohort. The
#' counts and bounds are computed as in [mbspls_participation()], with
#' [mbspls_partial_conjunction()], and the `print()` and
#' [autoplot()][autoplot.mbspls_participation] methods apply. `S` is the number
#' of centres in `centres` and `K` the number of blocks.
#'
#' @param centres (`data.frame`) One row per centre and block, with columns
#'   `site`, `block`, `n` (centre size), `r` (held-out correlation) and `p`
#'   (one-sided p-value), and optionally `train_r` (the training correlation of
#'   the block with its training composite). A block without a row at some
#'   centre is recorded there with `p = 1`.
#' @param alpha (`numeric(1)`) Level of the lower confidence bounds. Default
#'   `0.05`.
#' @param component (`integer(1)`) Component the results refer to. Default
#'   `1`.
#' @param n_perm (`integer(1)`) Permutations per centre, if known. Only
#'   printed.
#' @param strata (`logical(1)`) Were the permutations restricted to strata?
#'   Only printed. Default `FALSE`.
#'
#' @return An object of class `mbspls_participation`, as described in
#'   [mbspls_participation()].
#'
#' @seealso [mbspls_participation()], [autoplot.mbspls_participation()]
#' @examples
#' centres = data.frame(
#'   site = rep(c("a", "b", "c"), each = 2),
#'   block = rep(c("X", "Y"), 3),
#'   n = rep(c(50, 40, 30), each = 2),
#'   r = c(0.42, 0.05, 0.38, -0.02, 0.29, 0.11),
#'   p = c(0.0001, 0.36, 0.0008, 0.55, 0.06, 0.28)
#' )
#' as_mbspls_participation(centres)
#' @export
as_mbspls_participation = function(
  centres,
  alpha = 0.05,
  component = 1L,
  n_perm = NA_integer_,
  strata = FALSE
) {
  checkmate::assert_data_frame(centres, min.rows = 1L)
  missing_cols = setdiff(c("site", "block", "n", "r", "p"), names(centres))
  if (length(missing_cols)) {
    stop(sprintf(
      "`centres` lacks the column(s) %s.",
      paste0("'", missing_cols, "'", collapse = ", ")
    ), call. = FALSE)
  }
  checkmate::assert_number(alpha, lower = 0, upper = 1)
  component = checkmate::assert_int(component, lower = 1L, coerce = TRUE)
  n_perm = checkmate::assert_int(n_perm, lower = 1L, na.ok = TRUE, coerce = TRUE)
  checkmate::assert_flag(strata)

  centres = data.table::data.table(
    site = as.character(centres$site),
    n = as.integer(centres$n),
    block = as.character(centres$block),
    r = as.numeric(centres$r),
    train_r = if ("train_r" %in% names(centres)) as.numeric(centres$train_r) else NA_real_,
    p = as.numeric(centres$p)
  )
  if (anyNA(centres$site) || anyNA(centres$block)) {
    stop("`site` and `block` must not contain missing values.", call. = FALSE)
  }
  if (anyDuplicated(centres, by = c("site", "block"))) {
    stop("`centres` has more than one row for some centre and block.", call. = FALSE)
  }
  if (any(centres$p < 0 | centres$p > 1, na.rm = TRUE)) {
    stop("`p` must contain p-values in [0, 1].", call. = FALSE)
  }

  sites = unique(centres$site)
  block_names = unique(centres$block)
  # A block missing at a centre (for example, no usable columns in that
  # training set) is recorded with p = 1.
  site_n = centres[!is.na(n), list(n_site = max(n)), by = "site"]
  grid = data.table::CJ(site = sites, block = block_names, sorted = FALSE)
  centres = centres[grid, on = c("site", "block")]
  centres[is.na(n), n := site_n$n_site[match(site, site_n$site)]]
  centres[is.na(p), p := 1]

  S = length(sites)
  p_mat = matrix(NA_real_, S, length(block_names), dimnames = list(sites, block_names))
  p_mat[cbind(match(centres$site, sites), match(centres$block, block_names))] = centres$p
  pc = mbspls_partial_conjunction(p_mat, alpha = alpha)
  pc_sim = mbspls_partial_conjunction(p_mat, alpha = alpha, simultaneous = TRUE)

  blocks = centres[, list(
    S = .N,
    median_r = stats::median(r, na.rm = TRUE),
    n_significant = sum(p <= alpha)
  ), by = "block"]
  blocks[, `:=`(
    r_hat = unname(pc$r_hat[block]),
    r_hat_simultaneous = unname(pc_sim$r_hat[block]),
    r_count = unname(pc$r_count[block]),
    r_count_simultaneous = unname(pc_sim$r_count[block])
  )]

  structure(list(
    centres = centres[],
    blocks = blocks[],
    partial_conjunction = pc$p_pc,
    alpha = alpha,
    n_perm = n_perm,
    component = component,
    strata = strata
  ), class = "mbspls_participation")
}


#' Partial conjunction bounds from centre-level p-values
#'
#' @description
#' Combines the p-values of one or several blocks across `S` centres into
#' Bonferroni partial conjunction p-values and two lower confidence bounds for
#' the number of centres at which each null hypothesis is false: the partial
#' conjunction bound `r_hat` and the fixed-threshold count `r_count`.
#'
#' @details
#' For sorted p-values `p_(1) <= ... <= p_(S)`, the partial conjunction
#' p-value for "at least `r` of the `S` hypotheses are false" is
#' `min(1, (S - r + 1) p_(r))`, made monotone in `r` by a running maximum. It
#' is valid under arbitrary dependence between the p-values (Benjamini and
#' Heller, 2008), as needed for leave-one-centre-out p-values, whose training
#' sets overlap. `r_hat` is the largest `r` whose partial conjunction p-value
#' is at most the level, and `0` if there is none. With `simultaneous = TRUE`
#' the level is `alpha / K` for `K` blocks (columns), so that all bounds hold
#' jointly with probability at least `1 - alpha`.
#'
#' The validity of `r_hat` requires the set of true null hypotheses to be
#' fixed. In leave-one-centre-out validation the hypothesis tested at one
#' centre concerns weights learned at the other centres, so it depends on their
#' data, and `r_hat` can then exceed the number of false nulls with probability
#' above the level: with two centres whose data each decide whether the
#' other's null holds, the probability can reach `2 alpha (1 - alpha)`.
#' `r_count`, the number of p-values at most `level / S` (equivalently the
#' largest `r` with `S p_(r) <= level`), remains valid there: a false count
#' requires some true-null p-value at most the fixed threshold `level / S`,
#' and each p-value is valid conditionally on the data that decide whether its
#' null hypothesis holds. `r_count` never exceeds `r_hat`, and both agree on
#' whether at least one hypothesis is false.
#'
#' @param p (`numeric()` | `matrix`) P-values in `[0, 1]`: a vector (one per
#'   centre) or a matrix with centres in rows and blocks in columns.
#' @param alpha (`numeric(1)`) Level. Default `0.05`.
#' @param simultaneous (`logical(1)`) Adjust the level for the number of
#'   columns? Default `FALSE`.
#'
#' @return A list with the partial conjunction p-values `p_pc` (matrix,
#'   `r = 1, ..., S` by block), the bounds `r_hat` (named integer vector), the
#'   `level` used, and the counts `r_count` (named integer vector) of p-values
#'   at most `level / S`.
#'
#' @references
#' Benjamini, Y. and Heller, R. (2008). Screening for partial conjunction
#' hypotheses. *Biometrics*, 64(4), 1215-1222.
#'
#' @seealso [mbspls_participation()]
#' @examples
#' p = cbind(A = c(0.001, 0.004, 0.02, 0.3), B = c(0.2, 0.5, 0.04, 0.9))
#' mbspls_partial_conjunction(p)$r_hat
#' mbspls_partial_conjunction(p)$r_count
#' mbspls_partial_conjunction(p, simultaneous = TRUE)$r_hat
#' @export
mbspls_partial_conjunction = function(p, alpha = 0.05, simultaneous = FALSE) {
  if (is.null(dim(p))) {
    p = matrix(p, ncol = 1L, dimnames = list(names(p), "p"))
  }
  checkmate::assert_matrix(p, mode = "numeric", any.missing = FALSE, min.rows = 1L, min.cols = 1L)
  if (any(p < 0 | p > 1)) {
    stop("`p` must contain p-values in [0, 1].", call. = FALSE)
  }
  checkmate::assert_number(alpha, lower = 0, upper = 1)
  checkmate::assert_flag(simultaneous)
  if (is.null(colnames(p))) colnames(p) = paste0("block_", seq_len(ncol(p)))

  S = nrow(p)
  level = if (simultaneous) alpha / ncol(p) else alpha
  p_pc = vapply(seq_len(ncol(p)), function(j) {
    cummax(pmin(1, (S - seq_len(S) + 1) * sort(p[, j])))
  }, numeric(S))
  p_pc = matrix(p_pc, nrow = S, dimnames = list(paste0("r", seq_len(S)), colnames(p)))
  r_hat = vapply(seq_len(ncol(p)), function(j) sum(p_pc[, j] <= level), integer(1L))
  names(r_hat) = colnames(p)
  # Fixed threshold level / S: valid also when the hypotheses depend on the data.
  r_count = vapply(seq_len(ncol(p)), function(j) sum(p[, j] <= level / S), integer(1L))
  names(r_count) = colnames(p)
  list(p_pc = p_pc, r_hat = r_hat, level = level, r_count = r_count)
}


#' @export
print.mbspls_participation = function(x, ...) {
  cat(sprintf(
    "Held-out block participation, component %d: %d centres%s%s.\n",
    x$component, length(unique(x$centres$site)),
    if (is.na(x$n_perm)) "" else sprintf(", %d permutations per centre", x$n_perm),
    if (isTRUE(x$strata)) " (within strata)" else ""
  ))
  S = length(unique(x$centres$site))
  K = length(unique(x$centres$block))
  cat(sprintf(
    paste0(
      "r_count: centres with p <= %s (r_count_simultaneous: p <= %s), a lower %g%% confidence bound for the\n",
      "  number of centres at which the block participates in the rule fitted without that centre.\n",
      "r_hat: partial conjunction bound, valid when the null hypotheses do not depend on the training data\n",
      "  (a frozen signature tested at independent centres, or association of the raw blocks).\n"
    ),
    format(signif(x$alpha / S, 3L)), format(signif(x$alpha / (K * S), 3L)), 100 * (1 - x$alpha)
  ))
  print(x$blocks)
  invisible(x)
}


# Returns function(task, row_ids) -> trained GraphLearner.
.mbspls_participation_fitter = function(learner, mbspls_id) {
  if (is.function(learner)) {
    return(function(task, row_ids) {
      gl = learner(task, row_ids)
      if (!inherits(gl, "GraphLearner") || is.null(gl$model)) {
        stop("`learner(task, row_ids)` must return a trained GraphLearner.", call. = FALSE)
      }
      gl
    })
  }
  if (inherits(learner, "AutoTuner")) {
    inner = learner$learner
    if (!inherits(inner, "GraphLearner")) {
      stop("An AutoTuner `learner` must wrap a GraphLearner.", call. = FALSE)
    }
    .mbspls_pipeop_id(inner$graph, mbspls_id = mbspls_id, where = "learner$learner$graph")
    return(function(task, row_ids) {
      at = learner$clone(deep = TRUE)
      at$train(task, row_ids = row_ids)
      at$learner
    })
  }
  if (!inherits(learner, "GraphLearner")) {
    stop("`learner` must be a GraphLearner, an AutoTuner wrapping one, or a function.", call. = FALSE)
  }
  .mbspls_pipeop_id(learner$graph, mbspls_id = mbspls_id, where = "learner$graph")
  function(task, row_ids) {
    gl = learner$clone(deep = TRUE)
    gl$train(task, row_ids = row_ids)
    gl
  }
}

# One label per task row, from a backend column or a vector.
.mbspls_participation_labels = function(task, x, what) {
  if (is.character(x) && length(x) == 1L) {
    if (!x %in% task$backend$colnames) {
      stop(sprintf("Column '%s' (`%s`) is not in the task's backend.", x, what), call. = FALSE)
    }
    key = task$backend$primary_key
    d = task$backend$data(rows = task$row_ids, cols = c(key, x))
    return(as.character(d[[x]][match(task$row_ids, d[[key]])]))
  }
  if (length(x) != task$nrow) {
    stop(sprintf("`%s` must be a column name or have one value per task row.", what), call. = FALSE)
  }
  as.character(x)
}

# Block scores of one component for the given rows, with the training
# preprocessing and weights of the trained learner.
.mbspls_participation_scores = function(gl, task, rows, component, mbspls_id) {
  ev = mbspls_eval_new_data(gl, task$clone(deep = TRUE)$filter(rows), mbspls_id = mbspls_id)
  cols = paste0("LV", component, "_", ev$blocks)
  if (!all(cols %in% colnames(ev$T_mat))) {
    stop(sprintf("Component %d was not fitted in a training set.", component), call. = FALSE)
  }
  scores = ev$T_mat[, cols, drop = FALSE]
  colnames(scores) = ev$blocks
  scores
}

.mbspls_participation_std = function(v) {
  s = stats::sd(v)
  if (!is.finite(s) || s < 1e-12) {
    return(NULL)
  }
  (v - mean(v)) / s
}

# Tests at one held-out centre for every block.
.mbspls_participation_centre = function(t_train, t_test, n_perm, strata = NULL) {
  blocks = colnames(t_test)
  n = nrow(t_test)
  K = length(blocks)
  sg = suppressWarnings(sign(stats::cor(t_train)))
  sg[!is.finite(sg)] = 0
  z = lapply(seq_len(K), function(k) .mbspls_participation_std(t_test[, k]))
  z_train = scale(t_train)
  z_train[!is.finite(z_train)] = 0

  perms = if (n >= 3L) {
    if (is.null(strata)) {
      vapply(seq_len(n_perm), function(i) sample.int(n), integer(n))
    } else {
      idx = split(seq_len(n), strata)
      vapply(seq_len(n_perm), function(i) {
        o = seq_len(n)
        for (g in idx) o[g] = g[sample.int(length(g))]
        o
      }, integer(n))
    }
  }

  out = lapply(seq_len(K), function(k) {
    others = setdiff(seq_len(K), k)
    train_r = if (length(others)) {
      suppressWarnings(stats::cor(t_train[, k], drop(z_train[, others, drop = FALSE] %*% sg[k, others])))
    } else {
      NA_real_
    }
    res = data.table::data.table(block = blocks[k], r = NA_real_, train_r = as.numeric(train_r), p = 1)
    if (n < 3L || is.null(z[[k]]) || !length(others)) {
      return(res)
    }
    zo = vapply(others, function(j) if (is.null(z[[j]])) numeric(n) else z[[j]], numeric(n))
    comp = .mbspls_participation_std(drop(matrix(zo, n) %*% sg[k, others]))
    if (is.null(comp)) {
      return(res)
    }
    zk = z[[k]]
    observed = sum(zk * comp) / (n - 1)
    null = drop(crossprod(matrix(zk[perms], n, n_perm), comp)) / (n - 1)
    res[, `:=`(r = observed, p = mb_permutation_pvalue(observed, null, alternative = "greater"))]
    res
  })
  data.table::rbindlist(out)
}
