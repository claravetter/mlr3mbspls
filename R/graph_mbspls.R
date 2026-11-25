#' Create a k-NN Imputation Graph
#' @description k-NN imputation for numeric and factor columns with the Gower
#'   k-NN learners [LearnerRegrKNNGower] and [LearnerClassifKNNGower]. By
#'   default each feature is imputed from all other features. With `blocks`,
#'   each block is imputed from its own features only.
#' @details
#' With `blocks`, the imputed values of a block depend only on that block's
#' observed values (and on the training rows). Imputing a block from the other
#' blocks would create association between the block and the others, both in
#' the training data and in new data, wherever values are missing. Use
#' block-wise imputation whenever the association between blocks is the object
#' of inference, as in [mbspls_participation()].
#'
#' In the block-wise graph the learners use `no_neighbors = "marginal"`: a row
#' with too few observed features in the block to have an eligible neighbour
#' (see `min_feature_frac`), for example a participant missing the whole
#' block, receives the training mean or the most frequent training level.
#' Features outside all blocks are imputed from the other features outside all
#' blocks. Blocks are declared on the original feature names and resolved as in
#' [PipeOpMBsPLS], so factor columns replaced by an encoder still belong to
#' their block.
#' @param k Integer, number of neighbors (default 5).
#' @param blocks Optional named list of character vectors (feature names per
#'   block). `NULL` (default) imputes every feature from all other features.
#' @return [`Graph`]
#' @import mlr3 mlr3pipelines checkmate
#' @export
impute_knn_graph = function(k = 5, blocks = NULL) {
  checkmate::assert_integerish(k, lower = 1, len = 1)
  num_types = c("numeric", "integer")
  if (is.null(blocks)) {
    imp_num = po("imputelearner",
      learner = lrn("regr.knngower", k = k),
      affect_columns = selector_union(selector_type("numeric"), selector_type("integer"))
    )
    imp_num$id = "impute_num_knn"
    imp_fac = po("imputelearner",
      learner = lrn("classif.knngower", k = k),
      affect_columns = selector_type("factor")
    )
    imp_fac$id = "impute_fac_knn"
    return(imp_num %>>% imp_fac)
  }

  blocks = mb_normalize_blocks(blocks, .var.name = "blocks")
  # Columns of one group (a block, or all features outside the blocks),
  # resolved on the task the PipeOp receives.
  group_columns = function(group) {
    force(group)
    function(task) {
      resolved = mb_resolve_block_columns(task$feature_names, blocks)
      if (is.null(group)) {
        return(setdiff(task$feature_names, unlist(resolved, use.names = FALSE)))
      }
      resolved[[group]]
    }
  }
  typed_columns = function(columns, types) {
    force(columns)
    force(types)
    function(task) {
      cols = columns(task)
      ft = task$feature_types
      cols[ft$type[match(cols, ft$id)] %in% types]
    }
  }
  groups = c(lapply(names(blocks), identity), list(NULL))
  ops = list()
  for (group in groups) {
    columns = group_columns(group)
    suffix = if (is.null(group)) "unblocked" else paste0("block_", group)
    ops[[length(ops) + 1L]] = po("imputelearner",
      id = paste0("impute_num_knn_", suffix),
      learner = lrn("regr.knngower", k = k, no_neighbors = "marginal"),
      affect_columns = typed_columns(columns, num_types),
      context_columns = columns
    )
    ops[[length(ops) + 1L]] = po("imputelearner",
      id = paste0("impute_fac_knn_", suffix),
      learner = lrn("classif.knngower", k = k, no_neighbors = "marginal"),
      affect_columns = typed_columns(columns, "factor"),
      context_columns = columns
    )
  }
  Reduce(`%>>%`, ops)
}


#' Default MB-sPLS Preprocessing Graph
#' @description
#' Site correction -> encoding/imputation -> scaling. No MB-sPLS here.
#' @param blocks Named list of character vectors per block. If `NULL`, the
#'   mapping is taken from `task$blocks`.
#' @param task Optional [mlr3::Task] carrying multi-block metadata via
#'   [TaskMultiBlock()]. Used when `blocks = NULL`.
#' @param site_correction Named list of features used for site correction.
#'   Defaults to `list()` (no correction). When `task` is supplied, the columns
#'   read at prediction time (all `"partial_corr"` and `"dir"` columns and the
#'   ComBat `site`) must exist and must not be target columns; ComBat
#'   `covariates` may name the target (see [PipeOpSiteCorrection]).
#' @param site_correction_methods Named list of methods for site correction.
#'   Defaults to `list()`.
#' @param keep_site_col Keep site column after correction?
#' @param k k for kNN imputation (default 5).
#' @param id_suffix Optional suffix for PipeOp ids.
#' @param impute_within_blocks Logical; impute each block from its own
#'   features only (see [impute_knn_graph()]). Default `FALSE` imputes every
#'   feature from all other features. Use `TRUE` for [mbspls_participation()].
#'   Last in the signature, so that positional calls written for earlier
#'   versions keep their meaning.
#' @return [`Graph`]
#' @import mlr3 mlr3pipelines checkmate
#' @export
mbspls_preproc_graph = function(
  blocks = NULL,
  task = NULL,
  site_correction = list(),
  site_correction_methods = list(),
  keep_site_col = FALSE,
  k = 5,
  id_suffix = NULL,
  impute_within_blocks = FALSE
) {
  blocks = mb_graph_blocks(blocks = blocks, task = task, context = "mbspls_preproc_graph")
  assert_list(blocks, types = "character", names = "unique")
  assert_list(site_correction, types = c("character", "list"), names = "unique")
  assert_list(site_correction_methods, types = "character", names = "unique")
  mb_validate_site_correction(
    task = task,
    site_correction = site_correction,
    context = "mbspls_preproc_graph",
    methods = site_correction_methods
  )

  id_with_suffix = function(base) {
    if (is.null(id_suffix)) base else paste0(base, "_", id_suffix)
  }

  site_cols = unique(unlist(site_correction, recursive = TRUE, use.names = FALSE))
  affect_columns = if (length(site_cols)) {
    selector_invert(selector_name(site_cols))
  } else {
    selector_all()
  }

  checkmate::assert_flag(impute_within_blocks)
  ppl_convert_types = ppl("convert_types", "character", "factor")
  ppl_impute = ppl("imputeknn", k = k, blocks = if (impute_within_blocks) blocks else NULL)

  if (!is.null(id_suffix)) {
    ppl_convert_types$update_ids(postfix = paste0("_", id_suffix))
    ppl_impute$update_ids(postfix = paste0("_", id_suffix))
  }

  graph = ppl_convert_types %>>%
    po("encode",
      id = id_with_suffix("encode"),
      method = "treatment",
      affect_columns = affect_columns
    ) %>>%
    ppl_impute %>>%
    po("sitecorr",
      id = id_with_suffix("sitecorr"),
      blocks = blocks,
      site_correction = site_correction,
      method = site_correction_methods,
      keep_site_col = keep_site_col
    ) %>>%
    po("scale",
      id = id_with_suffix("scale")
    )
  graph
}


#' MB-sPLS Graph: preproc -> MB-sPLS -> bootstrap-select (two selection methods)
#'
#' @param blocks,task,site_correction,site_correction_methods,keep_site_col Preproc settings.
#' @param ncomp Number of MB-sPLS components.
#' @param k k-NN for imputation (preproc).
#' @param performance_metric "mac" or "frobenius".
#' @param correlation_method "pearson" or "spearman".
#' @param c_matrix Optional L1 constraint matrix for MB-sPLS.
#'
#' @param permutation_test,n_perm,perm_alpha Conditional component-wise
#'   train-time permutation diagnostic (MB-sPLS), not a full-pipeline test.
#' @param predict_weights character; one of "auto","raw","stable_ci","stable_frequency".
#'   Weights that [PipeOpMBsPLS] evaluates in its prediction-side payload and
#'   `val_test` diagnostics. The emitted LV features always use the training
#'   weights; with bootstrap selection active (and `stability_only = FALSE`),
#'   [PipeOpMBsPLSBootstrapSelect] replaces them by stable-weight LVs at train
#'   and predict time. With `stability_only = TRUE`, `"auto"` evaluates the raw
#'   weights and `"stable_ci"`/`"stable_frequency"` are rejected, because no
#'   feature of that graph is derived from the stable weights.
#' @param val_test,val_test_alpha,val_test_n,val_test_permute_all Prediction-side
#'   conditional permutation diagnostic or descriptive bootstrap uncertainty.
#' @param seed_validation Optional seed for prediction-side diagnostics; one
#'   deterministic stream is assigned per component.
#'
#' @param bootstrap Logical; run bootstrap selection (default TRUE).
#' @param store_train_blocks Logical; retain the fitted (training-centred)
#'   block matrices for bootstrap selection, post-fit summaries and plots. The
#'   default follows `bootstrap`. It must be `TRUE` when bootstrap selection is
#'   enabled.
#' @param bootstrap_selection Logical; whether to run bootstrap-based feature
#'   selection inside [PipeOpMBsPLSBootstrapSelect] (default TRUE).
#' @param stability_only Logical; only compute stability, no selection (default FALSE).
#' @param B Integer; bootstrap replicates (default 500).
#' @param alpha Numeric; CI alpha (default 0.05).
#' @param align Bootstrap sign alignment passed to [PipeOpMBsPLSBootstrapSelect]:
#'   `"block_sign"` (default) or `"score_correlation"`.
#' @param selection_method "ci" (default) or "frequency".
#' @param frequency_threshold Numeric in `[0,1]`; only if selection_method="frequency" (default 0.6).
#' @param stable_weight_source "training" (default) or "bootstrap_mean".
#' @param stratify_by_block Optional dummy block for stratified bootstrap (e.g., "Studygroup").
#' @param bootstrap_groups Optional exchangeability-group vector for cluster
#'   bootstrap sampling. Named vectors are aligned to task row IDs.
#' @param workers Integer; number of cross-platform bootstrap workers.
#' @param seed_train Optional seed for random draws during MB-sPLS training,
#'   i.e. the permutations of `permutation_test`. The fit itself is
#'   deterministic and does not depend on the seed.
#' @param seed_bootstrap Optional seed for the bootstrap replicates of
#'   [PipeOpMBsPLSBootstrapSelect]; when set, each replicate gets its own
#'   deterministic RNG stream. The default `NULL` matches the PipeOp and uses
#'   the ambient RNG, so results are reproducible with [set.seed()].
#' @param id_suffix Optional suffix for PipeOp ids.
#' @param log_env Shared environment (created if NULL).
#' @param impute_within_blocks Logical; impute each block from its own
#'   features only (see [impute_knn_graph()]). Default `FALSE`. Last in the
#'   signature, so that positional calls written for earlier versions keep
#'   their meaning.
#'
#' @return [mlr3pipelines::Graph]
#' @import mlr3 mlr3pipelines checkmate
#' @export
mbspls_graph = function(
  blocks = NULL,
  task = NULL,
  site_correction = list(),
  site_correction_methods = list(),
  keep_site_col = FALSE,
  ncomp,
  k = 5,
  performance_metric = c("mac", "frobenius"),
  correlation_method = c("pearson", "spearman"),
  c_matrix = NULL,

  permutation_test = FALSE,
  n_perm = 500L,
  perm_alpha = 0.05,
  predict_weights = c("auto", "raw", "stable_ci", "stable_frequency"),
  val_test = c("none", "permutation", "bootstrap"),
  val_test_alpha = 0.05,
  val_test_n = 1000L,
  val_test_permute_all = TRUE,
  seed_validation = NULL,

  bootstrap = TRUE,
  store_train_blocks = bootstrap,
  stability_only = FALSE,
  B = 500L,
  alpha = 0.05,
  align = c("block_sign", "score_correlation"),
  bootstrap_selection = TRUE,
  selection_method = c("ci", "frequency"),
  frequency_threshold = 0.60,
  stable_weight_source = c("training", "bootstrap_mean"),
  stratify_by_block = NULL,
  bootstrap_groups = NULL,
  seed_train = NULL,
  seed_bootstrap = NULL,
  workers = 1L,
  id_suffix = NULL,
  log_env = NULL,
  impute_within_blocks = FALSE
) {
  blocks = mb_graph_blocks(blocks = blocks, task = task, context = "mbspls_graph")
  checkmate::assert_list(blocks, types = "character", names = "unique")
  checkmate::assert_list(site_correction, types = c("character", "list"), names = "unique")
  checkmate::assert_list(site_correction_methods, types = "character", names = "unique")
  checkmate::assert_int(ncomp, lower = 1)
  performance_metric = match.arg(performance_metric)
  correlation_method = match.arg(correlation_method)
  predict_weights = match.arg(predict_weights)
  if (isTRUE(stability_only) && predict_weights %in% c("stable_ci", "stable_frequency")) {
    stop(
      sprintf(
        paste0(
          "`predict_weights = \"%s\"` cannot be combined with `stability_only = TRUE`: ",
          "that graph passes the LV features of the raw training weights to the learner, ",
          "so no feature is derived from the stable weights. Use \"raw\" or \"auto\", ",
          "or set `stability_only = FALSE`."
        ),
        predict_weights
      ),
      call. = FALSE
    )
  }
  if (isTRUE(stability_only) && identical(predict_weights, "auto")) {
    predict_weights = "raw"
  }
  val_test = match.arg(val_test)
  align = match.arg(align)
  selection_method = match.arg(selection_method)
  stable_weight_source = match.arg(stable_weight_source)
  checkmate::assert_flag(store_train_blocks)

  if (is.null(log_env)) {
    log_env = new.env(parent = emptyenv())
    # suppress overwrite warnings in typical resampling/tuning workflows
    log_env$warn_overwrite = FALSE
  }

  if (isTRUE(stability_only) && !isTRUE(bootstrap_selection)) {
    stop("`stability_only = TRUE` has no effect when `bootstrap_selection = FALSE` because there is no PipeOpMBsPLSBootstrapSelect in the graph. Set `bootstrap_selection = TRUE` or remove `stability_only`.", call. = FALSE)
  }
  if (isTRUE(bootstrap_selection) && isTRUE(bootstrap) &&
    !isTRUE(store_train_blocks)) {
    stop(
      "`store_train_blocks` must be TRUE when bootstrap selection is enabled.",
      call. = FALSE
    )
  }

  if (isTRUE(bootstrap_selection)) {
    po_bootstrap_select = po("mbspls_bootstrap_select",
      id = if (is.null(id_suffix)) "mbspls_bootstrap_select" else paste0("mbspls_bootstrap_select_", id_suffix),
      log_env = log_env,
      bootstrap = bootstrap,
      stability_only = stability_only,
      B = B,
      alpha = alpha,
      align = align,
      selection_method = selection_method,
      frequency_threshold = frequency_threshold,
      stable_weight_source = stable_weight_source,
      stratify_by_block = stratify_by_block,
      bootstrap_groups = bootstrap_groups,
      seed_bootstrap = seed_bootstrap,
      workers = workers
    )
  }

  graph = ppl("mbspls_preproc",
    blocks = blocks,
    task = task,
    site_correction = site_correction,
    site_correction_methods = site_correction_methods,
    keep_site_col = keep_site_col,
    k = k,
    impute_within_blocks = impute_within_blocks,
    id_suffix = id_suffix
  ) %>>%
    po("mbspls",
      id = if (is.null(id_suffix)) "mbspls" else paste0("mbspls_", id_suffix),
      blocks = blocks,
      ncomp = ncomp,
      performance_metric = performance_metric,
      correlation_method = correlation_method,
      c_matrix = c_matrix,

      # optional train-time permutation
      permutation_test = permutation_test,
      n_perm = n_perm,
      perm_alpha = perm_alpha,

      # which weights to use at predict/validation time?
      predict_weights = predict_weights,

      # optional prediction-side validation
      val_test = val_test,
      val_test_alpha = val_test_alpha,
      val_test_n = val_test_n,
      val_test_permute_all = val_test_permute_all,
      seed_validation = seed_validation,

      # expose training snapshot for selection
      store_train_blocks = store_train_blocks,
      append = isTRUE(bootstrap_selection) && isTRUE(bootstrap) && !isTRUE(stability_only),
      seed_train = seed_train,
      log_env = log_env
    )
  if (isTRUE(bootstrap_selection)) {
    graph = graph %>>% po_bootstrap_select
  }
  return(graph)
}




#' MB-sPLS GraphLearner: preproc -> MB-sPLS -> bootstrap-select (two selection methods) -> learner
#'
#' @param learner Downstream learner (default k-means with 1 center).
#' @inheritParams mbspls_graph
#'
#' @return [mlr3pipelines::GraphLearner]
#' @import mlr3 mlr3pipelines checkmate
#' @importFrom mlr3cluster LearnerClust
#' @export
mbspls_graph_learner = function(
  learner = lrn("clust.kmeans", centers = 1L),
  blocks = NULL,
  task = NULL,
  site_correction = list(),
  site_correction_methods = list(),
  keep_site_col = FALSE,
  ncomp,
  k = 5,
  performance_metric = c("mac", "frobenius"),
  correlation_method = c("pearson", "spearman"),
  c_matrix = NULL,

  permutation_test = FALSE,
  n_perm = 500L,
  perm_alpha = 0.05,
  predict_weights = c("auto", "raw", "stable_ci", "stable_frequency"),
  val_test = c("none", "permutation", "bootstrap"),
  val_test_alpha = 0.05,
  val_test_n = 1000L,
  val_test_permute_all = TRUE,
  seed_validation = NULL,

  bootstrap = TRUE,
  store_train_blocks = bootstrap,
  stability_only = FALSE,
  B = 500L,
  alpha = 0.05,
  align = c("block_sign", "score_correlation"),
  bootstrap_selection = TRUE,
  selection_method = c("ci", "frequency"),
  frequency_threshold = 0.60,
  stable_weight_source = c("training", "bootstrap_mean"),
  stratify_by_block = NULL,
  bootstrap_groups = NULL,
  seed_train = NULL,
  seed_bootstrap = NULL,
  workers = 1L,
  id_suffix = NULL,
  log_env = NULL,
  impute_within_blocks = FALSE
) {
  checkmate::assert_class(learner, "Learner")

  graph = mbspls_graph(
    blocks = blocks,
    task = task,
    site_correction = site_correction,
    site_correction_methods = site_correction_methods,
    keep_site_col = keep_site_col,
    ncomp = ncomp,
    k = k,
    impute_within_blocks = impute_within_blocks,
    performance_metric = performance_metric,
    correlation_method = correlation_method,
    c_matrix = c_matrix,
    permutation_test = permutation_test,
    n_perm = n_perm,
    perm_alpha = perm_alpha,
    predict_weights = predict_weights,
    val_test = val_test,
    val_test_alpha = val_test_alpha,
    val_test_n = val_test_n,
    val_test_permute_all = val_test_permute_all,
    seed_validation = seed_validation,
    bootstrap = bootstrap,
    store_train_blocks = store_train_blocks,
    stability_only = stability_only,
    B = B,
    alpha = alpha,
    align = align,
    bootstrap_selection = bootstrap_selection,
    selection_method = selection_method,
    frequency_threshold = frequency_threshold,
    stable_weight_source = stable_weight_source,
    stratify_by_block = stratify_by_block,
    bootstrap_groups = bootstrap_groups,
    seed_train = seed_train,
    seed_bootstrap = seed_bootstrap,
    workers = workers,
    id_suffix = id_suffix,
    log_env = log_env
  )
  as_learner(graph %>>%
    po("learner", learner))
}


#' Supervised MB-sPLS-XY graph: preproc -> MB-sPLS-XY
#'
#' @param blocks,task,site_correction,site_correction_methods,keep_site_col Preproc settings.
#' @param ncomp Number of MB-sPLS-XY components.
#' @param k k-NN for imputation (preproc).
#' @param performance_metric "mac" or "frobenius".
#' @param correlation_method "pearson" or "spearman".
#' @param c_matrix Optional L1 constraint matrix for MB-sPLS-XY.
#' @param permutation_test,n_perm,perm_alpha Conditional component-wise
#'   train-time permutation diagnostic for MB-sPLS-XY.
#' @param y_rep Integer replication count for the target block.
#' @param emit_y_scores Logical; store the training target-block scores in the
#'   state of [PipeOpMBsPLSXY] (`$state$scores_y`) for inspection. They are
#'   never added to the task features.
#' @param center_y,scale_y Logical target centering/scaling flags. The target
#'   block is always centred by its training means (`center_y = FALSE` is
#'   ignored with a warning); `scale_y` scales it to unit training standard
#'   deviation.
#' @param id_suffix Optional suffix for PipeOp ids.
#' @param log_env Shared environment (created if NULL).
#'
#' @return [mlr3pipelines::Graph]
#' @import mlr3 mlr3pipelines checkmate
#' @export
mbsplsxy_graph = function(
  blocks = NULL,
  task = NULL,
  site_correction = list(),
  site_correction_methods = list(),
  keep_site_col = FALSE,
  ncomp,
  k = 5,
  performance_metric = c("mac", "frobenius"),
  correlation_method = c("pearson", "spearman"),
  c_matrix = NULL,
  permutation_test = FALSE,
  n_perm = 500L,
  perm_alpha = 0.05,
  y_rep = 1L,
  emit_y_scores = FALSE,
  center_y = TRUE,
  scale_y = TRUE,
  id_suffix = NULL,
  log_env = NULL
) {
  blocks = mb_graph_blocks(blocks = blocks, task = task, context = "mbsplsxy_graph")
  checkmate::assert_list(blocks, types = "character", names = "unique")
  checkmate::assert_list(site_correction, types = c("character", "list"), names = "unique")
  checkmate::assert_list(site_correction_methods, types = "character", names = "unique")
  mb_validate_site_correction(
    task = task,
    site_correction = site_correction,
    context = "mbsplsxy_graph",
    methods = site_correction_methods
  )
  mb_validate_supervised_task(task = task, context = "mbsplsxy_graph")
  checkmate::assert_int(ncomp, lower = 1)
  performance_metric = match.arg(performance_metric)
  correlation_method = match.arg(correlation_method)

  log_env = if (is.null(log_env)) new.env(parent = emptyenv()) else log_env
  mbsplsxy_id = if (is.null(id_suffix)) "mbsplsxy" else paste0("mbsplsxy_", id_suffix)

  ppl("mbspls_preproc",
    blocks = blocks,
    task = task,
    site_correction = site_correction,
    site_correction_methods = site_correction_methods,
    keep_site_col = keep_site_col,
    k = k,
    id_suffix = id_suffix
  ) %>>%
    po("mbsplsxy",
      id = mbsplsxy_id,
      blocks = blocks,
      ncomp = ncomp,
      performance_metric = performance_metric,
      correlation_method = correlation_method,
      c_matrix = c_matrix,
      permutation_test = permutation_test,
      n_perm = n_perm,
      perm_alpha = perm_alpha,
      y_rep = y_rep,
      emit_y_scores = emit_y_scores,
      center_y = center_y,
      scale_y = scale_y,
      log_env = log_env
    )
}


#' MB-sPLS-XY GraphLearner: preproc -> MB-sPLS-XY -> learner
#'
#' @param learner Downstream learner. If `NULL`, a featureless learner matching
#'   the inferred supervised task type is used.
#' @param task_type One of `"classif"` or `"regr"`. Only used if `learner` and
#'   `task` do not jointly determine the target mode.
#' @inheritParams mbsplsxy_graph
#'
#' @return [mlr3pipelines::GraphLearner]
#' @import mlr3 mlr3pipelines checkmate
#' @export
mbsplsxy_graph_learner = function(
  learner = NULL,
  blocks = NULL,
  task = NULL,
  task_type = c("classif", "regr"),
  site_correction = list(),
  site_correction_methods = list(),
  keep_site_col = FALSE,
  ncomp,
  k = 5,
  performance_metric = c("mac", "frobenius"),
  correlation_method = c("pearson", "spearman"),
  c_matrix = NULL,
  permutation_test = FALSE,
  n_perm = 500L,
  perm_alpha = 0.05,
  y_rep = 1L,
  emit_y_scores = FALSE,
  center_y = TRUE,
  scale_y = TRUE,
  id_suffix = NULL,
  log_env = NULL
) {
  inferred_type = if (!is.null(task)) {
    mb_validate_supervised_task(task = task, context = "mbsplsxy_graph_learner")
  } else {
    match.arg(task_type)
  }

  if (is.null(learner)) {
    learner = if (identical(inferred_type, "classif")) {
      lrn("classif.featureless")
    } else {
      lrn("regr.featureless")
    }
  }

  checkmate::assert_class(learner, "Learner")
  mb_validate_supervised_learner(learner = learner, expected_type = inferred_type, context = "mbsplsxy_graph_learner")

  graph = mbsplsxy_graph(
    blocks = blocks,
    task = task,
    site_correction = site_correction,
    site_correction_methods = site_correction_methods,
    keep_site_col = keep_site_col,
    ncomp = ncomp,
    k = k,
    performance_metric = performance_metric,
    correlation_method = correlation_method,
    c_matrix = c_matrix,
    permutation_test = permutation_test,
    n_perm = n_perm,
    perm_alpha = perm_alpha,
    y_rep = y_rep,
    emit_y_scores = emit_y_scores,
    center_y = center_y,
    scale_y = scale_y,
    id_suffix = id_suffix,
    log_env = log_env
  )

  as_learner(graph %>>% po("learner", learner))
}
