#' k-Nearest Neighbours with Gower Distance - Regression
#'
#' @title Learner: `regr.knngower`
#' @name LearnerRegrKNNGower
#' @format R6 class inheriting from [mlr3::LearnerRegr]
#'
#' @description
#' A *k*-nearest neighbours regressor using **Gower distance**
#' (Rcpp/Armadillo backend). Supports numeric, unordered factors, and
#' ordered factors with the same encodings as the classification learner.
#' Predictions are (optionally inverse-distance) **weighted means** of the
#' neighbour responses; the (weighted) variance provides an approximate
#' `se` when requested. If a row has no eligible neighbour, prediction fails
#' with an error by default; with `no_neighbors = "marginal"` that row
#' receives the training mean (and the training variance for `se`).
#'
#' @section Parameters (in `param_set`):
#' \describe{
#'   \item{`k`}{`integer(1)`. Number of neighbours. Default: `5`.}
#'   \item{`weights`}{`character(1)`. `"uniform"` or `"inverse"`.
#'     Default: `"inverse"`.}
#'   \item{`min_feature_frac`}{`numeric(1)` in `[0,1]`. Minimum fraction of
#'     comparable features required for a neighbour. Default: `0.2`.}
#'   \item{`na_handling`}{`character(1)`. `"pairwise"` (skip missing values
#'     per feature) or `"fail"` to error if training or prediction features
#'     contain `NA`. Default: `"pairwise"`.}
#'   \item{`no_neighbors`}{`character(1)`. What to do with a prediction row
#'     that has no eligible neighbour (fewer than `min_feature_frac` of the
#'     features observed in common with every training row, e.g. a row with
#'     all features missing): `"error"` (default) stops; `"marginal"` predicts
#'     the training mean. Used by the block-wise [impute_knn_graph()].}
#' }
#'
#' @section Prediction:
#' Returns numeric `response`; if requested, `se` is the square root of the
#' (weighted) neighbour variance (non-unbiased).
#'
#' @param id `character(1)`. Learner identifier (default `"regr.knngower"`).
#' @param param_vals `list()`. Named list of initial hyperparameters.
#'
#' @return An object of class `LearnerRegrKNNGower`.
#'
#' @examples
#' if (requireNamespace("mlr3", quietly = TRUE)) {
#'   library(mlr3)
#'   set.seed(1)
#'   n = 100
#'   df = data.frame(x1 = runif(n), x2 = runif(n))
#'   df$y = sin(3 * df$x1) + df$x2 + rnorm(n, sd = 0.1)
#'   task = TaskRegr$new("toy", backend = df, target = "y")
#'   lrn = lrn("regr.knngower", k = 10, weights = "inverse")
#'   lrn$train(task)
#'   lrn$predict(task)$score(msr("regr.rmse"))
#' }
#'
#' @export
LearnerRegrKNNGower = R6::R6Class("LearnerRegrKNNGower",
  inherit = mlr3::LearnerRegr,
  public = list(
    #' @description Create a new `LearnerRegrKNNGower` instance.
    initialize = function() {
      ps = paradox::ps(
        k = paradox::p_int(lower = 1L, default = 5L, tags = c("train", "predict")),
        weights = paradox::p_fct(levels = c("uniform", "inverse"),
          default = "inverse", tags = "predict"),
        min_feature_frac = paradox::p_dbl(lower = 0, upper = 1,
          default = 0.2, tags = "predict"),
        na_handling = paradox::p_fct(levels = c("pairwise", "fail"),
          default = "pairwise", tags = c("train", "predict")),
        no_neighbors = paradox::p_fct(levels = c("error", "marginal"),
          default = "error", tags = "predict")
      )
      ps$set_values(k = 5L, weights = "inverse",
        min_feature_frac = 0.2, na_handling = "pairwise", no_neighbors = "error")

      super$initialize(
        id = "regr.knngower",
        feature_types = c("numeric", "integer", "factor", "ordered"),
        predict_types = c("response", "se"),
        param_set = ps,
        properties = c("missings"),
        packages = "mlr3mbspls"
      )
    }
  ),
  private = list(

    .train = function(task) {
      pv = self$param_set$get_values(tags = "train")
      y = task$truth()
      if (anyNA(y)) stop("Regression target contains missing values.")

      if (identical(pv$na_handling, "fail")) {
        x_all = task$data(cols = task$feature_names)
        if (anyNA(x_all)) stop("Training data contain missing features (na_handling = 'fail').")
      }

      df = task$data(cols = task$feature_names)
      if (pv$k > nrow(df)) {
        warning(sprintf(
          "regr.knngower: k (%d) exceeds the number of training observations (%d); all training points will be used as neighbors.",
          pv$k, nrow(df)))
      }

      types = task$feature_types
      num_cols = types[type %in% c("numeric", "integer"), id]
      cat_cols = types[type %in% c("factor"), id]
      ord_cols = types[type %in% c("ordered"), id]

      enc = knn_gower_encode_blocks(df, num_cols, cat_cols, ord_cols, ref = NULL)

      y_num = as.numeric(y)
      y_var = stats::var(y_num, na.rm = TRUE)
      if (!is.finite(y_var) || is.na(y_var)) {
        y_var = 0
      }

      list(
        num_cols = num_cols,
        cat_cols = cat_cols,
        ord_cols = ord_cols,
        Xnum = enc$Xnum,
        Xcat = enc$Xcat,
        Xord = enc$Xord,
        ranges_num = enc$ranges_num,
        cat_levels = enc$cat_levels,
        ord_levels = enc$ord_levels,
        y = y_num,
        y_mean = mean(y_num, na.rm = TRUE),
        y_var = y_var
      )
    },

    .predict = function(task) {
      st = self$model
      pv = self$param_set$get_values(tags = "predict")

      df = task$data(cols = task$feature_names)
      if (identical(pv$na_handling, "fail") && anyNA(df)) {
        stop("Prediction data contain missing features (na_handling = 'fail').")
      }
      enc_te = knn_gower_encode_blocks(
        df,
        num_cols = st$num_cols,
        cat_cols = st$cat_cols,
        ord_cols = st$ord_cols,
        ref = list(
          ranges_num = st$ranges_num,
          cat_levels = st$cat_levels,
          ord_levels = st$ord_levels
        )
      )

      # Rows without an eligible neighbour make the native kernel stop. With
      # no_neighbors = "marginal" they are excluded from the kernel call and
      # receive the training mean and variance instead.
      n_te = nrow(df)
      ok = if (identical(pv$no_neighbors, "marginal")) {
        knn_gower_has_neighbors(st, enc_te, pv$min_feature_frac)
      } else {
        rep(TRUE, n_te)
      }
      out = list(mean = rep(as.numeric(st$y_mean), n_te), var = rep(as.numeric(st$y_var), n_te))
      if (any(ok)) {
        te = if (all(ok)) enc_te else knn_gower_subset_rows(enc_te, ok)
        # NOTE: call positionally
        fit = knn_predict_regr_gower_cpp(
          st$Xnum, st$Xcat, st$Xord,
          te$Xnum, te$Xcat, te$Xord,
          st$ranges_num, st$y,
          as.integer(pv$k),
          pv$weights,
          as.numeric(pv$min_feature_frac),
          as.numeric(st$y_mean), # fallback_mean: passed for API compatibility; unused by C++ (it errors on no neighbors)
          as.numeric(st$y_var) # fallback_var:  passed for API compatibility; unused by C++ (it errors on no neighbors)
        )
        out$mean[ok] = as.numeric(fit$mean)
        out$var[ok] = as.numeric(fit$var)
      }

      response = as.numeric(out$mean)
      if ("se" %in% self$predict_type) {
        se = sqrt(pmax(0, as.numeric(out$var)))
        mlr3::PredictionRegr$new(task = task, response = response, se = se)
      } else {
        mlr3::PredictionRegr$new(task = task, response = response)
      }
    }
  )
)
