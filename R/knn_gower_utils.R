#' Encode features for the Gower k-NN learners
#'
#' Shared encoder of [LearnerClassifKNNGower] and [LearnerRegrKNNGower]. It
#' splits the features into the numeric, categorical and ordered blocks used by
#' the native Gower kernels. Without `ref`, the encoding (numeric ranges and
#' level sets) is fitted on `df`; with `ref`, the stored training encoding is
#' applied unchanged.
#'
#' - Numeric: training range (max - min) per column; non-positive or
#'   non-finite ranges are replaced by 1.
#' - Categorical: integer codes with `0 = NA`, `1..L` for known levels and
#'   `-1` for levels unseen during training (forces a mismatch).
#' - Ordered: codes mapped to `[0, 1]` as `(code - 1) / (L - 1)`; missing and
#'   unseen levels are `NA` and therefore skipped pairwise. A single-level
#'   ordered feature is encoded as `0`, keeping `NA` for missing values.
#'
#' @param df (`data.table`) Feature data.
#' @param num_cols,cat_cols,ord_cols (`character()`) Numeric, unordered factor
#'   and ordered factor columns.
#' @param ref (`list()` or `NULL`) Stored training encoding with elements
#'   `ranges_num`, `cat_levels` and `ord_levels`.
#' @return A list with the encoded matrices `Xnum`, `Xcat` and `Xord` and the
#'   encoding (`ranges_num`, `cat_levels`, `ord_levels`).
#' @noRd
knn_gower_encode_blocks = function(df, num_cols, cat_cols, ord_cols, ref = NULL) {
  n = nrow(df)

  # numeric
  if (length(num_cols)) {
    Xn = as.matrix(df[, num_cols, with = FALSE])
    storage.mode(Xn) = "double"
    if (is.null(ref)) {
      r_min = suppressWarnings(apply(Xn, 2, min, na.rm = TRUE))
      r_max = suppressWarnings(apply(Xn, 2, max, na.rm = TRUE))
      rng = r_max - r_min
      rng[!is.finite(rng) | rng <= 0] = 1.0
    } else {
      rng = ref$ranges_num
    }
  } else {
    Xn = matrix(numeric(0), nrow = n, ncol = 0)
    rng = numeric(0)
  }

  # categorical (unordered) -> integer codes (0=NA, 1..L known, -1 unseen in predict)
  if (length(cat_cols)) {
    if (is.null(ref)) {
      cat_levels = lapply(cat_cols, function(cn) levels(as.factor(df[[cn]])))
    } else {
      cat_levels = ref$cat_levels
    }
    Xc = matrix(0L, nrow = n, ncol = length(cat_cols))
    for (j in seq_along(cat_cols)) {
      x = df[[cat_cols[j]]]
      lv = cat_levels[[j]]
      if (is.null(ref)) {
        code = as.integer(factor(x, levels = lv))
        code[is.na(code)] = 0L
      } else {
        m = match(as.character(x), lv)
        code = ifelse(is.na(x), 0L, ifelse(is.na(m), -1L, as.integer(m)))
      }
      Xc[, j] = code
    }
    storage.mode(Xc) = "integer"
  } else {
    Xc = matrix(integer(0), nrow = n, ncol = 0)
    cat_levels = list()
  }

  # ordered -> [0,1] via (code-1)/(L-1); missing or unseen in predict => NA
  if (length(ord_cols)) {
    if (is.null(ref)) {
      ord_levels = lapply(ord_cols, function(cn) levels(as.ordered(df[[cn]])))
    } else {
      ord_levels = ref$ord_levels
    }
    Xo = matrix(NA_real_, nrow = n, ncol = length(ord_cols))
    for (j in seq_along(ord_cols)) {
      x = df[[ord_cols[j]]]
      lv = ord_levels[[j]]
      if (is.null(ref)) {
        code = as.integer(as.ordered(x))
      } else {
        m = match(as.character(x), lv)
        code = ifelse(is.na(x), NA_integer_, as.integer(m))
      }
      L = length(lv)
      if (L <= 1L) {
        Xo[, j] = ifelse(is.na(code), NA_real_, 0)
      } else {
        Xo[, j] = (as.numeric(code) - 1) / (L - 1)
      }
    }
    storage.mode(Xo) = "double"
  } else {
    Xo = matrix(numeric(0), nrow = n, ncol = 0)
    ord_levels = list()
  }

  list(
    Xnum = Xn, Xcat = Xc, Xord = Xo,
    ranges_num = as.numeric(rng),
    cat_levels = cat_levels,
    ord_levels = ord_levels
  )
}

#' Test rows with at least one eligible Gower neighbour
#'
#' Applies the eligibility rule of the native kernels: a training row is an
#' eligible neighbour of a test row when both rows have at least
#' `max(1, ceiling(min_feature_frac * p))` features observed in common, where
#' `p` is the total number of features. An unseen test level of an unordered
#' factor counts as observed (it is compared, as a mismatch).
#'
#' @param enc_train,enc_test (`list()`) Encoded matrices `Xnum`, `Xcat` and
#'   `Xord`, as returned by [knn_gower_encode_blocks()].
#' @param min_feature_frac (`numeric(1)`) Minimum fraction of comparable
#'   features.
#' @param chunk (`integer(1)`) Number of test rows processed at once.
#' @return A logical vector with one element per test row.
#' @noRd
knn_gower_has_neighbors = function(enc_train, enc_test, min_feature_frac, chunk = 1000L) {
  observed = function(e) {
    o = cbind(is.finite(e$Xnum), e$Xcat != 0L, is.finite(e$Xord))
    storage.mode(o) = "double"
    o
  }
  o_train = observed(enc_train)
  o_test = observed(enc_test)
  m = nrow(o_test)
  p = ncol(o_train)
  if (!m) {
    return(logical(0))
  }
  if (!p || !nrow(o_train)) {
    return(rep(FALSE, m))
  }
  min_cnt = max(1, ceiling(max(0, min(1, min_feature_frac)) * p))
  out = logical(m)
  for (start in seq.int(1L, m, by = chunk)) {
    idx = start:min(m, start + chunk - 1L)
    out[idx] = rowSums(tcrossprod(o_test[idx, , drop = FALSE], o_train) >= min_cnt) > 0
  }
  out
}

#' Subset the rows of encoded Gower matrices
#' @noRd
knn_gower_subset_rows = function(enc, rows) {
  list(
    Xnum = enc$Xnum[rows, , drop = FALSE],
    Xcat = enc$Xcat[rows, , drop = FALSE],
    Xord = enc$Xord[rows, , drop = FALSE]
  )
}
